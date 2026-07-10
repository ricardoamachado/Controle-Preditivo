import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt

class StateSpaceMPC:
    def __init__(self, A, B, C, N_p, N_c, Q, R, u_lim=None, du_lim=None):
        """
        State Space MPC with Integral Action (Velocity Formulation using augmented model)
        
        Parameters:
        A, B, C : State space matrices
        N_p     : Prediction horizon
        N_c     : Control horizon
        Q       : State/output tracking cost matrix (p x p)
        R       : Control effort cost matrix (m x m)
        u_lim   : Tuple (u_min, u_max) for absolute control constraints
        du_lim  : Tuple (du_min, du_max) for control rate constraints
        """
        self.A = np.array(A)
        self.B = np.array(B)
        self.C = np.array(C)
        
        self.n = self.A.shape[0]
        self.m = self.B.shape[1]
        self.p = self.C.shape[0]
        
        self.N_p = N_p
        self.N_c = N_c
        
        self.Q = np.array(Q)
        self.R = np.array(R)
        
        self.u_lim = u_lim
        self.du_lim = du_lim
        
        # Build augmented matrices as requested:
        # Aaug = [ A 0; C*A I ]
        # Baug = [B; C*B]
        # Caug = [0 I]
        self.A_aug = np.vstack([
            np.hstack([self.A, np.zeros((self.n, self.p))]),
            np.hstack([self.C @ self.A, np.eye(self.p)])
        ])
        self.B_aug = np.vstack([self.B, self.C @ self.B])
        self.C_aug = np.hstack([np.zeros((self.p, self.n)), np.eye(self.p)])
        
        # Precompute prediction matrices F and Phi
        self._build_prediction_matrices()
        
        # Cost weight matrices
        self.Q_bar = np.kron(np.eye(N_p), self.Q)
        self.R_bar = np.kron(np.eye(N_c), self.R)

    def _build_prediction_matrices(self):
        """Precompute the free response matrix (F) and forced response matrix (Phi)"""
        # F matrix: Y_free = F * x_aug
        F_rows = []
        for i in range(1, self.N_p + 1):
            F_rows.append(self.C_aug @ np.linalg.matrix_power(self.A_aug, i))
        self.F = np.vstack(F_rows)  # Shape: (p*N_p, n+p)
        
        # Phi matrix: Y_forced = Phi * Delta_U
        self.Phi = np.zeros((self.p * self.N_p, self.m * self.N_c))
        for i in range(self.N_p):
            for j in range(min(i + 1, self.N_c)):
                power = i - j
                if power == 0:
                    block = self.C_aug @ self.B_aug
                else:
                    block = self.C_aug @ np.linalg.matrix_power(self.A_aug, power) @ self.B_aug
                self.Phi[i*self.p:(i+1)*self.p, j*self.m:(j+1)*self.m] = block

    def solve(self, x_curr, x_prev, u_prev, ref, disturbance=0.0):
        """
        Solves the MPC optimization problem for the current step.
        
        Parameters:
        x_curr : Current state (n x 1)
        x_prev : Previous state (n x 1) (used to compute Delta x)
        u_prev : Previous control input (m x 1)
        ref    : Reference trajectory (p x 1) - assumed constant over horizon
        
        Returns:
        u_next : Next control input (m x 1)
        """
        x_curr = np.array(x_curr).reshape(-1, 1)
        x_prev = np.array(x_prev).reshape(-1, 1)
        u_prev = np.array(u_prev).reshape(-1, 1)
        
        # Compute Delta x and current output y
        delta_x = x_curr - x_prev
        y_curr = self.C @ x_curr + disturbance
        
        # Augmented state: X_aug = [Delta_x; y_curr]
        x_aug = np.vstack([delta_x, y_curr])
        
        # Reference over the prediction horizon
        ref_traj = np.tile(np.array(ref).reshape(-1, 1), (self.N_p, 1))
        
        # Define optimization variable (Change in control)
        dU = cp.Variable((self.m * self.N_c, 1))
        
        # Predicted output
        Y_pred = self.F @ x_aug + self.Phi @ dU
        
        # Cost function
        cost = cp.quad_form(Y_pred - ref_traj, cp.psd_wrap(self.Q_bar)) + \
               cp.quad_form(dU, cp.psd_wrap(self.R_bar))
               
        # Constraints
        constraints = []
        if self.du_lim is not None:
            constraints.append(dU >= self.du_lim[0])
            constraints.append(dU <= self.du_lim[1])
            
        if self.u_lim is not None:
            # Predict absolute control inputs: U = u_prev + cumsum(dU)
            U_pred = np.kron(np.ones((self.N_c, 1)), u_prev) + \
                     np.tril(np.ones((self.N_c, self.N_c))) @ dU
            constraints.append(U_pred >= self.u_lim[0])
            constraints.append(U_pred <= self.u_lim[1])
            
        # Solve QP
        prob = cp.Problem(cp.Minimize(cost), constraints)
        prob.solve(solver=cp.OSQP) # OSQP is fast for MPC QPs
        
        if prob.status not in ["optimal", "optimal_inaccurate"]:
            print(f"Warning: Solver status is {prob.status}")
            return u_prev # Fallback to previous input
            
        # Apply first control action
        dU_opt = dU.value
        u_next = u_prev + dU_opt[:self.m]
        
        return u_next

# ==========================================
# Simulation Example
# ==========================================
if __name__ == "__main__":
    # Discrete State Space Model (Sample time = 0.1s)
    # A simple 2nd order system (e.g., position and velocity)
    dt = 1/25000
    A = [[0.97672, -0.011917], [3.660966, 0.900450]]
    B = [[0.606945], [-0.713054]]
    C = [[0, 1]]
    
    # MPC Parameters
    N_p = 50  # Prediction horizon
    N_c = 5   # Control horizon
    Q = [[1.0]] # Output tracking weight
    R = [[500.0]]  # Control effort weight
    
    u_lim = (0, 1)  # Absolute input limits
    du_lim = (-np.inf, np.inf) # Input rate limits
    
    mpc = StateSpaceMPC(A, B, C, N_p, N_c, Q, R, u_lim, du_lim)
    
    # Simulation parameters
    sim_steps = 200
    x_prev = np.array([[0.25], [12]])    # Initial previous state
    x_curr = np.array([[0.25], [12]])    # Initial current state
    u_prev = np.array([[0.25]])           # Initial control input
    ref = 25.0                            # Setpoint (target)
    
    # Unmodeled constant disturbance (to prove integral action works)
    disturbance = 0.2
    
    history_u = []
    history_y = []
    history_delta_u = []

    for k in range(sim_steps):
        # Solve MPC
        if k > sim_steps //2:
            disturbance = 0.2 * ref
        else:
            disturbance = 0
        u = mpc.solve(x_curr, x_prev, u_prev, ref, disturbance)
        history_delta_u.append(u[0, 0] - u_prev[0,0])
        u_prev = u
        # Store current state before advancing simulation
        x_prev = x_curr
        x_next = A @ x_curr + B @ u
        # Advance state
        x_curr = x_next
        y_curr = C @ x_curr + disturbance
        history_y.append(y_curr[0, 0])
        history_u.append(u[0, 0])
        
    # Plotting results
    t = np.arange(sim_steps) * dt
    plt.figure(figsize=(14, 8))
    
    plt.subplot(3, 1, 1)
    plt.step(t, history_y, label='Output (y)')
    plt.axhline(ref, color='r', linestyle='--', label='Reference')
    plt.title('MPC with Integral Action (Augmented Model: [dx; y])')
    plt.ylabel('Voltage (V)')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(3, 1, 2)
    plt.step(t, history_u, label='Control Input (u)')
    plt.ylabel('Input')
    plt.xlabel('Time (s)')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(3, 1, 3)
    plt.step(t, history_delta_u, label='Control Incremenet (Δu)')
    plt.ylabel('Input')
    plt.xlabel('Time (s)')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()