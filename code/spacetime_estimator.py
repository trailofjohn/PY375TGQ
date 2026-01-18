
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from dataclasses import dataclass
from typing import List, Tuple, Optional
from scipy.linalg import block_diag

# --- Part 1: Physics & Math Core ---

# Constants
AU = 1.496e11          # Astronomical Unit (m)
C = 2.99792458e8       # Speed of light (m/s)
MU_SUN = 1.3271244e20  # Sun gravitational parameter (m^3/s^2)
DAY = 86400.0          # Seconds in a day

@dataclass
class SpaceTimeState:
    """
    9-State Vector:
    [x, y, z, vx, vy, vz, bias, drift, drift_rate]
    
    Units:
    Position: m
    Velocity: m/s
    Clock Bias: s
    Clock Drift: s/s (dimensionless)
    Clock Drift Rate: s/s^2
    """
    vec: np.ndarray # Shape (9,)

    @property
    def pos(self) -> np.ndarray: return self.vec[0:3]
    @property
    def vel(self) -> np.ndarray: return self.vec[3:6]
    @property
    def clock(self) -> np.ndarray: return self.vec[6:9] 
    
    @property
    def bias(self) -> float: return self.vec[6]
    @property
    def drift(self) -> float: return self.vec[7]
    @property
    def drift_rate(self) -> float: return self.vec[8]

class DynamicsModel:
    def __init__(self, process_noise_cov: np.ndarray):
        """
        Args:
            process_noise_cov (Q): Covariance matrix of process noise.
        """
        self.Q = process_noise_cov

    def gravity(self, pos: np.ndarray) -> np.ndarray:
        """Newtonian gravity acceleration."""
        r = np.linalg.norm(pos)
        return -MU_SUN * pos / (r**3)

    def propagate(self, state: SpaceTimeState, dt: float) -> SpaceTimeState:
        """RK4 Integration for Orbit, Analytic for Clock."""
        # Unpack
        r0 = state.pos
        v0 = state.vel
        clk0 = state.clock

        # 1. Integrate Orbit (RK4)
        k1_v = self.gravity(r0)
        k1_r = v0

        k2_v = self.gravity(r0 + 0.5 * dt * k1_r)
        k2_r = v0 + 0.5 * dt * k1_v

        k3_v = self.gravity(r0 + 0.5 * dt * k2_r)
        k3_r = v0 + 0.5 * dt * k2_v

        k4_v = self.gravity(r0 + dt * k3_r)
        k4_r = v0 + dt * k3_v

        pos_new = r0 + (dt / 6.0) * (k1_r + 2*k2_r + 2*k3_r + k4_r)
        vel_new = v0 + (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)

        # 2. Integrate Clock (Analytic Linear Matrix)
        # B_k+1 = B_k + D_k*dt + 0.5*DR_k*dt^2
        # D_k+1 = D_k + DR_k*dt
        # DR_k+1 = DR_k
        F_clk = np.array([
            [1, dt, 0.5*dt**2],
            [0, 1, dt],
            [0, 0, 1]
        ])
        clk_new = F_clk @ clk0

        # Inject Process Noise (simulated reality only)
        # Note: In a real filter we don't add noise during prediction, 
        # but for the "True" simulation track we do. 
        # This method is used for both Filter Prediction (no noise added here usually, Q applied to P)
        # and Truth Propagation (noise added externally or handled by a separate simulate method).
        # We will assume this is the Deterministic Propagator.
        
        return SpaceTimeState(np.concatenate([pos_new, vel_new, clk_new]))

    def jacobian(self, state: SpaceTimeState, dt: float) -> np.ndarray:
        """
        Computes the Jacobian Matrix F (Phi) for the state transition.
        F = I + A*dt (approx) or exact analytical blocks.
        """
        r = state.pos
        rn = np.linalg.norm(r)
        
        # Gravity Gradient Tensor G
        # G = -mu/r^3 * (I - 3*r*r.T/r^2)
        I3 = np.eye(3)
        rrT = np.outer(r, r)
        G = (-MU_SUN / rn**3) * (I3 - 3 * rrT / rn**2)

        # State Transition Matrix F (9x9)
        # Structure:
        # [ I  I*dt  0 ]
        # [ G*dt I   0 ]
        # [ 0  0   F_clk ]
        
        F = np.eye(9)
        
        # Position affected by Velocity
        F[0:3, 3:6] = I3 * dt
        
        # Velocity affected by Gravity (Position)
        F[3:6, 0:3] = G * dt # First order approx of integral G dt
        
        # Clock Block
        F[6:9, 6:9] = np.array([
            [1, dt, 0.5*dt**2],
            [0, 1, dt],
            [0, 0, 1]
        ])
        
        return F

class MeasurementModel:
    def __init__(self, measurement_noise_cov: np.ndarray):
        """
        Args:
            measurement_noise_cov (R): Covariance.
        """
        self.R = measurement_noise_cov

    def shapiro_delay(self, r_earth: np.ndarray, r_probe: np.ndarray) -> float:
        """
        Calculates One-Way Shapiro Delay (in meters).
        Delta_rho = (2*GM/c^2) * ln( (re + rp + rep) / (re + rp - rep) )
        """
        re = np.linalg.norm(r_earth)
        rp = np.linalg.norm(r_probe)
        rep = np.linalg.norm(r_earth - r_probe)
        
        # Schwarzschild radius of Sun approx 3km
        Rs = 2 * MU_SUN / (C**2)
        
        # Log term
        num = re + rp + rep
        den = re + rp - rep
        
        # Protect against singularity if probe is exactly on line behind sun
        if den < 1e-3: 
            return 0.0 # Occulted/Singularity
            
        val = (Rs / 2.0) * np.log(num / den) 
        # Note: The factor is gamma+1. GR gamma=1. So (1+1)/2 * Rs * ln(...) ?
        # Standard one-way form is 2GM/c^2 * ln(...).
        
        return val * C # Return distance equivalent (meters) to match range

    def predict(self, state: SpaceTimeState, r_earth: np.ndarray) -> float:
        """
        Predicts the measurement z (Range).
        z = GeometricRange + ShapiroDelay + ClockBiasDistance
        """
        r_probe = state.pos
        geo_range = np.linalg.norm(r_probe - r_earth)
        shapiro = self.shapiro_delay(r_earth, r_probe)
        clock_dist = state.bias * C
        
        return geo_range + shapiro + clock_dist

    def jacobian(self, state: SpaceTimeState, r_earth: np.ndarray) -> np.ndarray:
        """
        H matrix (1x9). Linearization of the measurement equation.
        """
        r_p = state.pos
        re_vec = r_earth
        
        re = np.linalg.norm(re_vec)
        rp = np.linalg.norm(r_p)
        rho_vec = r_p - re_vec
        rep = np.linalg.norm(rho_vec) # range rho
        
        # 1. Geometric Term Jacobian: d(rho)/d(r_p) = (r_p - r_e)^T / rho
        H_geo = (r_p - re_vec).T / rep
        
        # 2. Shapiro Term Jacobian (Analytical)
        # Gamma = 2GM/c^2
        Gamma = 2 * MU_SUN / (C**2) 
        # f = ln(N/D); N = re + rp + rep; D = re + rp - rep
        # d(Shapiro)/d(vec_rp) = (Gamma/2) * ( (1/N)*dN - (1/D)*dD )
        # dN/d(vec_rp) = d(rp)/d(vec_rp) + d(rep)/d(vec_rp) = r_p/rp + (r_p - r_e)/rep
        # dD/d(vec_rp) = r_p/rp - (r_p - r_e)/rep
        
        num = re + rp + rep
        den = re + rp - rep
        
        if den < 1e-3:
            H_shapiro = np.zeros(3)
        else:
            drp = r_p / rp
            drep = (r_p - re_vec) / rep
            
            dN = drp + drep
            dD = drp - drep
            
            # The Shapiro delay function returns METERS in my impl, so we multiply by C?
            # Wait, standard formula 2GM/c^2 has units of length.
            # So Gamma is length.
            
            term = (1.0/num)*dN - (1.0/den)*dD
            H_shapiro = (Gamma) * term * C # The factor C is because my predict returns distance? 
            # Re-checking predict: self.shapiro_delay returns meters. 
            # Equation: Delta t = (2GM/c^3) * ln(...).
            # Delta rho = c * Delta t = (2GM/c^2) * ln(...)
            # So the unit of Gamma is Meters. No extra C.
            # My predict function implemented: (Rs/2) * ln... * C.
            # Rs is 2GM/c^2 (meters).
            # So I should define my predict consistent with units.
            # Let's align:
            # Shapiro Delay in DISTANCE (meters) = (2GM/c^2) * ln(...)
            # So Gamma = (2GM/c^2) meters.
            # My H_shapiro should just be Gamma * term.
            
            H_shapiro = Gamma * term

        # 3. Clock Term Jacobian: d(c*b)/db = c
        H_clk = np.array([C, 0, 0])
        
        # Assemble H [1x9]
        # [ H_geo_x + H_shap_x, H_geo_y..., H_geo_z..., 0, 0, 0, C, 0, 0]
        
        H = np.zeros((1, 9))
        H[0, 0:3] = H_geo + H_shapiro
        H[0, 6:9] = H_clk
        
        return H

class ExtendedKalmanFilter:
    def __init__(self, x0: SpaceTimeState, P0: np.ndarray, dynamics: DynamicsModel, measurement: MeasurementModel):
        self.state = x0
        self.P = P0
        self.dynamics = dynamics
        self.measurement = measurement

    def predict(self, dt: float):
        # 1. State Prediction (Non-linear)
        self.state = self.dynamics.propagate(self.state, dt)
        
        # 2. Covariance Prediction (Linearized)
        F = self.dynamics.jacobian(self.state, dt)
        Q = self.dynamics.Q
        
        self.P = F @ self.P @ F.T + Q

    def update(self, z_measured: float, r_earth: np.ndarray):
        # 1. Predicted Measurement
        z_pred = self.measurement.predict(self.state, r_earth)
        
        # 2. Measurement Residual
        y = z_measured - z_pred
        
        # 3. Jacobian H
        H = self.measurement.jacobian(self.state, r_earth)
        
        # 4. Kalman Gain
        R = self.measurement.R
        S = H @ self.P @ H.T + R
        K = self.P @ H.T @ np.linalg.inv(S)
        
        # 5. State Update
        dx = K @ np.array([y]) # (9,)
        
        # Update raw vector
        self.state.vec = self.state.vec + dx.flatten()
        
        # 6. Covariance Update
        # P = (I - KH)P
        I = np.eye(9)
        self.P = (I - K @ H) @ self.P

# --- Part 2: Simulation Scenarios & Main ---

def run_simulation():
    # Setup Earth Orbit (Simplified: Circular at 1 AU, angular velocity omega)
    omega_earth = np.sqrt(MU_SUN / AU**3)
    
    def get_earth_pos(t):
        theta = omega_earth * t
        return np.array([AU * np.cos(theta), AU * np.sin(theta), 0.0])

    # --- Scenario Configuration ---
    # Select Scenario via string or uncomment
    SCENARIO = "C" # A=Mars, B=Jupiter, C=SolarGrazer
    
    print(f"Initializing Spacetime Estimator - Scenario {SCENARIO}")
    
    if SCENARIO == "A": # Mars Transfer
        name = "Martian Standard"
        tf = 250 * DAY 
        dt = 1 * DAY
        # Hohmann Transfer initial condition
        r1 = AU
        r2 = 1.52 * AU
        # Vis-Viva Eq at Perihelion: v = sqrt(mu * (2/r - 1/a))
        a_transfer = (r1 + r2) / 2
        v_dep = np.sqrt(MU_SUN * (2/r1 - 1/a_transfer))
        
        # Start at Earth's position (x=AU, y=0) ? No, Earth moves.
        # Let's start at x=AU, y=0. Earth starts there too.
        # Velocity is purely tangential y-direction.
        init_pos = np.array([r1, 0, 0])
        init_vel = np.array([0, v_dep, 0])
        
    elif SCENARIO == "B": # Jovian Deep-Dive
        name = "Jovian Deep-Dive"
        tf = 800 * DAY
        dt = 2 * DAY
        r1 = AU
        r2 = 5.2 * AU
        a_transfer = (r1 + r2) / 2
        v_dep = np.sqrt(MU_SUN * (2/r1 - 1/a_transfer))
        init_pos = np.array([r1, 0, 0])
        init_vel = np.array([0, v_dep, 0])
        
    elif SCENARIO == "C": # Solar Grazer
        name = "Solar Grazer"
        tf = 150 * DAY
        dt = 0.5 * DAY # finer steps for fast dynamics
        # Ellipse with perihelion very close to sun
        rp = 0.05 * AU # 0.05 AU perihelion (inside Mercury)
        ra = 1.0 * AU  # Aphelion at Earth
        a_orbit = (rp + ra) / 2
        # Start at Aphelion (Earth) going to Perihelion
        # Velocity at Aphelion
        v_ap = np.sqrt(MU_SUN * (2/ra - 1/a_orbit))
        
        init_pos = np.array([ra, 0, 0])
        init_vel = np.array([0, v_ap, 0]) # Tangential

    # Initial State (Truth)
    # Add random clock initialization
    true_bias = 1e-5 # 10 us
    true_drift = 1e-10 
    true_driftrate = 1e-15
    
    true_state = SpaceTimeState(np.concatenate([
        init_pos, 
        init_vel, 
        [true_bias, true_drift, true_driftrate]
    ]))
    
    # Uncertainty / Noise
    sig_pos = 1000.0 # 1km initial uncertainty
    sig_vel = 0.1 # 10 cm/s
    sig_clk_b = 100e-9 # 100 ns
    
    # Process Noise Q
    # Clock noise is random walk in driftrate (integral -> drift -> integral -> bias)
    # Simple diagonal approx for now
    q_pos = 1.0 # m
    q_clock_dr = 1e-18
    Q = np.diag([
        q_pos**2, q_pos**2, q_pos**2, 
        0.001**2, 0.001**2, 0.001**2,
        0, 1e-20, q_clock_dr**2 # only drift rate is driven noise usually, or drift
    ])
    
    # Measurement Noise R
    sig_range = 10.0 # 10m ranging error
    R = np.array([[sig_range**2]])
    
    dynamics = DynamicsModel(Q)
    meas_model = MeasurementModel(R)
    
    # Initialize Filter State (Perturbed)
    pert_pos = init_pos + np.random.randn(3)*sig_pos
    pert_vel = init_vel + np.random.randn(3)*sig_vel
    pert_clk = np.array([0, 0, 0]) # Start blind on clock
    
    est_state = SpaceTimeState(np.concatenate([pert_pos, pert_vel, pert_clk]))
    P0 = np.diag([
        sig_pos**2, sig_pos**2, sig_pos**2,
        sig_vel**2, sig_vel**2, sig_vel**2,
        1e-3**2, 1e-9**2, 1e-12**2 # High uncertainty on clock
    ])
    
    ekf = ExtendedKalmanFilter(est_state, P0, dynamics, meas_model)
    
    # History Storage
    history = {
        't': [],
        'true_pos': [],
        'est_pos': [],
        'dead_reckoning': [],
        'shapiro': [],
        'clock_err': [],
        'P_trace': []
    }
    
    # Dead Reckoning State (for visual comparison - no updates)
    dr_state = SpaceTimeState(est_state.vec.copy())
    
    current_time = 0.0
    steps = int(tf / dt)
    
    print(f"Simulating {steps} steps...")
    
    for i in range(steps):
        # 1. Propagate Truth
        # Add random walk noise to clock drift rate
        true_state.vec[8] += np.random.randn() * 1e-16 # process noise
        true_state = dynamics.propagate(true_state, dt)
        
        # 2. Propagate Dead Reckoning (No Update)
        dr_state = dynamics.propagate(dr_state, dt)
        
        # 3. Filter Prediction
        ekf.predict(dt)
        
        # 4. Filter Update (Measurement)
        # Generate noisy measurement
        r_earth = get_earth_pos(current_time)
        z_true = meas_model.predict(true_state, r_earth)
        z_noise = np.random.randn() * sig_range
        z_meas = z_true + z_noise
        
        ekf.update(z_meas, r_earth)
        
        # Store Data
        history['t'].append(current_time)
        history['true_pos'].append(true_state.pos)
        history['est_pos'].append(ekf.state.pos)
        history['dead_reckoning'].append(dr_state.pos)
        
        # Calc Shapiro for viz
        shap = meas_model.shapiro_delay(r_earth, true_state.pos)
        history['shapiro'].append(shap)
        
        # Clock residual (convert bias time to distance error for easier viz context or just ns)
        clk_err = (ekf.state.bias - true_state.bias) * C # in meters
        history['clock_err'].append(clk_err)
        history['P_trace'].append(np.trace(ekf.P))
        
        current_time += dt

    # Convert to arrays
    for k in history:
        history[k] = np.array(history[k])
        
    return history, name

def run_visualization(hist, title):
    print("Generating Animation...")
    fig = plt.figure(figsize=(16, 9))
    
    # Layout: 1 large 3D plot, 2 smaller subplots
    gs = fig.add_gridspec(2, 3)
    ax3d = fig.add_subplot(gs[:, 0:2], projection='3d')
    ax_clk = fig.add_subplot(gs[0, 2])
    ax_shap = fig.add_subplot(gs[1, 2])
    
    # --- 3D Plot Initial Setup ---
    ax3d.set_title(f"Mission: {title} (EKF Tracking)")
    ax3d.set_xlabel('X (AU)')
    ax3d.set_ylabel('Y (AU)')
    ax3d.set_zlabel('Z (AU)')
    
    # Plot Sun
    ax3d.scatter([0], [0], [0], color='yellow', s=200, label='Sun', edgecolors='orange')
    
    # Lines
    line_true, = ax3d.plot([], [], [], 'g-', label='True Traj', linewidth=1)
    line_dr, = ax3d.plot([], [], [], 'r--', label='Dead Reckoning', alpha=0.5)
    line_est, = ax3d.plot([], [], [], 'b:', label='EKF Est')
    point_est, = ax3d.plot([], [], [], 'bo', markersize=4)
    
    # Correction Spring/Line
    spring_line, = ax3d.plot([], [], [], 'm-', linewidth=1, alpha=0.7) # Magenta line connecting Est to True
    
    ax3d.legend()
    
    # Normalize scales for nice view
    max_range = np.max(np.abs(hist['true_pos'])) / AU
    ax3d.set_xlim(-max_range, max_range)
    ax3d.set_ylim(-max_range, max_range)
    ax3d.set_zlim(-max_range/2, max_range/2)
    
    # --- Dashboard Setup ---
    # Clock Error
    ax_clk.set_title("Clock Bias Residual (m)")
    ax_clk.set_xlabel("Time (days)")
    ax_clk.grid(True)
    line_clk, = ax_clk.plot([], [], 'k-', lw=1)
    
    # Shapiro
    ax_shap.set_title("Relativistic Shapiro Delay Magnitude (m)")
    ax_shap.set_xlabel("Time (days)")
    ax_shap.grid(True)
    line_shap, = ax_shap.plot([], [], 'purple', lw=1.5)
    
    # Animation Function
    days = hist['t'] / DAY
    
    def update(frame):
        # Downsample for speed if needed, but lets try 1:1 or 1:5
        idx = frame * 5
        if idx >= len(days): idx = len(days) - 1
        
        # Update 3D Paths (history up to idx)
        ax3d.plot(hist['true_pos'][:idx, 0]/AU, hist['true_pos'][:idx, 1]/AU, hist['true_pos'][:idx, 2]/AU, 'g-', lw=0.5)
        # Just update data for lines usually to avoid accumulations, but here we redraw full path?
        # Better: set_data
        
        line_true.set_data(hist['true_pos'][:idx, 0]/AU, hist['true_pos'][:idx, 1]/AU)
        line_true.set_3d_properties(hist['true_pos'][:idx, 2]/AU)
        
        line_dr.set_data(hist['dead_reckoning'][:idx, 0]/AU, hist['dead_reckoning'][:idx, 1]/AU)
        line_dr.set_3d_properties(hist['dead_reckoning'][:idx, 2]/AU)
        
        line_est.set_data(hist['est_pos'][:idx, 0]/AU, hist['est_pos'][:idx, 1]/AU)
        line_est.set_3d_properties(hist['est_pos'][:idx, 2]/AU)
        
        # Current Head
        point_est.set_data([hist['est_pos'][idx, 0]/AU], [hist['est_pos'][idx, 1]/AU])
        point_est.set_3d_properties([hist['est_pos'][idx, 2]/AU])
        
        # Spring (Error Line)
        # Connect Est[idx] to True[idx]
        sx = [hist['est_pos'][idx, 0]/AU, hist['true_pos'][idx, 0]/AU]
        sy = [hist['est_pos'][idx, 1]/AU, hist['true_pos'][idx, 1]/AU]
        sz = [hist['est_pos'][idx, 2]/AU, hist['true_pos'][idx, 2]/AU]
        spring_line.set_data(sx, sy)
        spring_line.set_3d_properties(sz)
        
        # Dashboard
        line_clk.set_data(days[:idx], hist['clock_err'][:idx])
        ax_clk.set_xlim(0, max(days[:idx]) + 1)
        ax_clk.set_ylim(min(hist['clock_err']), max(hist['clock_err']) + 1)
        
        line_shap.set_data(days[:idx], hist['shapiro'][:idx])
        ax_shap.set_xlim(0, max(days[:idx]) + 1)
        ax_shap.set_ylim(0, max(hist['shapiro']) * 1.1)
        
        return line_true, line_dr, line_est, point_est, spring_line, line_clk, line_shap

    ani = animation.FuncAnimation(fig, update, frames=len(days)//5, interval=20, blit=False)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    history_data, scenario_name = run_simulation()
    run_visualization(history_data, scenario_name)
