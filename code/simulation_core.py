
import numpy as np
from dataclasses import dataclass

# Constants
AU = 1.496e11          
C = 2.99792458e8       
MU_SUN = 1.3271244e20 
DAY = 86400.0          

@dataclass
class SpaceTimeState:
    """
    10-State Vector (V1):
    [x, y, z, vx, vy, vz, bias, drift, drift_rate, proper_time]
    
    New: proper_time (tau) accumulation.
    """
    vec: np.ndarray # Shape (10,)

    @property
    def pos(self) -> np.ndarray: return self.vec[0:3]
    @property
    def vel(self) -> np.ndarray: return self.vec[3:6]
    @property
    def clock(self) -> np.ndarray: return self.vec[6:9] 
    @property
    def tau(self) -> float: return self.vec[9] # Proper time
    
    @property
    def bias(self) -> float: return self.vec[6]
    @property
    def drift(self) -> float: return self.vec[7]
    @property
    def drift_rate(self) -> float: return self.vec[8]

class DynamicsModel:
    def __init__(self, process_noise_cov: np.ndarray):
        # Q is now 10x10
        self.Q = process_noise_cov

    def gravity(self, pos: np.ndarray) -> np.ndarray:
        r = np.linalg.norm(pos)
        return -MU_SUN * pos / (r**3)

    def proper_time_derivative(self, pos: np.ndarray, vel: np.ndarray) -> float:
        """
        Calculates d(tau)/dt based on GR/SR.
        dtau/dt = sqrt(1 - 2U/c^2 - v^2/c^2)
        Approx: 1 - U/c^2 - v^2/2c^2
        """
        r = np.linalg.norm(pos)
        v_sq = np.dot(vel, vel)
        U = MU_SUN / r # Potential
        
        # dtau/dt approx 1 - (Phi/c^2 + v^2/2c^2)
        # We model the NEGATIVE effect (slowing down)
        correction = -(U/(C**2) + v_sq/(2*C**2))
        return 1.0 + correction

    def propagate(self, state: SpaceTimeState, dt: float) -> SpaceTimeState:
        r0 = state.pos
        v0 = state.vel
        clk0 = state.clock
        tau0 = state.tau

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

        # 2. Integrate Clock (Analytic)
        F_clk = np.array([
            [1, dt, 0.5*dt**2],
            [0, 1, dt],
            [0, 0, 1]
        ])
        clk_new = F_clk @ clk0
        
        # 3. Integrate Proper Time
        # Euler approx sufficient for visualization of the accumulation, or RK4 if precise
        dtau_dt = self.proper_time_derivative(r0, v0)
        tau_new = tau0 + dtau_dt * dt
        
        return SpaceTimeState(np.concatenate([pos_new, vel_new, clk_new, [tau_new]]))

    def jacobian(self, state: SpaceTimeState, dt: float) -> np.ndarray:
        r = state.pos
        rn = np.linalg.norm(r)
        
        I3 = np.eye(3)
        rrT = np.outer(r, r)
        G = (-MU_SUN / rn**3) * (I3 - 3 * rrT / rn**2)

        F = np.eye(10)
        F[0:3, 3:6] = I3 * dt
        F[3:6, 0:3] = G * dt 
        F[6:9, 6:9] = np.array([
            [1, dt, 0.5*dt**2],
            [0, 1, dt],
            [0, 0, 1]
        ])
        # Tau does not affect X/V/Bias in this simplified model, 
        # but X/V affects Tau. We ignore d(Tau)/dX for Position Feedback
        # because we don't measure Proper Time directly usually, we infer it.
        # So leaving row 9 as identity [0...0 1] is fine for EKF stability unless we fuse Tau.
        
        return F

class MeasurementModel:
    def __init__(self, measurement_noise_cov: np.ndarray):
        # R is now 2x2 (Range, Doppler)
        self.R = measurement_noise_cov

    def shapiro_delay(self, r_earth: np.ndarray, r_probe: np.ndarray) -> float:
        """
        Calculates One-Way Shapiro Delay (in meters).
        Corrected coeff: 4GM/c^2
        """
        re = np.linalg.norm(r_earth)
        rp = np.linalg.norm(r_probe)
        rep = np.linalg.norm(r_earth - r_probe)
        
        # 4GM/c^2
        K = 4.0 * MU_SUN / (C**2)
        
        num = re + rp + rep
        den = re + rp - rep
        
        if den < 1e-3: 
            return 0.0
            
        val = K * np.log(num / den) 
        return val # meters

    def solve_light_time(self, r_earth_t: np.ndarray, r_probe_t: np.ndarray, vel_probe_t: np.ndarray) -> tuple:
        """
        Iterative Light-Time Solution.
        Finds the position of the probe at (t - light_time).
        Simplified: We assume Earth is at t_recv (fixed), we find probe at t_emit.
        """
        # Initial guess: Instantaneous distance
        dist = np.linalg.norm(r_probe_t - r_earth_t)
        light_time = dist / C
        
        # Iterate (usually 2-3 passes)
        for _ in range(3):
            # Probe pos at t_emit = r(t) - v(t)*lt (Linear approx backward)
            # In simulation_core we strictly only have state at t.
            # This linear back-propagation is standard for LT correction in filters.
            r_emit = r_probe_t - vel_probe_t * light_time
            dist = np.linalg.norm(r_emit - r_earth_t)
            light_time = dist / C
            
        return dist, light_time, r_emit

    def predict(self, state: SpaceTimeState, r_earth: np.ndarray) -> np.ndarray:
        """
        Returns [Range (m), Doppler (m/s)]
        """
        r_current = state.pos
        v_current = state.vel
        
        # 1. Light Time Corrected Range
        geo_range_lt, lt, r_emit = self.solve_light_time(r_earth, r_current, v_current)
        
        # 2. Shapiro
        shapiro = self.shapiro_delay(r_earth, r_emit)
        
        # 3. Clock Bias (Range equivalent)
        clock_dist = state.bias * C
        
        total_range = geo_range_lt + shapiro + clock_dist
        
        # 4. Doppler (Range Rate)
        # geometric range rate = dot(r_rel, v_rel) / |r_rel|
        # Note: Earth velocity is needed for strict Doppler!
        # We will assume Earth is static relative to simulation frame in this specific call
        # OR we need to pass v_earth. 
        # For V1 improvement, let's assume v_earth ~ approximated or pass 0 for relative calculation if not provided.
        # But we need v_earth for accurate doppler.
        # Let's approximate: Range Rate is projection of probe velocity onto Line of Sight
        # (ignoring Earth velocity contribution for this specific modular block if not passed, 
        # but realistically we must change signature. I will update Main to pass v_earth).
        # For now, simplistic:
        
        los_vec = (r_current - r_earth)
        range_val = np.linalg.norm(los_vec)
        if range_val > 0:
            los_unit = los_vec / range_val
            range_rate = np.dot(v_current, los_unit)
        else:
            range_rate = 0.0
            
        # Add clock drift to doppler (drift is dimensionless s/s, * C = m/s)
        doppler = range_rate + state.drift * C
        
        return np.array([total_range, doppler])

    def jacobian(self, state: SpaceTimeState, r_earth: np.ndarray) -> np.ndarray:
        """
        Returns H (2 x 10)
        Row 1: Range Derivatives
        Row 2: Doppler Derivatives
        """
        r_p = state.pos
        v_p = state.vel
        r_e = r_earth
        
        rho_vec = r_p - r_e
        rho = np.linalg.norm(rho_vec)
        if rho < 1e-3: rho = 1.0
        
        # -- Row 1: Range --
        # d(Range)/dPos = Unit Vector
        H_range_pos = rho_vec.T / rho
        
        # Shapiro partials (using current pos approx)
        # ... reusing logic from V0, simplified for conciseness
        re = np.linalg.norm(r_e); rp = np.linalg.norm(r_p); rep = rho
        K = 4.0 * MU_SUN / (C**2)
        num = re + rp + rep; den = re + rp - rep
        if den > 1e-3:
             drp = r_p / rp; drep = rho_vec / rho
             dN = drp + drep; dD = drp - drep
             term = (1.0/num)*dN - (1.0/den)*dD
             H_shap_pos = K * term
        else:
             H_shap_pos = np.zeros(3)

        H_range = np.zeros(10)
        H_range[0:3] = H_range_pos + H_shap_pos
        H_range[6] = C # Bias
        
        # -- Row 2: Doppler --
        # Range Rate = (r . v) / r
        # d(RR)/dv = r / r = Unit Vector (los)
        # d(RR)/dr = v/r - (r.v)r / r^3
        
        u = rho_vec / rho
        v_dot_u = np.dot(v_p, u)
        
        H_dop_vel = u.T
        H_dop_pos = (v_p.T / rho) - (u.T * v_dot_u / rho)
        
        H_doppler = np.zeros(10)
        H_doppler[0:3] = H_dop_pos
        H_doppler[3:6] = H_dop_vel
        H_doppler[7] = C # Drift
        
        # Combine
        H = np.vstack([H_range, H_doppler])
        return H

class ExtendedKalmanFilter:
    def __init__(self, x0: SpaceTimeState, P0: np.ndarray, dynamics: DynamicsModel, measurement: MeasurementModel):
        self.state = x0
        self.P = P0
        self.dynamics = dynamics
        self.measurement = measurement

    def predict(self, dt: float):
        self.state = self.dynamics.propagate(self.state, dt)
        F = self.dynamics.jacobian(self.state, dt)
        Q = self.dynamics.Q
        self.P = F @ self.P @ F.T + Q

    def update(self, z_measured: np.ndarray, r_earth: np.ndarray):
        # z_measured is now (2,) [Range, Doppler]
        z_pred = self.measurement.predict(self.state, r_earth)
        
        y = z_measured - z_pred
        
        H = self.measurement.jacobian(self.state, r_earth)
        R = self.measurement.R
        
        S = H @ self.P @ H.T + R
        K = self.P @ H.T @ np.linalg.inv(S)
        
        dx = K @ y # (10,)
        
        self.state.vec = self.state.vec + dx
        I = np.eye(10)
        self.P = (I - K @ H) @ self.P

