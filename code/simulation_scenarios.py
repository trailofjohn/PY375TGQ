
import numpy as np
from simulation_core import SpaceTimeState, MU_SUN, AU, DAY

def get_scenario_config(scenario_code: str):
    omega_earth = np.sqrt(MU_SUN / AU**3)
    
    config = {
        'omega_earth': omega_earth,
        'scenario_code': scenario_code
    }
    
    if scenario_code == "A": # Mars Transfer
        config['name'] = "Martian Standard"
        config['tf'] = 250 * DAY 
        config['dt'] = 1 * DAY
        r1 = AU
        r2 = 1.52 * AU
        a_transfer = (r1 + r2) / 2
        v_dep = np.sqrt(MU_SUN * (2/r1 - 1/a_transfer))
        config['init_pos'] = np.array([r1, 0, 0])
        config['init_vel'] = np.array([0, v_dep, 0])
        
    elif scenario_code == "B": # Jovian Deep-Dive
        config['name'] = "Jovian Deep-Dive"
        config['tf'] = 800 * DAY
        config['dt'] = 2 * DAY
        r1 = AU
        r2 = 5.2 * AU
        a_transfer = (r1 + r2) / 2
        v_dep = np.sqrt(MU_SUN * (2/r1 - 1/a_transfer))
        config['init_pos'] = np.array([r1, 0, 0])
        config['init_vel'] = np.array([0, v_dep, 0])
        
    elif scenario_code == "C": # Solar Grazer
        config['name'] = "Solar Grazer"
        config['tf'] = 150 * DAY
        config['dt'] = 0.5 * DAY 
        rp = 0.05 * AU 
        ra = 1.0 * AU 
        a_orbit = (rp + ra) / 2
        v_ap = np.sqrt(MU_SUN * (2/ra - 1/a_orbit))
        config['init_pos'] = np.array([ra, 0, 0])
        config['init_vel'] = np.array([0, v_ap, 0])

    return config

def get_initial_states(config):
    # Initial State (Truth)
    true_bias = 1e-5 
    true_drift = 1e-10 
    true_driftrate = 1e-15
    true_tau = 0.0 # Proper time starts at 0
    
    true_state = SpaceTimeState(np.concatenate([
        config['init_pos'], 
        config['init_vel'], 
        [true_bias, true_drift, true_driftrate, true_tau]
    ]))
    
    # Process Noise Q (10x10)
    # [Pos(3), Vel(3), Bias(1), Drift(1), DR(1), Tau(1)]
    q_pos = 1.0 
    q_clock_dr = 1e-18
    # Tau noise? Tau is integral of physics, usually low process noise unless modeling unknown G-potential variance.
    q_tau = 1e-20 
    
    diag_vals = [
        q_pos**2, q_pos**2, q_pos**2, 
        0.001**2, 0.001**2, 0.001**2,
        0, 1e-20, q_clock_dr**2, 
        q_tau**2
    ]
    Q = np.diag(diag_vals)
    
    # Measurement Noise R (2x2)
    # [Range (m), Doppler (m/s)]
    # RELAXED DOPPLER: 1.0 m/s (was 0.001 m/s) to prevent numerical instability with Clock coupling
    sig_range = 10.0 
    sig_doppler = 0.5 # Relaxed from 0.001 to 0.5 m/s
    
    R = np.diag([sig_range**2, sig_doppler**2])
    
    # Initialize Filter State (Perturbed)
    sig_pos = 1000.0 
    sig_vel = 0.1 
    
    pert_pos = config['init_pos'] + np.random.randn(3)*sig_pos
    pert_vel = config['init_vel'] + np.random.randn(3)*sig_vel
    pert_clk = np.array([0, 0, 0]) 
    pert_tau = 0.0
    
    est_state = SpaceTimeState(np.concatenate([pert_pos, pert_vel, pert_clk, [pert_tau]]))
    
    # P0 (10x10)
    P0_diag = [
        sig_pos**2, sig_pos**2, sig_pos**2,
        sig_vel**2, sig_vel**2, sig_vel**2,
        1e-5**2, 1e-11**2, 1e-14**2, # Tighter Clock Initial Covariance (was 1e-9 drift)
        1.0**2 # Low uncertainty on start time
    ]
    P0 = np.diag(P0_diag)
    
    # Return explicit signatures
    return true_state, est_state, P0, Q, R, sig_range, sig_doppler
