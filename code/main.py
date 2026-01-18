
import numpy as np
import sys
import multiprocessing
from simulation_core import DynamicsModel, MeasurementModel, ExtendedKalmanFilter, SpaceTimeState, C, DAY, AU
from simulation_scenarios import get_scenario_config, get_initial_states
from visualization import run_visualization

def run_simulation(scenario_code):
    """
    Runs the simulation logic and returns the history data.
    """
    config = get_scenario_config(scenario_code)
    print(f"[Process {scenario_code}] Initializing: {config['name']}")
    
    true_state, est_state, P0, Q, R, sig_range, sig_doppler = get_initial_states(config)
    
    dynamics = DynamicsModel(Q)
    meas_model = MeasurementModel(R)
    
    ekf = ExtendedKalmanFilter(est_state, P0, dynamics, meas_model)
    
    history = {
        't': [],
        'true_pos': [],
        'est_pos': [],
        'dead_reckoning': [],
        'shapiro': [],
        'clock_err': [],
        'P_trace': [],
        'proper_time': [],
        'doppler_err': []
    }
    
    dr_state = SpaceTimeState(est_state.vec.copy())
    
    omega_earth = config['omega_earth']
    def get_earth_pos(t):
        theta = omega_earth * t
        return np.array([AU * np.cos(theta), AU * np.sin(theta), 0.0])
    
    current_time = 0.0
    dt = config['dt']
    steps = int(config['tf'] / dt)
    
    print(f"[Process {scenario_code}] Simulating {steps} steps...")
    
    for i in range(steps):
        # 1. Propagate Truth
        # Index 8 is Drift Rate in new vector [0..9]
        # [x, y, z, vx, vy, vz, b, d, dr, tau]
        true_state.vec[8] += np.random.randn() * 1e-16 
        true_state = dynamics.propagate(true_state, dt)
        
        # 2. Propagate Dead Reckoning
        dr_state = dynamics.propagate(dr_state, dt)
        
        # 3. Filter Prediction
        ekf.predict(dt)
        
        # 4. Filter Update
        r_earth = get_earth_pos(current_time)
        z_true = meas_model.predict(true_state, r_earth) # Returns [Range, Doppler]
        
        # Noise (2D)
        # sig_range and sig_doppler are needed. Update unpacking.
        # R was returned as diag matrix, so we can use R diagonals or the scalar returned
        # Logic update needed in line 16.
        # Let's fix line 16 unpacking below this block first? No, replace contents.
        # We need sig_doppler in main scope.
        # I'll rely on z_noise construction.
        
        n_range = np.random.randn() * sig_range
        n_doppler = np.random.randn() * sig_doppler 
        
        z_meas = z_true + np.array([n_range, n_doppler])
        
        ekf.update(z_meas, r_earth)
        
        # Store Data
        history['t'].append(current_time)
        history['true_pos'].append(true_state.pos)
        history['est_pos'].append(ekf.state.pos)
        history['dead_reckoning'].append(dr_state.pos)
        history['shapiro'].append(meas_model.shapiro_delay(r_earth, true_state.pos))
        history['clock_err'].append((ekf.state.bias - true_state.bias) * C)
        history['P_trace'].append(np.trace(ekf.P))
        
        # New History?
        # Maybe store Proper Time difference?
        history['proper_time'].append((ekf.state.tau - true_state.tau) * C) # Store proper time offset error
        history['doppler_err'].append(z_meas[1] - meas_model.predict(ekf.state, r_earth)[1]) # Store Doppler measurement residual
        
        current_time += dt

    # Convert to arrays
    for k in history:
        history[k] = np.array(history[k])
        
    return history, config['name']

def simulation_worker(scenario_code):
    """
    Worker function for multiprocessing. 
    Runs simulation AND visualization in this process.
    """
    try:
        hist, name = run_simulation(scenario_code)
        run_visualization(hist, name)
    except Exception as e:
        print(f"Error in process {scenario_code}: {e}")

if __name__ == "__main__":
    # If arguments provided, run specific scenario
    if len(sys.argv) > 1:
        arg = sys.argv[1].upper()
        if arg in ["A", "B", "C"]:
            simulation_worker(arg)
        elif arg == "ALL":
             # Launch all 3
            processes = []
            for sc in ["A", "B", "C"]:
                p = multiprocessing.Process(target=simulation_worker, args=(sc,))
                p.start()
                processes.append(p)
            
            for p in processes:
                p.join()
        else:
            print("Usage: python main.py [A|B|C|ALL]")
    else:
        # Default behavior update: User asked for "Run main -> all 3 popup"
        print("Launching All Scenarios (A, B, C)...")
        processes = []
        for sc in ["A", "B", "C"]:
            p = multiprocessing.Process(target=simulation_worker, args=(sc,))
            p.start()
            processes.append(p)
        
        for p in processes:
            p.join()
