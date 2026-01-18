# V2 Roadmap: Improvements & Future Work

This document outlines planned enhancements for the Spacetime State Estimator moving from V1 to V2+.

---

## High Priority: JPL-Grade Ephemerides

### 1. Integrate JPL HORIZONS / SPICE Kernels
**Current State**: Earth's orbit is modeled as a simple circular orbit at 1 AU.

**Improvement**: Replace the analytical Earth position function with actual ephemeris data.

| Option | Description | Effort |
|:---|:---|:---|
| **SPICE + SpiceyPy** | Use NASA's SPICE toolkit via the `spiceypy` Python wrapper. Provides exact positions for all major bodies. | Medium |
| **JPL HORIZONS API** | Query the web API for pre-computed positions for specific mission dates. | Low |
| **DE440 Binary Kernels** | Load the latest Development Ephemeris files directly for maximum precision (sub-meter). | High |

**Impact**: Eliminates the ~0.5% error in Earth's position and allows simulation of real mission timelines (e.g., "Launch on 2028-10-01").

---

### 2. Add Multiple Gravitating Bodies
**Current State**: Only the Sun's gravity is modeled.

**Improvement**: Implement an N-body integrator or at least add Jupiter's perturbation for the Jovian scenario.

```python
# Example structure for multi-body gravity
def gravity_nbody(pos, body_positions, body_mus):
    accel = np.zeros(3)
    for r_body, mu_body in zip(body_positions, body_mus):
        r_rel = pos - r_body
        accel -= mu_body * r_rel / np.linalg.norm(r_rel)**3
    return accel
```

**Impact**: Critical for accurate trajectory prediction beyond Mars. Jupiter alone can perturb a trajectory by thousands of km.

---

## Medium Priority: Filter Enhancements

### 3. Implement Unscented Kalman Filter (UKF)
**Current State**: Extended Kalman Filter (EKF) uses first-order Taylor expansion (Jacobians).

**Improvement**: UKF captures nonlinearities better using sigma-point sampling.

**Impact**: More accurate uncertainty propagation, especially for the highly eccentric Solar Grazer scenario where the linearization assumption breaks down.

---

### 4. Add Ionospheric / Tropospheric Delay Model
**Current State**: Signal delay is purely geometric + Shapiro.

**Improvement**: Model the ~10-meter delay caused by Earth's ionosphere and troposphere (frequency-dependent).

**Impact**: Closer to reality for sub-meter accuracy. Use dual-frequency ranging to cancel ionospheric effects.

---

### 5. Range-Rate (Doppler) from Two-Way Tracking
**Current State**: Doppler is modeled as one-way (probe velocity only).

**Improvement**: Implement the full two-way Doppler equation:
$$\delta f = f_0 \frac{v_{probe} + v_{Earth}}{c}$$

This requires tracking both uplink and downlink frequencies.

---

## Low Priority: Visualization & UX

### 6. Save Animation as Video File
**Current State**: Animation only plays in a live Matplotlib window.

**Improvement**: Add a `--save` flag to export the animation as MP4/GIF.

```python
ani.save('martian_standard.mp4', writer='ffmpeg', fps=30, dpi=150)
```

---

### 7. Add Earth Tracking Station
**Current State**: Earth is a point mass.

**Improvement**: Visualize Earth's rotation and which tracking station (Goldstone, Madrid, Canberra) is in view. Add DSN downtime windows (no measurements during hand-off).

---

### 8. Configuration File for Scenarios
**Current State**: Scenarios are hardcoded in `simulation_scenarios.py`.

**Improvement**: Load scenarios from a YAML or JSON config file for easier experimentation.

```yaml
# scenarios.yaml
martian_standard:
  target: Mars
  launch_date: 2028-10-01
  transfer_type: hohmann
  dt_days: 1.0
```

---

### 9. Interactive Sliders for Noise Parameters
**Current State**: Measurement noise (`R`) and process noise (`Q`) are fixed.

**Improvement**: Add sliders in the dashboard to tune `sigma_range`, `sigma_doppler`, and process noise in real-time and observe EKF convergence behavior.

---

## Research / Stretch Goals

### 10. Add Solar Radiation Pressure
Model the non-gravitational acceleration from photons hitting the spacecraft (significant for Parker Solar Probe-class missions).

### 11. Add Relativistic Orbit Propagation
Currently, GR is used for time dilation only. For extreme missions (e.g., orbits near black holes), use the post-Newtonian equations of motion.

### 12. Monte Carlo Covariance Analysis
Run 1000+ simulations with varying initial conditions and noise realizations to validate the filter's covariance estimates.

---

## Summary: V2 Target Feature List

| # | Feature | Priority |
|:---|:---|:---|
| 1 | JPL SPICE/HORIZONS ephemerides | 🔴 High |
| 2 | Multi-body gravity (Sun + Jupiter) | 🔴 High |
| 3 | Unscented Kalman Filter option | 🟡 Medium |
| 4 | Ionospheric delay model | 🟡 Medium |
| 5 | Two-way Doppler | 🟡 Medium |
| 6 | Save animation to video | 🟢 Low |
| 7 | Tracking station visualization | 🟢 Low |
| 8 | YAML scenario configuration | 🟢 Low |
| 9 | Interactive noise parameter tuning | 🟢 Low |
| 10 | Solar radiation pressure | 🔵 Research |
| 11 | Relativistic orbit propagation | 🔵 Research |
| 12 | Monte Carlo validation | 🔵 Research |

---

*Last updated: January 2026*
