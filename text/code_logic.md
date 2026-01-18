# Code Logic Explainer: The Spacetime Estimator (V1)

This document explains the software architecture of the V1 "Flight-Ready" implementation.

## 1. Core Architecture Changes (vs V0)
While V0 was a standard MVC simulation, V1 introduces **coupled relativistic integration** and **multi-modal observable fusion**.
*   **State Space**: Increased from 9 to 10 dimensions.
*   **Observables**: Increased from 1 (Range) to 2 (Range + Doppler).
*   **Physics**: Non-linear light-time solution added to the measurement path.

## 2. The Augmented State Vector (`SpaceTimeState`)
We track 10 numbers (`v1/simulation_core.py`):
$$ X = [x, y, z, v_x, v_y, v_z, b, \dot{b}, \ddot{b}, \tau]^T $$
*   `0-2` ($r$): Position (Heliocentric Inertial).
*   `3-5` ($v$): Velocity.
*   `6-8`: Clock Bias, Drift, Drift Rate.
*   `9` ($\tau$): **Proper Time**. This is the time *experienced* by the probe, integrated along its world-line.

## 3. The Relativistic Physics Engine (`DynamicsModel`)
### Proper Time Integration
We integrate $\tau$ alongside position and velocity:
$$ \frac{d\tau}{dt} = \sqrt{1 - \frac{2U}{c^2} - \frac{v^2}{c^2}} \approx 1 - \frac{U}{c^2} - \frac{v^2}{2c^2} $$
where $U = \mu/r$. This captures gravitational time dilation and kinematic time dilation.

## 4. The Advanced Measurement Engine (`MeasurementModel`)
### A. The Light-Time Solution
In V0, we measured the probe where it *is*. In V1, we measure where it *was* when the signal left.
Method `solve_light_time` performs Fixed-Point Iteration:
1.  Guess $t_{emit} = t_{now} - |r|/c$.
2.  Update $r_{emit} \approx r(t_{now}) - v(t_{now}) \cdot \Delta t$.
3.  Recalculate distance and repeat.
This ensures the range measurement is physically causal.

### B. Doppler Observable
We emulate a "One-Way Doppler" link.
$$ z_{doppler} \approx \dot{\rho} + c \cdot \dot{b} $$
The geometric range rate $\dot{\rho}$ is calculated by projecting the velocity vector onto the line-of-sight unit vector.
$$ \dot{\rho} = \frac{\vec{v} \cdot \vec{r}}{|\vec{r}|} $$

### C. The Measurement Jacobian ($H$)
The $H$ matrix is now **2x10**.
*   **Row 1 (Range)**: Includes partials for Geometry, Shapiro, and Clock Bias ($c$).
*   **Row 2 (Doppler)**: Includes partials for Velocity ($H_{dop/vel} = \hat{u}$) and Position ($H_{dop/pos}$ accounting for changing angle).

## 5. The Estimator Tuning
The Kalman Filter matrices were retuned for stability:
*   **Doppler Noise ($R_{2,2}$)**: Set to $0.5 \text{ m/s}$. Tighter values caused numerical instability given the coupling with the speed of light ($c$).
*   **Clock Covariance ($P_0$)**: Tightened to prevent initial divergence.

## 6. Visualization
The `visualization.py` module was upgraded to handle the "Run All" capability using `multiprocessing` safe execution. It includes:
*   **Log-Scale Error**: To visualize the exponential convergence of the filter.
*   **Uncertainty Wireframes**: Visualizing the diagonal of $P$ as spatial volume.
