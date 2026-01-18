# Project Defense: Anticipated Evaluator Questions
**Spacetime Estimator (V0 & V1)**

---

## 1. Physics & Relativity
**Q: Why do you need General Relativity? Isn't Newton enough?**
**A:** For getting to Mars? Maybe. But for deep space tracking (like Parker Solar Probe), the **Shapiro Delay** adds ~30km of perceived distance error near the Sun. Precision navigation (meters) fails if you ignore a 30km error. Also, clocks on orbit drift by microseconds per day due to relativity; at light speed, 1 microsecond = 300 meters of error.

**Q: Explain the Shapiro Delay formula in your code.**
**A:** We use the Schwarzschild metric approx: $\Delta t = \frac{4GM}{c^3} \ln(\frac{4x_1 x_2}{b^2})$. The term `4GM/c^2` is the key coefficient. In V0, we initially had `2GM/c^2` (incorrect), which we fixed. It represents the "curvature" the signal has to traverse.

**Q: What is "Proper Time" and why is it a state in V1?**
**A:** Proper Time ($\tau$) is the time counted by the onboard clock. Orbit determination requires correlating *Onboard Time* to *Earth Time*. Since $\tau$ diverges from $t$ non-linearly (due to velocity $v^2$ and gravity $U$), we must estimate it effectively as a "dynamic bias."

## 2. The Kalman Filter (Estimation)
**Q: Why use an Extended Kalman Filter (EKF) and not a linear KF?**
**A:** The measurement function $z = \sqrt{x^2 + y^2 + z^2}$ is non-linear. A linear KF would fail. The EKF linearizes this around the current estimate using the Jacobian matrix (Taylor Series expansion).

**Q: Your filter "spiraled" in V1 initially. Why?**
**A:** Numerical Instability. We had Doppler noise set to `1 mm/s` ($10^{-3}$). The speed of light is $10^8$. The condition number of the matrix (ratio of largest to smallest eigenvalue) became too large for floating-point precision. Relaxing the Doppler noise to `0.5 m/s` stabilized the matrix inversion.

**Q: What does the "Trace of P" tell you?**
**A:** The Trace of the Covariance Matrix ($P$) is the sum of the variances. It represents the "total uncertainty" volume. In the simulation, you see it drop when measurements are good, and grow during prediction (coast) phases.

## 3. Software Architecture
**Q: How did you handle the visualization lag?**
**A:** In V0, everything ran in one thread. In V1 (and the optimized launcher), we used Python's `multiprocessing`. Each scenario runs in a separate OS process with its own memory space. This prevents the heavy Matplotlib rendering loop from blocking the numerical integration loop.

**Q: Why separate `DynamicsModel` and `MeasurementModel`?**
**A:** Separation of concerns.
*   `DynamicsModel` is the "Physics" (Gravity, Mechanics).
*   `MeasurementModel` is the "Sensor" (Radios, Telescopes).
This allows us to swap sensors (e.g., adding a Telescope/Angle measurement) without changing the physics of how the probe moves.

## 4. Advanced Topics
**Q: How would you make this better?**
**A:**
1.  **Iterative Smoother**: Use an RTS Smoother to process the data backwards after the mission.
2.  **Ephemerides**: Use real JPL planet positions (DE440) instead of circular orbits.
3.  **Unscented Kalman Filter (UKF)**: To handle the highly non-linear Solar Conjunctions without Jacobian linearization errors.

**Q: What is the "Light-Time Solution" in V1?**
**A:** A radio signal from Jupiter takes 40 minutes. We can't use the probe's *current* position to calculate range. We solved the implicit equation $t_{emit} = t_{recv} - |r(t_{emit})|/c$ iteratively to find where the probe *physically was* when it sent the ping.
