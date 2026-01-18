# Spacetime Estimator V1: Walkthrough

## 1. Launching the Simulation
```bash
cd v1
python main.py
```
This command spawns 3 independent processes (Martian, Jovian, Solar Grazer).

## 2. The Dashboard Interface
*   **Top**: 3D Trajectory. (Green=True, Cyan=Est, Magenta=Drift).
*   **Right Hand Side**:
    *   **Logic Scale Error**: This is your primary metric. Watch it drop from 1000km to <1km.
    *   **Clock Bias**: Shows the filter estimating the "Time Offset".
    *   **Shapiro Delay**: Shows how many meters "extra" the signal traveled due to gravity.

## 3. Key Observations
*   **The "Snap"**: In the first few frames, you will see the Cyan line snap instantly to the Green line. This is the **Doppler** measurement constraining the velocity.
*   **The "Spike" (Scenario B)**: Around Day 300, watch the Shapiro graph spike. This is the Solar Conjunction.
*   **The "Drift" (Scenario C)**: Watch the Clock Bias grow parabolically. This is the Relativistic Time Dilation accumulating.

## 4. Interaction
*   **Slider**: Drag to scrub time. It is calibrated (0.1x - 10x).
*   **Pause**: Click to stop execution.
