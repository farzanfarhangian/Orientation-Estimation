# Orientation Estimation

**Real-time 3D attitude estimation using Madgwick, Mahony, and Extended Kalman Filter — benchmarked against u-blox NAV-ATT ground truth.**

[**Live Demo →**](https://farzanfarhangian.github.io/Orientation-Estimation/)

---

## Overview

This project implements and compares three industry-standard orientation estimation algorithms entirely in C, evaluated on real IMU data collected from a u-blox EVK GNSS/IMU module. A Python pipeline feeds the C output into an interactive Three.js visualizer that lets you scrub through the full dataset and watch all three algorithms track the device's attitude in real time.

All algorithms are **6-DoF** (accelerometer + gyroscope only). Yaw observability is restored by fusing **GPS heading** from the u-blox NAV-ATT message — giving each algorithm a real heading reference rather than leaving yaw as an unobservable degree of freedom.

---

## Live Demo

> **[farzanfarhangian.github.io/Orientation-Estimation](https://farzanfarhangian.github.io/Orientation-Estimation/)**

The visualizer shows:
- Three synchronized 3D device boxes (one per algorithm), with a yellow wireframe ground-truth ghost overlay
- Real-time roll / pitch / yaw readouts per algorithm
- Roll, pitch, and yaw error time-series charts vs. NAV-ATT ground truth
- RMSE badges computed over the full run
- Playback controls: play / pause, speed (0.25× – 8×), scrub, reset

---

## Results

Evaluated on a ~10-minute u-blox EVK dataset at ~107 Hz IMU, downsampled to ~10 Hz output.  
Ground truth: u-blox NAV-ATT (GPS/INS fused heading + tilt).

| Algorithm | Roll RMSE | Pitch RMSE | Yaw RMSE |
|-----------|:---------:|:----------:|:--------:|
| Madgwick  | 2.6°      | 5.3°       | 10.9°    |
| Mahony    | 2.1°      | 5.4°       | 10.9°    |
| EKF       | 4.5°      | 5.1°       | 26.1°    |

- **Roll / Pitch**: Corrected by accelerometer. Madgwick and Mahony converge faster due to their complementary filter structure; the EKF is more conservative with tighter process noise.
- **Yaw**: All three algorithms use GPS heading injection at ~5 Hz. Madgwick/Mahony use a fixed-gain yaw blend (α = 0.3); the EKF uses a full 1-D Kalman update with a correctly derived atan2 Jacobian.

---

## Algorithm Details

### Madgwick Filter
A gradient-descent-based complementary filter that minimises the rotation between the measured and predicted gravity direction. Computationally cheap and highly effective for roll/pitch. Parameters: `β = 0.033` (convergence rate), `ζ = 0` (gyro bias correction).

### Mahony Filter
A nonlinear complementary filter using proportional-integral feedback of the cross-product error between measured and predicted gravity. Well-suited to embedded systems. Parameters: `Kp = 0.5`, `Ki = 0.1`.

### Extended Kalman Filter
A full quaternion-state EKF:
- **State**: unit quaternion `q = [w, x, y, z]`
- **Process model**: gyroscope integration with linearised `F` matrix and process noise `Q = W·σ²_gyro·Wᵀ`
- **Measurement model**: gravity direction in body frame (6-DoF), with Joseph-form covariance update for numerical stability
- **Heading correction**: 1-D scalar Kalman update using `h(q) = atan2(2(wz+xy), 1−2(y²+z²))` with analytically derived Jacobian `H = ∂h/∂q`
- Parameters: `σ²_gyro = 0.01`, `σ²_accel = 0.001`, `σ²_heading = 0.0005`

### GPS Heading Fusion
At each GPS heading sample (~5 Hz), the current yaw estimate is corrected toward the GPS-reported heading using a wrapped innovation `v = ψ_GPS − ψ_est ∈ (−π, π]`:
- **EKF**: standard Kalman gain `K = PH'(HPH'+R)⁻¹`, Joseph-form `P` update
- **Madgwick / Mahony**: yaw-only quaternion rotation `q ← q ⊗ [cos(αv/2), 0, 0, sin(αv/2)]`

---

## Repository Structure

```
.
├── main.c                        # Entry point — loads CSVs, runs all three filters, outputs CSV
├── Makefile
├── CMakeLists.txt
│
├── include/
│   └── algorithms/
│       ├── orientationDefines.h  # Shared types: QUATERNION, VECTOR_3D, MATRIX_*
│       ├── ekf.h
│       ├── madgwick.h
│       └── mahony.h
│
├── src/
│   └── algorithms/
│       ├── ekf.c                 # EKF prediction, accel correction, heading correction
│       ├── madgwick.c            # Madgwick 6-DoF / 9-DoF update
│       ├── mahony.c              # Mahony 6-DoF / 9-DoF update
│       └── mathTransform.c       # Quaternion / DCM / Euler math library
│
├── data/
│   ├── imu_raw.csv               # timestamp_s, ax, ay, az (m/s²), gx, gy, gz (deg/s)
│   ├── ground_truth.csv          # timestamp_s, roll, pitch, yaw (deg) — u-blox NAV-ATT
│   └── results.json              # Pre-computed output for the visualizer
│
├── precompute_results.py         # Runs C binary → parses output → writes results.json
├── index.html                    # Three.js visualizer (single-file, no build step)
└── ublox_util.py                 # u-blox UBX binary parser used to generate the CSVs
```

---

## Build & Run

### Prerequisites

- GCC (C11) and GNU Make
- Python 3.8+ with `numpy` and `pandas` (for the pipeline script)

### Steps

```bash
# 1. Clone
git clone https://github.com/farzanfarhangian/Orientation-Estimation.git
cd Orientation-Estimation

# 2. Build the C binary
make

# 3. Run the pipeline (C binary → results.json)
python3 precompute_results.py

# 4. Serve the visualizer locally
python3 -m http.server 8000
# Open http://localhost:8000 in your browser
```

> The C binary writes CSV to stdout and progress to stderr.  
> `precompute_results.py` captures stdout, interpolates ground truth, and writes `data/results.json`.

---

## Data Format

### `data/imu_raw.csv`
```
timestamp_s, accel_x, accel_y, accel_z, gyro_x, gyro_y, gyro_z
```
Accelerometer in m/s², gyroscope in deg/s. Sampled at ~107 Hz.

### `data/ground_truth.csv`
```
timestamp_s, roll, pitch, yaw
```
All angles in degrees. Yaw in [0°, 360°]. Source: u-blox NAV-ATT GNSS/INS fusion at ~5 Hz.

---

## Coordinate Convention

- **Body frame**: right-hand, X-forward, Y-right, Z-down (NED-aligned)
- **Quaternion**: Hamilton convention `[w, x, y, z]`, body → NED rotation
- **Euler angles**: ZYX (aerospace) — yaw → pitch → roll, output in degrees
- **Visualizer axis remap** (NED → Three.js Y-up): `mesh.quaternion.set(x, −z, −y, w)`

---

## License

MIT © [Farzan Farhangian](https://github.com/farzanfarhangian)
