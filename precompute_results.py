"""
Runs Madgwick, Mahony, and EKF (9-DoF: accel + GPS heading correction) on imu_raw.csv
and exports data/results.json for the Three.js visualization.

Output is downsampled to ~10 Hz to keep the JSON small.
Ground truth (NAV-ATT) is interpolated onto the same time grid.
GPS heading from NAV-ATT is fed into the EKF at ~5 Hz.
"""

import json
import math
import numpy as np
import pandas as pd

IMU_CSV = "data/imu_raw.csv"
GT_CSV  = "data/ground_truth.csv"
OUT_JSON = "data/results.json"

TARGET_HZ = 10          # output rate for visualizer
DEG2RAD   = math.pi / 180.0
RAD2DEG   = 180.0 / math.pi

# ── quaternion helpers ────────────────────────────────────────────────────────

def qnorm(q):
    n = math.sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3])
    if n == 0: return q
    return [x/n for x in q]

def qmul(a, b):
    aw,ax,ay,az = a
    bw,bx,by,bz = b
    return [
        aw*bw - ax*bx - ay*by - az*bz,
        aw*bx + ax*bw + ay*bz - az*by,
        aw*by - ax*bz + ay*bw + az*bx,
        aw*bz + ax*by - ay*bx + az*bw,
    ]

def quat_to_euler(q):
    """Returns (roll, pitch, yaw) in degrees."""
    w,x,y,z = q
    roll  = RAD2DEG * math.atan2(2*(w*x + y*z), 1 - 2*(x*x + y*y))
    sinp  = 2*(w*y - z*x)
    sinp  = max(-1.0, min(1.0, sinp))
    pitch = RAD2DEG * math.asin(sinp)
    yaw   = RAD2DEG * math.atan2(2*(w*z + x*y), 1 - 2*(y*y + z*z))
    return roll, pitch, yaw

def euler_to_quat(roll_deg, pitch_deg, yaw_deg):
    """Convert roll/pitch/yaw (degrees) to quaternion [w,x,y,z]."""
    r = roll_deg  * DEG2RAD
    p = pitch_deg * DEG2RAD
    y = yaw_deg   * DEG2RAD
    cr,sr = math.cos(r/2), math.sin(r/2)
    cp,sp = math.cos(p/2), math.sin(p/2)
    cy,sy = math.cos(y/2), math.sin(y/2)
    return qnorm([
        cr*cp*cy + sr*sp*sy,
        sr*cp*cy - cr*sp*sy,
        cr*sp*cy + sr*cp*sy,
        cr*cp*sy - sr*sp*cy,
    ])

def init_quat_from_accel(ax, ay, az):
    """Rough tilt-only initialisation from a single accelerometer reading."""
    n = math.sqrt(ax*ax + ay*ay + az*az)
    if n == 0: return [1,0,0,0]
    ax /= n; ay /= n; az /= n
    roll  = math.atan2(ay, az)
    pitch = math.atan2(-ax, math.sqrt(ay*ay + az*az))
    # yaw = 0 (no magnetometer)
    cr,sr = math.cos(roll/2),  math.sin(roll/2)
    cp,sp = math.cos(pitch/2), math.sin(pitch/2)
    cy,sy = 1.0, 0.0
    return qnorm([
        cr*cp*cy + sr*sp*sy,
        sr*cp*cy - cr*sp*sy,
        cr*sp*cy + sr*cp*sy,
        cr*cp*sy - sr*sp*cy,
    ])

# ── Madgwick 6-DoF (direct port of madgwick.c) ───────────────────────────────

def madgwick_update(q, ax, ay, az, gx, gy, gz, dt, beta=0.033):
    gx *= DEG2RAD; gy *= DEG2RAD; gz *= DEG2RAD

    qw,qx,qy,qz = q

    # quaternion derivative from gyro
    qdot_w = 0.5*(-qx*gx - qy*gy - qz*gz)
    qdot_x = 0.5*( qw*gx + qy*gz - qz*gy)
    qdot_y = 0.5*( qw*gy - qx*gz + qz*gx)
    qdot_z = 0.5*( qw*gz + qx*gy - qy*gx)

    n = math.sqrt(ax*ax + ay*ay + az*az)
    if n != 0:
        ax /= n; ay /= n; az /= n

        d2w = 2*qw; d2x = 2*qx; d2y = 2*qy; d2z = 2*qz
        d4w = 4*qw; d4x = 4*qx; d4y = 4*qy
        d8x = 8*qx; d8y = 8*qy
        ww = qw*qw; xx = qx*qx; yy = qy*qy; zz = qz*qz

        gW = d4w*yy + d2y*ax + d4w*xx - d2x*ay
        gX = d4x*zz - d2z*ax + 4*ww*qx - d2w*ay - d4x + d8x*xx + d8x*yy + d4x*az
        gY = 4*ww*qy + d2w*ax + d4y*zz - d2z*ay - d4y + d8y*xx + d8y*yy + d4y*az
        gZ = 4*xx*qz - d2x*ax + 4*yy*qz - d2y*ay

        gn = math.sqrt(gW*gW + gX*gX + gY*gY + gZ*gZ)
        if gn != 0:
            gW /= gn; gX /= gn; gY /= gn; gZ /= gn

        qdot_w -= beta*gW
        qdot_x -= beta*gX
        qdot_y -= beta*gY
        qdot_z -= beta*gZ

    q = qnorm([
        qw + qdot_w*dt,
        qx + qdot_x*dt,
        qy + qdot_y*dt,
        qz + qdot_z*dt,
    ])
    return q

# ── Mahony 6-DoF (direct port of mahony.c) ───────────────────────────────────

def mahony_update(q, ax, ay, az, gx, gy, gz, dt, kP=0.1, kI=1.0, bias=[0,0,0]):
    gx *= DEG2RAD; gy *= DEG2RAD; gz *= DEG2RAD

    n = math.sqrt(ax*ax + ay*ay + az*az)
    if n != 0:
        ax /= n; ay /= n; az /= n

        # DCM row 3 from quaternion
        qw,qx,qy,qz = q
        r31 = 2*(qx*qz - qw*qy)
        r32 = 2*(qy*qz + qw*qx)
        r33 = qw*qw - qx*qx - qy*qy + qz*qz

        # cross product: accel_norm x DCM_row3
        ex = ay*r33 - az*r32
        ey = az*r31 - ax*r33
        ez = ax*r32 - ay*r31

        # integrate bias
        bias[0] += -kI * ex * dt
        bias[1] += -kI * ey * dt
        bias[2] += -kI * ez * dt

        gx += kP*ex - bias[0]
        gy += kP*ey - bias[1]
        gz += kP*ez - bias[2]

    qw,qx,qy,qz = q
    dw = 0.5*(-qx*gx - qy*gy - qz*gz)
    dx = 0.5*( qw*gx + qy*gz - qz*gy)
    dy = 0.5*( qw*gy - qx*gz + qz*gx)
    dz = 0.5*( qw*gz + qx*gy - qy*gx)

    q = qnorm([qw+dw*dt, qx+dx*dt, qy+dy*dt, qz+dz*dt])
    return q, bias

# ── EKF matrix helpers (using numpy for clarity) ─────────────────────────────

class EKF9:
    """
    EKF with accelerometer correction every IMU step and GPS heading correction
    whenever a new heading measurement is available (~5 Hz from NAV-ATT).

    State: quaternion [w, x, y, z]  (Hamilton convention)
    Accel measurement model: normalized gravity in body frame (3-axis)
    Heading measurement model: yaw = atan2(2*(qx*qy + qw*qz), qw²+qx²-qy²-qz²)
    """

    def __init__(self, q, gyro_var=0.005, accel_var=0.1, heading_var=0.003):
        self.q = np.array(q, dtype=np.float64)
        self.q /= np.linalg.norm(self.q)
        # Initial covariance — start tight on tilt (known from accel init),
        # slightly looser on yaw (known from GPS init but GPS has noise)
        self.P = np.diag([0.001, 0.001, 0.001, 0.005])
        self.gv = gyro_var      # (rad/s)² per axis
        self.av = accel_var     # normalized accel noise variance per axis
        self.hv = heading_var   # (rad)² heading noise variance

    def predict(self, gx, gy, gz, dt):
        gx *= DEG2RAD; gy *= DEG2RAD; gz *= DEG2RAD
        h = dt / 2.0
        qw, qx, qy, qz = self.q

        # Quaternion kinematics: q_new = (I + h*Omega) * q
        self.q = np.array([
            qw + h*(-gx*qx - gy*qy - gz*qz),
            qx + h*( gx*qw + gz*qy - gy*qz),
            qy + h*( gy*qw - gz*qx + gx*qz),
            qz + h*( gz*qw + gy*qx - gx*qy),
        ])
        self.q /= np.linalg.norm(self.q)

        # State transition Jacobian F = I + h*Omega_matrix
        F = np.array([
            [1,      -gx*h, -gy*h, -gz*h],
            [gx*h,    1,     gz*h, -gy*h],
            [gy*h,   -gz*h,  1,     gx*h],
            [gz*h,    gy*h, -gx*h,  1   ],
        ])

        # Process noise: Q = h² * Xi * R_gyro * Xi^T
        qw, qx, qy, qz = self.q
        Xi = h * np.array([
            [-qx, -qy, -qz],
            [ qw, -qz,  qy],
            [ qz,  qw, -qx],
            [-qy,  qx,  qw],
        ])
        R_gyro = self.gv * np.eye(3)
        Q = Xi @ R_gyro @ Xi.T

        self.P = F @ self.P @ F.T + Q

    def correct_accel(self, ax, ay, az):
        n = math.sqrt(ax*ax + ay*ay + az*az)
        if n == 0:
            return
        ax /= n; ay /= n; az /= n

        qw, qx, qy, qz = self.q

        # Expected normalized gravity in body frame (third row of C_nb^T)
        hx = 2*(qx*qz - qw*qy)
        hy = 2*(qy*qz + qw*qx)
        hz = qw*qw - qx*qx - qy*qy + qz*qz

        innov = np.array([ax - hx, ay - hy, az - hz])

        # Jacobian: d[hx,hy,hz]/d[qw,qx,qy,qz]
        H = 2.0 * np.array([
            [-qy,  qz, -qw,  qx],
            [ qx,  qw,  qz,  qy],
            [ qw, -qx, -qy,  qz],
        ])

        R = self.av * np.eye(3)
        S = H @ self.P @ H.T + R
        K = self.P @ H.T @ np.linalg.inv(S)

        self.q = self.q + K @ innov
        self.q /= np.linalg.norm(self.q)

        # Joseph form for numerical stability
        I_KH = np.eye(4) - K @ H
        self.P = I_KH @ self.P @ I_KH.T + K @ R @ K.T

    def correct_heading(self, yaw_deg):
        """Update using GPS/NAV-ATT heading. yaw_deg in degrees."""
        yaw_meas = yaw_deg * DEG2RAD
        qw, qx, qy, qz = self.q

        # Expected yaw from quaternion: atan2(2*(qx*qy + qw*qz), qw²+qx²-qy²-qz²)
        num = 2.0*(qx*qy + qw*qz)
        den = qw*qw + qx*qx - qy*qy - qz*qz
        yaw_pred = math.atan2(num, den)

        # Innovation with wrap-around
        innov = math.atan2(math.sin(yaw_meas - yaw_pred),
                           math.cos(yaw_meas - yaw_pred))

        # Jacobian: d(yaw)/d(q)  via quotient rule on atan2(num, den)
        d = num*num + den*den
        if d < 1e-9:
            return  # near-singularity (pitch ≈ ±90°), skip
        H = np.array([[
            2*(-qz*den - qw*num) / d,   # d/dqw
            2*( qy*den - qx*num) / d,   # d/dqx — wait, sign check below
            2*( qx*den + qy*num) / d,   # d/dqy
            2*( qw*den - qz*num) / d,   # d/dqz
        ]])
        # Verify signs:
        #   d(num)/dq = [2qz, 2qy, 2qx, 2qw]
        #   d(den)/dq = [2qw, 2qx, -2qy, -2qz]
        #   H = (den*d_num - num*d_den) / (num²+den²)
        d_num = np.array([2*qz, 2*qy, 2*qx, 2*qw])
        d_den = np.array([2*qw, 2*qx, -2*qy, -2*qz])
        H = ((den * d_num - num * d_den) / d).reshape(1, 4)

        R = np.array([[self.hv]])
        S = H @ self.P @ H.T + R
        K = self.P @ H.T / float(S[0, 0])

        self.q = self.q + K.flatten() * innov
        self.q /= np.linalg.norm(self.q)

        I_KH = np.eye(4) - K @ H
        self.P = I_KH @ self.P @ I_KH.T + K @ R @ K.T


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading CSVs...")
    imu = pd.read_csv(IMU_CSV)
    gt  = pd.read_csv(GT_CSV)

    # Keep only the window where GPS heading is valid
    gt = gt[gt['roll'] != 0.0].reset_index(drop=True)
    T_START = float(gt['timestamp_s'].iloc[0])   # first valid GT sample (~50.4 s)
    T_END   = T_START + 600.0                     # 10 minutes of aligned data

    imu = imu[(imu['timestamp_s'] >= T_START) & (imu['timestamp_s'] <= T_END)].reset_index(drop=True)
    gt  = gt[ (gt['timestamp_s']  >= T_START) & (gt['timestamp_s']  <= T_END) ].reset_index(drop=True)
    print(f"Window: {T_START:.1f}s – {T_END:.1f}s  |  IMU={len(imu)} samples  GT={len(gt)} samples")

    # determine stride to reach TARGET_HZ
    dt_imu = float(imu['timestamp_s'].diff().median())
    stride  = max(1, round(1.0 / (TARGET_HZ * dt_imu)))
    print(f"IMU dt={dt_imu*1000:.2f}ms, stride={stride} → ~{1/(stride*dt_imu):.1f} Hz output")

    gt_ts_all  = gt['timestamp_s'].values
    gt_roll_all= gt['roll'].values
    gt_pit_all = gt['pitch'].values
    gt_yaw_all = gt['yaw'].values

    # Initialise all filters from the first GT sample (roll + pitch + yaw all known)
    q0 = euler_to_quat(gt_roll_all[0], gt_pit_all[0], gt_yaw_all[0])
    print(f"Init from GT: roll={gt_roll_all[0]:.2f}°  pitch={gt_pit_all[0]:.2f}°  yaw={gt_yaw_all[0]:.2f}°")

    q_madg  = list(q0)
    q_maho  = list(q0)
    bias_maho = [0.0, 0.0, 0.0]
    ekf = EKF9(q0, gyro_var=0.005, accel_var=0.1, heading_var=0.003)

    out_t, out_madg, out_maho, out_ekf, out_gt = [], [], [], [], []

    gt_ts  = gt_ts_all
    gt_roll = gt_roll_all
    gt_pit  = gt_pit_all
    gt_yaw  = gt_yaw_all

    # Track which GT sample was last used for heading correction
    last_gt_idx = -1

    print(f"Running algorithms on {len(imu)} samples...")
    for idx, row in enumerate(imu.itertuples(index=False)):
        dt  = dt_imu if idx > 0 else dt_imu
        ax,ay,az = row.accel_x, row.accel_y, row.accel_z
        gx,gy,gz = row.gyro_x,  row.gyro_y,  row.gyro_z
        t = row.timestamp_s

        q_madg = madgwick_update(q_madg, ax,ay,az, gx,gy,gz, dt, beta=0.033)
        q_maho, bias_maho = mahony_update(q_maho, ax,ay,az, gx,gy,gz, dt,
                                          kP=0.5, kI=0.1, bias=bias_maho)
        ekf.predict(gx,gy,gz, dt)
        ekf.correct_accel(ax,ay,az)

        # Feed GPS heading into EKF whenever a new GT sample is available
        gi = int(np.searchsorted(gt_ts, t, side='right')) - 1
        if gi >= 0 and gi < len(gt_ts) and gi != last_gt_idx:
            ekf.correct_heading(gt_yaw[gi])
            last_gt_idx = gi

        if idx % stride == 0:
            out_t.append(round(t, 4))
            out_madg.append([round(v,6) for v in q_madg])
            out_maho.append([round(v,6) for v in q_maho])
            out_ekf.append([round(v,6) for v in ekf.q])

            # interpolate ground truth
            gi = np.searchsorted(gt_ts, t)
            if gi == 0:
                out_gt.append([round(gt_roll[0],4), round(gt_pit[0],4), round(gt_yaw[0],4)])
            elif gi >= len(gt_ts):
                out_gt.append([round(gt_roll[-1],4), round(gt_pit[-1],4), round(gt_yaw[-1],4)])
            else:
                t0,t1 = gt_ts[gi-1], gt_ts[gi]
                alpha = (t-t0)/(t1-t0) if t1 != t0 else 0
                r = gt_roll[gi-1] + alpha*(gt_roll[gi]-gt_roll[gi-1])
                p = gt_pit[gi-1]  + alpha*(gt_pit[gi]-gt_pit[gi-1])
                y = gt_yaw[gi-1]  + alpha*(gt_yaw[gi]-gt_yaw[gi-1])
                out_gt.append([round(r,4), round(p,4), round(y,4)])

        if idx % 10000 == 0:
            print(f"  {idx}/{len(imu)} ({100*idx/len(imu):.0f}%)")

    # convert quaternions to euler for chart
    euler_madg = [quat_to_euler(q) for q in out_madg]
    euler_maho = [quat_to_euler(q) for q in out_maho]
    euler_ekf  = [quat_to_euler(q) for q in out_ekf]

    result = {
        "meta": {
            "hz": TARGET_HZ,
            "n": len(out_t),
            "duration_s": round(out_t[-1], 1),
        },
        "t":      out_t,
        "madgwick": {"q": out_madg, "euler": [[round(v,3) for v in e] for e in euler_madg]},
        "mahony":   {"q": out_maho, "euler": [[round(v,3) for v in e] for e in euler_maho]},
        "ekf":      {"q": out_ekf,  "euler": [[round(v,3) for v in e] for e in euler_ekf]},
        "gt":       {"euler": out_gt},
    }

    with open(OUT_JSON, 'w') as f:
        json.dump(result, f, separators=(',',':'))

    size_mb = len(json.dumps(result)) / 1e6
    print(f"\nSaved {OUT_JSON}  ({size_mb:.1f} MB, {len(out_t)} frames)")

if __name__ == "__main__":
    main()
