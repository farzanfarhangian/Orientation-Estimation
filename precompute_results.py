"""
Runs Madgwick, Mahony, and EKF (6-DoF, accel-only correction) on imu_raw.csv
and exports data/results.json for the Three.js visualization.

Output is downsampled to ~20 Hz to keep the JSON small (~2 MB).
Ground truth (NAV-ATT) is interpolated onto the same time grid.
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

# ── EKF 6-DoF accel-only correction (port of ekf.c, no-mag path) ─────────────

def mat_mul(A, B, ra, ca, cb):
    C = [[0.0]*cb for _ in range(ra)]
    for i in range(ra):
        for j in range(cb):
            for k in range(ca):
                C[i][j] += A[i][k]*B[k][j]
    return C

def mat_add(A, B, r, c):
    return [[A[i][j]+B[i][j] for j in range(c)] for i in range(r)]

def mat_scale(A, s, r, c):
    return [[A[i][j]*s for j in range(c)] for i in range(r)]

def mat_T(A, r, c):
    return [[A[i][j] for i in range(r)] for j in range(c)]

def mat_inv3(M):
    a,b,c = M[0]; d,e,f = M[1]; g,h,k = M[2]
    det = a*(e*k-f*h) - b*(d*k-f*g) + c*(d*h-e*g)
    if abs(det) < 1e-12: return [[1 if i==j else 0 for j in range(3)] for i in range(3)]
    inv = [
        [(e*k-f*h)/det, -(b*k-c*h)/det,  (b*f-c*e)/det],
        [-(d*k-f*g)/det,  (a*k-c*g)/det, -(a*f-c*d)/det],
        [(d*h-e*g)/det, -(a*h-b*g)/det,  (a*e-b*d)/det],
    ]
    return inv

class EKF6:
    def __init__(self, q, gyro_var=0.01, accel_var=0.01):
        self.q = list(q)
        self.P = [[1 if i==j else 0 for j in range(4)] for i in range(4)]
        self.gv = gyro_var
        self.av = accel_var

    def predict(self, gx, gy, gz, dt):
        gx *= DEG2RAD; gy *= DEG2RAD; gz *= DEG2RAD
        h = dt/2.0
        qw,qx,qy,qz = self.q
        self.q = qnorm([
            qw + h*(-gx*qx - gy*qy - gz*qz),
            qx + h*( gx*qw - gy*qz + gz*qy),
            qy + h*( gx*qz + gy*qw - gz*qx),
            qz + h*(-gx*qy + gy*qx + gz*qw),
        ])

        F = [
            [1,       -gx*h, -gy*h, -gz*h],
            [gx*h,    1,      gz*h, -gy*h],
            [gy*h,   -gz*h,   1,     gx*h],
            [gz*h,    gy*h,  -gx*h,  1   ],
        ]
        qw,qx,qy,qz = self.q
        W = [
            [-qx*h, -qy*h, -qz*h],
            [ qw*h, -qz*h,  qy*h],
            [ qz*h,  qw*h, -qx*h],
            [-qy*h,  qx*h,  qw*h],
        ]
        Q = mat_scale(mat_mul(W, mat_T(W,4,3), 4, 3, 4), self.gv, 4, 4)
        FP  = mat_mul(F, self.P, 4, 4, 4)
        FPFt= mat_mul(FP, mat_T(F,4,4), 4, 4, 4)
        self.P = mat_add(FPFt, Q, 4, 4)

    def correct(self, ax, ay, az):
        n = math.sqrt(ax*ax + ay*ay + az*az)
        if n == 0: return
        ax /= n; ay /= n; az /= n

        qw,qx,qy,qz = self.q
        # predicted gravity in body frame (DCM row 3 transposed)
        hx = 2*(qx*qz - qw*qy)
        hy = 2*(qy*qz + qw*qx)
        hz = qw*qw - qx*qx - qy*qy + qz*qz

        innov = [ax-hx, ay-hy, az-hz]

        # H = d(h)/d(q) * 2
        H = [
            [-qy,  qz, -qw, qx],
            [ qx,  qw,  qz, qy],
            [ qw, -qx, -qy, qz],
        ]
        H = mat_scale(H, 2.0, 3, 4)

        HP  = mat_mul(H, self.P, 3, 4, 4)
        HPHt= mat_mul(HP, mat_T(H,3,4), 3, 4, 3)
        R   = [[self.av if i==j else 0 for j in range(3)] for i in range(3)]
        S   = mat_add(HPHt, R, 3, 3)
        Sinv= mat_inv3(S)

        PHt = mat_mul(self.P, mat_T(H,3,4), 4, 4, 3)
        K   = mat_mul(PHt, Sinv, 4, 3, 3)

        # update q
        for i,name in enumerate(['q0','q1','q2','q3']):
            self.q[i] += K[i][0]*innov[0] + K[i][1]*innov[1] + K[i][2]*innov[2]
        self.q = qnorm(self.q)

        # update P
        KH  = mat_mul(K, H, 4, 3, 4)
        I   = [[1 if i==j else 0 for j in range(4)] for i in range(4)]
        IKH = mat_add(I, mat_scale(KH,-1,4,4), 4, 4)
        self.P = mat_mul(IKH, self.P, 4, 4, 4)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading CSVs...")
    imu = pd.read_csv(IMU_CSV)
    gt  = pd.read_csv(GT_CSV)

    # skip NAV-ATT init period (first 475 samples have roll=0)
    gt = gt[gt['roll'] != 0.0].reset_index(drop=True)

    # determine stride to reach TARGET_HZ
    dt_imu = float(imu['timestamp_s'].diff().median())
    stride  = max(1, round(1.0 / (TARGET_HZ * dt_imu)))
    print(f"IMU dt={dt_imu*1000:.2f}ms, stride={stride} → ~{1/(stride*dt_imu):.1f} Hz output")

    # initialise filters from first sample
    a0 = imu.iloc[0]
    q0 = init_quat_from_accel(a0.accel_x, a0.accel_y, a0.accel_z)

    q_madg  = list(q0)
    q_maho  = list(q0)
    bias_maho = [0.0, 0.0, 0.0]
    ekf = EKF6(q0)

    out_t, out_madg, out_maho, out_ekf, out_gt = [], [], [], [], []

    gt_ts  = gt['timestamp_s'].values
    gt_roll= gt['roll'].values
    gt_pit = gt['pitch'].values
    gt_yaw = gt['yaw'].values

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
        ekf.correct(ax,ay,az)

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
