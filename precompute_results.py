"""
Calls the compiled C binary (./my_prog) which runs Madgwick, Mahony, and EKF
on data/imu_raw.csv and outputs results as CSV to stdout.
This script reads that CSV, interpolates ground-truth angles, and writes
data/results.json for the Three.js visualizer.
"""

import json
import math
import subprocess
import numpy as np
import pandas as pd

GT_CSV   = "data/ground_truth.csv"
OUT_JSON = "data/results.json"
BINARY   = "./my_prog"

RAD2DEG = 180.0 / math.pi

def ang_wrap(deg):
    return (deg + 180.0) % 360.0 - 180.0

def main():
    print("Loading ground-truth CSV...")
    gt = pd.read_csv(GT_CSV)
    gt = gt[gt['roll'] != 0.0].reset_index(drop=True)
    T_START = float(gt['timestamp_s'].iloc[0])
    T_END   = T_START + 600.0
    gt = gt[(gt['timestamp_s'] >= T_START) & (gt['timestamp_s'] <= T_END)].reset_index(drop=True)
    print(f"GT window: {T_START:.1f}s – {T_END:.1f}s  ({len(gt)} samples)")

    gt_ts  = gt['timestamp_s'].values
    gt_roll = gt['roll'].values
    gt_pit  = gt['pitch'].values
    gt_yaw  = gt['yaw'].values

    print(f"Running C binary: {BINARY} ...")
    proc = subprocess.run([BINARY], capture_output=True, text=True)
    if proc.returncode != 0:
        print("ERROR: binary failed:\n", proc.stderr)
        raise SystemExit(1)

    # print progress lines from stderr
    for line in proc.stderr.strip().splitlines():
        print(" ", line)

    print("Parsing binary output...")
    from io import StringIO
    df = pd.read_csv(StringIO(proc.stdout))
    N = len(df)
    print(f"  {N} output frames")

    out_t     = df['t'].round(4).tolist()
    out_madg  = [[round(df['madg_w'][i],6), round(df['madg_x'][i],6),
                   round(df['madg_y'][i],6), round(df['madg_z'][i],6)] for i in range(N)]
    out_maho  = [[round(df['maho_w'][i],6), round(df['maho_x'][i],6),
                   round(df['maho_y'][i],6), round(df['maho_z'][i],6)] for i in range(N)]
    out_ekf   = [[round(df['ekf_w'][i],6),  round(df['ekf_x'][i],6),
                   round(df['ekf_y'][i],6),  round(df['ekf_z'][i],6)]  for i in range(N)]

    euler_madg = [[round(df['madg_roll'][i],3), round(df['madg_pitch'][i],3),
                    round(df['madg_yaw'][i],3)] for i in range(N)]
    euler_maho = [[round(df['maho_roll'][i],3), round(df['maho_pitch'][i],3),
                    round(df['maho_yaw'][i],3)] for i in range(N)]
    euler_ekf  = [[round(df['ekf_roll'][i],3),  round(df['ekf_pitch'][i],3),
                    round(df['ekf_yaw'][i],3)]  for i in range(N)]

    # interpolate ground truth onto the same time grid
    print("Interpolating ground-truth...")
    out_gt = []
    for t in df['t'].values:
        gi = np.searchsorted(gt_ts, t)
        if gi == 0:
            out_gt.append([round(gt_roll[0],4), round(gt_pit[0],4), round(gt_yaw[0],4)])
        elif gi >= len(gt_ts):
            out_gt.append([round(gt_roll[-1],4), round(gt_pit[-1],4), round(gt_yaw[-1],4)])
        else:
            t0, t1 = gt_ts[gi-1], gt_ts[gi]
            alpha = (t - t0) / (t1 - t0) if t1 != t0 else 0.0
            r = gt_roll[gi-1] + alpha * (gt_roll[gi] - gt_roll[gi-1])
            p = gt_pit[gi-1]  + alpha * (gt_pit[gi]  - gt_pit[gi-1])
            y = gt_yaw[gi-1]  + alpha * (gt_yaw[gi]  - gt_yaw[gi-1])
            out_gt.append([round(r,4), round(p,4), round(y,4)])

    # output rate (from actual sample times)
    dt_median = float(np.median(np.diff(df['t'].values)))
    hz = round(1.0 / dt_median, 1) if dt_median > 0 else 10.0

    result = {
        "meta": {
            "hz": hz,
            "n": N,
            "duration_s": round(out_t[-1], 1),
        },
        "t": out_t,
        "madgwick": {"q": out_madg, "euler": euler_madg},
        "mahony":   {"q": out_maho, "euler": euler_maho},
        "ekf":      {"q": out_ekf,  "euler": euler_ekf},
        "gt":       {"euler": out_gt},
    }

    with open(OUT_JSON, 'w') as f:
        json.dump(result, f, separators=(',', ':'))

    size_mb = len(json.dumps(result)) / 1e6
    print(f"\nSaved {OUT_JSON}  ({size_mb:.1f} MB, {N} frames, {hz:.1f} Hz)")

if __name__ == "__main__":
    main()
