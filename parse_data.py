"""
Parses the UBX file and exports:
  data/imu_raw.csv      - raw IMU at ~104Hz (timestamp_s, accel_x/y/z m/s², gyro_x/y/z deg/s)
  data/ground_truth.csv - NAV-ATT attitude at ~10Hz (timestamp_s, roll, pitch, yaw deg)
"""

import os
import sys
import numpy as np
import pandas as pd
from ublox_util import parse_ubx_file

UBX_FILE   = os.path.join("data", "2025-5-25_204743_serial-COM9.ubx")
OUTPUT_DIR = "data"

def main():
    print("Parsing UBX file (this may take a moment)...")
    (df_positions, df_satellite, df_esf, df_esf_cal,
     df_esf_alg, df_esf_raw, df_esf_ins, df_nav_att,
     df_nav_dop, raw_msgs) = parse_ubx_file(UBX_FILE)

    # ── IMU raw ──────────────────────────────────────────────────────────────
    required_imu = ["time", "raw_accel_x", "raw_accel_y", "raw_accel_z",
                    "raw_gyro_x",  "raw_gyro_y",  "raw_gyro_z"]

    if df_esf_raw.empty or not all(c in df_esf_raw.columns for c in required_imu):
        print("ERROR: ESF-RAW data missing or incomplete.")
        sys.exit(1)

    df_imu = df_esf_raw[required_imu].dropna().copy()
    df_imu["time"] = pd.to_datetime(df_imu["time"])
    df_imu.sort_values("time", inplace=True)

    t0 = df_imu["time"].iloc[0]
    df_imu["timestamp_s"] = (df_imu["time"] - t0).dt.total_seconds()

    imu_out = df_imu[["timestamp_s",
                       "raw_accel_x", "raw_accel_y", "raw_accel_z",
                       "raw_gyro_x",  "raw_gyro_y",  "raw_gyro_z"]].copy()
    imu_out.columns = ["timestamp_s",
                       "accel_x", "accel_y", "accel_z",
                       "gyro_x",  "gyro_y",  "gyro_z"]

    imu_path = os.path.join(OUTPUT_DIR, "imu_raw.csv")
    imu_out.to_csv(imu_path, index=False, float_format="%.6f")
    print(f"Saved {len(imu_out)} IMU samples → {imu_path}")
    print(f"  Duration : {imu_out['timestamp_s'].iloc[-1]:.1f} s")
    print(f"  Avg rate : {len(imu_out) / imu_out['timestamp_s'].iloc[-1]:.1f} Hz")

    # ── Ground truth (NAV-ATT) ────────────────────────────────────────────────
    required_gt = ["time", "roll", "pitch", "heading"]

    if df_nav_att.empty or not all(c in df_nav_att.columns for c in required_gt):
        print("WARNING: NAV-ATT ground truth not available; trying ESF-ALG...")
        if not df_esf_alg.empty and all(c in df_esf_alg.columns for c in ["time", "roll", "pitch", "yaw"]):
            df_gt = df_esf_alg[["time", "roll", "pitch", "yaw"]].dropna().copy()
            df_gt.rename(columns={"yaw": "heading"}, inplace=True)
        else:
            print("ERROR: No ground truth data found.")
            sys.exit(1)
    else:
        df_gt = df_nav_att[required_gt].dropna().copy()

    df_gt["time"] = pd.to_datetime(df_gt["time"])
    df_gt.sort_values("time", inplace=True)
    df_gt["timestamp_s"] = (df_gt["time"] - t0).dt.total_seconds()

    gt_out = df_gt[["timestamp_s", "roll", "pitch", "heading"]].copy()
    gt_out.columns = ["timestamp_s", "roll", "pitch", "yaw"]

    gt_path = os.path.join(OUTPUT_DIR, "ground_truth.csv")
    gt_out.to_csv(gt_path, index=False, float_format="%.6f")
    print(f"Saved {len(gt_out)} ground-truth samples → {gt_path}")

    print("\nDone. Ready to feed into orientation estimation algorithms.")

if __name__ == "__main__":
    main()
