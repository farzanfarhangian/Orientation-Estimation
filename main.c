/**
 * @author Farzan Farhangian
 * @file main.c
 *
 * Runs Madgwick, Mahony, and EKF on IMU data from data/imu_raw.csv and
 * ground-truth initialisation from data/ground_truth.csv.
 * Prints downsampled results as CSV to stdout:
 *   t,mw,mx,my,mz,mr,mp,my_deg,hw,hx,hy,hz,hr,hp,hy_deg,ew,ex,ey,ez,er,ep,ey_deg
 *
 * Gyro in CSV is deg/s → converted to rad/s before passing to algorithms.
 * All three algorithms are 6-DoF (accel only); no magnetometer is used.
 *
 * Build:  make
 * Run:    ./my_prog
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "algorithms/madgwick.h"
#include "algorithms/mahony.h"
#include "algorithms/ekf.h"

/* ── tuning parameters ─────────────────────────────────────────────────── */
#define MADGWICK_BETA   0.033f
#define MADGWICK_ZETA   0.0f
#define MAHONY_BETA     0.0f
#define MAHONY_KP       0.5f
#define MAHONY_KI       0.1f
#define EKF_GYRO_VAR    0.01f  /* sigma_gyro=0.1 rad/s → variance=0.1²; tuned for this dataset */
#define EKF_ACCEL_VAR   0.001f /* sigma_accel=0.0316 m/s²; tight trust in accel measurement */
#define EKF_MAG_VAR     1e6f   /* large → effectively disables magnetometer */
#define EKF_DIP         0.0f

#define TARGET_HZ       10
#define IMU_HZ_APPROX   107    /* used to compute stride; will be measured */

/* ── CSV paths ─────────────────────────────────────────────────────────── */
#define IMU_CSV  "data/imu_raw.csv"
#define GT_CSV   "data/ground_truth.csv"

/* ── simple dynamic array of doubles ───────────────────────────────────── */
typedef struct {
    double *d;
    size_t  n;
    size_t  cap;
} DblArr;

static void da_push(DblArr *a, double v) {
    if (a->n >= a->cap) {
        a->cap = a->cap ? a->cap * 2 : 4096;
        a->d = realloc(a->d, a->cap * sizeof(double));
    }
    a->d[a->n++] = v;
}

static void da_free(DblArr *a) { free(a->d); a->d = NULL; a->n = a->cap = 0; }

/* ── CSV parsing helpers ────────────────────────────────────────────────── */
static int skip_header(FILE *f) {
    char buf[512];
    return fgets(buf, sizeof(buf), f) != NULL;
}

/* ── main ──────────────────────────────────────────────────────────────── */
int main(void)
{
    /* ── 1. Load ground-truth CSV to find the valid window and init angles ── */
    FILE *fgt = fopen(GT_CSV, "r");
    if (!fgt) { fprintf(stderr, "Cannot open %s\n", GT_CSV); return 1; }
    skip_header(fgt);

    /* columns: timestamp_s, roll, pitch, yaw (degrees) */
    DblArr gt_t = {0}, gt_r = {0}, gt_p = {0}, gt_y = {0};
    double ts, roll, pitch, yaw_deg;
    char line[256];
    while (fgets(line, sizeof(line), fgt)) {
        if (sscanf(line, "%lf,%lf,%lf,%lf", &ts, &roll, &pitch, &yaw_deg) == 4) {
            /* skip rows where GPS heading is not yet valid (roll == 0 at start) */
            if (roll == 0.0 && gt_t.n == 0) continue;
            da_push(&gt_t, ts);
            da_push(&gt_r, roll);
            da_push(&gt_p, pitch);
            da_push(&gt_y, yaw_deg);
        }
    }
    fclose(fgt);

    if (gt_t.n == 0) { fprintf(stderr, "No valid ground-truth rows.\n"); return 1; }

    double T_START = gt_t.d[0];
    double T_END   = T_START + 600.0;   /* 10-minute window */

    fprintf(stderr, "GT window: %.1f s – %.1f s  (%zu samples)\n",
            T_START, T_END, gt_t.n);
    fprintf(stderr, "Init: roll=%.2f  pitch=%.2f  yaw=%.2f deg\n",
            gt_r.d[0], gt_p.d[0], gt_y.d[0]);

    /* ── 2. Load IMU CSV into memory (trimmed to window) ─────────────────── */
    FILE *fimu = fopen(IMU_CSV, "r");
    if (!fimu) { fprintf(stderr, "Cannot open %s\n", IMU_CSV); return 1; }
    skip_header(fimu);

    /* columns: timestamp_s, accel_x/y/z (m/s²), gyro_x/y/z (deg/s) */
    DblArr imu_t  = {0};
    DblArr imu_ax = {0}, imu_ay = {0}, imu_az = {0};
    DblArr imu_gx = {0}, imu_gy = {0}, imu_gz = {0};

    double ax, ay, az, gx, gy, gz;
    while (fgets(line, sizeof(line), fimu)) {
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &ts, &ax, &ay, &az, &gx, &gy, &gz) == 7) {
            if (ts < T_START || ts > T_END) continue;
            da_push(&imu_t,  ts);
            da_push(&imu_ax, ax);  da_push(&imu_ay, ay);  da_push(&imu_az, az);
            da_push(&imu_gx, gx);  da_push(&imu_gy, gy);  da_push(&imu_gz, gz);
        }
    }
    fclose(fimu);

    size_t N = imu_t.n;
    fprintf(stderr, "IMU samples in window: %zu\n", N);
    if (N < 2) { fprintf(stderr, "Too few IMU samples.\n"); return 1; }

    /* median dt from first 1000 diffs */
    double dt_imu;
    {
        size_t check = N < 1001 ? N - 1 : 1000;
        double dts[1001];
        for (size_t i = 0; i < check; i++)
            dts[i] = imu_t.d[i+1] - imu_t.d[i];
        /* simple selection for median */
        for (size_t i = 0; i < check - 1; i++)
            for (size_t j = i+1; j < check; j++)
                if (dts[j] < dts[i]) { double tmp=dts[i]; dts[i]=dts[j]; dts[j]=tmp; }
        dt_imu = dts[check/2];
    }
    int stride = (int)round(1.0 / (TARGET_HZ * dt_imu));
    if (stride < 1) stride = 1;
    fprintf(stderr, "dt_imu=%.3f ms  stride=%d  → ~%.1f Hz output\n",
            dt_imu*1000.0, stride, 1.0/(stride*dt_imu));

    /* ── 3. Initialise all three filters from first GT sample ───────────── */
    float init_heading_deg = (float)gt_y.d[0];
    VECTOR_3D init_accel = {
        (float)imu_ax.d[0],
        (float)imu_ay.d[0],
        (float)imu_az.d[0]
    };

    /* Use euler_to_quat logic via mathTransform_InitializeQuaternion,
       but that only sets tilt from accel and optionally yaw from heading.
       We want roll+pitch+yaw all from GT, so build the quaternion manually
       from GT Euler angles using mathTransform_EulerToQuat. */
    VECTOR_3D init_euler_rad = {
        (float)(gt_r.d[0] * PI / 180.0),   /* roll  */
        (float)(gt_p.d[0] * PI / 180.0),   /* pitch */
        (float)(gt_y.d[0] * PI / 180.0)    /* yaw   */
    };
    QUATERNION q_init = {1.0f, 0.0f, 0.0f, 0.0f};
    mathTransform_EulerToQuat(&q_init, &init_euler_rad);

    /* Madgwick */
    MADGWICK_DATA_FRAME madg = {0};
    {
        VECTOR_3D gyro_bias = {0.0f, 0.0f, 0.0f};
        madgwick_Init(&madg, &init_heading_deg, &init_accel, &gyro_bias,
                      MADGWICK_BETA, MADGWICK_ZETA);
        madg.qMadgwickQuat = q_init;  /* override with full GT init */
    }

    /* Mahony */
    MAHONY_DATA_FRAME maho = {0};
    mahony_Init(&maho, &init_heading_deg, &init_accel,
                MAHONY_BETA, MAHONY_KP, MAHONY_KI);
    maho.qMahonyQuat = q_init;

    /* EKF */
    EKF_DATA_FRAME ekf = {0};
    ekf_Init(&ekf, &init_heading_deg, &init_accel,
             EKF_GYRO_VAR, EKF_ACCEL_VAR, EKF_MAG_VAR, EKF_DIP);
    ekf.qEkfQuat = q_init;
    /* Tight initial covariance — trust the GT initialisation (mirrors attitude.py P0) */
    ekf.vCovarianceP.dROW1 = (QUATERNION){0.001f, 0.0f, 0.0f, 0.0f};
    ekf.vCovarianceP.dROW2 = (QUATERNION){0.0f, 0.001f, 0.0f, 0.0f};
    ekf.vCovarianceP.dROW3 = (QUATERNION){0.0f, 0.0f, 0.001f, 0.0f};
    ekf.vCovarianceP.dROW4 = (QUATERNION){0.0f, 0.0f, 0.0f, 0.001f};

    /* zero-magnitude mag → ekf_Update will call ekf_Correction_No_Mag */
    VECTOR_3D zero_mag = {0.0f, 0.0f, 0.0f};

    /* ── 4. Print CSV header ──────────────────────────────────────────────── */
    printf("t,"
           "madg_w,madg_x,madg_y,madg_z,madg_roll,madg_pitch,madg_yaw,"
           "maho_w,maho_x,maho_y,maho_z,maho_roll,maho_pitch,maho_yaw,"
           "ekf_w,ekf_x,ekf_y,ekf_z,ekf_roll,ekf_pitch,ekf_yaw\n");

    /* ── 5. Main loop ─────────────────────────────────────────────────────── */
    float deg = (float)(180.0 / PI);

    for (size_t idx = 0; idx < N; idx++) {
        float dt  = (idx == 0) ? (float)dt_imu : (float)(imu_t.d[idx] - imu_t.d[idx-1]);
        if (dt <= 0.0f) dt = (float)dt_imu;

        /* gyro: CSV is deg/s, algorithms expect rad/s */
        VECTOR_3D gyro = {
            (float)(imu_gx.d[idx] * PI / 180.0),
            (float)(imu_gy.d[idx] * PI / 180.0),
            (float)(imu_gz.d[idx] * PI / 180.0)
        };
        VECTOR_3D accel = {
            (float)imu_ax.d[idx],
            (float)imu_ay.d[idx],
            (float)imu_az.d[idx]
        };

        madgwick_6Dof_Update(&madg, &accel, &gyro, dt);
        mahony_6Dof_Update(&maho, &accel, &gyro, dt);
        ekf_Update(&ekf, &accel, &gyro, &zero_mag, dt);

        if ((int)idx % stride == 0) {
            QUATERNION qm = madgwick_GetQuat(&madg);
            QUATERNION qh = mahony_GetQuat(&maho);
            QUATERNION qe = ekf_GetQuat(&ekf);

            VECTOR_3D em, eh, ee;
            mathTransform_QuatToEuler(&em, &qm);
            mathTransform_QuatToEuler(&eh, &qh);
            mathTransform_QuatToEuler(&ee, &qe);

            printf("%.4f,"
                   "%.6f,%.6f,%.6f,%.6f,%.3f,%.3f,%.3f,"
                   "%.6f,%.6f,%.6f,%.6f,%.3f,%.3f,%.3f,"
                   "%.6f,%.6f,%.6f,%.6f,%.3f,%.3f,%.3f\n",
                   imu_t.d[idx],
                   qm.dW, qm.dX, qm.dY, qm.dZ,
                   em.dX*deg, em.dY*deg, em.dZ*deg,
                   qh.dW, qh.dX, qh.dY, qh.dZ,
                   eh.dX*deg, eh.dY*deg, eh.dZ*deg,
                   qe.dW, qe.dX, qe.dY, qe.dZ,
                   ee.dX*deg, ee.dY*deg, ee.dZ*deg);
        }

        if (idx % 10000 == 0)
            fprintf(stderr, "  %zu / %zu (%.0f%%)\n", idx, N, 100.0*idx/N);
    }

    da_free(&gt_t); da_free(&gt_r); da_free(&gt_p); da_free(&gt_y);
    da_free(&imu_t);
    da_free(&imu_ax); da_free(&imu_ay); da_free(&imu_az);
    da_free(&imu_gx); da_free(&imu_gy); da_free(&imu_gz);

    fprintf(stderr, "Done.\n");
    return 0;
}
