/**
 * @author Farzan Farhangian
 * @file ekf.h
 */

#include "orientationDefines.h"

typedef struct
{
	MATRIX_4X4 vCovarianceP;
	float dMagneticDip;
	QUATERNION qEkfQuat;
	float dGyroSigma2;
	float dAccelSigma2;
	float dMagnetSigma2;
	bool fInitialized;
	bool fGpsReceived;
} EKF_DATA_FRAME;

/* Initialize the EKF data frame and parameters */
void ekf_Init(
    EKF_DATA_FRAME * pDataFrame,
    const float * pdHeading,
    const VECTOR_3D * pvDeviceAccel,
    const float dGyroNoiseVar,
    const float dAccelNoiseVar,
    const float dMagnetNoiseVar,
    const float dDip);


/* EKF prediction step */
void ekf_Prediction(EKF_DATA_FRAME * pDataFrame, const VECTOR_3D * pvDeviceGyro, const float dSampleTime);

/* EKF correction step 9 degree of freedom */
void ekf_Correction_With_Mag(EKF_DATA_FRAME * pDataFrame, const VECTOR_3D * pvDeviceAccel, const VECTOR_3D * pvDeviceMagnet);


/**
 * @brief EKF correction step 6 degree of freedom
 */
void ekf_Correction_No_Mag(
    EKF_DATA_FRAME * pDataFrame,
    const VECTOR_3D * pvDeviceAccel);


/**
 * @brief EKF main update
 */
void ekf_Update(
    EKF_DATA_FRAME * ptEkf,
    const VECTOR_3D * vAccel,
    const VECTOR_3D * vGyro,
    const VECTOR_3D * vMag,
    const float dSampleTime);


/**
 * @brief Get quaternion from EKF dataframe
 */
QUATERNION ekf_GetQuat(EKF_DATA_FRAME * pDataFrame);

/**
 * @brief Multiply two matrices with maximum 6 dimension
 */
void ekf_MultiplyTwoMatrix(MATRIX_f * pInputA, MATRIX_f * pInputB, MATRIX_f * pOutput) ;

/**
 * @brief Create a matrix data frame
 */
void ekf_AllocateMatrix(MATRIX_f * pInput, uint16_t uiNumberOfRows, uint16_t uiNumberOfColumns);

/**
 * @brief Create a matrix data frame with double elements
 */
void ekf_AllocateMatrix_d(MATRIX_d * pInput, uint16_t uiNumberOfRows, uint16_t uiNumberOfColumns);

/**
 * @brief Create an identity matrix
 */
void ekf_eyeMatrix(MATRIX_f * pInput);

/**
 * @brief Multiply a scalar to a matrix
 */
void ekf_MultiplyScalarToMatrix(MATRIX_f * pInput, float dScalar, MATRIX_f * pOutput);

/**
 * @brief Add two matrices
 */
void ekf_AddTwoMatrix(MATRIX_f * pInputA, MATRIX_f * pInputB, MATRIX_f * pOutput) ;

/**
 * @brief Transpose of a matrix
 */
void ekf_TransposeMatrix(MATRIX_f * pInput, MATRIX_f * pOutput) ;

/**
 * @brief Obtain determinant of a matrix with special dimension gauss elimination method
 */
float ekf_MatrixDeterminant(MATRIX_f * pInput, uint16_t uiSize);

/**
 * @brief Obtain determinant of a matrix with special dimension with recursive method
 */
float ekf_RecursiveMatrixDeterminant(MATRIX_f * pInput, uint16_t uiSize);

/**
 * @brief Obtain deinverse of a matrix with special dimension
 */
void ekf_InverseMatrix(MATRIX_f * pInput, MATRIX_f * pOutput, uint16_t uiSize);

/**
 * @brief Copy a matrix to another matrix
 */
void ekf_CopyMatrix(MATRIX_f * pFromA, MATRIX_f * pToB);

/**
 * @brief Lower/Upper (LU) decomposition of a matrix
 */
float ekf_LUDeterminant(MATRIX_f * pInput, uint16_t uiSize);

/**
 * @brief Copy a matrix to another matrix
 */
void ekf_LUDecomposition(MATRIX_f * pInput, MATRIX_f * pU, MATRIX_f * pL, uint16_t uiSize);
