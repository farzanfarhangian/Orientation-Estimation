/**
 * @file madgwick.c
 * @author Farzan Farhangian
 *
 * @link
 * [Madgwick]: https://ahrs.readthedocs.io/en/latest/filters/madgwick.html
 */

#include <string.h>
#include <math.h>
#include "algorithms/madgwick.h"

void madgwick_Init(
    MADGWICK_DATA_FRAME * pDataFrame,
    const float * pdHeading,
    const VECTOR_3D * pvDeviceAccel,
    const VECTOR_3D * pvGyroscopeBias,
    const float dBeta,
    const float dZeta)
{
	memset(pDataFrame, 0, sizeof(MADGWICK_DATA_FRAME));

	//Initialize the quaternion using the accel data
	mathTransform_InitializeQuaternion(&pDataFrame->qMadgwickQuat, pvDeviceAccel, pdHeading);

	// Initialize the data frame
	pDataFrame->dMadgwickBeta = dBeta;
	pDataFrame->dMadgwickZeta = dZeta;

	pDataFrame->vGyroBias.dX = pvGyroscopeBias->dX;
	pDataFrame->vGyroBias.dY = pvGyroscopeBias->dY;
	pDataFrame->vGyroBias.dZ = pvGyroscopeBias->dZ;
	pDataFrame->fInitialized = true;
}


void madgwick_6Dof_Update(
    MADGWICK_DATA_FRAME * pDataFrame,
    const VECTOR_3D * pvDeviceAccel,
    const VECTOR_3D * pvDeviceGyro,
    const float dSampleTime)
{
	if (!pDataFrame->fInitialized)
	{
		float heading = 0.0f;
		VECTOR_3D vGyroBias = {0};
		madgwick_Init(pDataFrame, &heading, pvDeviceAccel, &vGyroBias, 0.001f, 0.001f);
	}

	// Main and auxilary variable definitions
	QUATERNION qGradient = {1, 0, 0, 0};
	QUATERNION qQuaternionDot = {1, 0, 0, 0};

	// Quaternion rate from angular rate data
	mathTransform_QuaternionDerivative(&qQuaternionDot, &pDataFrame->qMadgwickQuat, pvDeviceGyro);

	// Validate the accel data
	float dAccelMag = mathTransform_MagnitudeVector(pvDeviceAccel);
	if (dAccelMag != 0.0f)
	{
		// Normalize the accel
		VECTOR_3D accel;
		mathTransform_NormalizeVector(&accel, pvDeviceAccel);

		// Auxilary variables
		float d2Aux_W = 2.0f * pDataFrame->qMadgwickQuat.dW;
		float d2Aux_X = 2.0f * pDataFrame->qMadgwickQuat.dX;
		float d2Aux_Y = 2.0f * pDataFrame->qMadgwickQuat.dY;
		float d2Aux_Z = 2.0f * pDataFrame->qMadgwickQuat.dZ;

		float d4Aux_W = 4.0f * pDataFrame->qMadgwickQuat.dW;
		float d4Aux_X = 4.0f * pDataFrame->qMadgwickQuat.dX;
		float d4Aux_Y = 4.0f * pDataFrame->qMadgwickQuat.dY;

		float d8Aux_X = 8.0f * pDataFrame->qMadgwickQuat.dX;
		float d8Aux_Y = 8.0f * pDataFrame->qMadgwickQuat.dY;

		float dAux_WW = pDataFrame->qMadgwickQuat.dW * pDataFrame->qMadgwickQuat.dW;
		float dAux_XX = pDataFrame->qMadgwickQuat.dX * pDataFrame->qMadgwickQuat.dX;
		float dAux_YY = pDataFrame->qMadgwickQuat.dY * pDataFrame->qMadgwickQuat.dY;
		float dAux_ZZ = pDataFrame->qMadgwickQuat.dZ * pDataFrame->qMadgwickQuat.dZ;

		// Gradient Decent correction
		qGradient.dW = d4Aux_W * dAux_YY + d2Aux_Y * accel.dX + d4Aux_W * dAux_XX - d2Aux_X * accel.dY;
		qGradient.dX = d4Aux_X * dAux_ZZ - d2Aux_Z * accel.dX + 4.0f * dAux_WW * pDataFrame->qMadgwickQuat.dX - d2Aux_W * accel.dY - d4Aux_X + d8Aux_X * dAux_XX + d8Aux_X * dAux_YY + d4Aux_X * accel.dZ;
		qGradient.dY = 4.0f * dAux_WW * pDataFrame->qMadgwickQuat.dY + d2Aux_W * accel.dX + d4Aux_Y * dAux_ZZ - d2Aux_Z * accel.dY - d4Aux_Y + d8Aux_Y * dAux_XX + d8Aux_Y * dAux_YY + d4Aux_Y * accel.dZ;
		qGradient.dZ = 4.0f * dAux_XX * pDataFrame->qMadgwickQuat.dZ - d2Aux_X * accel.dX + 4.0f * dAux_YY * pDataFrame->qMadgwickQuat.dZ - d2Aux_Y * accel.dY;
		mathTransform_NormalizeQuat(&qGradient, &qGradient);

		// Apply zeta feedback gain
		qQuaternionDot.dW -= pDataFrame->dMadgwickBeta * qGradient.dW;
		qQuaternionDot.dX -= pDataFrame->dMadgwickBeta * qGradient.dX;
		qQuaternionDot.dY -= pDataFrame->dMadgwickBeta * qGradient.dY;
		qQuaternionDot.dZ -= pDataFrame->dMadgwickBeta * qGradient.dZ;

	}

	// Integrate to obtain the quaternion
	pDataFrame->qMadgwickQuat.dW += qQuaternionDot.dW * dSampleTime;
	pDataFrame->qMadgwickQuat.dX += qQuaternionDot.dX * dSampleTime;
	pDataFrame->qMadgwickQuat.dY += qQuaternionDot.dY * dSampleTime;
	pDataFrame->qMadgwickQuat.dZ += qQuaternionDot.dZ * dSampleTime;

	mathTransform_NormalizeQuat(&pDataFrame->qMadgwickQuat, &pDataFrame->qMadgwickQuat);
}



void madgwick_9Dof_Update(
    MADGWICK_DATA_FRAME * pDataFrame,
    const VECTOR_3D * pvDeviceAccel,
    const VECTOR_3D * pvDeviceGyro,
    const VECTOR_3D * pvDeviceMagnet,
    const float dSampleTime)
{
	if (!pDataFrame->fInitialized)
	{
		float heading = 0.0f;
		VECTOR_3D vGyroBias = {0};
		madgwick_Init(pDataFrame, &heading, pvDeviceAccel, &vGyroBias, 0.001f, 0.001f);
	}

	// Main and auxilary variable definitions
	QUATERNION qGradient = {1, 0, 0, 0};
	QUATERNION qQuaternionDot = {1, 0, 0, 0};

	if (!pvDeviceMagnet)
	{
		madgwick_6Dof_Update(pDataFrame,  pvDeviceAccel, pvDeviceGyro, dSampleTime);
		return;
	}

	// Validate the accel data
	float dAccelMag = mathTransform_MagnitudeVector(pvDeviceAccel);
	if (dAccelMag != 0.0f)
	{
		// Normalize the accel and magnetometer
		VECTOR_3D accel, mag;
		mathTransform_NormalizeVector(&accel, pvDeviceAccel);
		mathTransform_NormalizeVector(&mag, pvDeviceMagnet);

		// Auxilary variables
		float d2AuxMag_WX = 2.0f * pDataFrame->qMadgwickQuat.dW * mag.dX;
		float d2AuxMag_WY = 2.0f * pDataFrame->qMadgwickQuat.dW * mag.dY;
		float d2AuxMag_WZ = 2.0f * pDataFrame->qMadgwickQuat.dW * mag.dZ;
		float d2AuxMag_XX = 2.0f * pDataFrame->qMadgwickQuat.dX * mag.dX;

		float d2Aux_W = 2.0f * pDataFrame->qMadgwickQuat.dW;
		float d2Aux_X = 2.0f * pDataFrame->qMadgwickQuat.dX;
		float d2Aux_Y = 2.0f * pDataFrame->qMadgwickQuat.dY;
		float d2Aux_Z = 2.0f * pDataFrame->qMadgwickQuat.dZ;

		float d2Aux_WY = 2.0f * pDataFrame->qMadgwickQuat.dW * pDataFrame->qMadgwickQuat.dY;
		float d2Aux_YZ = 2.0f * pDataFrame->qMadgwickQuat.dY * pDataFrame->qMadgwickQuat.dZ;

		float dAux_WW = pDataFrame->qMadgwickQuat.dW * pDataFrame->qMadgwickQuat.dW;
		float dAux_WX = pDataFrame->qMadgwickQuat.dW * pDataFrame->qMadgwickQuat.dX;
		float dAux_WY = pDataFrame->qMadgwickQuat.dW * pDataFrame->qMadgwickQuat.dY;
		float dAux_WZ = pDataFrame->qMadgwickQuat.dW * pDataFrame->qMadgwickQuat.dZ;
		float dAux_XX = pDataFrame->qMadgwickQuat.dX * pDataFrame->qMadgwickQuat.dX;
		float dAux_XY = pDataFrame->qMadgwickQuat.dX * pDataFrame->qMadgwickQuat.dY;
		float dAux_XZ = pDataFrame->qMadgwickQuat.dX * pDataFrame->qMadgwickQuat.dZ;
		float dAux_YY = pDataFrame->qMadgwickQuat.dY * pDataFrame->qMadgwickQuat.dY;
		float dAux_YZ = pDataFrame->qMadgwickQuat.dY * pDataFrame->qMadgwickQuat.dZ;
		float dAux_ZZ = pDataFrame->qMadgwickQuat.dZ * pDataFrame->qMadgwickQuat.dZ;

		// Auxilary variables - Direction of Earth's magnetic field calculation
		float dAux_hX = mag.dX * dAux_WW - d2AuxMag_WY * pDataFrame->qMadgwickQuat.dZ + d2AuxMag_WZ * pDataFrame->qMadgwickQuat.dY + mag.dX * dAux_XX
		                + d2Aux_X * mag.dY * pDataFrame->qMadgwickQuat.dY + d2Aux_X * mag.dZ * pDataFrame->qMadgwickQuat.dZ - mag.dX * dAux_YY
		                - mag.dX * dAux_ZZ;

		float dAux_hY = d2AuxMag_WX * pDataFrame->qMadgwickQuat.dZ + mag.dY * dAux_WW - d2AuxMag_WZ * pDataFrame->qMadgwickQuat.dX + d2AuxMag_XX * pDataFrame->qMadgwickQuat.dY
		                - mag.dY * dAux_XX + mag.dY * dAux_YY + d2Aux_Y * mag.dZ * pDataFrame->qMadgwickQuat.dZ - mag.dY * dAux_ZZ;

		float d2Aux_bX = sqrtf(dAux_hX * dAux_hX + dAux_hY * dAux_hY);

		float d2Aux_bZ = -d2AuxMag_WX * pDataFrame->qMadgwickQuat.dY + d2AuxMag_WY * pDataFrame->qMadgwickQuat.dX + mag.dZ * dAux_WW
		                 + d2AuxMag_XX * pDataFrame->qMadgwickQuat.dZ - mag.dZ * dAux_XX + d2Aux_Y * mag.dY * pDataFrame->qMadgwickQuat.dZ
		                 - mag.dZ * dAux_YY + mag.dZ * dAux_ZZ;

		float d4Aux_bX = 2.0f * d2Aux_bX;
		float d4Aux_bZ = 2.0f * d2Aux_bZ;

		// Gradient Decent correction
		float dObjective_1 = 2.0f * dAux_XZ - d2Aux_WY - accel.dX;
		float dObjective_2 = 2.0f * dAux_WX + d2Aux_YZ - accel.dY;
		float dObjective_3 = 1 - 2.0f * dAux_XX - 2.0f * dAux_YY - accel.dZ;
		float dObjective_4 = d2Aux_bX * (0.5f - dAux_YY - dAux_ZZ) + d2Aux_bZ * (dAux_XZ - dAux_WY) - mag.dX;
		float dObjective_5 = d2Aux_bX * (dAux_XY - dAux_WZ) + d2Aux_bZ * (dAux_WX + dAux_YZ) - mag.dY;
		float dObjective_6 = d2Aux_bX * (dAux_WY + dAux_XZ) + d2Aux_bZ * (0.5f - dAux_XX - dAux_YY) - mag.dZ;

		float dJacobian_11 = -d2Aux_Y;
		float dJacobian_12 = d2Aux_Z;
		float dJacobian_13 = -d2Aux_W;
		float dJacobian_14 = d2Aux_X;

		float dJacobian_21 = d2Aux_X;
		float dJacobian_22 = d2Aux_W;
		float dJacobian_23 = d2Aux_Z;
		float dJacobian_24 = d2Aux_Y;

		float dJacobian_31 = 0.0f;
		float dJacobian_32 = - 4.0f * pDataFrame->qMadgwickQuat.dX;
		float dJacobian_33 = - 4.0f * pDataFrame->qMadgwickQuat.dY;
		float dJacobian_34 = 0.0f;

		float dJacobian_41 = - d2Aux_bZ * pDataFrame->qMadgwickQuat.dY;
		float dJacobian_42 = d2Aux_bZ * pDataFrame->qMadgwickQuat.dZ;
		float dJacobian_43 = -d4Aux_bX * pDataFrame->qMadgwickQuat.dY - d2Aux_bZ * pDataFrame->qMadgwickQuat.dW;
		float dJacobian_44 = -d4Aux_bX * pDataFrame->qMadgwickQuat.dZ + d2Aux_bZ * pDataFrame->qMadgwickQuat.dX;

		float dJacobian_51 = -d2Aux_bX * pDataFrame->qMadgwickQuat.dZ + d2Aux_bZ * pDataFrame->qMadgwickQuat.dX;
		float dJacobian_52 = d2Aux_bX * pDataFrame->qMadgwickQuat.dY + d2Aux_bZ * pDataFrame->qMadgwickQuat.dW;
		float dJacobian_53 = d2Aux_bX * pDataFrame->qMadgwickQuat.dX + d2Aux_bZ * pDataFrame->qMadgwickQuat.dZ;
		float dJacobian_54 = -d2Aux_bX * pDataFrame->qMadgwickQuat.dW + d2Aux_bZ * pDataFrame->qMadgwickQuat.dY;

		float dJacobian_61 = d2Aux_bX * pDataFrame->qMadgwickQuat.dY;
		float dJacobian_62 = d2Aux_bX * pDataFrame->qMadgwickQuat.dZ - d4Aux_bZ * pDataFrame->qMadgwickQuat.dX;
		float dJacobian_63 = d2Aux_bX * pDataFrame->qMadgwickQuat.dW - d4Aux_bZ * pDataFrame->qMadgwickQuat.dY;
		float dJacobian_64 = d2Aux_bX * pDataFrame->qMadgwickQuat.dX;

		qGradient.dW = dJacobian_11 * (dObjective_1) + dJacobian_21 * (dObjective_2) + dJacobian_31 * (dObjective_3) + dJacobian_41 * (dObjective_4) + dJacobian_51 * (dObjective_5) + dJacobian_61 * (dObjective_6);
		qGradient.dX = dJacobian_12 * (dObjective_1) + dJacobian_22 * (dObjective_2) + dJacobian_32 * (dObjective_3) + dJacobian_42 * (dObjective_4) + dJacobian_52 * (dObjective_5) + dJacobian_62 * (dObjective_6);
		qGradient.dY = dJacobian_13 * (dObjective_1) + dJacobian_23 * (dObjective_2) + dJacobian_33 * (dObjective_3) + dJacobian_43 * (dObjective_4) + dJacobian_53 * (dObjective_5) + dJacobian_63 * (dObjective_6);
		qGradient.dZ = dJacobian_14 * (dObjective_1) + dJacobian_24 * (dObjective_2) + dJacobian_34 * (dObjective_3) + dJacobian_44 * (dObjective_4) + dJacobian_54 * (dObjective_5) + dJacobian_64 * (dObjective_6);
		mathTransform_NormalizeQuat(&qGradient, &qGradient);

		// Compensate gyroscope drift
		float dBiasRateX = 2.0f * pDataFrame->qMadgwickQuat.dW * qGradient.dX - 2.0f * pDataFrame->qMadgwickQuat.dX * qGradient.dW - 2.0f * pDataFrame->qMadgwickQuat.dY * qGradient.dZ + 2.0f * pDataFrame->qMadgwickQuat.dZ * qGradient.dY;
		float dBiasRateY = 2.0f * pDataFrame->qMadgwickQuat.dW * qGradient.dY + 2.0f * pDataFrame->qMadgwickQuat.dX * qGradient.dZ - 2.0f * pDataFrame->qMadgwickQuat.dY * qGradient.dW - 2.0f * pDataFrame->qMadgwickQuat.dZ * qGradient.dX;
		float dBiasRateZ = 2.0f * pDataFrame->qMadgwickQuat.dW * qGradient.dZ - 2.0f * pDataFrame->qMadgwickQuat.dX * qGradient.dY + 2.0f * pDataFrame->qMadgwickQuat.dY * qGradient.dX - 2.0f * pDataFrame->qMadgwickQuat.dZ * qGradient.dW;

		float dBiasX = pDataFrame->vGyroBias.dX + dBiasRateX * dSampleTime * pDataFrame->dMadgwickZeta;
		float dBiasY = pDataFrame->vGyroBias.dY + dBiasRateY * dSampleTime * pDataFrame->dMadgwickZeta;
		float dBiasZ = pDataFrame->vGyroBias.dZ + dBiasRateZ * dSampleTime * pDataFrame->dMadgwickZeta;

		VECTOR_3D gyro;
		gyro.dX = pvDeviceGyro->dX - dBiasX;
		gyro.dY = pvDeviceGyro->dY - dBiasY;
		gyro.dZ = pvDeviceGyro->dZ - dBiasZ;

		mathTransform_QuaternionDerivative(&qQuaternionDot, &pDataFrame->qMadgwickQuat, &gyro);

		// Apply Beta feedback gain
		qQuaternionDot.dW -= pDataFrame->dMadgwickBeta * qGradient.dW;
		qQuaternionDot.dX -= pDataFrame->dMadgwickBeta * qGradient.dX;
		qQuaternionDot.dY -= pDataFrame->dMadgwickBeta * qGradient.dY;
		qQuaternionDot.dZ -= pDataFrame->dMadgwickBeta * qGradient.dZ;

	}
	else
	{
		// Quaternion rate from angular rate data
		mathTransform_QuaternionDerivative(&qQuaternionDot, &pDataFrame->qMadgwickQuat, pvDeviceGyro);
	}

	// Integrate to obtain the quaternion
	pDataFrame->qMadgwickQuat.dW += qQuaternionDot.dW * dSampleTime;
	pDataFrame->qMadgwickQuat.dX += qQuaternionDot.dX * dSampleTime;
	pDataFrame->qMadgwickQuat.dY += qQuaternionDot.dY * dSampleTime;
	pDataFrame->qMadgwickQuat.dZ += qQuaternionDot.dZ * dSampleTime;

	mathTransform_NormalizeQuat(&pDataFrame->qMadgwickQuat, &pDataFrame->qMadgwickQuat);

}


QUATERNION madgwick_GetQuat(MADGWICK_DATA_FRAME * pDataFrame)
{
	return pDataFrame->qMadgwickQuat;
}
