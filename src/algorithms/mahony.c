/**
 * @file mahony.c
 * @author Farzan Farhangian
 */

#include <string.h>
#include <math.h>
#include "algorithms/mahony.h"

void mahony_Init(
    MAHONY_DATA_FRAME * pDataFrame,
    const float * pdHeading,
    const VECTOR_3D * pvDeviceAccel,
    float dBeta,
    float dkP,
    float dkI)
{
	if (pDataFrame)
	{
		memset(pDataFrame, 0, sizeof(MAHONY_DATA_FRAME));

		//Initialize the quaternion using the accel data
		mathTransform_InitializeQuaternion(&(pDataFrame->qMahonyQuat), pvDeviceAccel, pdHeading);

		// Initialize the data frame
		pDataFrame->dBeta = dBeta;
		pDataFrame->dkP = dkP;
		pDataFrame->dkI = dkI;
		pDataFrame->fInitialized = true;

	}
}

void mahony_6Dof_Update(
    MAHONY_DATA_FRAME * pDataFrame,
    const VECTOR_3D * pvDeviceAccel,
    const VECTOR_3D * pvDeviceGyro,
    float dSampleTime)
{
	if (!pDataFrame->fInitialized)
	{
		float heading = 0.0f;
		mahony_Init(pDataFrame, &heading, pvDeviceAccel, 0.001f, 0.1f, 1.0f);
	}

	// Copy of Gyro Data
	VECTOR_3D vOmegaCopy;
	memcpy(&vOmegaCopy, pvDeviceGyro, sizeof(VECTOR_3D));

	// Validate Accel Data
	float dAccNorm = mathTransform_MagnitudeVector(pvDeviceAccel);
	if (dAccNorm != 0.0f)
	{
		// Initialize Variables
		MATRIX_3X3 vDcm;
		VECTOR_3D vAccNorm;
		VECTOR_3D vOmegaCost;
		VECTOR_3D vGyroBiasDelta;

		// Clear Variables
		mathTransform_ClearMatrix(&vDcm);
		mathTransform_ClearVector(&vAccNorm);
		mathTransform_ClearVector(&vOmegaCost);
		mathTransform_ClearVector(&vGyroBiasDelta);

		// Conversions and Normalizations
		mathTransform_QuatToDirectCosineMatrix(&vDcm, &(pDataFrame->qMahonyQuat));
		mathTransform_NormalizeVector(&vAccNorm, pvDeviceAccel);

		// Obtaining gyro bias
		mathTransform_CrossVectorWithVector(&vOmegaCost, &vAccNorm, &(vDcm.dROW3));
		mathTransform_MultiplyVectorByScalar(&vGyroBiasDelta, &vOmegaCost, (-pDataFrame->dkI) * dSampleTime);
		vGyroBiasDelta.dX += pDataFrame->dBeta;
		vGyroBiasDelta.dY += pDataFrame->dBeta;
		vGyroBiasDelta.dZ += pDataFrame->dBeta;

		// Correcting Gyro Bias
		mathTransform_SubtractVectorFromVector(&vOmegaCopy, &vOmegaCopy, &vGyroBiasDelta);
		mathTransform_MultiplyVectorByScalar(&vOmegaCost, &vOmegaCost, pDataFrame->dkP);
		mathTransform_AddVectorWithVector(&vOmegaCopy, &vOmegaCopy, &vOmegaCost);
	}

	// Obtain Rate of Change of Quaternion
	QUATERNION qTempQuat = {0, vOmegaCopy.dX, vOmegaCopy.dY, vOmegaCopy.dZ};
	mathTransform_MultiplyQuatByQuat(&qTempQuat, &(pDataFrame->qMahonyQuat), &qTempQuat);
	mathTransform_MultiplyQuatByScalar(&qTempQuat, &qTempQuat, 0.5f * dSampleTime);

	// Add Change to Quaternion and Normalize
	pDataFrame->qMahonyQuat.dW += qTempQuat.dW;
	pDataFrame->qMahonyQuat.dX += qTempQuat.dX;
	pDataFrame->qMahonyQuat.dY += qTempQuat.dY;
	pDataFrame->qMahonyQuat.dZ += qTempQuat.dZ;
	mathTransform_NormalizeQuat(&(pDataFrame->qMahonyQuat), &(pDataFrame->qMahonyQuat));
}

void mahony_9Dof_Update(
    MAHONY_DATA_FRAME * pDataFrame,
    const VECTOR_3D * pvDeviceAccel,
    const VECTOR_3D * pvDeviceGyro,
    const VECTOR_3D * pvDeviceMag,
    float dSampleTime)
{
	if (!pDataFrame->fInitialized)
	{
		float heading = 0.0f;
		mahony_Init(pDataFrame, &heading, pvDeviceAccel, 0.001f, 1.0f, 0.09f);
	}

	VECTOR_3D vOmegaCopy;
	memcpy(&vOmegaCopy, pvDeviceGyro, sizeof(VECTOR_3D));

	// Check if Mag Data Exists
	float dMagNorm = mathTransform_MagnitudeVector(pvDeviceMag);
	if (dMagNorm == 0.0f)
	{
		return mahony_6Dof_Update(pDataFrame, pvDeviceAccel, pvDeviceGyro, dSampleTime);
	}

	// Validate Accel Data
	float dAccNorm = mathTransform_MagnitudeVector(pvDeviceAccel);
	if (dAccNorm != 0.0f)
	{
		// Initialize Variables
		MATRIX_3X3 vDcm;
		VECTOR_3D vAccNorm;
		VECTOR_3D vMagNorm;
		VECTOR_3D vOmegaCost;
		VECTOR_3D vGyroBiasDelta;

		// Clear Variables
		mathTransform_ClearMatrix(&vDcm);
		mathTransform_ClearVector(&vAccNorm);
		mathTransform_ClearVector(&vMagNorm);
		mathTransform_ClearVector(&vOmegaCost);
		mathTransform_ClearVector(&vGyroBiasDelta);

		// Conversions and Normalizations
		mathTransform_QuatToDirectCosineMatrix(&vDcm, &(pDataFrame->qMahonyQuat));
		mathTransform_NormalizeVector(&vMagNorm, pvDeviceMag);
		mathTransform_NormalizeVector(&vAccNorm, pvDeviceAccel);

		// Rotate Magnetic Field to Inertial Frame
		VECTOR_3D vMagOrient =
		{
			.dX = vDcm.dROW1.dX * vMagNorm.dX + vDcm.dROW1.dY * vMagNorm.dY + vDcm.dROW1.dZ * vMagNorm.dZ,
			.dY = vDcm.dROW2.dX * vMagNorm.dX + vDcm.dROW2.dY * vMagNorm.dY + vDcm.dROW2.dZ * vMagNorm.dZ,
			.dZ = vDcm.dROW3.dX * vMagNorm.dX + vDcm.dROW3.dY * vMagNorm.dY + vDcm.dROW3.dZ * vMagNorm.dZ
		};
		vMagOrient = (VECTOR_3D)
		{
			0, sqrtf(vMagOrient.dX * vMagOrient.dX + vMagOrient.dY * vMagOrient.dY), vMagOrient.dZ
		};
		vMagOrient = (VECTOR_3D)
		{
			.dX = vDcm.dROW1.dX * vMagOrient.dX + vDcm.dROW2.dX * vMagOrient.dY + vDcm.dROW3.dX * vMagOrient.dZ,
			.dY = vDcm.dROW1.dY * vMagOrient.dX + vDcm.dROW2.dY * vMagOrient.dY + vDcm.dROW3.dY * vMagOrient.dZ,
			.dZ = vDcm.dROW1.dZ * vMagOrient.dX + vDcm.dROW2.dZ * vMagOrient.dY + vDcm.dROW3.dZ * vMagOrient.dZ
		};
		mathTransform_NormalizeVector(&vMagOrient, &vMagOrient);

		// Obtaining gyro bias
		mathTransform_CrossVectorWithVector(&vMagNorm, &vMagNorm, &vMagOrient);
		mathTransform_CrossVectorWithVector(&vOmegaCost, &vAccNorm, &(vDcm.dROW3));
		mathTransform_AddVectorWithVector(&vOmegaCost, &vOmegaCost, &vMagNorm);
		mathTransform_MultiplyVectorByScalar(&vGyroBiasDelta, &vOmegaCost, (-pDataFrame->dkI) * dSampleTime);
		vGyroBiasDelta.dX += pDataFrame->dBeta;
		vGyroBiasDelta.dY += pDataFrame->dBeta;
		vGyroBiasDelta.dZ += pDataFrame->dBeta;

		// Correcting Gyro Bias
		mathTransform_SubtractVectorFromVector(&vOmegaCopy, &vOmegaCopy, &vGyroBiasDelta);
		mathTransform_MultiplyVectorByScalar(&vOmegaCost, &vOmegaCost, pDataFrame->dkP);
		mathTransform_AddVectorWithVector(&vOmegaCopy, &vOmegaCopy, &vOmegaCost);
	}
	// Obtain Rate of Change of Quaternion
	QUATERNION qTempQuat = {0, vOmegaCopy.dX, vOmegaCopy.dY, vOmegaCopy.dZ};
	mathTransform_MultiplyQuatByQuat(&qTempQuat, &(pDataFrame->qMahonyQuat), &qTempQuat);
	mathTransform_MultiplyQuatByScalar(&qTempQuat, &qTempQuat, 0.5f * dSampleTime);

	// Add Change and Normalize
	pDataFrame->qMahonyQuat.dW += qTempQuat.dW;
	pDataFrame->qMahonyQuat.dX += qTempQuat.dX;
	pDataFrame->qMahonyQuat.dY += qTempQuat.dY;
	pDataFrame->qMahonyQuat.dZ += qTempQuat.dZ;
	mathTransform_NormalizeQuat(&(pDataFrame->qMahonyQuat), &(pDataFrame->qMahonyQuat));
}

QUATERNION mahony_GetQuat(MAHONY_DATA_FRAME * pDataFrame)
{
	return pDataFrame->qMahonyQuat;
}
