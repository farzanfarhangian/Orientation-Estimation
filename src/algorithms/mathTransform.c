/**
 * @author Farzan Farhangian
 * @file mathTransform.c
 */

#include "stdbool.h"
#include <string.h>
#include <math.h>
#include "algorithms/orientationDefines.h"

void mathTransform_ClearQuat(QUATERNION * pqInput)
{
	// zero-rotation quat is [1,0,0,0]
	pqInput->dW = 1;
	pqInput->dX = 0;
	pqInput->dY = 0;
	pqInput->dZ = 0;
}

void mathTransform_ClearVector(VECTOR_3D * pvInput)
{
	pvInput->dX = 0;
	pvInput->dY = 0;
	pvInput->dZ = 0;
}

void mathTransform_ClearMatrix(MATRIX_3X3 * pvInput)
{
	memset(pvInput, 0, sizeof(MATRIX_3X3));
}

void mathTransform_Conjugate(QUATERNION * pqOutput, const QUATERNION * pqInput)
{
	pqOutput->dW = pqInput->dW;
	pqOutput->dX = -pqInput->dX;
	pqOutput->dY = -pqInput->dY;
	pqOutput->dZ = -pqInput->dZ;
}

float mathTransform_SumOfSquaresQuat(const QUATERNION * pqInput)
{
	return ((pqInput->dW * pqInput->dW) + (pqInput->dX * pqInput->dX) + (pqInput->dY * pqInput->dY) + (pqInput->dZ * pqInput->dZ));
}

float mathTransform_SumOfSquaresVector(const VECTOR_3D * pvInput)
{
	return ((pvInput->dX * pvInput->dX) + (pvInput->dY * pvInput->dY) + (pvInput->dZ * pvInput->dZ));
}

float mathTransform_MagnitudeVector(const VECTOR_3D * pvInput)
{
	return sqrtf(mathTransform_SumOfSquaresVector(pvInput));
}

float mathTransform_MagnitudeQuat(const QUATERNION * pqInput)
{
	return sqrtf(mathTransform_SumOfSquaresQuat(pqInput));
}

void mathTransform_NormalizeQuat(QUATERNION * pqOutput, const QUATERNION * pqInput)
{
	float dMag = mathTransform_MagnitudeQuat(pqInput);
	if (dMag != 0.0f)
	{
		dMag = 1.0f / dMag; // use reciprocal so only need one divide
		pqOutput->dW = pqInput->dW * dMag;
		pqOutput->dX = pqInput->dX * dMag;
		pqOutput->dY = pqInput->dY * dMag;
		pqOutput->dZ = pqInput->dZ * dMag;
	}
}

void mathTransform_NormalizeVector(VECTOR_3D * pvOutput, const VECTOR_3D * pvInput)
{
	float dMag = mathTransform_MagnitudeVector(pvInput);
	if (dMag != 0.0f)
	{
		dMag = 1.0f / dMag; // use reciprocal so only need one divide
		pvOutput->dX = pvInput->dX * dMag;
		pvOutput->dY = pvInput->dY * dMag;
		pvOutput->dZ = pvInput->dZ * dMag;
	}
}

void mathTransform_MultiplyQuatByScalar(QUATERNION * pqOutput, const QUATERNION * pqInput, const float dScalar)
{
	// note: works with fractions if you ensure that dScalar is a float (ie 1.0/2.0), otherwise it returns all zeros
	pqOutput->dW = pqInput->dW * dScalar;
	pqOutput->dX = pqInput->dX * dScalar;
	pqOutput->dY = pqInput->dY * dScalar;
	pqOutput->dZ = pqInput->dZ * dScalar;
}

void mathTransform_Inverse(QUATERNION * pqOutput, const QUATERNION * pqInput)
{
	float dSumOfSquares = mathTransform_SumOfSquaresQuat(pqInput);
	mathTransform_Conjugate(pqOutput, pqInput);

	if (dSumOfSquares != 0.0f)
	{
		mathTransform_MultiplyQuatByScalar(pqOutput, pqOutput, (1.0f / dSumOfSquares));
	}
}

void mathTransform_MultiplyQuatByQuat(QUATERNION * pqOutput, const QUATERNION * pqInputA, const QUATERNION * pqInputB)
{
	// use a temporary quat so that `q1 = q1 * q2` is valid.
	QUATERNION qTemp = {1, 0, 0, 0};

	// find the Hamilton product (https://en.wikipedia.org/wiki/Quaternion#Hamilton_product)
	qTemp.dW = pqInputA->dW * pqInputB->dW - pqInputA->dX * pqInputB->dX - pqInputA->dY * pqInputB->dY - pqInputA->dZ * pqInputB->dZ;
	qTemp.dX = pqInputA->dW * pqInputB->dX + pqInputA->dX * pqInputB->dW + pqInputA->dY * pqInputB->dZ - pqInputA->dZ * pqInputB->dY;
	qTemp.dY = pqInputA->dW * pqInputB->dY + pqInputA->dY * pqInputB->dW + pqInputA->dZ * pqInputB->dX - pqInputA->dX * pqInputB->dZ;
	qTemp.dZ = pqInputA->dW * pqInputB->dZ + pqInputA->dZ * pqInputB->dW + pqInputA->dX * pqInputB->dY - pqInputA->dY * pqInputB->dX;

	memmove(pqOutput, &qTemp, sizeof(QUATERNION));
}

void mathTransform_MultiplyVectorByScalar(VECTOR_3D * pvOutput, const VECTOR_3D * pvInput, const float dScalar)
{
	// note: works with fractions if you ensure that dScalar is a float (ie 1.0/2.0), otherwise it returns all zeros
	pvOutput->dX = pvInput->dX * dScalar;
	pvOutput->dY = pvInput->dY * dScalar;
	pvOutput->dZ = pvInput->dZ * dScalar;
}

void mathTransform_CrossVectorWithVector(VECTOR_3D * pvOutput, const VECTOR_3D * pvInputA, const VECTOR_3D * pvInputB)
{
	// use a temporary vector so that `v1 = v1 x v2` is valid.
	VECTOR_3D vTemp = {0, 0, 0};

	// find the cross product (https://en.wikipedia.org/wiki/Cross_product#Computing_the_cross_product)
	vTemp.dX = pvInputA->dY * pvInputB->dZ - pvInputA->dZ * pvInputB->dY;
	vTemp.dY = pvInputA->dZ * pvInputB->dX - pvInputA->dX * pvInputB->dZ;
	vTemp.dZ = pvInputA->dX * pvInputB->dY - pvInputA->dY * pvInputB->dX;

	memmove(pvOutput, &vTemp, sizeof(VECTOR_3D));
}

void mathTransform_GetVectorComponent(VECTOR_3D * pvOutput, const QUATERNION * pqInput)
{
	pvOutput->dX = pqInput->dX;
	pvOutput->dY = pqInput->dY;
	pvOutput->dZ = pqInput->dZ;
}

void mathTransform_AddVectorWithVector(VECTOR_3D * pvOutput, const VECTOR_3D * pvInputA, const VECTOR_3D * pvInputB)
{
	pvOutput->dX = pvInputA->dX + pvInputB->dX;
	pvOutput->dY = pvInputA->dY + pvInputB->dY;
	pvOutput->dZ = pvInputA->dZ + pvInputB->dZ;
}

void mathTransform_SubtractVectorFromVector(VECTOR_3D * pvOutput, const VECTOR_3D * pvInputA, const VECTOR_3D * pvInputB)
{
	pvOutput->dX = pvInputA->dX - pvInputB->dX;
	pvOutput->dY = pvInputA->dY - pvInputB->dY;
	pvOutput->dZ = pvInputA->dZ - pvInputB->dZ;
}

/*
 * Spherical Linear intERPolation
 * Optimized formula based on Doom 3 implementation: https://fabiensanglard.net/doom3_documentation/37725-293747_293747.pdf
 */
void mathTransform_Slerp(QUATERNION * pqOutput, const QUATERNION * pqFrom, const QUATERNION * pqTo, const float dGain)
{
	float dCosOmega;
	float dAbsCosOmega;
	float dSinSqrOmega;
	float dCscOmega;
	float dOmega;
	float dScale0;
	float dScale1;

	// get cos(omega) using dot product
	dCosOmega = pqFrom->dW * pqTo->dW + pqFrom->dX * pqTo->dX + pqFrom->dY * pqTo->dY + pqFrom->dZ * pqTo->dZ;
	dAbsCosOmega = fabsf(dCosOmega);

	// use trig identity to avoid expensive acos(omega) call:
	// sinSqr(omega) + cosSqr(omega) = 1 therefore...
	dSinSqrOmega = 1.0f - dAbsCosOmega * dAbsCosOmega;

	// check for small angle and div0
	if ((1.0f - dAbsCosOmega) > 1e-6f && dSinSqrOmega != 0.0f)
	{
		dCscOmega = 1.0f / sqrtf(dSinSqrOmega); // recall: csc(omega) = 1 / sin(omega)
		dOmega = atan2f(dSinSqrOmega * dCscOmega, dAbsCosOmega); // acos(omega) avoided!

		// SLERP
		dScale0 = sinf((1.0f - dGain) * dOmega) * dCscOmega;
		dScale1 = sinf(dGain * dOmega) * dCscOmega;
	}
	else
	{
		// Use LERP for small angles or when dSinSqrOmega <= 0
		// note: LERP is not constant velocity so this will introduce error at large angles.
		dScale0 = 1.0f - dGain;
		dScale1 = dGain;
	}

	if (dCosOmega < 0.0f)
	{
		dScale1 *= -1.0f;
	}

	pqOutput->dW = dScale0 * pqFrom->dW + dScale1 * pqTo->dW;
	pqOutput->dX = dScale0 * pqFrom->dX + dScale1 * pqTo->dX;
	pqOutput->dY = dScale0 * pqFrom->dY + dScale1 * pqTo->dY;
	pqOutput->dZ = dScale0 * pqFrom->dZ + dScale1 * pqTo->dZ;

	mathTransform_NormalizeQuat(pqOutput, pqOutput);
}

float mathTransform_dotVectorWithVector(const VECTOR_3D * pvInputA, const VECTOR_3D * pvInputB)
{
	return (pvInputA->dX * pvInputB->dX + pvInputA->dY * pvInputB->dY + pvInputA->dZ * pvInputB->dZ);
}

void mathTransform_DirectCosineMatrixToQuat(QUATERNION * pqOutput, const MATRIX_3X3 * pvInput)
{
	QUATERNION qTempQuat = {1, 0, 0, 0};
	float dTrace = pvInput->dROW1.dX + pvInput->dROW2.dY + pvInput->dROW3.dZ;

	if (dTrace > 0)
	{
		float dTemp = 0.5f / sqrtf(dTrace + 1.0f);
		qTempQuat.dW = 0.25f / dTemp;
		qTempQuat.dX = (pvInput->dROW3.dY - pvInput->dROW2.dZ) * dTemp;
		qTempQuat.dY = (pvInput->dROW1.dZ - pvInput->dROW3.dX) * dTemp;
		qTempQuat.dZ = (pvInput->dROW2.dX - pvInput->dROW1.dY) * dTemp;
	}

	else
	{
		if (pvInput->dROW1.dX > pvInput->dROW2.dY && pvInput->dROW1.dX > pvInput->dROW3.dZ)
		{
			float dTemp = 2.0f * sqrtf(1.0f + pvInput->dROW1.dX - pvInput->dROW2.dY - pvInput->dROW3.dZ);
			qTempQuat.dW = (pvInput->dROW3.dY - pvInput->dROW2.dZ) / dTemp;
			qTempQuat.dX = 0.25f * dTemp;
			qTempQuat.dY = (pvInput->dROW1.dY + pvInput->dROW2.dX) / dTemp;
			qTempQuat.dZ = (pvInput->dROW1.dZ + pvInput->dROW3.dX) / dTemp;
		}

		else if (pvInput->dROW2.dY > pvInput->dROW3.dZ)
		{
			float dTemp = 2.0f * sqrtf(1.0f + pvInput->dROW2.dY - pvInput->dROW1.dX - pvInput->dROW3.dZ);
			qTempQuat.dW = (pvInput->dROW1.dZ - pvInput->dROW3.dX) / dTemp;
			qTempQuat.dX = (pvInput->dROW1.dY + pvInput->dROW2.dX) / dTemp;
			qTempQuat.dY = 0.25f * dTemp;
			qTempQuat.dZ = (pvInput->dROW2.dZ + pvInput->dROW3.dY) / dTemp;
		}

		else
		{
			float dTemp = 2.0f * sqrtf(1.0f + pvInput->dROW3.dZ - pvInput->dROW1.dX - pvInput->dROW2.dY);
			qTempQuat.dW = (pvInput->dROW2.dX - pvInput->dROW1.dY) / dTemp;
			qTempQuat.dX = (pvInput->dROW1.dZ + pvInput->dROW3.dX) / dTemp;
			qTempQuat.dY = (pvInput->dROW2.dZ + pvInput->dROW3.dY) / dTemp;
			qTempQuat.dZ = 0.25f * dTemp;
		}
	}

	mathTransform_NormalizeQuat(pqOutput, &qTempQuat);
}

void mathTransform_QuatToEuler(VECTOR_3D * pvOutput, const QUATERNION * pqInput)
{
	// calculate Roll angle (x-axis rotation)
	float dTemp_1 = 2.0f * (pqInput->dW * pqInput->dX + pqInput->dY * pqInput->dZ);
	float dTemp_2 = 1.0f - 2.0f * (pqInput->dX * pqInput->dX + pqInput->dY * pqInput->dY);

	pvOutput->dX = atan2f(dTemp_1, dTemp_2);

	// calculate Pitch angle (y-axis rotation)
	float sinp = 2.0f * (pqInput->dW * pqInput->dY - pqInput->dZ * pqInput->dX);

	if (fabsf(sinp) >= 1.0f)
	{
		pvOutput->dY = copysignf(PI / 2.0f, sinp);
	}
	else
	{
		pvOutput->dY = asinf(sinp);
	}

	// calculate yaw angle (z-axis rotation)
	dTemp_1 = 2.0f * (pqInput->dW * pqInput->dZ + pqInput->dX * pqInput->dY);
	dTemp_2 = 1.0f - 2.0f * (pqInput->dY * pqInput->dY + pqInput->dZ * pqInput->dZ);

	pvOutput->dZ = atan2f(dTemp_1, dTemp_2);
}

void mathTransform_QuatToDirectCosineMatrix(MATRIX_3X3 * pvOutput, const QUATERNION * pqInput)
{
	QUATERNION  pqNormalizedIn = {1, 0, 0, 0};
	mathTransform_NormalizeQuat(&pqNormalizedIn, pqInput);

	pvOutput->dROW1.dX = 1.0f - 2.0f * (powf(pqNormalizedIn.dY, 2) + powf(pqNormalizedIn.dZ, 2));
	pvOutput->dROW1.dY = 2.0f * (pqNormalizedIn.dX * pqNormalizedIn.dY - pqNormalizedIn.dW * pqNormalizedIn.dZ);
	pvOutput->dROW1.dZ = 2.0f * (pqNormalizedIn.dW * pqNormalizedIn.dY + pqNormalizedIn.dX * pqNormalizedIn.dZ);

	pvOutput->dROW2.dX = 2.0f * (pqNormalizedIn.dX * pqNormalizedIn.dY + pqNormalizedIn.dW * pqNormalizedIn.dZ);
	pvOutput->dROW2.dY = 1.0f - 2.0f * (powf(pqNormalizedIn.dX, 2) + powf(pqNormalizedIn.dZ, 2));
	pvOutput->dROW2.dZ = 2.0f * (pqNormalizedIn.dY * pqNormalizedIn.dZ - pqNormalizedIn.dW * pqNormalizedIn.dX);

	pvOutput->dROW3.dX = 2.0f * (pqNormalizedIn.dX * pqNormalizedIn.dZ - pqNormalizedIn.dW * pqNormalizedIn.dY);
	pvOutput->dROW3.dY = 2.0f * (pqNormalizedIn.dW * pqNormalizedIn.dX + pqNormalizedIn.dY * pqNormalizedIn.dZ);
	pvOutput->dROW3.dZ = 1.0f - 2.0f * (powf(pqNormalizedIn.dX, 2) + powf(pqNormalizedIn.dY, 2));
}

void mathTransform_EulerToQuat(QUATERNION * pqOutput, const VECTOR_3D * pvInput)
{
	float dCosYaw   = cosf(pvInput->dZ * 0.5f);
	float dSinYaw   = sinf(pvInput->dZ * 0.5f);
	float dCosPitch = cosf(pvInput->dY * 0.5f);
	float dSinPitch = sinf(pvInput->dY * 0.5f);
	float dCosRoll  = cosf(pvInput->dX * 0.5f);
	float dSinRoll  = sinf(pvInput->dX * 0.5f);

	pqOutput->dW = dCosRoll * dCosPitch * dCosYaw + dSinRoll * dSinPitch * dSinYaw;
	pqOutput->dX = dSinRoll * dCosPitch * dCosYaw - dCosRoll * dSinPitch * dSinYaw;
	pqOutput->dY = dCosRoll * dSinPitch * dCosYaw + dSinRoll * dCosPitch * dSinYaw;
	pqOutput->dZ = dCosRoll * dCosPitch * dSinYaw - dSinRoll * dSinPitch * dCosYaw;
}

void mathTransform_VectorToSkewSymmetric(MATRIX_3X3 * pvOutput, const VECTOR_3D * pvInput)
{
	pvOutput->dROW1.dX = 0.0;
	pvOutput->dROW1.dY = -pvInput->dZ;
	pvOutput->dROW1.dZ = pvInput->dY;

	pvOutput->dROW2.dX = pvInput->dZ;
	pvOutput->dROW2.dY = 0.0;
	pvOutput->dROW2.dZ = -pvInput->dX;

	pvOutput->dROW3.dX = -pvInput->dY;
	pvOutput->dROW3.dY = pvInput->dX;
	pvOutput->dROW3.dZ = 0.0;
}

void mathTransform_QuaternionDerivative(QUATERNION * pqQuaternionDot, const QUATERNION * pqQuaternion, const VECTOR_3D * pvDeviceGyro)
{
	pqQuaternionDot->dW = 0.5f * (-pqQuaternion->dX * pvDeviceGyro->dX - pqQuaternion->dY * pvDeviceGyro->dY - pqQuaternion->dZ * pvDeviceGyro->dZ);
	pqQuaternionDot->dX = 0.5f * (pqQuaternion->dW * pvDeviceGyro->dX + pqQuaternion->dY * pvDeviceGyro->dZ - pqQuaternion->dZ * pvDeviceGyro->dY);
	pqQuaternionDot->dY = 0.5f * (pqQuaternion->dW * pvDeviceGyro->dY - pqQuaternion->dX * pvDeviceGyro->dZ + pqQuaternion->dZ * pvDeviceGyro->dX);
	pqQuaternionDot->dZ = 0.5f * (pqQuaternion->dW * pvDeviceGyro->dZ + pqQuaternion->dX * pvDeviceGyro->dY - pqQuaternion->dY * pvDeviceGyro->dX);
}

void mathTransform_InitializeQuaternion(QUATERNION * pqOutput, const VECTOR_3D * pvDeviceAccel, const float * pdHeading)
{
	QUATERNION qTemp = {1, 0, 0, 0};
	VECTOR_3D vEulerAngles = {0, 0, 0};
	VECTOR_3D vAccelNormalized = {0, 0, 0};

	float dMag = mathTransform_MagnitudeVector(pvDeviceAccel);
	if (dMag != 0.0f)
	{
		mathTransform_NormalizeVector(&vAccelNormalized, pvDeviceAccel);

		vEulerAngles.dX = atan2f(vAccelNormalized.dY, vAccelNormalized.dZ);
		vEulerAngles.dY = atan2f(-vAccelNormalized.dX, sqrtf(powf(vAccelNormalized.dY, 2) + powf(vAccelNormalized.dZ, 2)));

		if (pdHeading)
		{
			vEulerAngles.dZ = (*pdHeading) * (PI / 180.0f);
		}

		mathTransform_EulerToQuat(&qTemp, &vEulerAngles);
	}

	memcpy(pqOutput, &qTemp, sizeof(QUATERNION));
	mathTransform_NormalizeQuat(pqOutput, pqOutput);
}

