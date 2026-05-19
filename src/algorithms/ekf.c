/**
 * @author Farzan Farhangian
 * @file ekf.c
 */


#include <string.h>
#include <assert.h>
#include <math.h>
#include "algorithms/ekf.h"

void ekf_Init(
    EKF_DATA_FRAME * pDataFrame,
    const float * pdHeading,
    const VECTOR_3D * pvDeviceAccel,
    const float dGyroNoiseVar,
    const float dAccelNoiseVar,
    const float dMagnetNoiseVar,
    const float dDip)
{
	memset(pDataFrame, 0, sizeof(EKF_DATA_FRAME));

	//Initialize the quaternion using the accel data
	mathTransform_InitializeQuaternion(&pDataFrame->qEkfQuat, pvDeviceAccel, pdHeading);

	// Initialize the data frame
	pDataFrame->dMagneticDip = dDip;

	pDataFrame->dAccelSigma2 = dAccelNoiseVar;
	pDataFrame->dGyroSigma2 = dGyroNoiseVar;
	pDataFrame->dMagnetSigma2 = dMagnetNoiseVar;

	pDataFrame->vCovarianceP.dROW1.dW = 1.0f;
	pDataFrame->vCovarianceP.dROW2.dX = 1.0f;
	pDataFrame->vCovarianceP.dROW3.dY = 1.0f;
	pDataFrame->vCovarianceP.dROW4.dZ = 1.0f;

	pDataFrame->fInitialized = true;
}


void ekf_Prediction(EKF_DATA_FRAME * pDataFrame, const VECTOR_3D * pvDeviceGyro, const float dSampleTime)
{
	float dTemp = dSampleTime / 2.0f;

	// 1. Calculate the quaternion in this iteration:
	pDataFrame->qEkfQuat.dW += (dTemp) * (-pvDeviceGyro->dX * pDataFrame->qEkfQuat.dX - pvDeviceGyro->dY * pDataFrame->qEkfQuat.dY - pvDeviceGyro->dZ * pDataFrame->qEkfQuat.dZ);
	pDataFrame->qEkfQuat.dX += (dTemp) * (pvDeviceGyro->dX * pDataFrame->qEkfQuat.dW - pvDeviceGyro->dY * pDataFrame->qEkfQuat.dZ + pvDeviceGyro->dZ * pDataFrame->qEkfQuat.dY);
	pDataFrame->qEkfQuat.dY += (dTemp) * (pvDeviceGyro->dX * pDataFrame->qEkfQuat.dZ + pvDeviceGyro->dY * pDataFrame->qEkfQuat.dW - pvDeviceGyro->dZ * pDataFrame->qEkfQuat.dX);
	pDataFrame->qEkfQuat.dZ += (dTemp) * (-pvDeviceGyro->dX * pDataFrame->qEkfQuat.dY + pvDeviceGyro->dY * pDataFrame->qEkfQuat.dX + pvDeviceGyro->dZ * pDataFrame->qEkfQuat.dW);

	// 2. Calculate the discretized F matrix
	MATRIX_f dF;
	ekf_AllocateMatrix(&dF, 4, 4);

	dF.matrix[0][0] =   1.0f;
	dF.matrix[0][1] = - pvDeviceGyro->dX * dTemp;
	dF.matrix[0][2] = - pvDeviceGyro->dY * dTemp;
	dF.matrix[0][3] = - pvDeviceGyro->dZ * dTemp;

	dF.matrix[1][0] =   pvDeviceGyro->dX * dTemp;
	dF.matrix[1][1] =   1.0f;
	dF.matrix[1][2] =   pvDeviceGyro->dZ * dTemp;
	dF.matrix[1][3] = - pvDeviceGyro->dY * dTemp;

	dF.matrix[2][0] =   pvDeviceGyro->dY * dTemp;
	dF.matrix[2][1] = - pvDeviceGyro->dZ * dTemp;
	dF.matrix[2][2] =   1.0f;
	dF.matrix[2][3] =   pvDeviceGyro->dX * dTemp;

	dF.matrix[3][0] =   pvDeviceGyro->dZ * dTemp;
	dF.matrix[3][1] =   pvDeviceGyro->dY * dTemp;
	dF.matrix[3][2] = - pvDeviceGyro->dX * dTemp;
	dF.matrix[3][3] =   1.0f;

	// 3. Calculate Q
	MATRIX_f dQ, dW, dW_Transpose;
	ekf_AllocateMatrix(&dQ, 4, 4);
	ekf_AllocateMatrix(&dW, 4, 3);
	ekf_AllocateMatrix(&dW_Transpose, 3, 4);

	dW.matrix[0][0] = - pDataFrame->qEkfQuat.dX * dTemp;
	dW.matrix[0][1] = - pDataFrame->qEkfQuat.dY * dTemp;
	dW.matrix[0][2] = - pDataFrame->qEkfQuat.dZ * dTemp;

	dW.matrix[1][0] =   pDataFrame->qEkfQuat.dW * dTemp;
	dW.matrix[1][1] = - pDataFrame->qEkfQuat.dZ * dTemp;
	dW.matrix[1][2] =   pDataFrame->qEkfQuat.dY * dTemp;

	dW.matrix[2][0] =   pDataFrame->qEkfQuat.dZ * dTemp;
	dW.matrix[2][1] =   pDataFrame->qEkfQuat.dW * dTemp;
	dW.matrix[2][2] = - pDataFrame->qEkfQuat.dX * dTemp;

	dW.matrix[3][0] = - pDataFrame->qEkfQuat.dY * dTemp;
	dW.matrix[3][1] =   pDataFrame->qEkfQuat.dX * dTemp;
	dW.matrix[3][2] =   pDataFrame->qEkfQuat.dW * dTemp;

	ekf_TransposeMatrix(&dW, &dW_Transpose);
	ekf_MultiplyTwoMatrix(&dW, &dW_Transpose, &dQ);
	ekf_MultiplyScalarToMatrix(&dQ, pDataFrame->dGyroSigma2, &dQ);

	// 4. Calculate the P primary matrix
	// P = FPF' + Q

	MATRIX_f dP, dF_Transpose;
	ekf_AllocateMatrix(&dP, 4, 4);
	ekf_AllocateMatrix(&dF_Transpose, 4, 4);

	dP.matrix[0][0] = pDataFrame->vCovarianceP.dROW1.dW;
	dP.matrix[0][1] = pDataFrame->vCovarianceP.dROW1.dX;
	dP.matrix[0][2] = pDataFrame->vCovarianceP.dROW1.dY;
	dP.matrix[0][3] = pDataFrame->vCovarianceP.dROW1.dZ;

	dP.matrix[1][0] = pDataFrame->vCovarianceP.dROW2.dW;
	dP.matrix[1][1] = pDataFrame->vCovarianceP.dROW2.dX;
	dP.matrix[1][2] = pDataFrame->vCovarianceP.dROW2.dY;
	dP.matrix[1][3] = pDataFrame->vCovarianceP.dROW2.dZ;

	dP.matrix[2][0] = pDataFrame->vCovarianceP.dROW3.dW;
	dP.matrix[2][1] = pDataFrame->vCovarianceP.dROW3.dX;
	dP.matrix[2][2] = pDataFrame->vCovarianceP.dROW3.dY;
	dP.matrix[2][3] = pDataFrame->vCovarianceP.dROW3.dZ;

	dP.matrix[3][0] = pDataFrame->vCovarianceP.dROW4.dW;
	dP.matrix[3][1] = pDataFrame->vCovarianceP.dROW4.dX;
	dP.matrix[3][2] = pDataFrame->vCovarianceP.dROW4.dY;
	dP.matrix[3][3] = pDataFrame->vCovarianceP.dROW4.dZ;

	ekf_TransposeMatrix(&dF, &dF_Transpose);   //F'
	ekf_MultiplyTwoMatrix(&dF, &dP, &dF);       //FP

	ekf_MultiplyTwoMatrix(&dF, &dF_Transpose, &dF);       //FPF'
	ekf_AddTwoMatrix(&dF, &dQ, &dP);

	pDataFrame->vCovarianceP.dROW1.dW = dP.matrix[0][0];
	pDataFrame->vCovarianceP.dROW1.dX = dP.matrix[0][1];
	pDataFrame->vCovarianceP.dROW1.dY = dP.matrix[0][2];
	pDataFrame->vCovarianceP.dROW1.dZ = dP.matrix[0][3];

	pDataFrame->vCovarianceP.dROW2.dW = dP.matrix[1][0];
	pDataFrame->vCovarianceP.dROW2.dX = dP.matrix[1][1];
	pDataFrame->vCovarianceP.dROW2.dY = dP.matrix[1][2];
	pDataFrame->vCovarianceP.dROW2.dZ = dP.matrix[1][3];

	pDataFrame->vCovarianceP.dROW3.dW = dP.matrix[2][0];
	pDataFrame->vCovarianceP.dROW3.dX = dP.matrix[2][1];
	pDataFrame->vCovarianceP.dROW3.dY = dP.matrix[2][2];
	pDataFrame->vCovarianceP.dROW3.dZ = dP.matrix[2][3];

	pDataFrame->vCovarianceP.dROW4.dW = dP.matrix[3][0];
	pDataFrame->vCovarianceP.dROW4.dX = dP.matrix[3][1];
	pDataFrame->vCovarianceP.dROW4.dY = dP.matrix[3][2];
	pDataFrame->vCovarianceP.dROW4.dZ = dP.matrix[3][3];
}


void ekf_Correction_With_Mag(
    EKF_DATA_FRAME * pDataFrame,
    const VECTOR_3D * pvDeviceAccel,
    const VECTOR_3D * pvDeviceMagnet)
{

	// 1. reference observation form magnetometer and accelerometer

	VECTOR_3D vAccelRef = {0.0, 0.0, 1.0};
	VECTOR_3D vMagnetRef = {0.0, 0.0, 0.0};

	vMagnetRef.dX = 0;
	vMagnetRef.dY = cosf(pDataFrame->dMagneticDip);// / sqrtf(cosf(pDataFrame->dMagneticDip)*cosf(pDataFrame->dMagneticDip) + sinf(pDataFrame->dMagneticDip)*sinf(pDataFrame->dMagneticDip));
	vMagnetRef.dZ = -sinf(pDataFrame->dMagneticDip);// / sqrtf(cosf(pDataFrame->dMagneticDip)*cosf(pDataFrame->dMagneticDip) + sinf(pDataFrame->dMagneticDip)*sinf(pDataFrame->dMagneticDip));

	mathTransform_NormalizeVector(&vMagnetRef, &vMagnetRef);
	mathTransform_NormalizeVector(&vAccelRef, &vAccelRef);

	// 2. Convert the reference observations to body frame and generate y vector

	MATRIX_3X3 vDCM;
	VECTOR_3D     vObservation1;
	VECTOR_3D     vObservation2;

	mathTransform_ClearMatrix(&vDCM);
	mathTransform_ClearVector(&vObservation1);
	mathTransform_ClearVector(&vObservation2);

	mathTransform_QuatToDirectCosineMatrix(&vDCM, &pDataFrame->qEkfQuat);

	vObservation1.dX = vAccelRef.dX * vDCM.dROW1.dX + vAccelRef.dY * vDCM.dROW2.dX + vAccelRef.dZ * vDCM.dROW3.dX;
	vObservation1.dY = vAccelRef.dX * vDCM.dROW1.dY + vAccelRef.dY * vDCM.dROW2.dY + vAccelRef.dZ * vDCM.dROW3.dY;
	vObservation1.dZ = vAccelRef.dX * vDCM.dROW1.dZ + vAccelRef.dY * vDCM.dROW2.dZ + vAccelRef.dZ * vDCM.dROW3.dZ;

	vObservation2.dX = vMagnetRef.dX * vDCM.dROW1.dX + vMagnetRef.dY * vDCM.dROW2.dX + vMagnetRef.dZ * vDCM.dROW3.dX;
	vObservation2.dY = vMagnetRef.dX * vDCM.dROW1.dY + vMagnetRef.dY * vDCM.dROW2.dY + vMagnetRef.dZ * vDCM.dROW3.dY;
	vObservation2.dZ = vMagnetRef.dX * vDCM.dROW1.dZ + vMagnetRef.dY * vDCM.dROW2.dZ + vMagnetRef.dZ * vDCM.dROW3.dZ;

	mathTransform_NormalizeVector(&vObservation1, &vObservation1);
	mathTransform_NormalizeVector(&vObservation2, &vObservation2);

	// 3. Obtain innovation vector and measurement vector (v = z - y)
	VECTOR_3D     vMeasurement1 = {0, 0, 0};
	VECTOR_3D     vMeasurement2 = {0, 0, 0};
	VECTOR_3D     vInnovation1 =  {0, 0, 0};
	VECTOR_3D     vInnovation2 =  {0, 0, 0};

	VECTOR_3D     pvDeviceAccel_norm = {0, 0, 0};
	mathTransform_NormalizeVector(&pvDeviceAccel_norm, pvDeviceAccel);
	vMeasurement1.dX = pvDeviceAccel_norm.dX;
	vMeasurement1.dY = pvDeviceAccel_norm.dY;
	vMeasurement1.dZ = pvDeviceAccel_norm.dZ;

	VECTOR_3D     pvDeviceMagnet_norm = {0, 0, 0};
	mathTransform_NormalizeVector(&pvDeviceMagnet_norm, pvDeviceMagnet);
	vMeasurement2.dX = pvDeviceMagnet_norm.dX;
	vMeasurement2.dY = pvDeviceMagnet_norm.dY;
	vMeasurement2.dZ = pvDeviceMagnet_norm.dZ;

	vInnovation1.dX = vMeasurement1.dX - vObservation1.dX;
	vInnovation1.dY = vMeasurement1.dY - vObservation1.dY;
	vInnovation1.dZ = vMeasurement1.dZ - vObservation1.dZ;

	vInnovation2.dX = vMeasurement2.dX - vObservation2.dX;
	vInnovation2.dY = vMeasurement2.dY - vObservation2.dY;
	vInnovation2.dZ = vMeasurement2.dZ - vObservation2.dZ;

	// 4. Obtain the measurement matrix H

	MATRIX_f dH;
	ekf_AllocateMatrix(&dH, 6, 4);

	dH.matrix[0][0] = -pDataFrame->qEkfQuat.dY * vAccelRef.dZ + pDataFrame->qEkfQuat.dZ * vAccelRef.dY; //checked
	dH.matrix[0][1] =  pDataFrame->qEkfQuat.dY * vAccelRef.dY + pDataFrame->qEkfQuat.dZ * vAccelRef.dZ; //checked
	dH.matrix[0][2] = -pDataFrame->qEkfQuat.dW * vAccelRef.dZ + pDataFrame->qEkfQuat.dX * vAccelRef.dY - 2.0f * pDataFrame->qEkfQuat.dY * vAccelRef.dX; //checked
	dH.matrix[0][3] =  pDataFrame->qEkfQuat.dW * vAccelRef.dY + pDataFrame->qEkfQuat.dX * vAccelRef.dZ - 2.0f * pDataFrame->qEkfQuat.dZ * vAccelRef.dX; //checked

	dH.matrix[1][0] =  pDataFrame->qEkfQuat.dX * vAccelRef.dZ - pDataFrame->qEkfQuat.dZ * vAccelRef.dX; //checked
	dH.matrix[1][1] =  pDataFrame->qEkfQuat.dW * vAccelRef.dZ - 2.0f * pDataFrame->qEkfQuat.dX * vAccelRef.dY + pDataFrame->qEkfQuat.dY * vAccelRef.dX; //checked
	dH.matrix[1][2] =  pDataFrame->qEkfQuat.dX * vAccelRef.dX + pDataFrame->qEkfQuat.dZ * vAccelRef.dZ; //checked
	dH.matrix[1][3] = -pDataFrame->qEkfQuat.dW * vAccelRef.dX + pDataFrame->qEkfQuat.dY * vAccelRef.dZ - 2.0f * pDataFrame->qEkfQuat.dZ * vAccelRef.dY; //checked

	dH.matrix[2][0] = -pDataFrame->qEkfQuat.dX * vAccelRef.dY + pDataFrame->qEkfQuat.dY * vAccelRef.dX; //checked
	dH.matrix[2][1] = -pDataFrame->qEkfQuat.dW * vAccelRef.dY - 2.0f * pDataFrame->qEkfQuat.dX * vAccelRef.dZ + pDataFrame->qEkfQuat.dZ * vAccelRef.dX; //checked
	dH.matrix[2][2] =  pDataFrame->qEkfQuat.dW * vAccelRef.dX - 2.0f * pDataFrame->qEkfQuat.dY * vAccelRef.dZ + pDataFrame->qEkfQuat.dZ * vAccelRef.dY; //checked
	dH.matrix[2][3] =  pDataFrame->qEkfQuat.dX * vAccelRef.dX + pDataFrame->qEkfQuat.dY * vAccelRef.dY; //checked

	dH.matrix[3][0] = -pDataFrame->qEkfQuat.dY * vMagnetRef.dZ + pDataFrame->qEkfQuat.dZ * vMagnetRef.dY; //checked
	dH.matrix[3][1] =  pDataFrame->qEkfQuat.dY * vMagnetRef.dY + pDataFrame->qEkfQuat.dZ * vMagnetRef.dZ; //checked
	dH.matrix[3][2] = -pDataFrame->qEkfQuat.dW * vMagnetRef.dZ + pDataFrame->qEkfQuat.dX * vMagnetRef.dY - 2.0f * pDataFrame->qEkfQuat.dY * vMagnetRef.dX; //checked
	dH.matrix[3][3] =  pDataFrame->qEkfQuat.dW * vMagnetRef.dY + pDataFrame->qEkfQuat.dX * vMagnetRef.dZ - 2.0f * pDataFrame->qEkfQuat.dZ * vMagnetRef.dX; //checked

	dH.matrix[4][0] =  pDataFrame->qEkfQuat.dX * vMagnetRef.dZ - pDataFrame->qEkfQuat.dZ * vMagnetRef.dX; //checked
	dH.matrix[4][1] =  pDataFrame->qEkfQuat.dW * vMagnetRef.dZ - 2.0f * pDataFrame->qEkfQuat.dX * vMagnetRef.dY + pDataFrame->qEkfQuat.dY * vMagnetRef.dX; //checked
	dH.matrix[4][2] =  pDataFrame->qEkfQuat.dX * vMagnetRef.dX + pDataFrame->qEkfQuat.dZ * vMagnetRef.dZ; //checked
	dH.matrix[4][3] = -pDataFrame->qEkfQuat.dW * vMagnetRef.dX + pDataFrame->qEkfQuat.dY * vMagnetRef.dZ - 2.0f * pDataFrame->qEkfQuat.dZ * vMagnetRef.dY; //checked

	dH.matrix[5][0] = -pDataFrame->qEkfQuat.dX * vMagnetRef.dY + pDataFrame->qEkfQuat.dY * vMagnetRef.dX; //checked
	dH.matrix[5][1] = -pDataFrame->qEkfQuat.dW * vMagnetRef.dY - 2.0f * pDataFrame->qEkfQuat.dX * vMagnetRef.dZ + pDataFrame->qEkfQuat.dZ * vMagnetRef.dX; //checked
	dH.matrix[5][2] =  pDataFrame->qEkfQuat.dW * vMagnetRef.dX - 2.0f * pDataFrame->qEkfQuat.dY * vMagnetRef.dZ + pDataFrame->qEkfQuat.dZ * vMagnetRef.dY; //checked
	dH.matrix[5][3] =  pDataFrame->qEkfQuat.dX * vMagnetRef.dX + pDataFrame->qEkfQuat.dY * vMagnetRef.dY; //checked

	ekf_MultiplyScalarToMatrix(&dH, 2.0, &dH);

	// 5. Measurement Prediction Covariance (S   = HPH' + R)

	MATRIX_f dS;
	MATRIX_f dR;
	MATRIX_f dP;
	MATRIX_f dH_Transpose;
	MATRIX_f dTempMatrix_6x4;

	ekf_AllocateMatrix(&dH_Transpose, 4, 6);
	ekf_AllocateMatrix(&dS, 6, 6);
	ekf_AllocateMatrix(&dR, 6, 6);
	ekf_AllocateMatrix(&dP, 4, 4);
	ekf_AllocateMatrix(&dTempMatrix_6x4, 6, 4);

	//dP.matrix[0][0] = pDataFrame->vCovarianceP.dROW1.dW;
	//dP.matrix[0][1] = pDataFrame->vCovarianceP.dROW1.dX;
	//dP.matrix[0][2] = pDataFrame->vCovarianceP.dROW1.dY;
	//dP.matrix[0][3] = pDataFrame->vCovarianceP.dROW1.dZ;

	dP.matrix[0][0] = pDataFrame->vCovarianceP.dROW1.dW;
	dP.matrix[0][1] = pDataFrame->vCovarianceP.dROW1.dX;
	dP.matrix[0][2] = pDataFrame->vCovarianceP.dROW1.dY;
	dP.matrix[0][3] = pDataFrame->vCovarianceP.dROW1.dZ;

	dP.matrix[1][0] = pDataFrame->vCovarianceP.dROW2.dW;
	dP.matrix[1][1] = pDataFrame->vCovarianceP.dROW2.dX;
	dP.matrix[1][2] = pDataFrame->vCovarianceP.dROW2.dY;
	dP.matrix[1][3] = pDataFrame->vCovarianceP.dROW2.dZ;

	dP.matrix[2][0] = pDataFrame->vCovarianceP.dROW3.dW;
	dP.matrix[2][1] = pDataFrame->vCovarianceP.dROW3.dX;
	dP.matrix[2][2] = pDataFrame->vCovarianceP.dROW3.dY;
	dP.matrix[2][3] = pDataFrame->vCovarianceP.dROW3.dZ;

	dP.matrix[3][0] = pDataFrame->vCovarianceP.dROW4.dW;
	dP.matrix[3][1] = pDataFrame->vCovarianceP.dROW4.dX;
	dP.matrix[3][2] = pDataFrame->vCovarianceP.dROW4.dY;
	dP.matrix[3][3] = pDataFrame->vCovarianceP.dROW4.dZ;

	dR.matrix[0][0] = pDataFrame->dAccelSigma2;
	dR.matrix[1][1] = pDataFrame->dAccelSigma2;
	dR.matrix[2][2] = pDataFrame->dAccelSigma2;
	dR.matrix[3][3] = pDataFrame->dMagnetSigma2;
	dR.matrix[4][4] = pDataFrame->dMagnetSigma2;
	dR.matrix[5][5] = pDataFrame->dMagnetSigma2;

	ekf_MultiplyTwoMatrix(&dH, &dP, &dTempMatrix_6x4);                          // HP
	ekf_TransposeMatrix(&dH, &dH_Transpose);                                   //H'
	ekf_MultiplyTwoMatrix(&dTempMatrix_6x4, &dH_Transpose, &dS);                // HPH'
	ekf_AddTwoMatrix(&dS, &dR, &dS);                                            // HPH'+R

	//6. Calculate Kalman gain (K = P.H'.inv(S))

	MATRIX_f dK;
	ekf_AllocateMatrix(&dK, 4, 6);

	ekf_InverseMatrix(&dS, &dS, 6);
	ekf_MultiplyTwoMatrix(&dP, &dH_Transpose, &dH_Transpose);               //P.H'
	ekf_MultiplyTwoMatrix(&dH_Transpose, &dS, &dK);                 // K = P.H'.inv(S)


	// 7. Update the quaternion
	pDataFrame->qEkfQuat.dW += dK.matrix[0][0] * vInnovation1.dX + dK.matrix[0][1] * vInnovation1.dY + dK.matrix[0][2] * vInnovation1.dZ +
	                           dK.matrix[0][3] * vInnovation2.dX + dK.matrix[0][4] * vInnovation2.dY + dK.matrix[0][5] * vInnovation2.dZ;

	pDataFrame->qEkfQuat.dX += dK.matrix[1][0] * vInnovation1.dX + dK.matrix[1][1] * vInnovation1.dY + dK.matrix[1][2] * vInnovation1.dZ +
	                           dK.matrix[1][3] * vInnovation2.dX + dK.matrix[1][4] * vInnovation2.dY + dK.matrix[1][5] * vInnovation2.dZ;

	pDataFrame->qEkfQuat.dY += dK.matrix[2][0] * vInnovation1.dX + dK.matrix[2][1] * vInnovation1.dY + dK.matrix[2][2] * vInnovation1.dZ +
	                           dK.matrix[2][3] * vInnovation2.dX + dK.matrix[2][4] * vInnovation2.dY + dK.matrix[2][5] * vInnovation2.dZ;

	pDataFrame->qEkfQuat.dZ += dK.matrix[3][0] * vInnovation1.dX + dK.matrix[3][1] * vInnovation1.dY + dK.matrix[3][2] * vInnovation1.dZ +
	                           dK.matrix[3][3] * vInnovation2.dX + dK.matrix[3][4] * vInnovation2.dY + dK.matrix[3][5] * vInnovation2.dZ;

	mathTransform_NormalizeQuat(&pDataFrame->qEkfQuat, &pDataFrame->qEkfQuat);


	// 8. Update the covariance matrix   P = (I-KH)P
	MATRIX_f dTempMatrix_4x4, dI;
	ekf_AllocateMatrix(&dTempMatrix_4x4, 4, 4);
	ekf_AllocateMatrix(&dI, 4, 4);

	ekf_MultiplyTwoMatrix(&dK, &dH, &dTempMatrix_4x4);
	ekf_MultiplyScalarToMatrix(&dTempMatrix_4x4, -1.0f, &dTempMatrix_4x4);     //-KH

	dI.matrix[0][0] = 1.0f;
	dI.matrix[1][1] = 1.0f;
	dI.matrix[2][2] = 1.0f;
	dI.matrix[3][3] = 1.0f;

	ekf_AddTwoMatrix(&dI, &dTempMatrix_4x4, &dTempMatrix_4x4);     //I-KH
	ekf_MultiplyTwoMatrix(&dTempMatrix_4x4, &dP, &dP);

	//ekf_TransposeMatrix(&dP, &dI);
	//ekf_AddTwoMatrix(&dP, &dI, &dP);
	//ekf_MultiplyScalarToMatrix(&dP, 1.5f, &dP);

	pDataFrame->vCovarianceP.dROW1.dW = dP.matrix[0][0];
	pDataFrame->vCovarianceP.dROW1.dX = dP.matrix[0][1];
	pDataFrame->vCovarianceP.dROW1.dY = dP.matrix[0][2];
	pDataFrame->vCovarianceP.dROW1.dZ = dP.matrix[0][3];

	pDataFrame->vCovarianceP.dROW2.dW = dP.matrix[1][0];
	pDataFrame->vCovarianceP.dROW2.dX = dP.matrix[1][1];
	pDataFrame->vCovarianceP.dROW2.dY = dP.matrix[1][2];
	pDataFrame->vCovarianceP.dROW2.dZ = dP.matrix[1][3];

	pDataFrame->vCovarianceP.dROW3.dW = dP.matrix[2][0];
	pDataFrame->vCovarianceP.dROW3.dX = dP.matrix[2][1];
	pDataFrame->vCovarianceP.dROW3.dY = dP.matrix[2][2];
	pDataFrame->vCovarianceP.dROW3.dZ = dP.matrix[2][3];

	pDataFrame->vCovarianceP.dROW4.dW = dP.matrix[3][0];
	pDataFrame->vCovarianceP.dROW4.dX = dP.matrix[3][1];
	pDataFrame->vCovarianceP.dROW4.dY = dP.matrix[3][2];
	pDataFrame->vCovarianceP.dROW4.dZ = dP.matrix[3][3];
}

void ekf_Correction_No_Mag(
    EKF_DATA_FRAME * pDataFrame,
    const VECTOR_3D * pvDeviceAccel)
{

	// 1. reference observation form magnetometer and accelerometer
	VECTOR_3D vAccelRef = {0.0, 0.0, 1.0};
	mathTransform_NormalizeVector(&vAccelRef, &vAccelRef);

	// 2. Convert the reference observations to body frame and generate y vector
	MATRIX_3X3 vDCM;
	VECTOR_3D     vObservation;

	mathTransform_ClearMatrix(&vDCM);
	mathTransform_ClearVector(&vObservation);

	mathTransform_QuatToDirectCosineMatrix(&vDCM, &pDataFrame->qEkfQuat);

	vObservation.dX = vAccelRef.dX * vDCM.dROW1.dX + vAccelRef.dY * vDCM.dROW2.dX + vAccelRef.dZ * vDCM.dROW3.dX;
	vObservation.dY = vAccelRef.dX * vDCM.dROW1.dY + vAccelRef.dY * vDCM.dROW2.dY + vAccelRef.dZ * vDCM.dROW3.dY;
	vObservation.dZ = vAccelRef.dX * vDCM.dROW1.dZ + vAccelRef.dY * vDCM.dROW2.dZ + vAccelRef.dZ * vDCM.dROW3.dZ;

	mathTransform_NormalizeVector(&vObservation, &vObservation);

	// 3. Obtain innovation vector and measurement vector (v = z - y)
	VECTOR_3D     vMeasurement = {0, 0, 0};
	VECTOR_3D     vInnovation =  {0, 0, 0};

	VECTOR_3D     pvDeviceAccel_norm = {0, 0, 0};

	mathTransform_NormalizeVector(&pvDeviceAccel_norm, pvDeviceAccel);

	vMeasurement.dX = pvDeviceAccel_norm.dX;
	vMeasurement.dY = pvDeviceAccel_norm.dY;
	vMeasurement.dZ = pvDeviceAccel_norm.dZ;

	vInnovation.dX = vMeasurement.dX - vObservation.dX;
	vInnovation.dY = vMeasurement.dY - vObservation.dY;
	vInnovation.dZ = vMeasurement.dZ - vObservation.dZ;

	// 4. Obtain the measurement matrix H
	MATRIX_f dH;
	ekf_AllocateMatrix(&dH, 3, 4);

	dH.matrix[0][0] = -pDataFrame->qEkfQuat.dY * vAccelRef.dZ + pDataFrame->qEkfQuat.dZ * vAccelRef.dY; //checked
	dH.matrix[0][1] =  pDataFrame->qEkfQuat.dY * vAccelRef.dY + pDataFrame->qEkfQuat.dZ * vAccelRef.dZ; //checked
	dH.matrix[0][2] = -pDataFrame->qEkfQuat.dW * vAccelRef.dZ + pDataFrame->qEkfQuat.dX * vAccelRef.dY - 2.0f * pDataFrame->qEkfQuat.dY * vAccelRef.dX; //checked
	dH.matrix[0][3] =  pDataFrame->qEkfQuat.dW * vAccelRef.dY + pDataFrame->qEkfQuat.dX * vAccelRef.dZ - 2.0f * pDataFrame->qEkfQuat.dZ * vAccelRef.dX; //checked

	dH.matrix[1][0] =  pDataFrame->qEkfQuat.dX * vAccelRef.dZ - pDataFrame->qEkfQuat.dZ * vAccelRef.dX; //checked
	dH.matrix[1][1] =  pDataFrame->qEkfQuat.dW * vAccelRef.dZ - 2.0f * pDataFrame->qEkfQuat.dX * vAccelRef.dY + pDataFrame->qEkfQuat.dY * vAccelRef.dX; //checked
	dH.matrix[1][2] =  pDataFrame->qEkfQuat.dX * vAccelRef.dX + pDataFrame->qEkfQuat.dZ * vAccelRef.dZ; //checked
	dH.matrix[1][3] = -pDataFrame->qEkfQuat.dW * vAccelRef.dX + pDataFrame->qEkfQuat.dY * vAccelRef.dZ - 2.0f * pDataFrame->qEkfQuat.dZ * vAccelRef.dY; //checked

	dH.matrix[2][0] = -pDataFrame->qEkfQuat.dX * vAccelRef.dY + pDataFrame->qEkfQuat.dY * vAccelRef.dX; //checked
	dH.matrix[2][1] = -pDataFrame->qEkfQuat.dW * vAccelRef.dY - 2.0f * pDataFrame->qEkfQuat.dX * vAccelRef.dZ + pDataFrame->qEkfQuat.dZ * vAccelRef.dX; //checked
	dH.matrix[2][2] =  pDataFrame->qEkfQuat.dW * vAccelRef.dX - 2.0f * pDataFrame->qEkfQuat.dY * vAccelRef.dZ + pDataFrame->qEkfQuat.dZ * vAccelRef.dY; //checked
	dH.matrix[2][3] =  pDataFrame->qEkfQuat.dX * vAccelRef.dX + pDataFrame->qEkfQuat.dY * vAccelRef.dY; //checked

	ekf_MultiplyScalarToMatrix(&dH, 2.0, &dH);

	// 5. Measurement Prediction Covariance (S   = HPH' + R)
	MATRIX_f dS;
	MATRIX_f dR;
	MATRIX_f dP;
	MATRIX_f dH_Transpose;
	MATRIX_f dTempMatrix_3x4;

	ekf_AllocateMatrix(&dH_Transpose, 4, 3);
	ekf_AllocateMatrix(&dS, 3, 3);
	ekf_AllocateMatrix(&dR, 3, 3);
	ekf_AllocateMatrix(&dP, 4, 4);
	ekf_AllocateMatrix(&dTempMatrix_3x4, 3, 4);

	dP.matrix[0][0] = pDataFrame->vCovarianceP.dROW1.dW;
	dP.matrix[0][1] = pDataFrame->vCovarianceP.dROW1.dX;
	dP.matrix[0][2] = pDataFrame->vCovarianceP.dROW1.dY;
	dP.matrix[0][3] = pDataFrame->vCovarianceP.dROW1.dZ;

	dP.matrix[1][0] = pDataFrame->vCovarianceP.dROW2.dW;
	dP.matrix[1][1] = pDataFrame->vCovarianceP.dROW2.dX;
	dP.matrix[1][2] = pDataFrame->vCovarianceP.dROW2.dY;
	dP.matrix[1][3] = pDataFrame->vCovarianceP.dROW2.dZ;

	dP.matrix[2][0] = pDataFrame->vCovarianceP.dROW3.dW;
	dP.matrix[2][1] = pDataFrame->vCovarianceP.dROW3.dX;
	dP.matrix[2][2] = pDataFrame->vCovarianceP.dROW3.dY;
	dP.matrix[2][3] = pDataFrame->vCovarianceP.dROW3.dZ;

	dP.matrix[3][0] = pDataFrame->vCovarianceP.dROW4.dW;
	dP.matrix[3][1] = pDataFrame->vCovarianceP.dROW4.dX;
	dP.matrix[3][2] = pDataFrame->vCovarianceP.dROW4.dY;
	dP.matrix[3][3] = pDataFrame->vCovarianceP.dROW4.dZ;

	dR.matrix[0][0] = pDataFrame->dAccelSigma2;
	dR.matrix[1][1] = pDataFrame->dAccelSigma2;
	dR.matrix[2][2] = pDataFrame->dAccelSigma2;

	ekf_MultiplyTwoMatrix(&dH, &dP, &dTempMatrix_3x4);                          // HP
	ekf_TransposeMatrix(&dH, &dH_Transpose);                                   //H'
	ekf_MultiplyTwoMatrix(&dTempMatrix_3x4, &dH_Transpose, &dS);                // HPH'
	ekf_AddTwoMatrix(&dS, &dR, &dS);                                            // HPH'+R

	//6. Calculate Kalman gain (K = P.H'.inv(S))
	MATRIX_f dK;
	ekf_AllocateMatrix(&dK, 4, 3);

	ekf_InverseMatrix(&dS, &dS, 3);
	ekf_MultiplyTwoMatrix(&dP, &dH_Transpose, &dH_Transpose);               //P.H'
	ekf_MultiplyTwoMatrix(&dH_Transpose, &dS, &dK);                 // K = P.H'.inv(S)


	// 7. Update the quaternion
	pDataFrame->qEkfQuat.dW += dK.matrix[0][0] * vInnovation.dX + dK.matrix[0][1] * vInnovation.dY + dK.matrix[0][2] * vInnovation.dZ;
	pDataFrame->qEkfQuat.dX += dK.matrix[1][0] * vInnovation.dX + dK.matrix[1][1] * vInnovation.dY + dK.matrix[1][2] * vInnovation.dZ;
	pDataFrame->qEkfQuat.dY += dK.matrix[2][0] * vInnovation.dX + dK.matrix[2][1] * vInnovation.dY + dK.matrix[2][2] * vInnovation.dZ;
	pDataFrame->qEkfQuat.dZ += dK.matrix[3][0] * vInnovation.dX + dK.matrix[3][1] * vInnovation.dY + dK.matrix[3][2] * vInnovation.dZ;

	mathTransform_NormalizeQuat(&pDataFrame->qEkfQuat, &pDataFrame->qEkfQuat);


	// 8. Update the covariance matrix — Joseph form: P = (I-KH)P(I-KH)' + KRK'
	MATRIX_f dTempMatrix_4x4, dI, dIKH, dIKH_T, dKR, dKR_KT;
	ekf_AllocateMatrix(&dTempMatrix_4x4, 4, 4);
	ekf_AllocateMatrix(&dI,     4, 4);
	ekf_AllocateMatrix(&dIKH,   4, 4);
	ekf_AllocateMatrix(&dIKH_T, 4, 4);
	ekf_AllocateMatrix(&dKR,    4, 3);
	ekf_AllocateMatrix(&dKR_KT, 4, 4);

	// Build I - KH
	ekf_MultiplyTwoMatrix(&dK, &dH, &dTempMatrix_4x4);              // KH
	ekf_MultiplyScalarToMatrix(&dTempMatrix_4x4, -1.0f, &dTempMatrix_4x4); // -KH
	dI.matrix[0][0] = 1.0f;
	dI.matrix[1][1] = 1.0f;
	dI.matrix[2][2] = 1.0f;
	dI.matrix[3][3] = 1.0f;
	ekf_AddTwoMatrix(&dI, &dTempMatrix_4x4, &dIKH);                 // I - KH

	// (I-KH) P (I-KH)'
	ekf_TransposeMatrix(&dIKH, &dIKH_T);
	ekf_MultiplyTwoMatrix(&dIKH, &dP, &dTempMatrix_4x4);            // (I-KH)P
	ekf_MultiplyTwoMatrix(&dTempMatrix_4x4, &dIKH_T, &dP);          // (I-KH)P(I-KH)'

	// K R K'
	ekf_MultiplyTwoMatrix(&dK, &dR, &dKR);                          // KR
	MATRIX_f dK_T;
	ekf_AllocateMatrix(&dK_T, 3, 4);
	ekf_TransposeMatrix(&dK, &dK_T);
	ekf_MultiplyTwoMatrix(&dKR, &dK_T, &dKR_KT);                    // KRK'

	ekf_AddTwoMatrix(&dP, &dKR_KT, &dP);                            // (I-KH)P(I-KH)' + KRK'

	pDataFrame->vCovarianceP.dROW1.dW = dP.matrix[0][0];
	pDataFrame->vCovarianceP.dROW1.dX = dP.matrix[0][1];
	pDataFrame->vCovarianceP.dROW1.dY = dP.matrix[0][2];
	pDataFrame->vCovarianceP.dROW1.dZ = dP.matrix[0][3];

	pDataFrame->vCovarianceP.dROW2.dW = dP.matrix[1][0];
	pDataFrame->vCovarianceP.dROW2.dX = dP.matrix[1][1];
	pDataFrame->vCovarianceP.dROW2.dY = dP.matrix[1][2];
	pDataFrame->vCovarianceP.dROW2.dZ = dP.matrix[1][3];

	pDataFrame->vCovarianceP.dROW3.dW = dP.matrix[2][0];
	pDataFrame->vCovarianceP.dROW3.dX = dP.matrix[2][1];
	pDataFrame->vCovarianceP.dROW3.dY = dP.matrix[2][2];
	pDataFrame->vCovarianceP.dROW3.dZ = dP.matrix[2][3];

	pDataFrame->vCovarianceP.dROW4.dW = dP.matrix[3][0];
	pDataFrame->vCovarianceP.dROW4.dX = dP.matrix[3][1];
	pDataFrame->vCovarianceP.dROW4.dY = dP.matrix[3][2];
	pDataFrame->vCovarianceP.dROW4.dZ = dP.matrix[3][3];
}

void ekf_Update(
    EKF_DATA_FRAME * ptEkf,
    const VECTOR_3D * vAccel,
    const VECTOR_3D * vGyro,
    const VECTOR_3D * vMag,
    const float dSampleTime)
{
	if (!ptEkf->fInitialized)
	{
		ekf_Init(ptEkf, NULL, vAccel, 1e-1f * 1e-1f, 1e-1f * 1e-1f, 1.0f * 1.0f, -10.14f * DEGREE_TO_RAD);
	}

	ekf_Prediction(ptEkf, vGyro, dSampleTime);

	if (mathTransform_MagnitudeVector(vMag) != 0)
	{
		ekf_Correction_With_Mag(ptEkf, vAccel, vMag);
	}
	else
	{
		ekf_Correction_No_Mag(ptEkf, vAccel);
	}

}

QUATERNION ekf_GetQuat(EKF_DATA_FRAME * pDataFrame)
{
	return pDataFrame->qEkfQuat;
}


// ********************************************************
// **************  some math functions  *******************
// ********************************************************

void ekf_MultiplyTwoMatrix(MATRIX_f * pInputA, MATRIX_f * pInputB, MATRIX_f * pOutput)
{
	int i = 0;
	int j = 0;
	int k = 0;
	float sum = 0.0;

	MATRIX_f mTemp;
	ekf_AllocateMatrix(&mTemp, pInputA->uiNumRows, pInputB->uiNumCols);

	for (i = 0; i < pInputA->uiNumRows; i++)
	{
		for (j = 0; j < pInputB->uiNumCols; j++)
		{
			for (k = 0; k < pInputA->uiNumCols; k++)
			{
				sum += pInputA->matrix[i][k] * pInputB->matrix[k][j];
			}
			mTemp.matrix[i][j] = sum;
			sum = 0;
		}
	}

	ekf_CopyMatrix(&mTemp, pOutput);
}

void ekf_AddTwoMatrix(MATRIX_f * pInputA, MATRIX_f * pInputB, MATRIX_f * pOutput)
{
	int i = 0;
	int j = 0;

	for (i = 0; i < pInputA->uiNumRows; i++)
	{
		for (j = 0; j < pInputA->uiNumCols; j++)
		{
			pOutput->matrix[i][j] = pInputB->matrix[i][j] + pInputA->matrix[i][j];
		}
	}

}

void ekf_AllocateMatrix(MATRIX_f * pInput, uint16_t uiNumberOfRows, uint16_t uiNumberOfColumns)
{
	memset(pInput, 0, sizeof(MATRIX_f));

	pInput->uiNumRows = uiNumberOfRows;
	pInput->uiNumCols = uiNumberOfColumns;
}

void ekf_AllocateMatrix_d(MATRIX_d * pInput, uint16_t uiNumberOfRows, uint16_t uiNumberOfColumns)
{
	memset(pInput, 0, sizeof(MATRIX_d));

	if (pInput != NULL)
	{
		pInput->uiNumRows = uiNumberOfRows;
		pInput->uiNumCols = uiNumberOfColumns;
	}
}

void ekf_eyeMatrix(MATRIX_f * pInput)
{
	int i = 0;
	int j = 0;

	if (pInput != NULL)
	{
		for (i = 0; i < pInput->uiNumCols; i++)
		{
			for (j = 0; j < pInput->uiNumCols; j++)
			{
				if (i == j)
				{
					pInput->matrix[i][j] = 1;
				}
				else
				{
					pInput->matrix[i][j] = 0;
				}
			}
		}
	}

}

void ekf_MultiplyScalarToMatrix(MATRIX_f * pInput, float dScalar, MATRIX_f * pOutput)
{
	int i = 0;
	int j = 0;

	for (i = 0; i < pInput->uiNumRows; i++)
	{
		for (j = 0; j < pInput->uiNumCols; j++)
		{
			pOutput->matrix[i][j] = dScalar * (pInput->matrix[i][j]);
		}
	}
}

void ekf_TransposeMatrix(MATRIX_f * pInput, MATRIX_f * pOutput)
{
	int i = 0;
	int j = 0;

	MATRIX_f mTemp;
	ekf_AllocateMatrix(&mTemp, pInput->uiNumCols, pInput->uiNumRows);

	for (i = 0; i < pInput->uiNumRows; i++)
	{
		for (j = 0; j < pInput->uiNumCols; j++)
		{
			mTemp.matrix[j][i] = pInput->matrix[i][j];
		}
	}

	ekf_CopyMatrix(&mTemp, pOutput);
}

float ekf_RecursiveMatrixDeterminant(MATRIX_f * pInput, uint16_t uiSize)
{
	float s = 1;
	float det = 0;

	MATRIX_f b;
	ekf_AllocateMatrix(&b, uiSize, uiSize);

	int i = 0;
	int j = 0;
	int m = 0;
	int n = 0;
	int c = 0;
	if (uiSize == 1)
	{
		return (pInput->matrix[0][0]);
	}
	else
	{
		det = 0;
		for (c = 0; c < uiSize; c++)
		{
			m = 0;
			n = 0;
			for (i = 0; i < uiSize; i++)
			{
				for (j = 0 ; j < uiSize; j++)
				{
					b.matrix[i][j] = 0;
					if (i != 0 && j != c)
					{
						b.matrix[m][n] = pInput->matrix[i][j];
						if (n < (uiSize - 2))
						{
							n++;
						}
						else
						{
							n = 0;
							m++;
						}
					}
				}
			}
			det = det + s * (pInput->matrix[0][c] * ekf_RecursiveMatrixDeterminant(&b, uiSize - 1));
			s = -1 * s;
		}
	}
	return (det);
}

float ekf_MatrixDeterminant(MATRIX_f * pInput, uint16_t uiSize)
{
	int i, j, a;
	float mult;
	float divider;
	float determinate = 1;
	float swap[uiSize];
	int sign = 1;

	MATRIX_f pInput_Copy;
	ekf_AllocateMatrix(&pInput_Copy, uiSize, uiSize);
	ekf_CopyMatrix(pInput, &pInput_Copy);

	for (a = 0; a < uiSize; a++) //for changing the value of the row
	{
		divider = pInput_Copy.matrix[a][a];
		if (divider == 0) //if diagonal element became 0 then swap the above rows.
		{
			if (a == uiSize - 1 && pInput_Copy.matrix[uiSize - 1][uiSize - 1] == 0)
			{
				determinate = 0;
				break;
			}
			for (j = 0; j < uiSize; j++)
			{
				swap[j] = pInput_Copy.matrix[a + 1][j];
				pInput_Copy.matrix[a + 1][j] = pInput_Copy.matrix[a][j];
				pInput_Copy.matrix[a][j] = swap[j];
			}
			sign = -1 * sign;
		}
		divider = pInput_Copy.matrix[a][a];
		for (j = a; j < uiSize; j++) //dividing to make 1 (R1---->R1/a)
		{
			pInput_Copy.matrix[a][j] = pInput_Copy.matrix[a][j] / divider;
		}
		determinate = determinate * divider;

		for (i = a + 1; i < uiSize; i++) //for changing the value of multiplying factor(a) and making 0(R2--->R2-aR1)
		{
			mult = pInput_Copy.matrix[i][a];
			for (j = 0; j < uiSize; j++) //effect of (R2--->R2-aR1)
			{
				pInput_Copy.matrix[i][j] = pInput_Copy.matrix[i][j] - mult * pInput_Copy.matrix[a][j];
			}
		}
	}

	return determinate;
}

void ekf_LUDecomposition(MATRIX_f * pInput, MATRIX_f * pU, MATRIX_f * pL, uint16_t uiSize)
{

	float sum;
	int i, j, t;

	//Puts 1 in the diagonal of L//
	for (i = 0; i < uiSize; i++)
	{
		pL->matrix[i][i] = 1;
	}

	//First row of U and First column of L//
	for (j = 0; j < uiSize; j++)
	{
		pU->matrix[0][j] = pInput->matrix[0][j];
		pL->matrix[j][0] = pInput->matrix[j][0] / pInput->matrix[0][0];
	}


	for (i = 1; i < (uiSize - 1); i++)
	{
		sum = 0;
		for (t = 0; t <= (i - 1); t++)
		{
			sum = sum + (pL->matrix[i][t] * pU->matrix[t][i]);
		}

		pU->matrix[i][i] = pInput->matrix[i][i] - sum; //The Diagonal of U

		for (j = i; j < uiSize; j++)
		{
			sum = 0;
			for (t = 0; t <= (i - 1); t++)
			{
				sum = sum + (pL->matrix[i][t] * pU->matrix[t][j]);
			}
			pU->matrix[i][j] = pInput->matrix[i][j] - sum; //ith row of U

			sum = 0;
			for (t = 0; t <= (i - 1); t++)
			{
				sum = sum + (pL->matrix[j][t] * pU->matrix[t][i]);
			}
			pL->matrix[j][i] = (pInput->matrix[j][i] - sum) / pU->matrix[i][i]; // ith column of L
		}
	}

	sum = 0;
	for (t = 0; t < uiSize - 1; t++)
	{
		sum = sum + (pL->matrix[uiSize - 1][t] * pU->matrix[t][uiSize - 1]);
	}

	pU->matrix[uiSize - 1][uiSize - 1] = pInput->matrix[uiSize - 1][uiSize - 1] - sum; //n,n element of U

}

float ekf_LUDeterminant(MATRIX_f * pInput, uint16_t uiSize)
{

	float result = 1;
	int i;
	MATRIX_f L, U;
	ekf_AllocateMatrix(&L, uiSize, uiSize);
	ekf_AllocateMatrix(&U, uiSize, uiSize);

	ekf_LUDecomposition(pInput, &U, &L, uiSize);

	for (i = 0; i < uiSize; i++)
	{
		result *= L.matrix[i][i];
		result *= U.matrix[i][i];
	}

	return result;
}


void ekf_InverseMatrix(MATRIX_f * pInput, MATRIX_f * pOutput, uint16_t uiSize)
{
	MATRIX_f b, fac;
	int p, q, m, n, i, j;

	ekf_AllocateMatrix(&b, uiSize, uiSize);
	ekf_AllocateMatrix(&fac, uiSize, uiSize);

	for (q = 0; q < uiSize; q++)
	{
		for (p = 0; p < uiSize; p++)
		{
			m = 0;
			n = 0;
			for (i = 0; i < uiSize; i++)
			{
				for (j = 0; j < uiSize; j++)
				{
					if (i != q && j != p)
					{
						b.matrix[m][n] = pInput->matrix[i][j];
						if (n < (uiSize - 2))
						{
							n++;
						}
						else
						{
							n = 0;
							m++;
						}
					}
				}
			}

			fac.matrix[q][p] = powf(-1, q + p) * ekf_RecursiveMatrixDeterminant(&b, uiSize - 1);
		}
	}

	// do transpose and inverse
	ekf_TransposeMatrix(&fac, &b);
	float dDet =  ekf_LUDeterminant(pInput, uiSize);

	if (dDet)
	{
		ekf_MultiplyScalarToMatrix(&b, 1.0f / dDet, pOutput);
	}
}

void ekf_CopyMatrix(MATRIX_f * pFromA, MATRIX_f * pToB)
{
	int i = 0;
	int j = 0;

	for (i = 0; i < pFromA->uiNumRows; i++)
	{
		for (j = 0; j < pFromA->uiNumCols; j++)
		{
			pToB->matrix[i][j] = pFromA->matrix[i][j];
		}
	}

}
