/**
 * @author Farzan Farhangian
 * @file mahony.h
 */

#ifndef MAHONY_H
#define MAHONY_H

#include "orientationDefines.h"

/**
 * @brief Mahony Struct
 */
typedef struct
{
	float dBeta;
	float dkP;
	float dkI;
	QUATERNION qMahonyQuat;
	bool fInitialized;
} MAHONY_DATA_FRAME;

/**
 * @brief Initialize Mahony Filter Data Structure
 */
void mahony_Init(
    MAHONY_DATA_FRAME * pDataFrame,
    const float * pdHeading,
    const VECTOR_3D * pvDeviceAccel,
    float dBeta,
    float dkP,
    float dkI);

/**
 * @brief Process Mahony Filter with 6 DOF using IMU data
 */
void mahony_6Dof_Update(
    MAHONY_DATA_FRAME * pDataFrame,
    const VECTOR_3D * pvDeviceAccel,
    const VECTOR_3D * pvDeviceGyro,
    float dSampleTime);

/**
 * @brief Process Mahony Filter with 9 DOF using MIMU data
 */
void mahony_9Dof_Update(
    MAHONY_DATA_FRAME * pDataFrame,
    const VECTOR_3D * pvDeviceAccel,
    const VECTOR_3D * pvDeviceGyro,
    const VECTOR_3D * pvDeviceMag,
    float dSampleTime);

/**
 * @brief Extract Mahony Quaternion from Data Frame
 */
QUATERNION mahony_GetQuat(MAHONY_DATA_FRAME * pDataFrame);

#endif /* MAHONY_H */