/**
 * @author Farzan Farhangian
 * @file madgwick.h
 */

#ifndef MADGWICK_H
#define MADGWICK_H

#include "stdbool.h"
#include "orientationDefines.h"

typedef struct
{
	float dMadgwickBeta;
	float dMadgwickZeta;
	QUATERNION qMadgwickQuat;
	VECTOR_3D vGyroBias;
	bool fInitialized;

} MADGWICK_DATA_FRAME;

/**
 * @brief Initialize the madgwick filter data structure
 */
void madgwick_Init(
    MADGWICK_DATA_FRAME * pDataFrame,
    const float * pdHeading,
    const VECTOR_3D * pvDeviceAccel,
    const VECTOR_3D * pvGyroscopeBias,
    const float dBeta,
    const float dZeta);

/**
 * @brief Madgwick 6 Degree of Freedom update
 */
void madgwick_6Dof_Update(
    MADGWICK_DATA_FRAME * pDataFrame,
    const VECTOR_3D * pvDeviceGyro,
    const VECTOR_3D * pvDeviceAccel,
    const float dSampleTime);

/**
 * @brief Madgwick 9 Degree of Freedom update
 */
void madgwick_9Dof_Update(
    MADGWICK_DATA_FRAME * pDataFrame,
    const VECTOR_3D * pvDeviceAccel,
    const VECTOR_3D * pvDeviceGyro,
    const VECTOR_3D * pvDeviceMagnet,
    const float dSampleTime);

/**
 * @brief Get quaternion for Madgwick filter
 */
QUATERNION madgwick_GetQuat(MADGWICK_DATA_FRAME * pDataFrame);

#endif /* MADGWICK_H */
