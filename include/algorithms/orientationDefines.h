/**
 * @author Farzan Farhangian
 * @file orientationDefines.h
 */

#ifndef ORIENTATION_DEFINES_H
#define ORIENTATION_DEFINES_H

#include <stdbool.h>

typedef unsigned char uint8_t;
typedef unsigned short int uint16_t;
typedef unsigned int uint32_t;
typedef unsigned  long int uint64_t;

#define PI 3.14159265358979323846
#define DEGREE_TO_RAD PI/180

typedef struct
{
    float dW;
    float dX;
    float dY;
    float dZ;
} QUATERNION;

typedef struct
{
    float dX;
    float dY;
    float dZ;
} VECTOR_3D;

typedef struct
{
    float matrix[6][6];
    uint16_t uiNumRows;
    uint16_t uiNumCols;
} MATRIX_f;

typedef struct
{
    double matrix[6][6];
    uint16_t uiNumRows;
    uint16_t uiNumCols;
} MATRIX_d;

typedef struct
{
    QUATERNION dROW1;
    QUATERNION dROW2;
    QUATERNION dROW3;
    QUATERNION dROW4;
} MATRIX_4X4;

typedef struct
{
    VECTOR_3D dROW1;
    VECTOR_3D dROW2;
    VECTOR_3D dROW3;
} MATRIX_3X3;


void mathTransform_ClearQuat(QUATERNION * pqInput);

void mathTransform_ClearVector(VECTOR_3D * pvInput);

void mathTransform_ClearMatrix(MATRIX_3X3 * pvInput);

void mathTransform_Conjugate(QUATERNION * pqOutput, const QUATERNION * pqInput);

float mathTransform_SumOfSquaresQuat(const QUATERNION * pqInput);

float mathTransform_SumOfSquaresVector(const VECTOR_3D * pvInput);

float mathTransform_MagnitudeVector(const VECTOR_3D * pvInput);

float mathTransform_MagnitudeQuat(const QUATERNION * pqInput);

void mathTransform_NormalizeQuat(QUATERNION * pqOutput, const QUATERNION * pqInput);

void mathTransform_NormalizeVector(VECTOR_3D * pvOutput, const VECTOR_3D * pvInput);

void mathTransform_MultiplyQuatByScalar(QUATERNION * pqOutput, const QUATERNION * pqInput, const float dScalar);

void mathTransform_Inverse(QUATERNION * pqOutput, const QUATERNION * pqInput);

void mathTransform_MultiplyQuatByQuat(QUATERNION * pqOutput, const QUATERNION * pqInputA, const QUATERNION * pqInputB);

void mathTransform_MultiplyVectorByScalar(VECTOR_3D * pvOutput, const VECTOR_3D * pvInput, const float dScalar);

void mathTransform_CrossVectorWithVector(VECTOR_3D * pvOutput, const VECTOR_3D * pvInputA, const VECTOR_3D * pvInputB);

void mathTransform_GetVectorComponent(VECTOR_3D * pvOutput, const QUATERNION * pqInput);

void mathTransform_AddVectorWithVector(VECTOR_3D * pvOutput, const VECTOR_3D * pvInputA, const VECTOR_3D * pvInputB);

void mathTransform_SubtractVectorFromVector(VECTOR_3D * pvOutput, const VECTOR_3D * pvInputA, const VECTOR_3D * pvInputB);

void mathTransform_Slerp(QUATERNION * pqOutput, const QUATERNION * pqFrom, const QUATERNION * pqTo, const float dGain);

float mathTransform_dotVectorWithVector(const VECTOR_3D * pvInputA, const VECTOR_3D * pvInputB);

void mathTransform_DirectCosineMatrixToQuat(QUATERNION * pqOutput, const MATRIX_3X3 * pvInput);

void mathTransform_QuatToEuler(VECTOR_3D * pvOutput, const QUATERNION * pqInput);

void mathTransform_QuatToDirectCosineMatrix(MATRIX_3X3 * pvOutput, const QUATERNION * pqInput);

void mathTransform_EulerToQuat(QUATERNION * pqOutput, const VECTOR_3D * pvInput);

void mathTransform_VectorToSkewSymmetric(MATRIX_3X3 * pvOutput, const VECTOR_3D * pvInput);

void mathTransform_QuaternionDerivative(QUATERNION * pqQuaternionDot, const QUATERNION * pqQuaternion, const VECTOR_3D * pvDeviceGyro);

void mathTransform_InitializeQuaternion(QUATERNION * pqOutput, const VECTOR_3D * pvDeviceAccel, const float * pdHeading);

#endif /* ORIENTATION_DEFINES_H */
