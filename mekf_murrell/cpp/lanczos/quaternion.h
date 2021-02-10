// Rishav (2020/12/7)

#ifndef _LANCZOS_QUATERNION_H_
#define _LANCZOS_QUATERNION_H_

#include "vector.h"
#include "matrix.h"

class Quaternion : public Vector<4>
{
public:
    Quaternion(const float (&array)[4]);
    Quaternion(const Quaternion &q);
    Quaternion(const Vector<4> &q);
    Quaternion();
    ~Quaternion() {}

    void operator=(const float (&array)[4]);
    void operator=(const Quaternion &q);

    Matrix<3, 3> getDCM();
    Vector<3> getYPR();
};

Quaternion::Quaternion(const float (&array)[4])
{
    for (uint8_t n = 0; n < 4; n++)
    {
        V[n] = array[n];
    }
}

Quaternion::Quaternion(const Quaternion &q)
{
    for (uint8_t n = 0; n < 4; n++)
    {
        V[n] = q.V[n];
    }
}

Quaternion::Quaternion(const Vector<4> &v)
{
    for (uint8_t n = 0; n < 4; n++)
    {
        V[n] = v.V[n];
    }
}

Quaternion::Quaternion()
{
    V[0] = 1;
    for (uint8_t n = 1; n < 4; n++)
    {
        V[n] = 0.0;
    }
}

void Quaternion::operator=(const float (&array)[4])
{
    for (uint8_t n = 0; n < 4; n++)
    {
        V[n] = array[n];
    }
}

void Quaternion::operator=(const Quaternion &q)
{
    for (uint8_t n = 0; n < 4; n++)
    {
        V[n] = q.V[n];
    }
}

Matrix<3, 3> Quaternion::getDCM()
{
    Matrix<3, 3> DCM;
    float q0 = V[3];
    float q1 = V[0];
    float q2 = V[1];
    float q3 = V[2];

    DCM.set(0, 0, q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3);
    DCM.set(0, 1, 2.0 * (q1 * q2 + q0 * q3));
    DCM.set(0, 2, 2.0 * (q1 * q3 - q0 * q2));
    DCM.set(1, 0, 2.0 * (q1 * q2 - q0 * q3));
    DCM.set(1, 1, q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3);
    DCM.set(1, 2, 2.0 * (q2 * q3 + q0 * q1));
    DCM.set(2, 0, 2.0 * (q1 * q3 + q0 * q2));
    DCM.set(2, 1, 2.0 * (q2 * q3 - q0 * q1));
    DCM.set(2, 2, q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3);

    return DCM;
}

Vector<3> Quaternion::getYPR()
{
    float q0 = V[0];
    float q1 = V[1];
    float q2 = V[2];
    float q3 = V[3];

    // radians
    float roll = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2)); // phi
    float pitch = asin(2 * (q0 * q2 - q3 * q1));                              // theta
    float yaw = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3));  // psi

    Vector<3> euler_angles({yaw, pitch, roll});
    return euler_angles;
}

#endif // quaternion.h
