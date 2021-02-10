// Rishav (2020/11/15)

#ifndef _LANCZOS_VECTOR_H_
#define _LANCZOS_VECTOR_H_

#include <cstdint>
#include <cmath>

#if DISPLAY_MATH == 1
#include <iostream>
#include <iomanip>
#endif

template <size_t N>
class Vector
{
public:
    float V[N];

    Vector(const float (&array)[N]);
    Vector(const Vector &m);
    Vector();
    ~Vector() {}

    void operator=(const float (&array)[N]);
    void operator=(const Vector &m);

    template <size_t U>
    friend float sum(const Vector<U> &v);

    template <size_t U>
    friend float var(const Vector<U> &v);

    template <size_t U>
    friend Vector<U> unit(const Vector<U> &v);

    // Vector-Vector operations
    template <size_t U>
    friend float dot(const Vector<U> &a, const Vector<U> &b);

    template <size_t U>
    friend Vector<U> cross(const Vector<U> &a, const Vector<U> &b);

    Vector operator+(const Vector &m) const;
    Vector operator-(const Vector &m) const;

    // Scalar-Vector operations
    Vector operator*(float s);
    Vector operator+(float s);
    Vector operator-(float s);
    Vector operator/(float s);

    // Vector-Scalar operations
    template <size_t U>
    friend Vector<U> operator*(float s, const Vector<U> &v);
    template <size_t U>
    friend Vector<U> operator+(float s, const Vector<U> &v);
    template <size_t U>
    friend Vector<U> operator-(float s, const Vector<U> &v);

    // Vector elements interface
    void set(uint8_t n, float value);
    float operator()(size_t n) const;

#if DISPLAY_MATH == 1
    template <size_t U>
    friend void operator<<(std::ostream &oc, const Vector<U> &v);
#endif
};

template <size_t N>
Vector<N>::Vector(const float (&array)[N])
{
    for (uint8_t n = 0; n < N; n++)
    {
        V[n] = array[n];
    }
}

template <size_t N>
Vector<N>::Vector(const Vector<N> &v)
{
    for (uint8_t n = 0; n < N; n++)
    {
        V[n] = v.V[n];
    }
}

template <size_t N>
Vector<N>::Vector()
{
    for (uint8_t n = 0; n < N; n++)
    {
        V[n] = 0.0;
    }
}

template <size_t N>
void Vector<N>::operator=(const float (&array)[N])
{
    for (uint8_t n = 0; n < N; n++)
    {
        V[n] = array[n];
    }
}

template <size_t N>
void Vector<N>::operator=(const Vector<N> &v)
{
    for (uint8_t n = 0; n < N; n++)
    {
        V[n] = v.V[n];
    }
}

template <size_t N>
float sum(const Vector<N> &v)
{
    float sum = 0.0;
    for (uint8_t n = 0; n < N; n++)
    {
        sum += v.V[n];
    }
    return sum;
}

template <size_t N>
float var(const Vector<N> &v)
{
    float mean = 0.0;
    for (uint8_t n = 0; n < N; n++)
    {
        mean += v.V[n];
    }
    mean = mean / float(N);

    float variance = 0.0;
    for (uint8_t n = 0; n < N; n++)
    {
        variance += (v.V[n] - mean) * (v.V[n] - mean);
    }
    variance = variance / (float(N) - 1);

    return variance;
}

template <size_t N>
float norm(const Vector<N> &v)
{
    float norm = 0.0;
    for (uint8_t n = 0; n < N; n++)
    {
        norm += v.V[n] * v.V[n];
    }
    return sqrt(norm);
}

template <size_t N>
Vector<N> unit(const Vector<N> &v)
{
    Vector<N> unit;
    float norm = 0.0;

    for (uint8_t n = 0; n < N; n++)
    {
        norm += v.V[n] * v.V[n];
    }
    norm = sqrt(norm);

    for (uint8_t n = 0; n < N; n++)
    {
        unit.V[n] = v.V[n] / norm;
    }
    return unit;
}

template <size_t N>
float dot(const Vector<N> &a, const Vector<N> &b)
{
    float output = 0;
    for (uint8_t n = 0; n < N; n++)
    {
        output += a.V[n] * b.V[n];
    }
    return output;
}

template <size_t N = 3>
Vector<N> cross(const Vector<N> &a, const Vector<N> &b)
{
    float i = a.V[1] * b.V[2] - a.V[2] * b.V[1];
    float j = a.V[0] * b.V[2] - a.V[2] * b.V[0];
    float k = a.V[0] * b.V[1] - a.V[1] * b.V[0];
    float out[3] = {i, j, k};

    Vector<3> output(out);
    return output;
}

template <size_t N>
Vector<N> Vector<N>::operator+(const Vector<N> &v) const
{
    Vector<N> output;
    for (uint8_t n = 0; n < N; n++)
    {
        output.V[n] = V[n] + v.V[n];
    }
    return output;
}

template <size_t N>
Vector<N> Vector<N>::operator-(const Vector<N> &v) const
{
    Vector<N> output;
    for (uint8_t n = 0; n < N; n++)
    {
        output.V[n] = V[n] - v.V[n];
    }
    return output;
}

template <size_t N>
Vector<N> Vector<N>::operator*(float s)
{
    Vector<N> output;
    for (uint8_t n = 0; n < N; n++)
    {
        output.V[n] = V[n] * s;
    }
    return output;
}

template <size_t N>
Vector<N> Vector<N>::operator+(float s)
{
    Vector<N> output;
    for (uint8_t n = 0; n < N; n++)
    {
        output.V[n] = V[n] + s;
    }
    return output;
}

template <size_t N>
Vector<N> Vector<N>::operator-(float s)
{
    Vector<N> output;
    for (uint8_t n = 0; n < N; n++)
    {
        output.V[n] = V[n] - s;
    }
    return output;
}

template <size_t N>
Vector<N> Vector<N>::operator/(float s)
{
    Vector<N> output;
    for (uint8_t n = 0; n < N; n++)
    {
        output.V[n] = V[n] / s;
    }
    return output;
}

template <size_t N>
Vector<N> operator*(float s, const Vector<N> &v)
{
    Vector<N> output;
    for (uint8_t n = 0; n < N; n++)
    {
        output.V[n] = v.V[n] * s;
    }
    return output;
}

template <size_t N>
Vector<N> operator+(float s, const Vector<N> &v)
{
    Vector<N> output;
    for (uint8_t n = 0; n < N; n++)
    {
        output.V[n] = v.V[n] + s;
    }
    return output;
}

template <size_t N>
Vector<N> operator-(float s, const Vector<N> &v)
{
    Vector<N> output;
    for (uint8_t n = 0; n < N; n++)
    {
        output.V[n] = v.V[n] - s;
    }
    return output;
}

template <size_t N>
void Vector<N>::set(uint8_t n, float value)
{
    V[n] = value;
}

template <size_t N>
float Vector<N>::operator()(size_t n) const
{
    return V[n];
}

#if DISPLAY_MATH == 1
template <size_t N>
void operator<<(std::ostream &os, const Vector<N> &v)
{
    os << '\n'
       << std::setw(10)
       << std::setfill(' ')
       << std::setprecision(10)
       << v.V[0];
    for (int n = 1; n < N; n++)
    {
        std::cout << std::setw(18)
                  << std::setfill(' ')
                  << std::setprecision(10)
                  << static_cast<double>(v.V[n]);
    }
    os << std::endl;
}
#endif

#endif // vector.h