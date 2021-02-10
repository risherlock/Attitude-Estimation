// Rishav (2020/11/11)


#ifndef _LANCZOS_MATRIX_H_
#define _LANCZOS_MATRIX_H_

#include "vector.h"
#include <cstdint>
#include <cmath>

#if DISPLAY_MATH == 1
#include <iostream>
#include <iomanip>
#endif

template <size_t R, size_t C>
class Matrix
{
private:
    uint8_t rows = R;
    uint8_t cols = C;
    
public:
    Matrix(const float (&array)[R][C]);
    Matrix(const Matrix &m);
    Matrix();
    ~Matrix() {}

    float M[R][C];
    
    void operator=(const float (&array)[R][C]);
    void operator=(const Matrix &m);

    Matrix compute_inv() const;
    float compute_det() const;

    template <size_t U, size_t V>
    friend float norm(const Matrix<U, V> &m);

    template <size_t U>
    friend float det(const Matrix<U, U> &m);

    template <size_t U>
    friend Matrix<U, U> inv(const Matrix<U, U> &m);

    template <size_t U, size_t V>
    friend Matrix<V, U> trans(const Matrix<U, V> &m);

    template <size_t U>
    friend float trace(const Matrix<U, U> &m);

    // Matrix-Matrix opertations
    template <size_t T>
    Matrix<R, T> operator*(const Matrix<C, T> &m) const;

    Matrix operator+(const Matrix &m) const;
    Matrix operator-(const Matrix &m) const;
    Matrix operator/(const Matrix &m) const;

    // Matrix-Vector operations
    Vector<R> operator*(const Vector<C> &v) const;
    template <size_t U>
    friend Matrix<U, U> vvT(const Vector<U> &v1, const Vector<U> &v2);

    // Scalar-Matrix operations
    Matrix operator*(float s);
    Matrix operator+(float s);
    Matrix operator-(float s);

    // Matrix-Scalar operations
    template <size_t U, size_t V>
    friend Matrix<U, V> operator*(float s, const Matrix<U, V> &m);
    template <size_t U, size_t V>
    friend Matrix<U, V> operator+(float s, const Matrix<U, V> &m);
    template <size_t U, size_t V>
    friend Matrix<U, V> operator-(float s, const Matrix<U, V> &m);

    // Matrix elements interface
    void set(uint8_t r, uint8_t c, float value);
    float operator()(size_t r, size_t c) const;
    Vector<C> get_row(uint8_t r) const;
    Vector<R> get_col(uint8_t c) const;

#if DISPLAY_MATH == 1
    template <size_t U, size_t V>
    friend void operator<<(std::ostream &oc, const Matrix<U, V> &m);
#endif
};

template <size_t R, size_t C>
Matrix<R, C>::Matrix(const float (&array)[R][C])
{
    for (uint8_t r = 0; r < rows; r++)
    {
        for (uint8_t c = 0; c < cols; c++)
        {
            M[r][c] = array[r][c];
        }
    }
}

template <size_t R, size_t C>
Matrix<R, C>::Matrix(const Matrix<R, C> &m)
{
    for (uint8_t r = 0; r < rows; r++)
    {
        for (uint8_t c = 0; c < cols; c++)
        {
            M[r][c] = m.M[r][c];
        }
    }
}

template <size_t R, size_t C>
Matrix<R, C>::Matrix()
{
    for (uint8_t r = 0; r < rows; r++)
    {
        for (uint8_t c = 0; c < cols; c++)
        {
            M[r][c] = 0.0;
        }
    }
}

template <size_t R, size_t C>
void Matrix<R, C>::operator=(const float (&array)[R][C])
{
    for (uint8_t r = 0; r < rows; r++)
    {
        for (uint8_t c = 0; c < cols; c++)
        {
            M[r][c] = array[r][c];
        }
    }
}

template <size_t R, size_t C>
void Matrix<R, C>::operator=(const Matrix &m)
{
    for (uint8_t r = 0; r < rows; r++)
    {
        for (uint8_t c = 0; c < cols; c++)
        {
            M[r][c] = m.M[r][c];
        }
    }
}

// Frobenius norm
template <size_t R, size_t C>
float norm(const Matrix<R, C> &m)
{
    float norm = 0.0;
    for (int r = 0; r < m.rows; r++)
    {
        for (int c = 0; c < m.cols; c++)
        {

            norm += m.M[r][c] * m.M[r][c];
        }
    }
    return sqrt(norm);
}

template <size_t U>
float trace(const Matrix<U, U> &a)
{
    float trace = 0.0;
    for (int i_iters = 0; i_iters < U; i_iters++)
    {
        trace += a.M[i_iters][i_iters];
    }
    return trace;
}

template <>
float Matrix<1, 1>::compute_det() const
{
    return M[0][0];
}

template <>
float Matrix<2, 2>::compute_det() const
{
    const float a = M[0][0];
    const float b = M[0][1];
    const float c = M[1][0];
    const float d = M[1][1];
    return a * d - b * c;
}

template <>
float Matrix<3, 3>::compute_det() const
{
    float d1 = M[1][1] * M[2][2] - M[2][1] * M[1][2];
    float d2 = M[1][0] * M[2][2] - M[2][0] * M[1][2];
    float d3 = M[1][0] * M[2][1] - M[2][0] * M[1][1];

    return M[0][0] * d1 - M[0][1] * d2 + M[0][2] * d3;
}

template <>
float Matrix<4, 4>::compute_det() const
{
    float m00 = M[0][0];
    float m10 = M[1][0];
    float m20 = M[2][0];
    float m30 = M[3][0];
    float m01 = M[0][1];
    float m11 = M[1][1];
    float m21 = M[2][1];
    float m31 = M[3][1];
    float m02 = M[0][2];
    float m12 = M[1][2];
    float m22 = M[2][2];
    float m32 = M[3][2];
    float m03 = M[0][3];
    float m13 = M[1][3];
    float m23 = M[2][3];
    float m33 = M[3][3];

    return m03 * m12 * m21 * m30 - m02 * m13 * m21 * m30 -
           m03 * m11 * m22 * m30 + m01 * m13 * m22 * m30 +
           m02 * m11 * m23 * m30 - m01 * m12 * m23 * m30 -
           m03 * m12 * m20 * m31 + m02 * m13 * m20 * m31 +
           m03 * m10 * m22 * m31 - m00 * m13 * m22 * m31 -
           m02 * m10 * m23 * m31 + m00 * m12 * m23 * m31 +
           m03 * m11 * m20 * m32 - m01 * m13 * m20 * m32 -
           m03 * m10 * m21 * m32 + m00 * m13 * m21 * m32 +
           m01 * m10 * m23 * m32 - m00 * m11 * m23 * m32 -
           m02 * m11 * m20 * m33 + m01 * m12 * m20 * m33 +
           m02 * m10 * m21 * m33 - m00 * m12 * m21 * m33 -
           m01 * m10 * m22 * m33 + m00 * m11 * m22 * m33;
}

template <size_t R>
float det(const Matrix<R, R> &m)
{
    return m.compute_det();
}

template <>
Matrix<3, 3> Matrix<3, 3>::compute_inv() const
{
    Matrix<3, 3> inv;

    float inv00 = M[1][1] * M[2][2] - M[1][2] * M[2][1];
    float inv01 = M[0][2] * M[2][1] - M[0][1] * M[2][2];
    float inv02 = M[0][1] * M[1][2] - M[0][2] * M[1][1];
    float inv10 = M[1][2] * M[2][0] - M[1][0] * M[2][2];
    float inv11 = M[0][0] * M[2][2] - M[0][2] * M[2][0];
    float inv12 = M[0][2] * M[1][0] - M[0][0] * M[1][2];
    float inv20 = M[1][0] * M[2][1] - M[1][1] * M[2][0];
    float inv21 = M[0][1] * M[2][0] - M[0][0] * M[2][1];
    float inv22 = M[0][0] * M[1][1] - M[0][1] * M[1][0];

    const float determinant =
        M[0][0] * inv00 + M[0][1] * inv10 + M[0][2] * inv20;

    const float detinv = (1.0) / determinant;

    inv.M[0][0] = inv00 *= detinv;
    inv.M[0][1] = inv01 *= detinv;
    inv.M[0][2] = inv02 *= detinv;
    inv.M[1][0] = inv10 *= detinv;
    inv.M[1][1] = inv11 *= detinv;
    inv.M[1][2] = inv12 *= detinv;
    inv.M[2][0] = inv20 *= detinv;
    inv.M[2][1] = inv21 *= detinv;
    inv.M[2][2] = inv22 *= detinv;

    return inv;
}

template <size_t R>
Matrix<R, R> inv(const Matrix<R, R> &m)
{
    return m.compute_inv();
}

template <size_t R, size_t C>
Matrix<C, R> trans(const Matrix<R, C> &m)
{
    Matrix<C, R> transpose;
    for (uint8_t r = 0; r < m.rows; r++)
    {
        for (uint8_t c = 0; c < m.cols; c++)
        {
            transpose.M[c][r] = m.M[r][c];
        }
    }
    return transpose;
}

template <size_t R, size_t C>
void Matrix<R, C>::set(uint8_t r, uint8_t c, float value)
{
    M[r][c] = value;
}

template <size_t R, size_t C>
float Matrix<R, C>::operator()(size_t r, size_t c) const
{
    return M[r][c];
}

template <size_t R, size_t C>
template <size_t T>
Matrix<R, T> Matrix<R, C>::operator*(const Matrix<C, T> &m) const
{
    Matrix<R, T> output;
    for (uint8_t r = 0;  r < rows; r++)
    {
        for (uint8_t c = 0; c < cols; c++)
        {
            float elem = 0.0;    
            for (uint8_t t = 0; t < cols; t++)
            {
               elem += this->M[r][t] * m.M[t][c];
            }
            output.M[r][c] = elem;
        }
    }
    return output;
}

template <size_t R, size_t C>
Matrix<R, C> Matrix<R, C>::operator+(const Matrix<R, C> &m) const
{
    Matrix<R, C> output;
    for (uint8_t r = 0; r < R; r++)
    {
        for (uint8_t c = 0; c < C; c++)
        {
            output.M[r][c] = M[r][c] + m.M[r][c];
        }
    }
    return output;
}

template <size_t R, size_t C>
Matrix<R, C> Matrix<R, C>::operator-(const Matrix<R, C> &m) const
{
    Matrix<R, C> output;
    for (uint8_t r = 0; r < R; r++)
    {
        for (uint8_t c = 0; c < C; c++)
        {
            output.M[r][c] = M[r][c] - m.M[r][c];
        }
    }
    return output;
}

template <size_t R, size_t C>
Matrix<R, C> Matrix<R, C>::operator/(const Matrix<R, C> &m) const
{
    Matrix<R, C> output;
    for (uint8_t r = 0; r < R; r++)
    {
        for (uint8_t c = 0; c < C; c++)
        {
            output.M[r][c] = M[r][c] / m.M[r][c];
        }
    }
    return output;
}

template <size_t R, size_t C>
Vector<R> Matrix<R, C>::operator*(const Vector<C> &v) const
{
    Vector<R> output;
    for (uint8_t r = 0; r < R; r++)
    {
        float elem = 0.0;
        for (size_t c = 0; c < C; c++)
        {
            elem += M[r][c] * v.V[c];
        }
        output.V[r] = elem;
    }
    return output;
}

template <size_t U>
Matrix<U, U> vvT(const Vector<U> &v1, const Vector<U> &v2)
{
    Matrix<U, U> output;
    for (uint8_t r = 0; r < U; r++)
    {
        for (uint8_t c = 0; c < U; c++)
        {
            output.M[r][c] = v1(r) * v2(c);
        }
    }
    return output;
}

template <size_t R, size_t C>
Matrix<R, C> Matrix<R, C>::operator*(float s)
{
    Matrix<R, C> output;
    for (uint8_t r = 0; r < R; r++)
    {
        for (uint8_t c = 0; c < C; c++)
        {
            output.M[r][c] = M[r][c] * s;
        }
    }
    return output;
}

template <size_t R, size_t C>
Matrix<R, C> Matrix<R, C>::operator+(float s)
{
    Matrix<R, C> output;
    for (uint8_t r = 0; r < R; r++)
    {
        for (uint8_t c = 0; c < C; c++)
        {
            output.M[r][c] = M[r][c] + s;
        }
    }
    return output;
}

template <size_t R, size_t C>
Matrix<R, C> Matrix<R, C>::operator-(float s)
{
    Matrix<R, C> output;
    for (uint8_t r = 0; r < R; r++)
    {
        for (uint8_t c = 0; c < C; c++)
        {
            output.M[r][c] = M[r][c] - s;
        }
    }
    return output;
}

template <size_t R, size_t C>
Matrix<R, C> operator*(float s, const Matrix<R, C> &m)
{
    Matrix<R, C> output;
    for (uint8_t r = 0; r < R; r++)
    {
        for (uint8_t c = 0; c < C; c++)
        {
            output.M[r][c] = m.M[r][c] * s;
        }
    }
    return output;
}

template <size_t R, size_t C>
Matrix<R, C> operator+(float s, const Matrix<R, C> &m)
{
    Matrix<R, C> output;
    for (uint8_t r = 0; r < R; r++)
    {
        for (uint8_t c = 0; c < C; c++)
        {
            output.M[r][c] = m.M[r][c] + s;
        }
    }
    return output;
}

template <size_t R, size_t C>
Matrix<R, C> operator-(float s, const Matrix<R, C> &m)
{
    Matrix<R, C> output;
    for (uint8_t r = 0; r < R; r++)
    {
        for (uint8_t c = 0; c < C; c++)
        {
            output.M[r][c] = m.M[r][c] - s;
        }
    }
    return output;
}

template <size_t R, size_t C>
Vector<C> Matrix<R, C>::get_row(uint8_t r) const
{
    Vector<C> row;
    for (uint8_t c = 0; c < C; c++)
    {
        row.V[c] = M[r][c];
    }
    return row;
}

template <size_t R, size_t C>
Vector<R> Matrix<R, C>::get_col(uint8_t c) const
{
    Vector<R> col;
    for (uint8_t r = 0; r < R; r++)
    {
        col.V[r] = M[r][c];
    }
    return col;
}

#if DISPLAY_MATH == 1
template <size_t R, size_t C>
void operator<<(std::ostream &os, const Matrix<R, C> &m)
{
    os << '\n';
    for (uint8_t r = 0; r < R; r++)
    {
        os << std::setw(10)
           << std::setfill(' ')
           << std::setprecision(10)
           << static_cast<double>(m.M[r][0]);
        for (int c = 1; c < C; c++)
        {
            std::cout << std::setw(18)
                      << std::setfill(' ')
                      << std::setprecision(10)
                      << static_cast<double>(m.M[r][c]);
        }
        os << '\n';
    }
    os << std::endl;
}
#endif

#endif // matrix.h
