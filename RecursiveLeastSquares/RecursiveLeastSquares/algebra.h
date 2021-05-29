#pragma once
#include <array>
template <size_t N> using Matrix = std::array<std::array<double, N>, N>;

template <int N>
std::array<double, N> operator*(double c, const std::array<double, N>& A)
{
    std::array<double, N> retv{};
    for (int i = 0; i < N; i++) {
        retv[i] = c * A[i];
    }
    return retv;
}

template <int N>
std::array<double, N> operator*(const std::array<double, N>& A, double c)
{
    std::array<double, N> retv{};
    for (int i = 0; i < N; i++) {
        retv[i] = c * A[i];
    }
    return retv;
}

template <int N>
Matrix<N> operator*(const Matrix<N>& L, const Matrix<N>& R)
{
    Matrix<N> retv{};
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                retv[i][j] += L[i][k] * R[k][j];
            }
        }
    }
    return retv;
}

template <int N>
Matrix<N> operator*(double c, const Matrix<N>& M)
{
    Matrix<N> retv;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            retv[i][j] = M[i][j] * c;
        }
    }
    return retv;
}

template <int N>
Matrix<N> operator*(const Matrix<N>& M, double c)
{
    Matrix<N> retv;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            retv[i][j] = M[i][j] * c;
        }
    }
    return retv;
}

template <int N>
Matrix<N> operator-(const Matrix<N>& L, const Matrix<N>& R)
{
    Matrix<N> retv;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            retv[i][j] = L[i][j] - R[i][j];
        }
    }
    return retv;
}

template <int N>
Matrix<N> operator+(const Matrix<N>& L, const Matrix<N>& R)
{
    Matrix<N> retv;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            retv[i][j] = L[i][j] + R[i][j];
        }
    }
    return retv;
}
