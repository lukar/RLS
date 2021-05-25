#pragma once
#include <iostream>
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
std::array<double, N> operator*( const std::array<double, N>& A,double c)
{
    std::array<double, N> retv{};
    for (int i = 0; i < N; i++) {
        retv[i] = c * A[i];
    }
    return retv;
}

template <int N>
Matrix<N> operator*(const Matrix<N> &L, const Matrix<N> &R)
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
    };
    return retv;
}

template <int N>
Matrix<N> operator*( const Matrix<N>& M, double c)
{
    Matrix<N> retv;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            retv[i][j] = M[i][j] * c;
        }
    };
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
    };
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
    };
    return retv;
}

//

template <int N>
class Rls
{
    Matrix<N> P_old{}, I{};
    std::array<double, N> theta_old{}, gamma{}, x{}, tempv{};
    double err=0.0, y=0.0, alpha=0.0;
    double lambda = 0.99;

    double dot(std::array<double, N> x, std::array<double, N> y);
    Matrix<N> outer_prod(std::array<double, N> x, std::array<double, N> y);
public:
    Rls(double alpha = 100, double lambda = 0.99); // alpha is initial variance for matrix P
    Rls(std::array<double, N> theta_initial,double alpha = 10, double lambda = 1); // alpha is initial variance for matrix P
    void update(std::array<double, N> x_new, double y_new);
    void print(); 
};


template<int N>
inline double Rls<N>::dot(std::array<double, N> x, std::array<double, N> y)
{
    double retv = 0;
    for (int i = 0; i < N; i++) {
        retv += x[i] * y[i]; 
    }
    return retv;
}

template<int N>
inline Matrix<N> Rls<N>::outer_prod(std::array<double, N> x, std::array<double, N> y)
{
    std::array<std::array<double, N>, N> retv;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            retv[i][j] = x[i] * y[j];
        }
    };
    return retv;
}

template<int N>
inline Rls<N>::Rls(double alpha,double lambda)
{
    for (int i = 0; i < N; i++) { P_old[i][i] = alpha; I[i][i] = 1.0; }
}

template<int N>
inline Rls<N>::Rls(std::array<double, N> theta_initial,double alpha, double lambda)
{
    theta_old = theta_initial;
    for (int i = 0; i < N; i++) {
        P_old[i][i] = alpha; I[i] = 1.0;
    };
}

template<int N>
void inline Rls<N>::update(std::array<double, N> x_new, double y_new)
{
    x = x_new; y = y_new;
    err = y - dot(x, theta_old);
    // new gamma
    for (int i = 0; i < N; i++) { tempv[i] = dot(P_old[i], x); }
    double c = 1.0 / (dot(x, tempv) + lambda);
    gamma = c * tempv;
    // new theta
    for (int i = 0; i < N; i++) { theta_old[i] += err * gamma[i]; }
    // new P
    P_old = (1.0 / lambda) * (I-outer_prod(gamma,x)) * P_old;    
}

template<int N>
inline void Rls<N>::print()
{
    std::cout << "Theta is: ";
    for (int i = 0; i < N; i++) { std::cout << theta_old[i] << " "; }
    std::cout << " error: " << err << "\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) std::cout << P_old[i][j] << "\t";
        std::cout << "\n";
    }
    std::cout << "\n";
}
