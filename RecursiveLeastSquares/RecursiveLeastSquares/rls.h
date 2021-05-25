#pragma once
#include <iostream>
#include <array> 

template <typename T,size_t Row, size_t Col> using Matrix = std::array<std::array<T, Col>, Row>;

template <int N>
class Rls
{
    Matrix<double, N, N> P_old{ { 0.0 } }, I{ { 0.0 } };
    std::array<double, N> theta_old{ 0.0 }, gamma{ 0.0 }, x{ 0.0 }, temp{ 0.0 };
    double err=0.0, y=0.0, alpha=0.0;
    double lambda = 0.99;

    double scalar(std::array<double, N> x, std::array<double, N> y);
    Matrix<double, N, N> outer_prod(std::array<double, N> x, std::array<double, N> y);
    Matrix<double, N, N> subtract_mat(Matrix<double, N, N> X, Matrix<double, N, N> Y);
    Matrix<double, N, N> product_mat(Matrix<double, N, N> X, Matrix<double, N, N> Y);
    Matrix<double, N, N> scalar_mat(double x, Matrix<double, N, N> Y);
public:
    Rls(double alpha = 100, double lambda = 0.99); // alpha is initial variance for matrix P
    Rls(std::array<double, N> theta_initial,double alpha = 10, double lambda = 1); // alpha is initial variance for matrix P
    void update(std::array<double, N> x_new, double y_new);
    void print(); 
};


template<int N>
inline double Rls<N>::scalar(std::array<double, N> x, std::array<double, N> y)
{
    double retv = 0;
    for (int i = 0; i < N; i++) {
        retv += x[i] * y[i]; 
    }
    return retv;
}

template<int N>
inline Matrix<double, N, N> Rls<N>::outer_prod(std::array<double, N> x, std::array<double, N> y)
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
inline Matrix<double, N, N> Rls<N>::subtract_mat(Matrix<double, N, N> X, Matrix<double, N, N> Y)
{
    Matrix<double, N, N> retv;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            retv[i][j] = X[i][j] - Y[i][j];
        }
    };
    return retv;
}

template<int N>
inline Matrix<double, N, N> Rls<N>::product_mat(Matrix<double, N, N> X, Matrix<double, N, N> Y)
{
    std::array<std::array<double, N>, N> retv{ {0.0} };
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                retv[i][j] += X[i][k] * Y[k][j];
            }
        }
    };
    return retv;
}

template<int N>
inline Matrix<double, N, N> Rls<N>::scalar_mat(double x, Matrix<double, N, N> Y)
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Y[i][j] *=x;
        }
    };
    return Y;
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

    err = y - scalar(x, theta_old);
    // new gamma
    for (int i = 0; i < N; i++) { temp[i] = scalar(P_old[i], x); }
    double c = 1.0 / (scalar(x, temp) + lambda);
    for (int i = 0; i < N; i++) { gamma[i] = c * temp[i]; }
    // new theta
    for (int i = 0; i < N; i++) { theta_old[i] += err * gamma[i]; }
    // new P
    P_old= scalar_mat(1.0 / lambda, product_mat(subtract_mat(I,outer_prod(gamma,x)), P_old));    
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
