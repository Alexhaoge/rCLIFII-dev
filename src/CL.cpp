#include <Rcpp.h>
#include <cstring>
#include <cmath>

using namespace Rcpp;

//' @export
//' @rdname CL
// [[Rcpp::export(LIR.CL)]]
double LIR_CL(NumericVector theta, Function model, NumericMatrix data, NumericVector tp, List model_args = List(), double mtau = R_PosInf) {
    int N = data.nrow(), T = data.ncol();
    int* n = (int*)malloc(sizeof(int) * T);
    for (int j = 0, i; j < T; j++)
        for (n[j] = 0, i = 0; i < N; i++)
            if (data(i, j) > 0)
                n[j]++;
    double cl = 0.0, tau;
    for (int i = 0, m; i < T; i++)
        for (int j = i + 1; j < T; j++) {
            tau = tp[j] - tp[i];
            if (tau <= mtau) {
                m = 0;
                for (int k = 0; k < N; k++)
                    if (data(i, k) == data(j, k) && data(i, k) > 0)
                        m += 1;
                double nR_tau = n[j] * (*REAL(model(theta, tau, Named("model_args", model_args))));
                if (nR_tau <= 0 || nR_tau >= 1) return -9e12;
                cl += m * log(nR_tau) + (n[i] - m) * log(1 - nR_tau);
            }
        }
    return cl;
}

//' @export
//' @rdname CLgrad
// [[Rcpp::export(LIR.CLgrad)]]
NumericVector LIR_Grad(NumericVector theta, Function model, Function grad, NumericMatrix data, NumericVector tp, List model_args = List(), double mtau = R_PosInf) {
    int N = data.nrow(), T = data.ncol();
    int* n = (int*)malloc(sizeof(int) * T);
    for (int j = 0, i; j < T; j++)
        for (n[j] = 0, i = 0; i < N; i++)
            if (data(i, j) > 0)
                n[j]++;
    double tau, R_tau;
    NumericVector cl_grad(theta.size()), R_tau_grad;
    for (int i = 0, m; i < T; i++)
        for (int j = i + 1; j < T; j++) {
            tau = tp[j] - tp[i];
            if (tau <= mtau) {
                m = 0;
                for (int k = 0; k < N; k++)
                    if (data(i, k) == data(j, k) && data(i, k) > 0)
                        m += 1;
                R_tau = *REAL(model(theta, Named("model_args", model_args)));
                R_tau_grad = as<NumericVector>(grad(theta, tau, Named("model_args", model_args)));
                cl_grad += R_tau_grad * (n[j] * (n[i] - m) / (1 - R_tau) - m / R_tau);
            }
        }
    return cl_grad;
}

//' @export
//' @rdname CLhessian
// [[Rcpp::export(LIR.CLhessian)]]
NumericMatrix LIR_Hessian(NumericVector theta, Function model, Function grad, Function hessian, NumericMatrix data, NumericVector tp, List model_args = List(), double mtau = R_PosInf) {
    int N = data.nrow(), T = data.ncol();
    int* n = (int*)malloc(sizeof(int) * T);
    for (int j = 0, i; j < T; j++)
        for (n[j] = 0, i = 0; i < N; i++)
            if (data(i, j) > 0)
                n[j]++;
    double tau, R_tau;
    NumericVector R_tau_grad;
    NumericMatrix cl_hessian(theta.size()), R_tau_hessian;
    for (int i = 0, m; i < T; i++)
        for (int j = i + 1; j < T; j++) {
            tau = tp[j] - tp[i];
            if (tau <= mtau) {
                m = 0;
                for (int k = 0; k < N; k++)
                    if (data(i, k) == data(j, k) && data(i, k) > 0)
                        m += 1;
                R_tau = *REAL(model(theta, Named("model_args", model_args)));
                R_tau_grad = as<NumericVector>(grad(theta, tau, Named("model_args", model_args)));
                R_tau_hessian = as<NumericMatrix>(hessian(theta, tau, Named("model_args", model_args)));
                for (int r = 0; r < theta.size(); r++)
                    for (int c = 0; c < theta.size(); c++)
                        cl_hessian(r, c) += R_tau_hessian(r, c) * (n[j] * (n[i] - m) / (1 - R_tau) - m / R_tau)\
                         + R_tau_grad[r] * R_tau_grad[c] * (n[j] * n[j] * (n[i] - m) / (1 - R_tau) / (1 - R_tau) + m / R_tau / R_tau);
            }
        }
    return cl_hessian;
}