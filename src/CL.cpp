#include <Rcpp.h>
#include <cstring>
#include <cmath>
#define CL_INF 9e12

//' @export
//' @rdname CL
// [[Rcpp::export(LIR.CL)]]
double LIR_CL(
    const Rcpp::NumericVector &theta, const Rcpp::Function &model,
    const Rcpp::NumericMatrix &data, const Rcpp::NumericVector &tp,
    const Rcpp::List &model_args, const double &mtau = -1.0)
{
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
            if (mtau < 0 || tau <= mtau) {
                m = 0;
                for (int k = 0; k < N; k++)
                    if (data(i, k) == data(j, k) && data(i, k) > 0)
                        m += 1;
                double nR_tau = n[j] * (*REAL(model(theta, tau, Rcpp::Named("model_args", model_args))));
                if (nR_tau <= 0 || nR_tau >= 1) return -CL_INF;
                cl += m * log(nR_tau) + (n[i] - m) * log(1 - nR_tau);
            }
        }
    return cl;
}

//' @export
//' @rdname CLgrad
// [[Rcpp::export(LIR.CLgrad)]]
Rcpp::NumericVector LIR_Grad(
    const Rcpp::NumericVector &theta,
    const Rcpp::Function &model, const Rcpp::Function &grad,
    const Rcpp::NumericMatrix &data, const Rcpp::NumericVector &tp,
    const Rcpp::List &model_args, const double &mtau = -1.0)
{
    int N = data.nrow(), T = data.ncol();
    int* n = (int*)malloc(sizeof(int) * T);
    for (int j = 0, i; j < T; j++)
        for (n[j] = 0, i = 0; i < N; i++)
            if (data(i, j) > 0)
                n[j]++;
    double tau, R_tau;
    Rcpp::NumericVector cl_grad(theta.size()), R_tau_grad;
    for (int i = 0, m; i < T; i++)
        for (int j = i + 1; j < T; j++) {
            tau = tp[j] - tp[i];
            if (mtau < 0 || tau <= mtau) {
                m = 0;
                for (int k = 0; k < N; k++)
                    if (data(i, k) == data(j, k) && data(i, k) > 0)
                        m += 1;
                R_tau = *REAL(model(theta, Rcpp::Named("model_args", model_args)));
                R_tau_grad = Rcpp::as<Rcpp::NumericVector>(grad(theta, tau, Rcpp::Named("model_args", model_args)));
                cl_grad += R_tau_grad * (n[j] * (n[i] - m) / (1 - R_tau) - m / R_tau);
            }
        }
    return cl_grad;
}

//' @export
//' @rdname CLhessian
// [[Rcpp::export(LIR.CLhessian)]]
Rcpp::NumericMatrix LIR_Hessian(
    const Rcpp::NumericVector &theta,
    const Rcpp::Function &model, const Rcpp::Function &grad, const Rcpp::Function &hessian,
    const Rcpp::NumericMatrix &data, const Rcpp::NumericVector &tp,
    const Rcpp::List &model_args, const double &mtau = -1.0)
{
    int N = data.nrow(), T = data.ncol();
    int* n = (int*)malloc(sizeof(int) * T);
    for (int j = 0, i; j < T; j++)
        for (n[j] = 0, i = 0; i < N; i++)
            if (data(i, j) > 0)
                n[j]++;
    double tau, R_tau;
    Rcpp::NumericVector R_tau_grad;
    Rcpp::NumericMatrix cl_hessian(theta.size()), R_tau_hessian;
    for (int i = 0, m; i < T; i++)
        for (int j = i + 1; j < T; j++) {
            tau = tp[j] - tp[i];
            if (mtau < 0 || tau <= mtau) {
                m = 0;
                for (int k = 0; k < N; k++)
                    if (data(i, k) == data(j, k) && data(i, k) > 0)
                        m += 1;
                R_tau = *REAL(model(theta, Named("model_args", model_args)));
                R_tau_grad = Rcpp::as<Rcpp::NumericVector>(grad(theta, tau, Rcpp::Named("model_args", model_args)));
                R_tau_hessian = Rcpp::as<Rcpp::NumericMatrix>(hessian(theta, tau, Rcpp::Named("model_args", model_args)));
                for (int r = 0; r < theta.size(); r++)
                    for (int c = 0; c < theta.size(); c++)
                        cl_hessian(r, c) += R_tau_hessian(r, c) * (n[j] * (n[i] - m) / (1 - R_tau) - m / R_tau)\
                        + R_tau_grad[r] * R_tau_grad[c] * (n[j] * n[j] * (n[i] - m) / (1 - R_tau) / (1 - R_tau) + m / R_tau / R_tau);
            }
        }
    return cl_hessian;
}

/********************************
 * CL with built-in models
 *******************************/

double LIR_model_builtin(const char &model, const Rcpp::NumericVector &theta, const double &tau) {
    switch (model) {
    case 'A':   return theta[0]; break;
    case 'B':   return theta[0] * exp(-theta[1] * tau); break;
    case 'C':   return theta[0] * exp(-theta[1] * tau) + theta[2]; break;
    default:    throw std::invalid_argument("Unsupported model, only\
                                            \"A\"/\"B\"/\"C\" are allowed.");
    }
}

void LIR_model_grad_builtin(const char &model, const Rcpp::NumericVector &theta, const double &tau, Rcpp::NumericVector &grad) {
    switch (model) {
    case 'A':   grad[0] = 1; break;
    case 'B':   
        grad[0] = exp(-theta[1]*tau), grad[1] = -theta[0] * tau * grad[0];
        break;
    case 'C':   
        grad[0] = exp(-theta[1]*tau), grad[1] = -theta[0] * tau * grad[0], grad[2] = 1;
        break;
    default:    throw std::invalid_argument("Unsupported model, only\
                                            \"A\"/\"B\"/\"C\" are allowed.");
    }
}

void LIR_model_hessian_builtin(const char &model, const Rcpp::NumericVector &theta, const double &tau, Rcpp::NumericMatrix &hessian) {
    switch (model) {
    case 'A':   hessian(0, 0) = 1; break;
    case 'B':
        hessian(0, 1) = hessian(1, 0) = -tau * exp(-theta[1]*tau);
        hessian(1, 1) = -tau * theta[0] * hessian(0, 1);
        break;
    case 'C':
        hessian(0, 1) = hessian(1, 0) = -tau * exp(-theta[1]*tau);
        hessian(1, 1) = -tau * theta[0] * hessian(0, 1);
        break;
    default:    throw std::invalid_argument("Unsupported model, only\
                                            \"A\"/\"B\"/\"C\" are allowed.");
    }
}


// Built-in function for CL, use Cpp function for model directly
double LIR_CL_builtin(const Rcpp::NumericVector &theta, const Rcpp::String &model,
    const Rcpp::NumericMatrix &data, const Rcpp::NumericVector &tp, const double &mtau = -1.0)
{
    char _model = model.get_cstring()[0];
    int N = data.nrow(), T = data.ncol();
    int* n = (int*)malloc(sizeof(int) * T);
    for (int j = 0, i; j < T; j++)
        for (n[j] = 0, i = 0; i < N; i++)
            if (data(i, j) > 0)
                n[j]++;
    double cl = 0.0, tau, nR_tau;
    for (int i = 0, m; i < T; i++)
        for (int j = i + 1; j < T; j++) {
            tau = tp[j] - tp[i];
            if (mtau < 0 || tau <= mtau) {
                m = 0;
                for (int k = 0; k < N; k++)
                    if (data(i, k) == data(j, k) && data(i, k) > 0)
                        m += 1;
                nR_tau = LIR_model_builtin(_model, theta, tau);
                if (nR_tau <= 0 || nR_tau >= 1) return -CL_INF;
                cl += m * log(nR_tau) + (n[i] - m) * log(1 - nR_tau);
            }
        }
    return cl;
}

Rcpp::NumericVector LIR_grad_builtin(
    Rcpp::NumericVector theta, Rcpp::String model,
    Rcpp::NumericMatrix data, Rcpp::NumericVector tp,
    Rcpp::List model_args, double mtau = -1.0)
{
    char _model = model.get_cstring()[0];
    int N = data.nrow(), T = data.ncol();
    int* n = (int*)malloc(sizeof(int) * T);
    for (int j = 0, i; j < T; j++)
        for (n[j] = 0, i = 0; i < N; i++)
            if (data(i, j) > 0)
                n[j]++;
    double tau, R_tau;
    Rcpp::NumericVector cl_grad(theta.size()), R_tau_grad(theta.size());
    for (int i = 0, m; i < T; i++)
        for (int j = i + 1; j < T; j++) {
            tau = tp[j] - tp[i];
            if (mtau < 0 || tau <= mtau) {
                m = 0;
                for (int k = 0; k < N; k++)
                    if (data(i, k) == data(j, k) && data(i, k) > 0)
                        m += 1;
                R_tau = LIR_model_builtin(_model, theta, tau);
                LIR_model_grad_builtin(_model, theta, tau, R_tau_grad);
                cl_grad += R_tau_grad * (n[j] * (n[i] - m) / (1 - R_tau) - m / R_tau);
            }
        }
    return cl_grad;
}

Rcpp::NumericMatrix LIR_Hessian_builtin(
    const Rcpp::NumericVector &theta, const Rcpp::String &model,
    const Rcpp::NumericMatrix &data, const Rcpp::NumericVector &tp,
    const Rcpp::List &model_args, const double &mtau = -1.0)
{
    char _model = model.get_cstring()[0];
    int N = data.nrow(), T = data.ncol();
    int* n = (int*)malloc(sizeof(int) * T);
    for (int j = 0, i; j < T; j++)
        for (n[j] = 0, i = 0; i < N; i++)
            if (data(i, j) > 0)
                n[j]++;
    double tau, R_tau;
    Rcpp::NumericVector R_tau_grad(theta.size());
    Rcpp::NumericMatrix cl_hessian(theta.size()), R_tau_hessian(theta.size());
    for (int i = 0, m; i < T; i++)
        for (int j = i + 1; j < T; j++) {
            tau = tp[j] - tp[i];
            if (mtau < 0 || tau <= mtau) {
                m = 0;
                for (int k = 0; k < N; k++)
                    if (data(i, k) == data(j, k) && data(i, k) > 0)
                        m += 1;
                R_tau = LIR_model_builtin(_model, theta, tau);
                LIR_model_grad_builtin(_model, theta, tau, R_tau_grad);
                LIR_model_hessian_builtin(_model, theta, tau, R_tau_hessian);
                for (int r = 0; r < theta.size(); r++)
                    for (int c = 0; c < theta.size(); c++)
                        cl_hessian(r, c) += R_tau_hessian(r, c) * (n[j] * (n[i] - m) / (1 - R_tau) - m / R_tau)\
                        + R_tau_grad[r] * R_tau_grad[c] * (n[j] * n[j] * (n[i] - m) / (1 - R_tau) / (1 - R_tau) + m / R_tau / R_tau);
            }
        }
    return cl_hessian;
}
