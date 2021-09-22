#include <Rcpp.h>
#include <algorithm>

//' Function for VIF(Variance Inflation Factor) in QAIC
//'
//' @description
//' \deqn{\hat{c}=\chi^2/df=\sum_{i,\tau}\frac{(m_{t_i,t_i+\tau}-E(m_{t_i,t_i+\tau}))^2}{E(m_{t_i,t_i+\tau})}/(unique(\tau)-k-1)}
//'
//' @note This function has memory complexity of O(N^2) and may suffer from out-of-memory
//' error if observation number is very large.
//'
//' @param data Observation matrix
//' @param tp List-like observation time(1d vector)
//' @param k the number of variable in the most general model
//' @param mtau The maximum allowable lag time. If a lagged pair has time \eqn{\tau}
//'   greater than `mtau`, it will not be used to calculate composite likelihood.
//'   If `mtau` is less than zero, all pairs will be used. Default -1.0.
//' @rdname chat
//' @export
//'
// [[Rcpp::export(LIR.chat)]]
double LIR_chat(Rcpp::NumericMatrix data, Rcpp::NumericVector tp, double k, double mtau = -1.0) {
    // This implementation uses three loops and looks clumsy
    // It is designed to avoid record all pairwise `m` into different `tau` categories,
    // which is O(N^2) in space. Still using tau_map for bin count of tau would be
    // O(N^2) in the worst case. Is there any better solution with O(N) space comlexity?
    // Maybe the best
    std::map<double, int> tau_map, cat_map;
    tau_map.clear(), cat_map.clear();
    int N = data.nrow(), T = data.ncol();

    double tau;
    // count
    for (int i = 0; i < T; i++)
        for (int j = 0; j < T; j++) {
            tau = tp[i] - tp[j];
            if (mtau < 0 || tau <= mtau) {
                if (tau_map.find(tau) == tau_map.end()) tau_map[tau]++;
                else    tau_map[tau] = 1;
            }
        }
    // lump the catgories
    std::map<double, int>::iterator it, it2, end=tau_map.end();
    std::map<double, int>::reverse_iterator last=tau_map.rbegin();
    int tmp = 0, ncat = 0;
    for (it = it2 = tau_map.begin(); it != end; it++) {
        tmp += it->second;
        if (tmp >= 6 || it->first == last->first) {
            ncat++;
            tmp = 0;
            while(true) {
                cat_map[it2->first] = ncat;
                if (it2->first == it->first) break;
                    else it2++;
            }
        }
    }
    // Use recurrence formula of variance to calculate Chi^2
    int n_m = 0;
    std::map<double, double> var, mean;
    var.clear(), mean.clear();
    n_m++;
    for (int i = 0, m; i < T; i++)
        for (int j = 0; j < T; j++) {
            tau = tp[i] - tp[j];
            if (mtau < 0 || tau <= mtau) {
                m = 0;
                for (int k = 0; k < N; k++)
                    if (data(k, i) == data(k, j) && data(k, i) > 0)
                        m += 1;
                if (mean.find(tau) == mean.end()) {
                    mean[tau] = m, var[tau] = 0.0;
                } else {
                    var[tau] += (m - mean[tau]) * (n_m - 1) / n_m;
                    mean[tau] += (m - mean[tau]) / n_m;
                }
            }
        }
    std::map<double, double>::iterator it3, end2 = mean.end();
    double X2 = 0.0;
    for (it3 = mean.begin(); it3 != end2; it++)
        X2 += var[it3->first] / it3->second;
    return X2 / (ncat - k - 1.0);
}
