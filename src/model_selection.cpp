#include <Rcpp.h>
#include <algorithm>

//' Function for VIF(Variance Inflation Factor) in QAIC
//' 
//' @description
//' \deqn{\hat{c}=\chi^2/df=\sum_{i,\tau}\frac{(m_{t_i,t_i+\tau}-E(m_{t_i,t_i+\tau}))^2}{E(m_{t_i,t_i+\tau})}/(unique(\tau)-k-1)}
//' 
//' @param data Observation matrix
//' @param tp List-like observation time(1d vector)
//' @param k the number 
//' @rdname chat
//' @export
//'
// [[Rcpp::export(LIR.chat)]]
double LIR_chat(Rcpp::NumericMatrix data, Rcpp::NumericVector tp, double k, double mtau = R_PosInf) {
    std::map<double, int> tau_map, cat_map;
    tau_map.clear(), cat_map.clear();
    int N = data.nrow(), T = data.ncol();

    double tau;
    // count
    for (int i = 0, m; i < T; i++)
        for (int j = 0; j < T; j++) {
            tau = tp[i] - tp[j];
            if (tau <= mtau) {
                if (tau_map.find(tau) == tau_map.end()) tau_map[tau]++;
                else    tau_map[tau] = 1;
                // m = 0;
                // for (int k = 0; k < N; k++)
                //     if (data(i, k) == data(j, k) && data(i, k) > 0)
                //         m += 1;
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
            bool loop = true;
            while(loop) {
                cat_map[it2->first] = ncat;
                it2++;
            }
        }
    }
    // Use recurrence formula of variance to calculate Chi^2
    int n_m = 0;
    double mean = 0, var = 0.0;
    n_m++;
        var += (m - mean) * (n_m - 1) / n_m;
        mean += (m - mean) / n_m;
}