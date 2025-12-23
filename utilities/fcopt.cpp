#include "fcopt.h"
#include "filters.h"
#include <cmath>
#include <algorithm>
#include <numeric>

static double psa95(const std::vector<double>& x, double fs) {
    int N = x.size();
    std::vector<double> psd(N/2);

    for (int k = 0; k < (int)psd.size(); ++k) {
        double re = 0, im = 0;
        for (int n = 0; n < N; ++n) {
            double w = 2 * M_PI * k * n / N;
            re += x[n] * cos(w);
            im -= x[n] * sin(w);
        }
        psd[k] = re*re + im*im;
    }

    double total = std::accumulate(psd.begin(), psd.end(), 0.0);
    double cum = 0;

    for (int k = 0; k < (int)psd.size(); ++k) {
        cum += psd[k];
        if (cum / total >= 0.95)
            return (double)k * fs / N;
    }
    return fs / 2.0;
}

static double residual_analysis(
    const std::vector<double>& x,
    double fs
) {
    std::vector<double> rms_residuals;
    std::vector<double> cutoffs;

    for (double fc = 1.0; fc <= 15.0; fc += 0.5) {
        ButterworthLowpass lp(fc, fs, 2);
        auto xf = lp.apply(x);

        double rms = 0;
        for (size_t i = 0; i < x.size(); ++i)
            rms += (x[i] - xf[i]) * (x[i] - xf[i]);
        rms = sqrt(rms / x.size());

        rms_residuals.push_back(rms);
        cutoffs.push_back(fc);
    }

    int idx = std::distance(
        rms_residuals.begin(),
        std::min_element(rms_residuals.begin(), rms_residuals.end())
    );

    return cutoffs[idx];
}

static double autocorr_method(
    const std::vector<double>& x,
    double fs
) {
    int N = x.size();
    std::vector<double> ac(N, 0.0);

    for (int lag = 0; lag < N; ++lag) {
        for (int i = 0; i + lag < N; ++i)
            ac[lag] += x[i] * x[i + lag];
    }

    int first_zero = 1;
    while (first_zero < N && ac[first_zero] > 0)
        first_zero++;

    double freq = fs / first_zero;
    return std::min(freq * 2.0, fs / 2.0);
}

static std::vector<double> downsample(
    const std::vector<double>& x,
    int factor
) {
    std::vector<double> y;
    for (size_t i = 0; i < x.size(); i += factor)
        y.push_back(x[i]);
    return y;
}

static double dominant_cluster_center(
    const std::vector<double>& values
) {
    std::vector<double> v = values;
    std::sort(v.begin(), v.end());

    double c1 = v[v.size()/3];
    double c2 = v[2*v.size()/3];

    for (int iter = 0; iter < 20; ++iter) {
        double s1=0,c1n=0, s2=0,c2n=0;
        for (double x : v) {
            if (fabs(x - c1) < fabs(x - c2)) {
                s1 += x; c1n++;
            } else {
                s2 += x; c2n++;
            }
        }
        if (c1n) c1 = s1 / c1n;
        if (c2n) c2 = s2 / c2n;
    }

    int n1=0,n2=0;
    for (double x : v)
        fabs(x - c1) < fabs(x - c2) ? n1++ : n2++;

    return (n1 > n2) ? c1 : c2;
}

double fcopt_cutoff_frequency(
    const std::vector<double>& signal,
    double fs
) {
    std::vector<double> candidates;

    for (int ds : {1, 2, 3}) {
        auto x = (ds == 1) ? signal : downsample(signal, ds);
        double fds = fs / ds;

        candidates.push_back(psa95(x, fds));
        candidates.push_back(residual_analysis(x, fds));
        candidates.push_back(autocorr_method(x, fds));
    }

    return dominant_cluster_center(candidates);
}
