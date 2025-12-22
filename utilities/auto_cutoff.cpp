#include "auto_cutoff.h"
#include <complex>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace AutoCutoff {

using cd = std::complex<double>;

// recursive FFT (inplace conceptual; returns vector)
static std::vector<cd> fft_rec(const std::vector<cd>& x) {
    size_t N = x.size();
    if (N == 1) return { x[0] };
    if (N % 2 != 0) {
        // fallback to DFT if not power of two
        std::vector<cd> X(N);
        for (size_t k = 0; k < N; ++k) {
            cd sum = 0;
            for (size_t n = 0; n < N; ++n)
                sum += x[n] * std::polar(1.0, -2.0 * M_PI * k * n / N);
            X[k] = sum;
        }
        return X;
    }
    std::vector<cd> even(N/2), odd(N/2);
    for (size_t i = 0; i < N/2; ++i) {
        even[i] = x[2*i];
        odd[i]  = x[2*i+1];
    }
    auto Fe = fft_rec(even);
    auto Fo = fft_rec(odd);
    std::vector<cd> X(N);
    for (size_t k = 0; k < N/2; ++k) {
        cd t = std::polar(1.0, -2.0 * M_PI * k / N) * Fo[k];
        X[k] = Fe[k] + t;
        X[k + N/2] = Fe[k] - t;
    }
    return X;
}

static size_t next_pow2(size_t n) {
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

static std::vector<double> smooth_vec(const std::vector<double>& v, int win=3) {
    int n = (int)v.size();
    if (n==0) return {};
    std::vector<double> out(n);
    int half = win/2;
    for (int i=0;i<n;++i) {
        int a = std::max(0, i-half);
        int b = std::min(n-1, i+half);
        double s=0; int cnt=0;
        for (int j=a;j<=b;++j){ s += v[j]; ++cnt; }
        out[i] = s / std::max(1,cnt);
    }
    return out;
}

double estimateCutoffKneePoint(const std::vector<double>& signal, double fs,
                               double f_min, double f_max, double fallback) {
    if (signal.empty() || fs <= 0) return fallback;

    // require at least 0.5 s
    if (signal.size() < (size_t)std::max(16UL, (size_t)std::ceil(0.5*fs))) return fallback;

    // prepare zero-padded signal length (power of two)
    size_t N0 = signal.size();
    size_t N = next_pow2(N0);
    // pad to at least N (if equal fine), but ensure N >= 32
    if (N < 32) N = 32;

    std::vector<cd> x(N, cd(0.0, 0.0));
    for (size_t i = 0; i < N0; ++i) x[i] = cd(signal[i], 0.0);

    auto X = fft_rec(x);

    // single-sided PSD magnitude^2, frequencies up to Nyquist (N/2)
    size_t M = N/2;
    std::vector<double> psd(M);
    for (size_t k = 0; k < M; ++k) {
        double mag2 = std::norm(X[k]) / (double)N; // scale by N
        // for single-sided, double non-DC bins
        if (k != 0 && k != M) mag2 *= 2.0;
        psd[k] = mag2;
    }

    // optional smoothing to reduce spiky PSD
    psd = smooth_vec(psd, 5);

    // cumulative energy normalized
    std::vector<double> cum(M);
    cum[0] = psd[0];
    for (size_t i = 1; i < M; ++i) cum[i] = cum[i-1] + psd[i];
    double total = cum.empty() ? 0.0 : cum.back();
    if (!(total > 0.0)) return fallback; // all zeros?

    for (size_t i = 0; i < M; ++i) cum[i] /= total; // normalize to [0,1]

    // elbow detection: find index with max perpendicular distance from line joining (0, cum0=0) to (M-1, cum_last=1)
    // line in (i, cum[i]) coordinates from (0,0) to (M-1,1)
    double x1 = 0.0, y1 = 0.0;
    double x2 = (double)(M-1), y2 = 1.0;
    double denom = std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    if (denom <= 0.0) return fallback;

    double maxDist = -1.0;
    size_t maxIdx = 0;
    for (size_t i = 0; i < M; ++i) {
        double xi = (double)i;
        double yi = cum[i];
        // perpendicular distance formula
        double num = std::abs( (y2-y1)*xi - (x2-x1)*yi + x2*y1 - y2*x1 );
        double dist = num / denom;
        if (dist > maxDist) { maxDist = dist; maxIdx = i; }
    }

    // convert index to frequency
    double freq_res = fs / (double)N; // bin width
    double cutoff = maxIdx * freq_res;

    // fallback if too small or NaN
    if (!std::isfinite(cutoff) || cutoff <= 0.0) cutoff = fallback;

    // clamp to [f_min, f_max]
    cutoff = std::max(f_min, std::min(f_max, cutoff));

    return cutoff;
}

} // namespace AutoCutoff