#include "StrokeSegmenter.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>

namespace Canoe {

static double compute_mac(const std::vector<double>& v) {
    if (v.size() <= 1) return 0.0;
    double diff_sum = 0.0;
    for (size_t i = 0; i < v.size() - 1; ++i) {
        diff_sum += std::abs(v[i + 1] - v[i]);
    }
    return diff_sum / (v.size() - 1);
}

static double compute_integral(const std::vector<double>& y, double dt) {
    if (y.size() <= 1) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < y.size() - 1; ++i) {
        sum += (y[i] + y[i + 1]);
    }
    return sum * 0.5 * dt;
}

StrokeSegmenter::StrokeSegmenter(double hysteresis, double min_phase_duration_s, double fs)
    : hysteresis_(hysteresis),
      min_phase_duration_s_(min_phase_duration_s),
      min_phase_samples_(static_cast<int>(min_phase_duration_s * fs)),
      fs_(fs),
      dt_(1.0 / fs),
      current_phase_(SegmenterPhase::ACCEL),
      has_last_accel_(false),
      initialized_(false),
      stroke_counter_(0) {
    current_buffer_.reserve(512);
    last_accel_buffer_.reserve(512);
}

std::optional<StrokeResult> StrokeSegmenter::process_sample(const ImuSample& sample) {
    if (!initialized_) {
        current_phase_ = (sample.acc_y >= 0.0) ? SegmenterPhase::ACCEL : SegmenterPhase::DECEL;
        initialized_ = true;
    }

    if (current_phase_ == SegmenterPhase::ACCEL && sample.acc_y < -hysteresis_) {
        if ((int)current_buffer_.size() >= min_phase_samples_) {
            // Accel phase finished. Store buffer as last_accel.
            last_accel_buffer_ = current_buffer_;
            has_last_accel_ = true;
            
            // Set state to decel and start decel phase with the current sample
            current_buffer_.clear();
            current_buffer_.push_back(sample);
            current_phase_ = SegmenterPhase::DECEL;
        } else {
            current_buffer_.push_back(sample);
        }
    } 
    else if (current_phase_ == SegmenterPhase::DECEL && sample.acc_y > hysteresis_) {
        if ((int)current_buffer_.size() >= min_phase_samples_) {
            // Decel phase finished. 
            std::optional<StrokeResult> result;
            if (has_last_accel_) {
                result = compute_stroke(last_accel_buffer_, current_buffer_);
            }
            
            // Set state to accel and start accel phase with the current sample
            last_accel_buffer_.clear();
            has_last_accel_ = false;
            
            current_buffer_.clear();
            current_buffer_.push_back(sample);
            current_phase_ = SegmenterPhase::ACCEL;
            
            return result;
        } else {
            current_buffer_.push_back(sample);
        }
    }
    else {
        current_buffer_.push_back(sample);
    }

    return std::nullopt;
}

std::optional<StrokeResult> StrokeSegmenter::flush() {
    if (current_phase_ == SegmenterPhase::DECEL && has_last_accel_ && (int)current_buffer_.size() >= min_phase_samples_) {
        StrokeResult result = compute_stroke(last_accel_buffer_, current_buffer_);
        last_accel_buffer_.clear();
        has_last_accel_ = false;
        current_buffer_.clear();
        return result;
    }
    return std::nullopt;
}

StrokeResult StrokeSegmenter::compute_stroke(const std::vector<ImuSample>& accel, const std::vector<ImuSample>& decel) {
    StrokeResult res;
    stroke_counter_++;
    res.stroke_id = stroke_counter_;

    // 1. Accel Phase Metrics
    res.accel = compute_phase_metrics(accel, true);
    res.time_to_pos_peak_s = res.accel.t_to_peak_pct * 0.01 * res.accel.duration_s;
    if (res.time_to_pos_peak_s > 0.0) {
        res.rad = res.accel.peak / res.time_to_pos_peak_s;
    } else {
        res.rad = 0.0;
    }

    // 2. Decel Phase Metrics
    res.decel = compute_phase_metrics(decel, false);
    res.time_to_neg_peak_s = res.decel.t_to_peak_pct * 0.01 * res.decel.duration_s;

    // 3. Whole Stroke Metrics (Accel + Decel combined)
    std::vector<double> stroke_acc_x;
    std::vector<double> stroke_gyr_x;
    std::vector<double> stroke_gyr_y;
    std::vector<double> stroke_gyr_z;

    stroke_acc_x.reserve(accel.size() + decel.size());
    stroke_gyr_x.reserve(accel.size() + decel.size());
    stroke_gyr_y.reserve(accel.size() + decel.size());
    stroke_gyr_z.reserve(accel.size() + decel.size());

    double g_x_max = -1e9;
    double g_x_min = 1e9;

    auto collect_samples = [&](const std::vector<ImuSample>& buf) {
        for (const auto& s : buf) {
            stroke_acc_x.push_back(s.acc_x);
            stroke_gyr_x.push_back(s.gyr_x);
            stroke_gyr_y.push_back(s.gyr_y);
            stroke_gyr_z.push_back(s.gyr_z);
            g_x_max = std::max(g_x_max, s.gyr_x);
            g_x_min = std::min(g_x_min, s.gyr_x);
        }
    };

    collect_samples(accel);
    collect_samples(decel);

    res.gyr_x_pos_peak = g_x_max;
    res.gyr_x_neg_peak = g_x_min;
    res.stroke_mac_acc_x = compute_mac(stroke_acc_x);
    res.stroke_mac_gyr_x = compute_mac(stroke_gyr_x);
    res.stroke_mac_gyr_y = compute_mac(stroke_gyr_y);
    res.stroke_mac_gyr_z = compute_mac(stroke_gyr_z);

    // 4. Consolidated Summary
    res.total_time_s = res.accel.duration_s + res.decel.duration_s;
    res.area_positive = res.accel.velocity_integral;
    res.area_negative = res.decel.velocity_integral;

    // Metadata indices (non-inclusive end matching Python)
    if (!accel.empty()) {
        res.accel_start_idx = accel.front().sample_idx;
        res.accel_end_idx = accel.back().sample_idx + 1;
    }
    if (!decel.empty()) {
        res.decel_start_idx = decel.front().sample_idx;
        res.decel_end_idx = decel.back().sample_idx + 1;
    }
    res.pos_peak_idx = res.accel.peak_sample_idx;
    res.neg_peak_idx = res.decel.peak_sample_idx;

    return res;
}

PhaseMetrics StrokeSegmenter::compute_phase_metrics(const std::vector<ImuSample>& phase_samples, bool is_accel) {
    PhaseMetrics m;
    int n = phase_samples.size();
    m.duration_s = n * dt_;

    std::vector<double> p_acc_x(n), p_acc_y(n), p_acc_z(n);
    std::vector<double> p_gyr_x(n), p_gyr_y(n), p_gyr_z(n);

    for (int i = 0; i < n; ++i) {
        p_acc_x[i] = phase_samples[i].acc_x;
        p_acc_y[i] = phase_samples[i].acc_y;
        p_acc_z[i] = phase_samples[i].acc_z;
        p_gyr_x[i] = phase_samples[i].gyr_x;
        p_gyr_y[i] = phase_samples[i].gyr_y;
        p_gyr_z[i] = phase_samples[i].gyr_z;
    }

    m.mac_acc_x = compute_mac(p_acc_x);
    m.mac_gyr_x = compute_mac(p_gyr_x);
    m.mac_gyr_y = compute_mac(p_gyr_y);
    m.mac_gyr_z = compute_mac(p_gyr_z);
    m.velocity_integral = compute_integral(p_acc_y, dt_);

    // Peak and time to peak
    if (n > 0) {
        if (is_accel) {
            auto max_it = std::max_element(p_acc_y.begin(), p_acc_y.end());
            m.peak = *max_it;
            int peak_idx = std::distance(p_acc_y.begin(), max_it);
            m.t_to_peak_pct = 100.0 * (peak_idx * dt_) / m.duration_s;
            m.peak_sample_idx = phase_samples[peak_idx].sample_idx;
        } else {
            auto min_it = std::min_element(p_acc_y.begin(), p_acc_y.end());
            m.peak = *min_it;
            int peak_idx = std::distance(p_acc_y.begin(), min_it);
            m.t_to_peak_pct = 100.0 * (peak_idx * dt_) / m.duration_s;
            m.peak_sample_idx = phase_samples[peak_idx].sample_idx;
        }
    }

    return m;
}

} // namespace Canoe
