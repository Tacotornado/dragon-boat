#pragma once

#include <vector>
#include <optional>
#include <string>

namespace Canoe {

struct ImuSample {
    uint64_t timestamp_us = 0;
    int sample_idx = 0;
    double acc_x = 0.0, acc_y = 0.0, acc_z = 0.0;
    double gyr_x = 0.0, gyr_y = 0.0, gyr_z = 0.0;
};

struct PhaseMetrics {
    double duration_s = 0.0;
    double peak = 0.0;
    double t_to_peak_pct = 0.0;
    double velocity_integral = 0.0;
    double mac_acc_x = 0.0;
    double mac_gyr_x = 0.0;
    double mac_gyr_y = 0.0;
    double mac_gyr_z = 0.0;
    int peak_sample_idx = 0;
};

struct StrokeResult {
    int stroke_id = 0;
    
    // Accel phase
    PhaseMetrics accel;
    double rad = 0.0;

    // Decel phase
    PhaseMetrics decel;

    // Whole stroke
    double gyr_x_pos_peak = 0.0;
    double gyr_x_neg_peak = 0.0;
    double stroke_mac_acc_x = 0.0;
    double stroke_mac_gyr_x = 0.0;
    double stroke_mac_gyr_y = 0.0;
    double stroke_mac_gyr_z = 0.0;

    // Consolidated metrics
    double total_time_s = 0.0;
    double time_to_pos_peak_s = 0.0;
    double time_to_neg_peak_s = 0.0;
    double area_positive = 0.0;
    double area_negative = 0.0;

    // Metadata for index comparison
    int accel_start_idx = 0;
    int accel_end_idx = 0;
    int decel_start_idx = 0;
    int decel_end_idx = 0;
    int pos_peak_idx = 0;
    int neg_peak_idx = 0;
};

enum class SegmenterPhase { ACCEL, DECEL };

class StrokeSegmenter {
public:
    StrokeSegmenter(double hysteresis = 0.15, double min_phase_duration_s = 0.08, double fs = 120.0);

    // Processes a single filtered sample and returns a completed stroke if available
    std::optional<StrokeResult> process_sample(const ImuSample& sample);
    std::optional<StrokeResult> flush();

private:
    double hysteresis_;
    double min_phase_duration_s_;
    int min_phase_samples_;
    double fs_;
    double dt_;

    SegmenterPhase current_phase_;
    std::vector<ImuSample> current_buffer_;
    std::vector<ImuSample> last_accel_buffer_;
    bool has_last_accel_;
    bool initialized_;
    int stroke_counter_;

    StrokeResult compute_stroke(const std::vector<ImuSample>& accel, const std::vector<ImuSample>& decel);
    PhaseMetrics compute_phase_metrics(const std::vector<ImuSample>& phase_samples, bool is_accel);
};

} // namespace Canoe
