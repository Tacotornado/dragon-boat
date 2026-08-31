import type { StrokeCycle, AccPhaseMetrics, DecPhaseMetrics, StrokePhaseMetrics, IMUDataPoint } from "./types"

// Realistic mock data based on the C++ stroke analysis output
// Simulates 120Hz IMU data from a kayak/dragon boat paddle sensor
export const mockStrokes: StrokeCycle[] = [
  { StrokeID: 1, StartZero: 45, PosPeak: 78, MidZero: 112, NegPeak: 138, EndZero: 165, TotalTime: 1.00, TimeToPosPeak: 0.275, TimeToNegPeak: 0.217, Area_Positive: 2.34, Area_Negative: 1.89 },
  { StrokeID: 2, StartZero: 165, PosPeak: 198, MidZero: 230, NegPeak: 256, EndZero: 282, TotalTime: 0.975, TimeToPosPeak: 0.275, TimeToNegPeak: 0.217, Area_Positive: 2.51, Area_Negative: 2.03 },
  { StrokeID: 3, StartZero: 282, PosPeak: 314, MidZero: 348, NegPeak: 375, EndZero: 401, TotalTime: 0.992, TimeToPosPeak: 0.267, TimeToNegPeak: 0.225, Area_Positive: 2.68, Area_Negative: 2.15 },
  { StrokeID: 4, StartZero: 401, PosPeak: 435, MidZero: 467, NegPeak: 494, EndZero: 520, TotalTime: 0.992, TimeToPosPeak: 0.283, TimeToNegPeak: 0.225, Area_Positive: 2.73, Area_Negative: 2.22 },
  { StrokeID: 5, StartZero: 520, PosPeak: 552, MidZero: 585, NegPeak: 611, EndZero: 638, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.217, Area_Positive: 2.81, Area_Negative: 2.35 },
  { StrokeID: 6, StartZero: 638, PosPeak: 670, MidZero: 702, NegPeak: 730, EndZero: 756, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.233, Area_Positive: 2.92, Area_Negative: 2.41 },
  { StrokeID: 7, StartZero: 756, PosPeak: 790, MidZero: 822, NegPeak: 848, EndZero: 875, TotalTime: 0.992, TimeToPosPeak: 0.283, TimeToNegPeak: 0.217, Area_Positive: 2.85, Area_Negative: 2.30 },
  { StrokeID: 8, StartZero: 875, PosPeak: 908, MidZero: 940, NegPeak: 966, EndZero: 993, TotalTime: 0.983, TimeToPosPeak: 0.275, TimeToNegPeak: 0.217, Area_Positive: 2.77, Area_Negative: 2.28 },
  { StrokeID: 9, StartZero: 993, PosPeak: 1025, MidZero: 1058, NegPeak: 1085, EndZero: 1112, TotalTime: 0.992, TimeToPosPeak: 0.267, TimeToNegPeak: 0.225, Area_Positive: 2.65, Area_Negative: 2.18 },
  { StrokeID: 10, StartZero: 1112, PosPeak: 1144, MidZero: 1177, NegPeak: 1203, EndZero: 1230, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.217, Area_Positive: 2.58, Area_Negative: 2.12 },
  { StrokeID: 11, StartZero: 1230, PosPeak: 1264, MidZero: 1296, NegPeak: 1322, EndZero: 1348, TotalTime: 0.983, TimeToPosPeak: 0.283, TimeToNegPeak: 0.217, Area_Positive: 2.70, Area_Negative: 2.25 },
  { StrokeID: 12, StartZero: 1348, PosPeak: 1380, MidZero: 1413, NegPeak: 1440, EndZero: 1466, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.225, Area_Positive: 2.88, Area_Negative: 2.38 },
  { StrokeID: 13, StartZero: 1466, PosPeak: 1498, MidZero: 1530, NegPeak: 1558, EndZero: 1584, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.233, Area_Positive: 2.95, Area_Negative: 2.44 },
  { StrokeID: 14, StartZero: 1584, PosPeak: 1618, MidZero: 1650, NegPeak: 1676, EndZero: 1702, TotalTime: 0.983, TimeToPosPeak: 0.283, TimeToNegPeak: 0.217, Area_Positive: 2.82, Area_Negative: 2.33 },
  { StrokeID: 15, StartZero: 1702, PosPeak: 1734, MidZero: 1766, NegPeak: 1793, EndZero: 1820, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.225, Area_Positive: 2.74, Area_Negative: 2.20 },
  { StrokeID: 16, StartZero: 1820, PosPeak: 1853, MidZero: 1885, NegPeak: 1912, EndZero: 1938, TotalTime: 0.983, TimeToPosPeak: 0.275, TimeToNegPeak: 0.225, Area_Positive: 2.61, Area_Negative: 2.10 },
  { StrokeID: 17, StartZero: 1938, PosPeak: 1970, MidZero: 2002, NegPeak: 2028, EndZero: 2056, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.217, Area_Positive: 2.55, Area_Negative: 2.05 },
  { StrokeID: 18, StartZero: 2056, PosPeak: 2088, MidZero: 2120, NegPeak: 2148, EndZero: 2174, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.233, Area_Positive: 2.49, Area_Negative: 1.98 },
  { StrokeID: 19, StartZero: 2174, PosPeak: 2206, MidZero: 2238, NegPeak: 2264, EndZero: 2292, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.217, Area_Positive: 2.42, Area_Negative: 1.92 },
  { StrokeID: 20, StartZero: 2292, PosPeak: 2324, MidZero: 2358, NegPeak: 2384, EndZero: 2410, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.217, Area_Positive: 2.38, Area_Negative: 1.88 },
  { StrokeID: 21, StartZero: 2410, PosPeak: 2442, MidZero: 2474, NegPeak: 2500, EndZero: 2528, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.217, Area_Positive: 2.90, Area_Negative: 2.40 },
  { StrokeID: 22, StartZero: 2528, PosPeak: 2560, MidZero: 2592, NegPeak: 2618, EndZero: 2644, TotalTime: 0.967, TimeToPosPeak: 0.267, TimeToNegPeak: 0.217, Area_Positive: 2.96, Area_Negative: 2.45 },
  { StrokeID: 23, StartZero: 2644, PosPeak: 2676, MidZero: 2708, NegPeak: 2735, EndZero: 2762, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.225, Area_Positive: 3.01, Area_Negative: 2.50 },
  { StrokeID: 24, StartZero: 2762, PosPeak: 2794, MidZero: 2826, NegPeak: 2854, EndZero: 2880, TotalTime: 0.983, TimeToPosPeak: 0.267, TimeToNegPeak: 0.233, Area_Positive: 2.87, Area_Negative: 2.36 },
]

// Mock acceleration phase metrics
export const mockAccMetrics: AccPhaseMetrics[] = mockStrokes.map(s => ({
  StrokeID: s.StrokeID,
  Duration: (s.MidZero - s.StartZero) / 120,
  AccY_PosPeak: 1.5 + Math.random() * 0.5,
  TimeToPosPeakPercent: 45 + Math.random() * 10,
  RAD: 3.5 + Math.random() * 1.5,
  VelocityChange: 0.15 + Math.random() * 0.08,
  MeanAbsAccX: 0.08 + Math.random() * 0.04,
  MeanAbsGyrX: 15 + Math.random() * 10,
  MeanAbsGyrY: 8 + Math.random() * 6,
  MeanAbsGyrZ: 12 + Math.random() * 8,
}))

// Mock deceleration phase metrics
export const mockDecMetrics: DecPhaseMetrics[] = mockStrokes.map(s => ({
  StrokeID: s.StrokeID,
  Duration: (s.EndZero - s.MidZero) / 120,
  AccY_NegPeak: -1.2 - Math.random() * 0.4,
  TimeToNegPeakPercent: 40 + Math.random() * 10,
  VelocityChange: -0.12 - Math.random() * 0.06,
  MeanAbsAccX: 0.06 + Math.random() * 0.03,
  MeanAbsGyrX: 12 + Math.random() * 8,
  MeanAbsGyrY: 6 + Math.random() * 5,
  MeanAbsGyrZ: 10 + Math.random() * 7,
}))

// Mock stroke phase metrics
export const mockStrokePhaseMetrics: StrokePhaseMetrics[] = mockStrokes.map(s => ({
  StrokeID: s.StrokeID,
  GyrX_PosPeak: 45 + Math.random() * 20,
  GyrX_NegPeak: -35 - Math.random() * 15,
  MeanAbsAccX: 0.07 + Math.random() * 0.035,
  MeanAbsGyrX: 13 + Math.random() * 9,
  MeanAbsGyrY: 7 + Math.random() * 5.5,
  MeanAbsGyrZ: 11 + Math.random() * 7.5,
}))

// Generate synthetic multi-axis waveform data for the overview chart
export function generateWaveformData(): IMUDataPoint[] {
  const data: IMUDataPoint[] = []
  const fs = 120
  const totalSamples = 2900

  for (let i = 0; i < totalSamples; i++) {
    const t = i / fs

    // Acc_Y: Main paddle stroke signal
    let accY = 1.8 * Math.sin(2 * Math.PI * 1.02 * t + 0.3)
      + 0.6 * Math.sin(2 * Math.PI * 2.04 * t + 1.1)
      + 0.3 * Math.sin(2 * Math.PI * 3.06 * t + 0.5)

    // Amplitude envelope
    const envelope = Math.min(1, i / (fs * 2)) * Math.min(1, (totalSamples - i) / (fs * 3))
    accY *= envelope
    accY += (Math.random() - 0.5) * 0.08

    // Acc_X: Forward/backward acceleration (smaller, phase-shifted)
    let accX = 0.6 * Math.sin(2 * Math.PI * 1.02 * t + 1.5)
      + 0.2 * Math.sin(2 * Math.PI * 2.04 * t + 2.2)
    accX *= envelope
    accX += (Math.random() - 0.5) * 0.04

    // Gyr_X: Roll rotation
    let gyrX = 35 * Math.sin(2 * Math.PI * 1.02 * t + 0.8)
      + 12 * Math.sin(2 * Math.PI * 2.04 * t + 1.6)
    gyrX *= envelope
    gyrX += (Math.random() - 0.5) * 3

    // Gyr_Y: Pitch rotation
    let gyrY = 25 * Math.sin(2 * Math.PI * 1.02 * t + 1.2)
      + 8 * Math.sin(2 * Math.PI * 2.04 * t + 2.0)
    gyrY *= envelope
    gyrY += (Math.random() - 0.5) * 2

    // Gyr_Z: Yaw rotation
    let gyrZ = 20 * Math.sin(2 * Math.PI * 1.02 * t + 0.5)
      + 7 * Math.sin(2 * Math.PI * 2.04 * t + 1.3)
    gyrZ *= envelope
    gyrZ += (Math.random() - 0.5) * 2.5

    data.push({
      sample: i,
      Acc_X: parseFloat(accX.toFixed(4)),
      Acc_Y: parseFloat(accY.toFixed(4)),
      Gyr_X: parseFloat(gyrX.toFixed(2)),
      Gyr_Y: parseFloat(gyrY.toFixed(2)),
      Gyr_Z: parseFloat(gyrZ.toFixed(2)),
    })
  }

  return data
}
