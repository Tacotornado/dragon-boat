export interface StrokeCycle {
  StrokeID: number
  StartZero: number
  PosPeak: number
  MidZero: number
  NegPeak: number
  EndZero: number
  TotalTime: number
  TimeToPosPeak: number
  TimeToNegPeak: number
  Area_Positive: number
  Area_Negative: number
}

export interface AccPhaseMetrics {
  StrokeID: number
  Duration: number
  AccY_PosPeak: number
  TimeToPosPeakPercent: number
  RAD: number
  VelocityChange: number
  MeanAbsAccX: number
  MeanAbsGyrX: number
  MeanAbsGyrY: number
  MeanAbsGyrZ: number
}

export interface DecPhaseMetrics {
  StrokeID: number
  Duration: number
  AccY_NegPeak: number
  TimeToNegPeakPercent: number
  VelocityChange: number
  MeanAbsAccX: number
  MeanAbsGyrX: number
  MeanAbsGyrY: number
  MeanAbsGyrZ: number
}

export interface StrokePhaseMetrics {
  StrokeID: number
  GyrX_PosPeak: number
  GyrX_NegPeak: number
  MeanAbsAccX: number
  MeanAbsGyrX: number
  MeanAbsGyrY: number
  MeanAbsGyrZ: number
}

// Multi-axis IMU data point
export interface IMUDataPoint {
  sample: number
  Acc_X: number
  Acc_Y: number
  Gyr_X: number
  Gyr_Y: number
  Gyr_Z: number
}

// Available axis types for visualization
export type AxisType = "Acc_X" | "Acc_Y" | "Gyr_X" | "Gyr_Y" | "Gyr_Z"

// Combined stroke data with all phase metrics
export interface FullStrokeData {
  cycle: StrokeCycle
  accPhase?: AccPhaseMetrics
  decPhase?: DecPhaseMetrics
  strokePhase?: StrokePhaseMetrics
}

// Summary statistics for dashboard
export interface SummaryStats {
  totalStrokes: number
  avgSPM: number // Strokes per minute
  avgVelocityAcc: number // Avg delta velocity in acceleration phase (m/s)
  avgVelocityDec: number // Avg delta velocity in deceleration phase (m/s)
  avgMACAccX: number // Mean Absolute Change in Acc-X (m/s²)
  avgMACGyrX: number // Mean Absolute Change in Gyr-X (deg/s)
  avgMACGyrY: number // Mean Absolute Change in Gyr-Y (deg/s)
  avgMACGyrZ: number // Mean Absolute Change in Gyr-Z (deg/s)
}

// Parameter definitions for configurable display
export interface ParameterDef {
  key: string
  label: string
  unit: string
  category: "acceleration" | "deceleration" | "stroke"
  format?: (v: number) => string
}

export const AVAILABLE_PARAMETERS: ParameterDef[] = [
  // Acceleration Phase
  { key: "acc_Duration", label: "Duration", unit: "s", category: "acceleration" },
  { key: "acc_AccY_PosPeak", label: "Acc-Y+ Peak", unit: "m/s²", category: "acceleration" },
  { key: "acc_TimeToPosPeakPercent", label: "Time to Peak", unit: "%", category: "acceleration" },
  { key: "acc_RAD", label: "RAD", unit: "m/s³", category: "acceleration" },
  { key: "acc_VelocityChange", label: "ΔVelocity", unit: "m/s", category: "acceleration" },
  { key: "acc_MeanAbsAccX", label: "|Acc-X|", unit: "m/s²", category: "acceleration" },
  { key: "acc_MeanAbsGyrX", label: "|Gyr-X|", unit: "deg/s", category: "acceleration" },
  { key: "acc_MeanAbsGyrY", label: "|Gyr-Y|", unit: "deg/s", category: "acceleration" },
  { key: "acc_MeanAbsGyrZ", label: "|Gyr-Z|", unit: "deg/s", category: "acceleration" },
  // Deceleration Phase
  { key: "dec_Duration", label: "Duration", unit: "s", category: "deceleration" },
  { key: "dec_AccY_NegPeak", label: "Acc-Y- Peak", unit: "m/s²", category: "deceleration" },
  { key: "dec_TimeToNegPeakPercent", label: "Time to Neg Peak", unit: "%", category: "deceleration" },
  { key: "dec_VelocityChange", label: "ΔVelocity", unit: "m/s", category: "deceleration" },
  { key: "dec_MeanAbsAccX", label: "|Acc-X|", unit: "m/s²", category: "deceleration" },
  { key: "dec_MeanAbsGyrX", label: "|Gyr-X|", unit: "deg/s", category: "deceleration" },
  { key: "dec_MeanAbsGyrY", label: "|Gyr-Y|", unit: "deg/s", category: "deceleration" },
  { key: "dec_MeanAbsGyrZ", label: "|Gyr-Z|", unit: "deg/s", category: "deceleration" },
  // Stroke Phase
  { key: "stroke_GyrX_PosPeak", label: "Gyr-X+ Peak", unit: "deg/s", category: "stroke" },
  { key: "stroke_GyrX_NegPeak", label: "Gyr-X- Peak", unit: "deg/s", category: "stroke" },
  { key: "stroke_MeanAbsAccX", label: "|Acc-X|", unit: "m/s²", category: "stroke" },
  { key: "stroke_MeanAbsGyrX", label: "|Gyr-X|", unit: "deg/s", category: "stroke" },
  { key: "stroke_MeanAbsGyrY", label: "|Gyr-Y|", unit: "deg/s", category: "stroke" },
  { key: "stroke_MeanAbsGyrZ", label: "|Gyr-Z|", unit: "deg/s", category: "stroke" },
  // Basic cycle metrics
  { key: "TotalTime", label: "Total Time", unit: "s", category: "stroke" },
  { key: "TimeToPosPeak", label: "T+ Peak", unit: "s", category: "stroke" },
  { key: "TimeToNegPeak", label: "T- Peak", unit: "s", category: "stroke" },
  { key: "Area_Positive", label: "Area+", unit: "", category: "stroke" },
  { key: "Area_Negative", label: "Area-", unit: "", category: "stroke" },
]

// ============================================
// SESSION & REAL-TIME TYPES
// ============================================

export type SessionStatus = "idle" | "recording" | "processing" | "completed" | "error"

export interface Session {
  id: string
  name: string
  athleteName?: string
  createdAt: string
  updatedAt: string
  status: SessionStatus
  duration: number // seconds
  sampleCount: number
  strokeCount: number
  notes?: string
}

export interface SessionWithData extends Session {
  rawData: IMUDataPoint[]
  strokes: StrokeCycle[]
  accMetrics: AccPhaseMetrics[]
  decMetrics: DecPhaseMetrics[]
  strokePhaseMetrics: StrokePhaseMetrics[]
}

// Real-time message types for WebSocket
export type WSMessageType = 
  | "session_start"
  | "session_stop"
  | "imu_data"
  | "stroke_detected"
  | "metrics_update"
  | "error"

export interface WSMessage {
  type: WSMessageType
  sessionId: string
  timestamp: number
  payload: unknown
}

export interface WSIMUData {
  samples: IMUDataPoint[]
}

export interface WSStrokeDetected {
  stroke: StrokeCycle
  accPhase?: AccPhaseMetrics
  decPhase?: DecPhaseMetrics
  strokePhase?: StrokePhaseMetrics
}

export interface WSMetricsUpdate {
  stats: SummaryStats
  latestStrokes: StrokeCycle[]
}

// Xsens API configuration
export interface XsensConfig {
  deviceId: string
  sampleRate: number // Hz (e.g., 120)
  apiEndpoint?: string
}
