import { pgTable, text, integer, doublePrecision, timestamp, date, bigint } from "drizzle-orm/pg-core"

export const strokeRuns = pgTable("stroke_runs", {
  id: text("id").primaryKey(),
  name: text("name").notNull(),
  racerName: text("racer_name").notNull(),
  uploaderName: text("uploader_name"),
  note: text("note"),
  runDate: date("run_date").notNull(),
  raceDate: date("race_date"),
  uploadedAt: timestamp("uploaded_at", { withTimezone: true }).notNull().defaultNow(),
  sourceFilename: text("source_filename").notNull(),
  sampleCount: integer("sample_count").notNull().default(0),
  strokeCount: integer("stroke_count").notNull().default(0),
  durationSeconds: doublePrecision("duration_seconds").notNull().default(0),
  rawCsv: text("raw_csv").notNull(),
  createdAt: timestamp("created_at", { withTimezone: true }).notNull().defaultNow(),
})

export const strokeImuSamples = pgTable("stroke_imu_samples", {
  id: bigint("id", { mode: "number" }).generatedAlwaysAsIdentity().primaryKey(),
  runId: text("run_id").notNull(),
  sampleIndex: integer("sample_index").notNull(),
  timestampSeconds: doublePrecision("timestamp_seconds"),
  accX: doublePrecision("acc_x"),
  accY: doublePrecision("acc_y"),
  gyrX: doublePrecision("gyr_x"),
  gyrY: doublePrecision("gyr_y"),
  gyrZ: doublePrecision("gyr_z"),
})

export const strokeCycles = pgTable("stroke_cycles", {
  id: bigint("id", { mode: "number" }).generatedAlwaysAsIdentity().primaryKey(),
  runId: text("run_id").notNull(),
  strokeId: integer("stroke_id").notNull(),
  startZero: integer("start_zero").notNull(),
  posPeak: integer("pos_peak").notNull(),
  midZero: integer("mid_zero").notNull(),
  negPeak: integer("neg_peak").notNull(),
  endZero: integer("end_zero").notNull(),
  totalTime: doublePrecision("total_time").notNull(),
  timeToPosPeak: doublePrecision("time_to_pos_peak").notNull(),
  timeToNegPeak: doublePrecision("time_to_neg_peak").notNull(),
  areaPositive: doublePrecision("area_positive").notNull(),
  areaNegative: doublePrecision("area_negative").notNull(),
})

export const strokePhaseMetrics = pgTable("stroke_phase_metrics", {
  id: bigint("id", { mode: "number" }).generatedAlwaysAsIdentity().primaryKey(),
  runId: text("run_id").notNull(),
  strokeId: integer("stroke_id").notNull(),
  phase: text("phase").notNull(),
  metricName: text("metric_name").notNull(),
  metricValue: doublePrecision("metric_value").notNull(),
  unit: text("unit"),
})

export const strokeRunArtifacts = pgTable("stroke_run_artifacts", {
  id: bigint("id", { mode: "number" }).generatedAlwaysAsIdentity().primaryKey(),
  runId: text("run_id").notNull(),
  kind: text("kind").notNull(),
  filename: text("filename").notNull(),
  relativePath: text("relative_path").notNull(),
  byteSize: integer("byte_size").notNull().default(0),
  strokeId: integer("stroke_id"),
  createdAt: timestamp("created_at", { withTimezone: true }).notNull().defaultNow(),
})

export type StrokeRunRow = typeof strokeRuns.$inferSelect
export type StrokeImuSampleRow = typeof strokeImuSamples.$inferSelect
export type StrokeCycleRow = typeof strokeCycles.$inferSelect
export type StrokePhaseMetricRow = typeof strokePhaseMetrics.$inferSelect
export type StrokeRunArtifactRow = typeof strokeRunArtifacts.$inferSelect
