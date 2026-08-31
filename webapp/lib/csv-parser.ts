import type { 
  StrokeCycle, 
  AccPhaseMetrics, 
  DecPhaseMetrics, 
  StrokePhaseMetrics, 
  IMUDataPoint,
  SummaryStats 
} from "./types"

/**
 * Parses the stroke_full_cycles.csv output from stroke_analysis.cpp
 */
export function parseStrokeCyclesCSV(csvText: string): StrokeCycle[] {
  const lines = csvText.trim().split(/\r?\n/)
  if (lines.length < 2) return []

  const header = lines[0].split(",").map((h) => h.trim())

  const colIndex = (name: string): number => {
    const idx = header.indexOf(name)
    if (idx === -1) throw new Error(`Missing column "${name}" in CSV header. Found: ${header.join(", ")}`)
    return idx
  }

  const idxStrokeID = colIndex("StrokeID")
  const idxStartZero = colIndex("StartZero")
  const idxPosPeak = colIndex("PosPeak")
  const idxMidZero = colIndex("MidZero")
  const idxNegPeak = colIndex("NegPeak")
  const idxEndZero = colIndex("EndZero")
  const idxTotalTime = colIndex("TotalTime")
  const idxTimeToPosPeak = colIndex("TimeToPosPeak")
  const idxTimeToNegPeak = colIndex("TimeToNegPeak")
  const idxAreaPositive = colIndex("Area_Positive")
  const idxAreaNegative = colIndex("Area_Negative")

  const strokes: StrokeCycle[] = []

  for (let i = 1; i < lines.length; i++) {
    const row = lines[i].split(",").map((c) => c.trim())
    if (row.length < header.length) continue

    const stroke: StrokeCycle = {
      StrokeID: parseInt(row[idxStrokeID], 10),
      StartZero: parseInt(row[idxStartZero], 10),
      PosPeak: parseInt(row[idxPosPeak], 10),
      MidZero: parseInt(row[idxMidZero], 10),
      NegPeak: parseInt(row[idxNegPeak], 10),
      EndZero: parseInt(row[idxEndZero], 10),
      TotalTime: parseFloat(row[idxTotalTime]),
      TimeToPosPeak: parseFloat(row[idxTimeToPosPeak]),
      TimeToNegPeak: parseFloat(row[idxTimeToNegPeak]),
      Area_Positive: parseFloat(row[idxAreaPositive]),
      Area_Negative: parseFloat(row[idxAreaNegative]),
    }

    if (Object.values(stroke).some((v) => typeof v === "number" && isNaN(v))) {
      continue
    }

    strokes.push(stroke)
  }

  return strokes
}

/**
 * Parses acceleration phase metrics from CSV (exported from Excel sheet 1)
 * Header: StrokeID,Duration,ACC-Y+Peak,TimeToPeak(%),RAD,ΔVelocity,|ACC-X|,|Gyr-X|,|Gyr-Y|,|Gyr-Z|
 */
export function parseAccPhaseCSV(csvText: string): AccPhaseMetrics[] {
  const lines = csvText.trim().split(/\r?\n/)
  if (lines.length < 2) return []

  const header = lines[0].split(",").map((h) => h.trim())
  
  // Flexible column matching
  const findCol = (names: string[]): number => {
    for (const name of names) {
      const idx = header.findIndex(h => 
        h.toLowerCase().replace(/[^a-z0-9]/g, '').includes(name.toLowerCase().replace(/[^a-z0-9]/g, ''))
      )
      if (idx !== -1) return idx
    }
    return -1
  }

  const metrics: AccPhaseMetrics[] = []

  for (let i = 1; i < lines.length; i++) {
    const row = lines[i].split(",").map((c) => c.trim())
    if (row.length < 2) continue

    const m: AccPhaseMetrics = {
      StrokeID: parseInt(row[findCol(["StrokeID"])] || "0", 10),
      Duration: parseFloat(row[findCol(["Duration"])] || "0"),
      AccY_PosPeak: parseFloat(row[findCol(["ACC-Y+Peak", "AccYPosPeak", "PosPeak"])] || "0"),
      TimeToPosPeakPercent: parseFloat(row[findCol(["TimeToPeak", "TimeToPosPeak"])] || "0"),
      RAD: parseFloat(row[findCol(["RAD"])] || "0"),
      VelocityChange: parseFloat(row[findCol(["Velocity", "ΔVelocity"])] || "0"),
      MeanAbsAccX: parseFloat(row[findCol(["ACC-X", "AccX"])] || "0"),
      MeanAbsGyrX: parseFloat(row[findCol(["Gyr-X", "GyrX"])] || "0"),
      MeanAbsGyrY: parseFloat(row[findCol(["Gyr-Y", "GyrY"])] || "0"),
      MeanAbsGyrZ: parseFloat(row[findCol(["Gyr-Z", "GyrZ"])] || "0"),
    }

    if (m.StrokeID > 0) metrics.push(m)
  }

  return metrics
}

/**
 * Parses deceleration phase metrics from CSV (exported from Excel sheet 2)
 */
export function parseDecPhaseCSV(csvText: string): DecPhaseMetrics[] {
  const lines = csvText.trim().split(/\r?\n/)
  if (lines.length < 2) return []

  const header = lines[0].split(",").map((h) => h.trim())
  
  const findCol = (names: string[]): number => {
    for (const name of names) {
      const idx = header.findIndex(h => 
        h.toLowerCase().replace(/[^a-z0-9]/g, '').includes(name.toLowerCase().replace(/[^a-z0-9]/g, ''))
      )
      if (idx !== -1) return idx
    }
    return -1
  }

  const metrics: DecPhaseMetrics[] = []

  for (let i = 1; i < lines.length; i++) {
    const row = lines[i].split(",").map((c) => c.trim())
    if (row.length < 2) continue

    const m: DecPhaseMetrics = {
      StrokeID: parseInt(row[findCol(["StrokeID"])] || "0", 10),
      Duration: parseFloat(row[findCol(["Duration"])] || "0"),
      AccY_NegPeak: parseFloat(row[findCol(["ACC-Y-Peak", "AccYNegPeak", "NegPeak"])] || "0"),
      TimeToNegPeakPercent: parseFloat(row[findCol(["TimeToNegPeak", "TimeToNeg"])] || "0"),
      VelocityChange: parseFloat(row[findCol(["Velocity", "ΔVelocity"])] || "0"),
      MeanAbsAccX: parseFloat(row[findCol(["ACC-X", "AccX"])] || "0"),
      MeanAbsGyrX: parseFloat(row[findCol(["Gyr-X", "GyrX"])] || "0"),
      MeanAbsGyrY: parseFloat(row[findCol(["Gyr-Y", "GyrY"])] || "0"),
      MeanAbsGyrZ: parseFloat(row[findCol(["Gyr-Z", "GyrZ"])] || "0"),
    }

    if (m.StrokeID > 0) metrics.push(m)
  }

  return metrics
}

/**
 * Parses stroke phase metrics from CSV (exported from Excel sheet 3)
 */
export function parseStrokePhaseCSV(csvText: string): StrokePhaseMetrics[] {
  const lines = csvText.trim().split(/\r?\n/)
  if (lines.length < 2) return []

  const header = lines[0].split(",").map((h) => h.trim())
  
  const findCol = (names: string[]): number => {
    for (const name of names) {
      const idx = header.findIndex(h => 
        h.toLowerCase().replace(/[^a-z0-9]/g, '').includes(name.toLowerCase().replace(/[^a-z0-9]/g, ''))
      )
      if (idx !== -1) return idx
    }
    return -1
  }

  const metrics: StrokePhaseMetrics[] = []

  for (let i = 1; i < lines.length; i++) {
    const row = lines[i].split(",").map((c) => c.trim())
    if (row.length < 2) continue

    const m: StrokePhaseMetrics = {
      StrokeID: parseInt(row[findCol(["StrokeID"])] || "0", 10),
      GyrX_PosPeak: parseFloat(row[findCol(["GyrX+Peak", "GyrXPosPeak"])] || "0"),
      GyrX_NegPeak: parseFloat(row[findCol(["GyrX-Peak", "GyrXNegPeak"])] || "0"),
      MeanAbsAccX: parseFloat(row[findCol(["ACC-X", "AccX"])] || "0"),
      MeanAbsGyrX: parseFloat(row[findCol(["Gyr-X", "GyrX"])] || "0"),
      MeanAbsGyrY: parseFloat(row[findCol(["Gyr-Y", "GyrY"])] || "0"),
      MeanAbsGyrZ: parseFloat(row[findCol(["Gyr-Z", "GyrZ"])] || "0"),
    }

    if (m.StrokeID > 0) metrics.push(m)
  }

  return metrics
}

/**
 * Parses the raw IMU CSV (Canoe_xsens_dot.csv) with all axes
 * Returns multi-axis IMU data points
 */
export function parseRawIMUCSV(csvText: string): IMUDataPoint[] {
  const lines = csvText
    .replace(/^\uFEFF/, "")
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter(Boolean)
  if (lines.length < 2) return []

  // Movella DOT exports include metadata rows before the actual header. They
  // may also contain a UTF-8 BOM, quoted headings, or tab/semicolon delimiters.
  const cleanHeading = (value: string) => value
    .replace(/^\uFEFF/, "")
    .replace(/^\s*[\"']|[\"']\s*$/g, "")
    .replace(/[^a-zA-Z0-9]/g, "")
    .toLowerCase()

  const delimiterFor = (line: string) => {
    const candidates = [",", "\t", ";"]
    return candidates.sort((a, b) => line.split(b).length - line.split(a).length)[0]
  }

  const isAxis = (column: string, family: "acc" | "gyr") => {
    const normalized = cleanHeading(column)
    if (family === "acc") return /^(acc|acceleration|linearacceleration)(x|y|z)$/.test(normalized)
    return /^(gyr|gyro|gyroscope|angularvelocity)(x|y|z)$/.test(normalized)
  }

  const headerIndex = lines.findIndex((line) => {
    const delimiter = delimiterFor(line)
    const columns = line.split(delimiter).map(cleanHeading)
    const axisColumns = columns.filter((column) => isAxis(column, "acc") || isAxis(column, "gyr"))
    return axisColumns.length > 0
  })

  if (headerIndex === -1) {
    throw new Error('Could not find an IMU header. Expected Acc_X/Acc_Y/Acc_Z or Gyr_X/Gyr_Y/Gyr_Z columns.')
  }

  const delimiter = delimiterFor(lines[headerIndex])
  const header = lines[headerIndex].split(delimiter).map((h) => h.trim())
  const normalizedHeader = header.map(cleanHeading)
  const findColumn = (...names: string[]) => {
    const normalizedNames = names.map(cleanHeading)
    return normalizedHeader.findIndex((heading) => normalizedNames.includes(heading))
  }

  const packetIdx = findColumn("PacketCounter", "Sample", "SampleCounter")
  const accXIdx = findColumn("Acc_X")
  const accYIdx = findColumn("Acc_Y")
  const gyrXIdx = findColumn("Gyr_X")
  const gyrYIdx = findColumn("Gyr_Y")
  const gyrZIdx = findColumn("Gyr_Z")

  const data: IMUDataPoint[] = []
  for (let i = headerIndex + 1; i < lines.length; i++) {
    const row = lines[i].split(delimiter).map((cell) => cell.trim().replace(/^\s*["']|["']\s*$/g, ""))
    if (row.length < header.length || row.every((cell) => cell === "")) continue

    const value = (index: number) => index >= 0 ? Number.parseFloat(row[index]) || 0 : 0
    const parsedSample = packetIdx >= 0 ? Number.parseInt(row[packetIdx], 10) : Number.NaN

    data.push({
      sample: Number.isFinite(parsedSample) ? parsedSample : data.length,
      Acc_X: value(accXIdx),
      Acc_Y: value(accYIdx),
      Gyr_X: value(gyrXIdx),
      Gyr_Y: value(gyrYIdx),
      Gyr_Z: value(gyrZIdx),
    })
  }

  return data
}

/**
 * Legacy function for backward compatibility - returns only accY
 */
export function parseRawIMUCSVLegacy(csvText: string): { sample: number; accY: number }[] {
  const data = parseRawIMUCSV(csvText)
  return data.map(d => ({ sample: d.sample, accY: d.Acc_Y }))
}

/**
 * Calculate summary statistics from stroke data
 */
export function calculateSummaryStats(
  strokes: StrokeCycle[],
  accMetrics?: AccPhaseMetrics[],
  decMetrics?: DecPhaseMetrics[],
  strokePhaseMetrics?: StrokePhaseMetrics[]
): SummaryStats {
  const n = strokes.length
  if (n === 0) {
    return {
      totalStrokes: 0,
      avgSPM: 0,
      avgVelocityAcc: 0,
      avgVelocityDec: 0,
      avgMACAccX: 0,
      avgMACGyrX: 0,
      avgMACGyrY: 0,
      avgMACGyrZ: 0,
    }
  }

  // Calculate total time span for SPM
  const totalTime = strokes.reduce((sum, s) => sum + s.TotalTime, 0)
  const avgSPM = totalTime > 0 ? (n / totalTime) * 60 : 0

  // Velocity changes from phase metrics
  const avgVelocityAcc = accMetrics && accMetrics.length > 0
    ? accMetrics.reduce((sum, m) => sum + m.VelocityChange, 0) / accMetrics.length
    : 0

  const avgVelocityDec = decMetrics && decMetrics.length > 0
    ? decMetrics.reduce((sum, m) => sum + m.VelocityChange, 0) / decMetrics.length
    : 0

  // MAC values from stroke phase metrics
  const avgMACAccX = strokePhaseMetrics && strokePhaseMetrics.length > 0
    ? strokePhaseMetrics.reduce((sum, m) => sum + m.MeanAbsAccX, 0) / strokePhaseMetrics.length
    : 0

  const avgMACGyrX = strokePhaseMetrics && strokePhaseMetrics.length > 0
    ? strokePhaseMetrics.reduce((sum, m) => sum + m.MeanAbsGyrX, 0) / strokePhaseMetrics.length
    : 0

  const avgMACGyrY = strokePhaseMetrics && strokePhaseMetrics.length > 0
    ? strokePhaseMetrics.reduce((sum, m) => sum + m.MeanAbsGyrY, 0) / strokePhaseMetrics.length
    : 0

  const avgMACGyrZ = strokePhaseMetrics && strokePhaseMetrics.length > 0
    ? strokePhaseMetrics.reduce((sum, m) => sum + m.MeanAbsGyrZ, 0) / strokePhaseMetrics.length
    : 0

  return {
    totalStrokes: n,
    avgSPM,
    avgVelocityAcc,
    avgVelocityDec,
    avgMACAccX,
    avgMACGyrX,
    avgMACGyrY,
    avgMACGyrZ,
  }
}

/**
 * Validates that the uploaded stroke CSV looks reasonable
 */
export function validateStrokeData(strokes: StrokeCycle[]): { valid: boolean; message: string } {
  if (strokes.length === 0) {
    return { valid: false, message: "No stroke cycles found in the CSV file." }
  }

  const ids = strokes.map((s) => s.StrokeID)
  const isSequential = ids.every((id, i) => i === 0 || id > ids[i - 1])
  if (!isSequential) {
    return { valid: true, message: `Loaded ${strokes.length} strokes (non-sequential IDs).` }
  }

  return { valid: true, message: `Loaded ${strokes.length} stroke cycles successfully.` }
}
