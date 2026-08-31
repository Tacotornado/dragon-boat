/**
 * Server-side normalizer for raw Xsens / Movella DOT exports.
 *
 * Movella writes 5-6 metadata rows (DeviceTag, FirmwareVersion, AppVersion,
 * StartTime, copyright) plus a blank row before the real column header. The C++
 * engine (stroke_analysis.cpp) reads line 1 as the header and matches columns by
 * the exact strings "Acc_X", "Acc_Y", "Gyr_X", "Gyr_Y", "Gyr_Z", so we hand it a
 * clean, canonical CSV instead of the raw file.
 */

export const ENGINE_HEADER = ["SampleTimeFine", "Acc_X", "Acc_Y", "Acc_Z", "Gyr_X", "Gyr_Y", "Gyr_Z"] as const

const ALIASES = {
  timestamp: ["sampletimefine", "sampletime", "timestamp", "times", "time", "packetcounter", "samplecounter", "sample"],
  accX: ["accx", "accelerationx", "linearaccelerationx", "freeaccx"],
  accY: ["accy", "accelerationy", "linearaccelerationy", "freeaccy"],
  accZ: ["accz", "accelerationz", "linearaccelerationz", "freeaccz"],
  gyrX: ["gyrx", "gyrox", "gyroscopex", "angularvelocityx", "ratex"],
  gyrY: ["gyry", "gyroy", "gyroscopey", "angularvelocityy", "ratey"],
  gyrZ: ["gyrz", "gyroz", "gyroscopez", "angularvelocityz", "ratez"],
} as const

const AXIS_ALIASES = [
  ...ALIASES.accX, ...ALIASES.accY, ...ALIASES.accZ,
  ...ALIASES.gyrX, ...ALIASES.gyrY, ...ALIASES.gyrZ,
]

export type AxisIndices = { [K in keyof typeof ALIASES]: number }

export interface NormalizedIMU {
  /** Index of the real column header inside the raw file. */
  headerIndex: number
  /** Detected field delimiter (comma, tab or semicolon). */
  delimiter: string
  /** Key/value pairs harvested from the metadata rows above the header. */
  metadata: Map<string, string>
  /** Resolved column positions for each axis (-1 when absent). */
  axis: AxisIndices
  /** Parsed numeric rows, ordered by sample index. */
  rows: IMURow[]
  /** Canonical CSV text ready to feed to the C++ binary. */
  engineCsv: string
  /** Run date derived from StartTime metadata, when present. */
  metadataDate: string | null
}

export interface IMURow {
  sampleIndex: number
  timestampSeconds: number | null
  accX: number | null
  accY: number | null
  accZ: number | null
  gyrX: number | null
  gyrY: number | null
  gyrZ: number | null
}

export function normalizeHeading(value: string) {
  return value
    .replace(/^\uFEFF/, "")
    .replace(/^\s*["']|["']\s*$/g, "")
    .replace(/[^a-zA-Z0-9]/g, "")
    .toLowerCase()
}

function delimiterFor(line: string) {
  return [",", "\t", ";"].sort((a, b) => line.split(b).length - line.split(a).length)[0]
}

export function splitCsvLine(line: string, delimiter = ",") {
  const values: string[] = []
  let value = ""
  let quoted = false
  for (const char of line) {
    if (char === '"') quoted = !quoted
    else if (char === delimiter && !quoted) {
      values.push(value.trim())
      value = ""
    } else value += char
  }
  values.push(value.trim())
  return values.map((v) => v.replace(/^"|"$/g, ""))
}

function toNumber(value: string | undefined) {
  if (value === undefined || value === "") return null
  const parsed = Number(value)
  return Number.isFinite(parsed) ? parsed : null
}

/**
 * Strips metadata, resolves aliased axis columns and emits a canonical CSV.
 * Throws when no recognizable Acc/Gyr header exists anywhere in the file.
 */
export function normalizeIMUCsv(rawCsv: string): NormalizedIMU {
  const lines = rawCsv
    .replace(/^\uFEFF/, "")
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter(Boolean)

  if (lines.length < 2) throw new Error("The CSV contains no data rows.")

  const headerIndex = lines.findIndex((line) =>
    [",", "\t", ";"].some((delimiter) => {
      const headings = line.split(delimiter).map(normalizeHeading)
      return headings.filter((heading) => AXIS_ALIASES.includes(heading)).length >= 2
    }),
  )

  if (headerIndex < 0) {
    throw new Error("Could not find an IMU header. Expected Acc_X/Acc_Y/Acc_Z or Gyr_X/Gyr_Y/Gyr_Z columns.")
  }

  const delimiter = delimiterFor(lines[headerIndex])
  const headings = lines[headerIndex].split(delimiter).map(normalizeHeading)
  // Alias order is priority order: SampleTimeFine wins over PacketCounter.
  const locate = (names: readonly string[]) => {
    for (const name of names) {
      const found = headings.indexOf(name)
      if (found >= 0) return found
    }
    return -1
  }

  const axis: AxisIndices = {
    timestamp: locate(ALIASES.timestamp),
    accX: locate(ALIASES.accX),
    accY: locate(ALIASES.accY),
    accZ: locate(ALIASES.accZ),
    gyrX: locate(ALIASES.gyrX),
    gyrY: locate(ALIASES.gyrY),
    gyrZ: locate(ALIASES.gyrZ),
  }

  if (axis.accY < 0) {
    throw new Error("The engine needs an Acc_Y column to detect stroke cycles, but none was found.")
  }

  const metadata = new Map<string, string>()
  for (const line of lines.slice(0, headerIndex)) {
    const cells = splitCsvLine(line, delimiter)
    const key = normalizeHeading(cells[0] ?? "")
    const inlineValue = cells.slice(1).filter(Boolean).join(" ").trim()
    if (key && inlineValue) {
      metadata.set(key.replace(/:$/, ""), inlineValue)
      continue
    }
    const separator = line.indexOf(":")
    if (separator > 0) metadata.set(normalizeHeading(line.slice(0, separator)), line.slice(separator + 1).trim())
  }

  const startTime = metadata.get("starttime") || metadata.get("rundate") || ""
  const dateMatch = startTime.match(/(20\d{2})[-_/](\d{1,2})[-_/](\d{1,2})/)
  const metadataDate = dateMatch
    ? `${dateMatch[1]}-${dateMatch[2].padStart(2, "0")}-${dateMatch[3].padStart(2, "0")}`
    : null

  const rows: IMURow[] = []
  for (const line of lines.slice(headerIndex + 1)) {
    const cells = splitCsvLine(line, delimiter)
    const row: IMURow = {
      sampleIndex: rows.length,
      timestampSeconds: toNumber(cells[axis.timestamp]),
      accX: toNumber(cells[axis.accX]),
      accY: toNumber(cells[axis.accY]),
      accZ: toNumber(cells[axis.accZ]),
      gyrX: toNumber(cells[axis.gyrX]),
      gyrY: toNumber(cells[axis.gyrY]),
      gyrZ: toNumber(cells[axis.gyrZ]),
    }
    const hasSignal = [row.accX, row.accY, row.accZ, row.gyrX, row.gyrY, row.gyrZ].some((v) => v !== null)
    if (hasSignal) rows.push({ ...row, sampleIndex: rows.length })
  }

  if (!rows.length) throw new Error("No numeric IMU samples were found below the header row.")

  const engineCsv = [
    ENGINE_HEADER.join(","),
    ...rows.map((row) =>
      [
        row.timestampSeconds ?? row.sampleIndex,
        row.accX ?? 0,
        row.accY ?? 0,
        row.accZ ?? 0,
        row.gyrX ?? 0,
        row.gyrY ?? 0,
        row.gyrZ ?? 0,
      ].join(","),
    ),
  ].join("\n")

  return { headerIndex, delimiter, metadata, axis, rows, engineCsv, metadataDate }
}
