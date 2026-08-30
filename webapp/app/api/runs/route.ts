import { NextResponse } from "next/server"
import { and, desc, eq, sql } from "drizzle-orm"
import { db } from "@/lib/db"
import { strokeRuns, strokeImuSamples, strokeCycles, strokeRunArtifacts } from "@/lib/db/schema"
import { normalizeIMUCsv } from "@/lib/imu-csv"
import { ENGINE_ENABLED, runStrokeEngine } from "@/lib/cpp-engine"

export const runtime = "nodejs"
export const maxDuration = 300

const SAMPLE_RATE_HZ = Number(process.env.SAMPLE_RATE_HZ || 120)

function fallbackDate(filename: string) {
  const match = filename.match(/(20\d{2})[-_](\d{1,2})[-_](\d{1,2})/)
  if (!match) return new Date().toISOString().slice(0, 10)
  return `${match[1]}-${match[2].padStart(2, "0")}-${match[3].padStart(2, "0")}`
}

export async function GET() {
  const runs = await db
    .select({
      id: strokeRuns.id,
      name: strokeRuns.name,
      racerName: strokeRuns.racerName,
      uploaderName: strokeRuns.uploaderName,
      note: strokeRuns.note,
      runDate: strokeRuns.runDate,
      raceDate: strokeRuns.raceDate,
      uploadedAt: strokeRuns.uploadedAt,
      sourceFilename: strokeRuns.sourceFilename,
      sampleCount: strokeRuns.sampleCount,
      strokeCount: strokeRuns.strokeCount,
      durationSeconds: strokeRuns.durationSeconds,
    })
    .from(strokeRuns)
    .orderBy(desc(strokeRuns.uploadedAt))
  return NextResponse.json(runs)
}

export async function POST(request: Request) {
  try {
    const form = await request.formData()
    const file = form.get("file")
    const racerName = String(form.get("racerName") || "").trim()

    if (!(file instanceof File) || !file.name.toLowerCase().endsWith(".csv")) {
      return NextResponse.json({ error: "A CSV file is required." }, { status: 400 })
    }
    if (!racerName) return NextResponse.json({ error: "Racer name is required." }, { status: 400 })
    if (file.size > 50 * 1024 * 1024) {
      return NextResponse.json({ error: "File exceeds the 50 MB limit." }, { status: 413 })
    }

    const rawCsv = await file.text()

    let normalized
    try {
      normalized = normalizeIMUCsv(rawCsv)
    } catch (error) {
      return NextResponse.json(
        { error: error instanceof Error ? error.message : "Could not parse the CSV." },
        { status: 400 },
      )
    }

    const runDate = String(form.get("runDate") || normalized.metadataDate || fallbackDate(file.name))
    const prefix = racerName.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "") || "racer"
    const countResult = await db
      .select({ count: sql<number>`count(*)` })
      .from(strokeRuns)
      .where(and(eq(strokeRuns.racerName, racerName), eq(strokeRuns.runDate, runDate)))
    const sequence = Number(countResult[0]?.count || 0) + 1
    const id = crypto.randomUUID()
    const name = `${prefix}-${runDate}-${sequence}`

    const sampleRows = normalized.rows.map((row) => ({
      runId: id,
      sampleIndex: row.sampleIndex,
      timestampSeconds: row.timestampSeconds,
      accX: row.accX,
      accY: row.accY,
      gyrX: row.gyrX,
      gyrY: row.gyrY,
      gyrZ: row.gyrZ,
    }))

    await db.insert(strokeRuns).values({
      id,
      name,
      racerName,
      uploaderName: String(form.get("uploaderName") || "").trim() || null,
      note: String(form.get("note") || "").trim() || null,
      runDate,
      raceDate: String(form.get("raceDate") || "") || null,
      sourceFilename: file.name,
      sampleCount: sampleRows.length,
      durationSeconds: sampleRows.length / SAMPLE_RATE_HZ,
      rawCsv,
    })

    for (let offset = 0; offset < sampleRows.length; offset += 1000) {
      await db.insert(strokeImuSamples).values(sampleRows.slice(offset, offset + 1000))
    }

    // Hand the canonical CSV to the C++ engine and persist what it produces.
    const engine = await runStrokeEngine(id, normalized.engineCsv)

    if (engine.cycles.length) {
      const cycleRows = engine.cycles.map((cycle) => ({
        runId: id,
        strokeId: cycle.strokeId,
        startZero: cycle.startZero,
        posPeak: cycle.posPeak,
        midZero: cycle.midZero,
        negPeak: cycle.negPeak,
        endZero: cycle.endZero,
        totalTime: cycle.totalTime,
        timeToPosPeak: cycle.timeToPosPeak,
        timeToNegPeak: cycle.timeToNegPeak,
        areaPositive: cycle.areaPositive,
        areaNegative: cycle.areaNegative,
      }))
      for (let offset = 0; offset < cycleRows.length; offset += 500) {
        await db.insert(strokeCycles).values(cycleRows.slice(offset, offset + 500))
      }
      await db.update(strokeRuns).set({ strokeCount: cycleRows.length }).where(eq(strokeRuns.id, id))
    }

    if (engine.artifacts.length) {
      await db.insert(strokeRunArtifacts).values(
        engine.artifacts.map((artifact) => ({
          runId: id,
          kind: artifact.kind,
          filename: artifact.filename,
          relativePath: artifact.relativePath,
          byteSize: artifact.byteSize,
          strokeId: artifact.strokeId,
        })),
      )
    }

    if (engine.error) console.error("[v0] stroke engine:", engine.error, engine.stderr)

    return NextResponse.json(
      {
        id,
        name,
        sampleCount: sampleRows.length,
        strokeCount: engine.cycles.length,
        artifactCount: engine.artifacts.length,
        engine: {
          enabled: ENGINE_ENABLED,
          ran: engine.ran,
          exitCode: engine.exitCode,
          error: engine.error ?? null,
        },
      },
      { status: 201 },
    )
  } catch (error) {
    console.error("[v0] Run upload failed", error)
    return NextResponse.json(
      { error: error instanceof Error ? error.message : "Run upload failed." },
      { status: 500 },
    )
  }
}
