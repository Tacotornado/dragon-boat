import { NextResponse } from "next/server"
import { asc, eq } from "drizzle-orm"
import { db } from "@/lib/db"
import { rm } from "fs/promises"
import { join } from "path"
import { strokeRuns, strokeImuSamples, strokeCycles, strokePhaseMetrics, strokeRunArtifacts } from "@/lib/db/schema"
import { ARTIFACT_ROOT } from "@/lib/cpp-engine"

export const runtime = "nodejs"

export async function GET(_request: Request, { params }: { params: Promise<{ id: string }> }) {
  const { id } = await params
  const run = await db.select().from(strokeRuns).where(eq(strokeRuns.id, id)).limit(1)
  if (!run[0]) return NextResponse.json({ error: "Run not found" }, { status: 404 })
  const [samples, strokes, metrics] = await Promise.all([
    db.select().from(strokeImuSamples).where(eq(strokeImuSamples.runId, id)).orderBy(asc(strokeImuSamples.sampleIndex)),
    db.select().from(strokeCycles).where(eq(strokeCycles.runId, id)).orderBy(asc(strokeCycles.strokeId)),
    db.select().from(strokePhaseMetrics).where(eq(strokePhaseMetrics.runId, id)).orderBy(asc(strokePhaseMetrics.strokeId)),
  ])
  return NextResponse.json({ run: run[0], samples, strokes, metrics })
}

export async function DELETE(_request: Request, { params }: { params: Promise<{ id: string }> }) {
  const { id } = await params
  await rm(join(ARTIFACT_ROOT, id), { recursive: true, force: true }).catch(() => {})
  await db.delete(strokeRunArtifacts).where(eq(strokeRunArtifacts.runId, id))
  await db.delete(strokePhaseMetrics).where(eq(strokePhaseMetrics.runId, id))
  await db.delete(strokeCycles).where(eq(strokeCycles.runId, id))
  await db.delete(strokeImuSamples).where(eq(strokeImuSamples.runId, id))
  const deleted = await db.delete(strokeRuns).where(eq(strokeRuns.id, id)).returning({ id: strokeRuns.id })
  if (!deleted[0]) return NextResponse.json({ error: "Run not found" }, { status: 404 })
  return NextResponse.json({ success: true })
}
