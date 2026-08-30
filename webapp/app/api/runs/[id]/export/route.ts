import { NextResponse } from "next/server"
import { eq } from "drizzle-orm"
import { db } from "@/lib/db"
import { strokeRuns } from "@/lib/db/schema"

export const runtime = "nodejs"

function escapeCsv(value: unknown) {
  const text = String(value ?? "")
  return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text
}

export async function GET(_request: Request, { params }: { params: Promise<{ id: string }> }) {
  const { id } = await params
  const result = await db.select().from(strokeRuns).where(eq(strokeRuns.id, id)).limit(1)
  const run = result[0]
  if (!run) return NextResponse.json({ error: "Run not found" }, { status: 404 })
  const metadata = [
    ["Run", run.name], ["Racer", run.racerName], ["Uploader", run.uploaderName],
    ["Run date", run.runDate], ["Race date", run.raceDate], ["Uploaded at", run.uploadedAt.toISOString()],
    ["Note", run.note], ["Source file", run.sourceFilename],
  ].map((row) => row.map(escapeCsv).join(",")).join("\n")
  const body = `${metadata}\n\n${run.rawCsv}`
  return new NextResponse(body, {
    headers: {
      "Content-Type": "text/csv; charset=utf-8",
      "Content-Disposition": `attachment; filename="${run.name}.csv"`,
    },
  })
}
