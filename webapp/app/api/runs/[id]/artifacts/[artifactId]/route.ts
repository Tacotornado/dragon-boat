import { NextResponse } from "next/server"
import { and, eq } from "drizzle-orm"
import { readFile } from "fs/promises"
import { join, normalize, resolve } from "path"
import { db } from "@/lib/db"
import { strokeRunArtifacts } from "@/lib/db/schema"
import { ARTIFACT_ROOT } from "@/lib/cpp-engine"

export const runtime = "nodejs"

const CONTENT_TYPES: Record<string, string> = {
  ".csv": "text/csv; charset=utf-8",
  ".png": "image/png",
  ".xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
}

export async function GET(
  request: Request,
  { params }: { params: Promise<{ id: string; artifactId: string }> },
) {
  const { id, artifactId } = await params
  const numericId = Number(artifactId)
  if (!Number.isInteger(numericId)) return NextResponse.json({ error: "Invalid artifact id" }, { status: 400 })

  const rows = await db
    .select()
    .from(strokeRunArtifacts)
    .where(and(eq(strokeRunArtifacts.id, numericId), eq(strokeRunArtifacts.runId, id)))
    .limit(1)
  const artifact = rows[0]
  if (!artifact) return NextResponse.json({ error: "Artifact not found" }, { status: 404 })

  // Guard against path traversal: the resolved file must stay inside ARTIFACT_ROOT.
  const root = resolve(ARTIFACT_ROOT)
  const target = resolve(join(root, normalize(artifact.relativePath)))
  if (!target.startsWith(root + "/")) return NextResponse.json({ error: "Invalid artifact path" }, { status: 400 })

  let file: Buffer
  try {
    file = await readFile(target)
  } catch {
    return NextResponse.json({ error: "Artifact file is no longer on disk" }, { status: 410 })
  }

  const extension = artifact.filename.slice(artifact.filename.lastIndexOf(".")).toLowerCase()
  const inline = new URL(request.url).searchParams.get("inline") === "1"
  const disposition = `${inline ? "inline" : "attachment"}; filename="${artifact.filename}"`

  return new NextResponse(new Uint8Array(file), {
    headers: {
      "Content-Type": CONTENT_TYPES[extension] || "application/octet-stream",
      "Content-Length": String(file.byteLength),
      "Content-Disposition": disposition,
      "Cache-Control": "private, max-age=3600",
      "X-Content-Type-Options": "nosniff",
    },
  })
}
