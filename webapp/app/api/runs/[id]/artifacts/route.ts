import { NextResponse } from "next/server"
import { asc, eq } from "drizzle-orm"
import { db } from "@/lib/db"
import { strokeRunArtifacts } from "@/lib/db/schema"

export const runtime = "nodejs"

export async function GET(_request: Request, { params }: { params: Promise<{ id: string }> }) {
  const { id } = await params
  const artifacts = await db
    .select()
    .from(strokeRunArtifacts)
    .where(eq(strokeRunArtifacts.runId, id))
    .orderBy(asc(strokeRunArtifacts.strokeId), asc(strokeRunArtifacts.id))
  return NextResponse.json(artifacts)
}
