import { NextResponse } from "next/server"
import { getSession, deleteSession, updateSessionStatus, replaceSessionAnalysis } from "@/lib/session-store"
import { runStrokeEngine } from "@/lib/cpp-engine"

// GET /api/sessions/[id] - Get session details with data
export async function GET(
  request: Request,
  { params }: { params: Promise<{ id: string }> }
) {
  const { id } = await params
  const session = await getSession(id)
  
  if (!session) {
    return NextResponse.json(
      { error: "Session not found" },
      { status: 404 }
    )
  }
  
  return NextResponse.json(session)
}

// PATCH /api/sessions/[id] - Update session status
export async function PATCH(
  request: Request,
  { params }: { params: Promise<{ id: string }> }
) {
  const { id } = await params
  
  try {
    const body = await request.json()
    const { status } = body
    
    if (!status) {
      return NextResponse.json(
        { error: "Status is required" },
        { status: 400 }
      )
    }
    
    const current = await getSession(id)
    if (!current) {
      return NextResponse.json(
        { error: "Session not found" },
        { status: 404 }
      )
    }

    let engine: { ran: boolean; cycles: unknown[]; artifacts: unknown[]; error?: string } | null = null
    if (status === "completed" && current.rawData.length > 0) {
      const engineCsv = [
        "SampleTimeFine,Acc_X,Acc_Y,Acc_Z,Gyr_X,Gyr_Y,Gyr_Z",
        ...current.rawData.map((sample) => [
          sample.sample,
          sample.Acc_X,
          sample.Acc_Y,
          sample.Acc_Z ?? 0,
          sample.Gyr_X,
          sample.Gyr_Y,
          sample.Gyr_Z,
        ].join(",")),
      ].join("\\n")
      const result = await runStrokeEngine(id, engineCsv)
      engine = { ran: result.ran, cycles: result.cycles, artifacts: result.artifacts, error: result.error }
      if (result.cycles.length > 0) await replaceSessionAnalysis(id, result.cycles.map((cycle) => ({
        StrokeID: cycle.strokeId,
        StartZero: cycle.startZero,
        PosPeak: cycle.posPeak,
        MidZero: cycle.midZero,
        NegPeak: cycle.negPeak,
        EndZero: cycle.endZero,
        TotalTime: cycle.totalTime,
        TimeToPosPeak: cycle.timeToPosPeak,
        TimeToNegPeak: cycle.timeToNegPeak,
        Area_Positive: cycle.areaPositive,
        Area_Negative: cycle.areaNegative,
      })))
    }

    const session = await updateSessionStatus(id, status)
    return NextResponse.json({ session, engine })
  } catch {
    return NextResponse.json(
      { error: "Failed to update session" },
      { status: 500 }
    )
  }
}

// DELETE /api/sessions/[id] - Delete a session
export async function DELETE(
  request: Request,
  { params }: { params: Promise<{ id: string }> }
) {
  const { id } = await params
  const deleted = await deleteSession(id)
  
  if (!deleted) {
    return NextResponse.json(
      { error: "Session not found" },
      { status: 404 }
    )
  }
  
  return NextResponse.json({ success: true })
}
