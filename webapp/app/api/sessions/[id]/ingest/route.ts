import { NextResponse } from "next/server"
import { appendIMUData, getSession, updateSessionStatus } from "@/lib/session-store"
import { IMUDataPoint } from "@/lib/types"
import { processIMUChunk } from "@/lib/stroke-processor"
import { broadcastIMUData, broadcastStrokeDetected } from "../stream/route"

// POST /api/sessions/[id]/ingest - Receive IMU data from Xsens
// This endpoint receives batches of IMU samples and triggers C++ processing
export async function POST(
  request: Request,
  { params }: { params: Promise<{ id: string }> }
) {
  const { id } = await params
  
  try {
    const session = await getSession(id)
    
    if (!session) {
      return NextResponse.json(
        { error: "Session not found" },
        { status: 404 }
      )
    }
    
    // Auto-start recording if idle
    if (session.status === "idle") {
      await updateSessionStatus(id, "recording")
    }
    
    const body = await request.json()
    const samples: IMUDataPoint[] = body.samples
    
    if (!samples || !Array.isArray(samples)) {
      return NextResponse.json(
        { error: "Invalid samples data" },
        { status: 400 }
      )
    }
    
    // Store raw IMU data
    const stored = await appendIMUData(id, samples)
    
    if (!stored) {
      return NextResponse.json(
        { error: "Failed to store data" },
        { status: 500 }
      )
    }
    
    // Broadcast raw IMU data to live viewers
    broadcastIMUData(id, samples)
    
    // Process the chunk through C++ engine (returns any newly detected strokes)
    const result = await processIMUChunk(id, samples)
    
    // Broadcast stroke detection events
    if (result.newStrokes.length > 0) {
      await broadcastStrokeDetected(id, result.newStrokes, result.totalStrokes)
    }
    
    return NextResponse.json({
      success: true,
      samplesReceived: samples.length,
      newStrokes: result.newStrokes,
      totalStrokes: result.totalStrokes,
    })
  } catch (error) {
    console.error("Ingest error:", error)
    return NextResponse.json(
      { error: "Failed to process data" },
      { status: 500 }
    )
  }
}
