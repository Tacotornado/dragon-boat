import { getSession, getSessionStats } from "@/lib/session-store"

// Store active connections for each session
const connections = new Map<string, Set<ReadableStreamDefaultController>>()

// GET /api/sessions/[id]/stream - SSE endpoint for real-time updates
export async function GET(
  request: Request,
  { params }: { params: Promise<{ id: string }> }
) {
  const { id } = await params
  
  const session = await getSession(id)
  if (!session) {
    return new Response("Session not found", { status: 404 })
  }
  
  const stream = new ReadableStream({
    start(controller) {
      // Register this connection
      if (!connections.has(id)) {
        connections.set(id, new Set())
      }
      connections.get(id)!.add(controller)
      
      // Send initial state
      const initialData = JSON.stringify({
        type: "connected",
        sessionId: id,
        status: session.status,
        strokeCount: session.strokeCount,
        sampleCount: session.sampleCount,
      })
      controller.enqueue(`data: ${initialData}\n\n`)
      
      // Cleanup on disconnect
      request.signal.addEventListener("abort", () => {
        connections.get(id)?.delete(controller)
        if (connections.get(id)?.size === 0) {
          connections.delete(id)
        }
      })
    },
  })
  
  return new Response(stream, {
    headers: {
      "Content-Type": "text/event-stream",
      "Cache-Control": "no-cache",
      "Connection": "keep-alive",
    },
  })
}

/**
 * Broadcast an event to all clients watching a session.
 * Call this from the ingest endpoint when new data arrives.
 */
export function broadcastToSession(
  sessionId: string,
  eventType: string,
  data: unknown
) {
  const sessionConnections = connections.get(sessionId)
  if (!sessionConnections) return
  
  const message = JSON.stringify({
    type: eventType,
    timestamp: Date.now(),
    ...data,
  })
  
  for (const controller of sessionConnections) {
    try {
      controller.enqueue(`data: ${message}\n\n`)
    } catch {
      // Connection closed, remove it
      sessionConnections.delete(controller)
    }
  }
}

/**
 * Broadcast stroke detected event
 */
export async function broadcastStrokeDetected(
  sessionId: string,
  newStrokes: unknown[],
  totalStrokes: number
) {
  const stats = await getSessionStats(sessionId)
  
  broadcastToSession(sessionId, "stroke_detected", {
    newStrokes,
    totalStrokes,
    stats,
  })
}

/**
 * Broadcast new IMU samples (for live waveform)
 */
export function broadcastIMUData(
  sessionId: string,
  samples: unknown[]
) {
  broadcastToSession(sessionId, "imu_data", {
    samples,
    sampleCount: samples.length,
  })
}
