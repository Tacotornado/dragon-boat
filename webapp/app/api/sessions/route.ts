import { NextResponse } from "next/server"
import { createSession, listSessions } from "@/lib/session-store"

// GET /api/sessions - List all sessions
export async function GET() {
  const sessions = await listSessions()
  return NextResponse.json(sessions)
}

// POST /api/sessions - Create a new session
export async function POST(request: Request) {
  try {
    const body = await request.json()
    const { name, athleteName } = body
    
    if (!name) {
      return NextResponse.json(
        { error: "Session name is required" },
        { status: 400 }
      )
    }
    
    const session = await createSession(name, athleteName)
    return NextResponse.json(session, { status: 201 })
  } catch {
    return NextResponse.json(
      { error: "Failed to create session" },
      { status: 500 }
    )
  }
}
