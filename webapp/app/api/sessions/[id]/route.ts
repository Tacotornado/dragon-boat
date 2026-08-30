import { NextResponse } from "next/server"
import { getSession, deleteSession, updateSessionStatus } from "@/lib/session-store"

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
    
    const session = await updateSessionStatus(id, status)
    
    if (!session) {
      return NextResponse.json(
        { error: "Session not found" },
        { status: 404 }
      )
    }
    
    return NextResponse.json(session)
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
