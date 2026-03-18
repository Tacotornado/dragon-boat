import { NextResponse } from "next/server"
import { mockStrokes } from "@/lib/mock-data"

export async function GET() {
  return NextResponse.json(mockStrokes)
}
