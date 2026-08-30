import { NextResponse } from "next/server"
import { generateWaveformData } from "@/lib/mock-data"

export async function GET() {
  const waveform = generateWaveformData()
  return NextResponse.json(waveform)
}
