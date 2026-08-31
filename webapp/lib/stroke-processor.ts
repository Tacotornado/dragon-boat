/**
 * Stroke Processor - Bridges Node.js and C++ analysis engine
 * 
 * For real-time processing, you have two options:
 * 
 * OPTION 1: Modify your C++ to be a streaming service
 *   - Run as a daemon that reads from stdin, writes to stdout
 *   - Use spawn() to communicate via pipes
 * 
 * OPTION 2: Call C++ binary per chunk (simpler, shown below)
 *   - Write chunk to temp file, run analysis, read results
 *   - Higher latency but easier to implement
 * 
 * This file demonstrates Option 2 with a TypeScript fallback
 * for environments where the C++ binary isn't available.
 */

import { spawn } from "child_process"
import { writeFile, readFile, unlink, mkdir } from "fs/promises"
import { join } from "path"
import { tmpdir } from "os"
import {
  IMUDataPoint,
  StrokeCycle,
  AccPhaseMetrics,
  DecPhaseMetrics,
  StrokePhaseMetrics,
} from "./types"
import { addStroke, getSession } from "./session-store"

// Buffer to accumulate samples between stroke boundaries
const sessionBuffers = new Map<string, {
  samples: IMUDataPoint[]
  lastProcessedIndex: number
  strokeIdCounter: number
}>()

interface ProcessResult {
  newStrokes: StrokeCycle[]
  totalStrokes: number
}

/**
 * Process a chunk of IMU data for a session.
 * In production, this calls your C++ binary.
 * For development/demo, it uses a TypeScript implementation.
 */
export async function processIMUChunk(
  sessionId: string,
  newSamples: IMUDataPoint[]
): Promise<ProcessResult> {
  // Initialize buffer for this session if needed
  if (!sessionBuffers.has(sessionId)) {
    sessionBuffers.set(sessionId, {
      samples: [],
      lastProcessedIndex: 0,
      strokeIdCounter: 1,
    })
  }
  
  const buffer = sessionBuffers.get(sessionId)!
  buffer.samples.push(...newSamples)
  
  // Try to use C++ binary first
  const useCpp = process.env.USE_CPP_ENGINE === "true"
  
  if (useCpp) {
    return await processWithCpp(sessionId, buffer)
  } else {
    return await processWithTypeScript(sessionId, buffer)
  }
}

/**
 * Process using the C++ binary
 */
async function processWithCpp(
  sessionId: string,
  buffer: { samples: IMUDataPoint[]; lastProcessedIndex: number; strokeIdCounter: number }
): Promise<ProcessResult> {
  const newStrokes: StrokeCycle[] = []
  
  // Need minimum samples to detect a stroke (roughly 1 stroke worth at 120Hz)
  const MIN_SAMPLES = 120 // ~1 second
  
  if (buffer.samples.length - buffer.lastProcessedIndex < MIN_SAMPLES) {
    const session = await getSession(sessionId)
    return { newStrokes: [], totalStrokes: session?.strokeCount || 0 }
  }
  
  try {
    // Create temp directory for this session
    const tempDir = join(tmpdir(), "stroke-analysis", sessionId)
    await mkdir(tempDir, { recursive: true })
    
    // Write samples to CSV (format expected by C++)
    const csvPath = join(tempDir, "chunk.csv")
    const csvContent = [
      "SampleTimeFine,Acc_X,Acc_Y,Gyr_X,Gyr_Y,Gyr_Z",
      ...buffer.samples.slice(buffer.lastProcessedIndex).map((s, i) =>
        `${buffer.lastProcessedIndex + i},${s.Acc_X},${s.Acc_Y},${s.Gyr_X},${s.Gyr_Y},${s.Gyr_Z}`
      ),
    ].join("\n")
    
    await writeFile(csvPath, csvContent)
    
    // Call C++ binary
    const cppPath = process.env.CPP_ENGINE_PATH || "./stroke_analysis_realtime"
    const result = await runCppProcess(cppPath, [csvPath, tempDir])
    
    if (result.success) {
      // Read output CSV
      const outputPath = join(tempDir, "stroke_full_cycles.csv")
      const outputContent = await readFile(outputPath, "utf-8")
      const parsedStrokes = parseStrokeCsv(outputContent, buffer.strokeIdCounter)
      
      for (const stroke of parsedStrokes) {
        newStrokes.push(stroke)
        await addStroke(sessionId, stroke)
        buffer.strokeIdCounter++
      }
      
      buffer.lastProcessedIndex = buffer.samples.length
    }
    
    // Cleanup temp files
    await unlink(csvPath).catch(() => {})
    
  } catch (error) {
    console.error("C++ processing error:", error)
  }
  
  const session = await getSession(sessionId)
  return { newStrokes, totalStrokes: session?.strokeCount || 0 }
}

/**
 * TypeScript fallback for when C++ binary isn't available.
 * Implements basic stroke detection (simplified version of your C++ algorithm).
 */
async function processWithTypeScript(
  sessionId: string,
  buffer: { samples: IMUDataPoint[]; lastProcessedIndex: number; strokeIdCounter: number }
): Promise<ProcessResult> {
  const newStrokes: StrokeCycle[] = []
  const samples = buffer.samples
  
  // Simple zero-crossing detection on Acc_Y
  // This is a simplified version - your C++ has better filtering
  const zeroCrossings: number[] = []
  
  for (let i = buffer.lastProcessedIndex + 1; i < samples.length; i++) {
    const prev = samples[i - 1].Acc_Y
    const curr = samples[i].Acc_Y
    
    // Detect zero crossing (with small hysteresis)
    if ((prev < -0.1 && curr >= 0) || (prev > 0.1 && curr <= 0)) {
      zeroCrossings.push(i)
    }
  }
  
  // Need at least 3 zero crossings for a full stroke cycle
  // Pattern: start(0) -> mid(crossing) -> end(0)
  while (zeroCrossings.length >= 3) {
    const startZero = zeroCrossings.shift()!
    const midZero = zeroCrossings.shift()!
    const endZero = zeroCrossings[0] // peek, don't remove
    
    // Find peaks between crossings
    let posPeak = startZero
    let negPeak = midZero
    let maxVal = -Infinity
    let minVal = Infinity
    
    for (let i = startZero; i < midZero; i++) {
      if (samples[i].Acc_Y > maxVal) {
        maxVal = samples[i].Acc_Y
        posPeak = i
      }
    }
    
    for (let i = midZero; i < endZero; i++) {
      if (samples[i].Acc_Y < minVal) {
        minVal = samples[i].Acc_Y
        negPeak = i
      }
    }
    
    // Calculate metrics (simplified - your C++ has trapz integration)
    const sampleRate = 120 // Hz
    const totalTime = (endZero - startZero) / sampleRate
    const timeToPosPeak = (posPeak - startZero) / sampleRate
    const timeToNegPeak = (negPeak - startZero) / sampleRate
    
    // Rough area calculation (your C++ uses proper trapz)
    let areaPos = 0
    let areaNeg = 0
    for (let i = startZero; i < midZero; i++) {
      areaPos += Math.max(0, samples[i].Acc_Y) / sampleRate
    }
    for (let i = midZero; i < endZero; i++) {
      areaNeg += Math.abs(Math.min(0, samples[i].Acc_Y)) / sampleRate
    }
    
    const stroke: StrokeCycle = {
      StrokeID: buffer.strokeIdCounter,
      StartZero: startZero,
      PosPeak: posPeak,
      MidZero: midZero,
      NegPeak: negPeak,
      EndZero: endZero,
      TotalTime: totalTime,
      TimeToPosPeak: timeToPosPeak,
      TimeToNegPeak: timeToNegPeak,
      Area_Positive: areaPos,
      Area_Negative: areaNeg,
    }
    
    newStrokes.push(stroke)
    await addStroke(sessionId, stroke)
    buffer.strokeIdCounter++
    
    buffer.lastProcessedIndex = endZero
  }
  
  const session = await getSession(sessionId)
  return { newStrokes, totalStrokes: session?.strokeCount || 0 }
}

/**
 * Run the C++ binary as a child process
 */
function runCppProcess(
  binaryPath: string,
  args: string[]
): Promise<{ success: boolean; stdout: string; stderr: string }> {
  return new Promise((resolve) => {
    const proc = spawn(binaryPath, args)
    let stdout = ""
    let stderr = ""
    
    proc.stdout.on("data", (data) => {
      stdout += data.toString()
    })
    
    proc.stderr.on("data", (data) => {
      stderr += data.toString()
    })
    
    proc.on("close", (code) => {
      resolve({
        success: code === 0,
        stdout,
        stderr,
      })
    })
    
    proc.on("error", (error) => {
      resolve({
        success: false,
        stdout,
        stderr: error.message,
      })
    })
  })
}

/**
 * Parse the C++ output CSV
 */
function parseStrokeCsv(content: string, startId: number): StrokeCycle[] {
  const lines = content.trim().split("\n")
  const strokes: StrokeCycle[] = []
  
  // Skip header
  for (let i = 1; i < lines.length; i++) {
    const cols = lines[i].split(",")
    if (cols.length >= 11) {
      strokes.push({
        StrokeID: startId + i - 1,
        StartZero: parseInt(cols[1]),
        PosPeak: parseInt(cols[2]),
        MidZero: parseInt(cols[3]),
        NegPeak: parseInt(cols[4]),
        EndZero: parseInt(cols[5]),
        TotalTime: parseFloat(cols[6]),
        TimeToPosPeak: parseFloat(cols[7]),
        TimeToNegPeak: parseFloat(cols[8]),
        Area_Positive: parseFloat(cols[9]),
        Area_Negative: parseFloat(cols[10]),
      })
    }
  }
  
  return strokes
}

/**
 * Clear session buffer (call when session ends)
 */
export function clearSessionBuffer(sessionId: string): void {
  sessionBuffers.delete(sessionId)
}
