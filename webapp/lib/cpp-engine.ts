/**
 * Bridge between the Next.js server and the compiled C++ stroke_analysis binary.
 *
 * The binary's contract (from stroke_analysis.cpp / cxxopts):
 *   stroke_analysis -i <input.csv> -o <output_dir>
 *
 * It writes into <output_dir>:
 *   stroke_full_cycles.csv      StrokeID,StartZero,PosPeak,MidZero,NegPeak,EndZero,
 *                               TotalTime,TimeToPosPeak,TimeToNegPeak,Area_Positive,Area_Negative
 *   stroke_phase_metrics.xlsx   per-stroke phase metrics workbook
 *   strokes_overview.png        full-session Acc_Y overview plot
 *   stroke_plots/stroke_N.png   one plot per detected stroke
 */

import { spawn } from "child_process"
import { mkdir, writeFile, readFile, readdir, rm, stat, cp } from "fs/promises"
import { join } from "path"
import { tmpdir } from "os"

export const ENGINE_BINARY = process.env.CPP_ENGINE_PATH || "/app/engine/stroke_analysis"
export const ENGINE_ENABLED = process.env.USE_CPP_ENGINE === "true"
export const ARTIFACT_ROOT = process.env.ARTIFACT_DIR || "/app/data/runs"
const ENGINE_TIMEOUT_MS = Number(process.env.CPP_ENGINE_TIMEOUT_MS || 120_000)

export interface EngineCycle {
  strokeId: number
  startZero: number
  posPeak: number
  midZero: number
  negPeak: number
  endZero: number
  totalTime: number
  timeToPosPeak: number
  timeToNegPeak: number
  areaPositive: number
  areaNegative: number
}

export interface EngineArtifact {
  kind: "cycles_csv" | "phase_metrics_xlsx" | "overview_plot" | "stroke_plot"
  filename: string
  /** Path relative to ARTIFACT_ROOT, e.g. "<runId>/stroke_plots/stroke_3.png" */
  relativePath: string
  byteSize: number
  strokeId: number | null
}

export interface EngineResult {
  ran: boolean
  exitCode: number | null
  cycles: EngineCycle[]
  artifacts: EngineArtifact[]
  stdout: string
  stderr: string
  error?: string
}

function runBinary(binary: string, args: string[], cwd: string) {
  return new Promise<{ code: number | null; stdout: string; stderr: string; error?: string }>((resolve) => {
    const child = spawn(binary, args, {
      cwd,
      env: {
        ...process.env,
        // matplotlib-cpp needs a headless backend inside the container.
        MPLBACKEND: process.env.MPLBACKEND || "Agg",
        // Point matplotlib at a writable cache dir (the container home is not
        // writable), avoiding the per-run temp-dir warning and slow imports.
        MPLCONFIGDIR: process.env.MPLCONFIGDIR || join(cwd, ".mplconfig"),
      },
    })
    let stdout = ""
    let stderr = ""
    const timer = setTimeout(() => {
      child.kill("SIGKILL")
      resolve({ code: null, stdout, stderr, error: `Engine timed out after ${ENGINE_TIMEOUT_MS} ms` })
    }, ENGINE_TIMEOUT_MS)

    child.stdout.on("data", (chunk) => {
      stdout += chunk.toString()
    })
    child.stderr.on("data", (chunk) => {
      stderr += chunk.toString()
    })
    child.on("error", (error) => {
      clearTimeout(timer)
      resolve({ code: null, stdout, stderr, error: error.message })
    })
    child.on("close", (code) => {
      clearTimeout(timer)
      resolve({ code, stdout, stderr })
    })
  })
}

export function parseCyclesCsv(content: string): EngineCycle[] {
  const lines = content.trim().split(/\r?\n/).filter(Boolean)
  const cycles: EngineCycle[] = []
  for (const line of lines.slice(1)) {
    const cells = line.split(",").map((cell) => cell.trim())
    if (cells.length < 11) continue
    const numbers = cells.map(Number)
    if (numbers.some((value) => !Number.isFinite(value))) continue
    cycles.push({
      strokeId: Math.round(numbers[0]),
      startZero: Math.round(numbers[1]),
      posPeak: Math.round(numbers[2]),
      midZero: Math.round(numbers[3]),
      negPeak: Math.round(numbers[4]),
      endZero: Math.round(numbers[5]),
      totalTime: numbers[6],
      timeToPosPeak: numbers[7],
      timeToNegPeak: numbers[8],
      areaPositive: numbers[9],
      areaNegative: numbers[10],
    })
  }
  return cycles
}

async function collectArtifacts(runId: string, resultDir: string): Promise<EngineArtifact[]> {
  const destination = join(ARTIFACT_ROOT, runId)
  await mkdir(join(destination, "stroke_plots"), { recursive: true })
  const artifacts: EngineArtifact[] = []

  const record = async (
    kind: EngineArtifact["kind"],
    source: string,
    filename: string,
    relativePath: string,
    strokeId: number | null,
  ) => {
    try {
      const target = join(ARTIFACT_ROOT, relativePath)
      await cp(source, target)
      const info = await stat(target)
      artifacts.push({ kind, filename, relativePath, byteSize: info.size, strokeId })
    } catch {
      // Missing artifact (e.g. plotting disabled) is not fatal.
    }
  }

  await record("cycles_csv", join(resultDir, "stroke_full_cycles.csv"), "stroke_full_cycles.csv", `${runId}/stroke_full_cycles.csv`, null)
  await record("phase_metrics_xlsx", join(resultDir, "stroke_phase_metrics.xlsx"), "stroke_phase_metrics.xlsx", `${runId}/stroke_phase_metrics.xlsx`, null)
  await record("overview_plot", join(resultDir, "strokes_overview.png"), "strokes_overview.png", `${runId}/strokes_overview.png`, null)

  try {
    const plots = await readdir(join(resultDir, "stroke_plots"))
    for (const plot of plots.filter((name) => name.endsWith(".png")).sort()) {
      const strokeId = Number(plot.match(/stroke_(\d+)\.png$/)?.[1] ?? NaN)
      await record(
        "stroke_plot",
        join(resultDir, "stroke_plots", plot),
        plot,
        `${runId}/stroke_plots/${plot}`,
        Number.isFinite(strokeId) ? strokeId : null,
      )
    }
  } catch {
    // No per-stroke plots produced.
  }

  return artifacts
}

/**
 * Writes the normalized CSV to a scratch dir, runs the binary, then parses the
 * cycles CSV and persists plots/xlsx under ARTIFACT_ROOT/<runId>/.
 */
export async function runStrokeEngine(runId: string, engineCsv: string): Promise<EngineResult> {
  if (!ENGINE_ENABLED) {
    return { ran: false, exitCode: null, cycles: [], artifacts: [], stdout: "", stderr: "", error: "USE_CPP_ENGINE is not enabled" }
  }

  const workDir = join(tmpdir(), "stroke-engine", runId)
  const resultDir = join(workDir, "results")
  const inputPath = join(workDir, "input.csv")

  try {
    await mkdir(resultDir, { recursive: true })
    await writeFile(inputPath, engineCsv, "utf-8")

    const run = await runBinary(ENGINE_BINARY, ["-i", inputPath, "-o", resultDir], workDir)

    let cycles: EngineCycle[] = []
    try {
      cycles = parseCyclesCsv(await readFile(join(resultDir, "stroke_full_cycles.csv"), "utf-8"))
    } catch {
      cycles = []
    }

    if (!cycles.length && run.code !== 0) {
      return {
        ran: true,
        exitCode: run.code,
        cycles: [],
        artifacts: [],
        stdout: run.stdout,
        stderr: run.stderr,
        error: run.error || `Engine exited with code ${run.code}`,
      }
    }

    const artifacts = await collectArtifacts(runId, resultDir)
    return { ran: true, exitCode: run.code, cycles, artifacts, stdout: run.stdout, stderr: run.stderr, error: run.error }
  } catch (error) {
    return {
      ran: true,
      exitCode: null,
      cycles: [],
      artifacts: [],
      stdout: "",
      stderr: "",
      error: error instanceof Error ? error.message : "Engine invocation failed",
    }
  } finally {
    await rm(workDir, { recursive: true, force: true }).catch(() => {})
  }
}
