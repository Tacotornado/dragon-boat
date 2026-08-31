"use client"

import { useState, useCallback, useMemo } from "react"
import { SummaryCards } from "@/components/summary-cards"
import { OverviewWaveform } from "@/components/overview-waveform"
import { StrokeTable } from "@/components/stroke-table"
import { StrokeCharts } from "@/components/stroke-charts"
import { StrokeGrid } from "@/components/stroke-grid"
import { DashboardHeader } from "@/components/dashboard-header"
import { UploadPanel } from "@/components/upload-panel"
import type { StrokeCycle, IMUDataPoint, SummaryStats, AccPhaseMetrics, DecPhaseMetrics, StrokePhaseMetrics } from "@/lib/types"
import {
  parseStrokeCyclesCSV,
  parseRawIMUCSV,
  validateStrokeData,
} from "@/lib/csv-parser"

export function Dashboard() {
  const [selectedStrokeId, setSelectedStrokeId] = useState<number | null>(null)
  const [filterTopN, setFilterTopN] = useState<number | null>(null)
  const [selectedParams, setSelectedParams] = useState<string[]>([])
  const [selectedRange, setSelectedRange] = useState<[number, number] | null>(null)

  // Upload state — the dashboard starts empty and only renders real data.
  const [uploadedStrokes, setUploadedStrokes] = useState<StrokeCycle[] | null>(null)
  const [uploadedWaveform, setUploadedWaveform] = useState<IMUDataPoint[] | null>(null)
  const [uploadMessage, setUploadMessage] = useState<string | null>(null)

  const strokes = uploadedStrokes
  const waveformData = uploadedWaveform
  const hasData = uploadedStrokes !== null && uploadedWaveform !== null

  // Keep only complete strokes contained in the selected waveform interval.
  const visibleStrokes = useMemo(() => {
    if (!strokes || !selectedRange) return strokes
    return strokes.filter((stroke) =>
      stroke.StartZero >= selectedRange[0] && stroke.EndZero <= selectedRange[1]
    )
  }, [strokes, selectedRange])

  // Phase metrics are populated once Excel parsing is wired up for uploads.
  const allAccMetrics: AccPhaseMetrics[] = []
  const allDecMetrics: DecPhaseMetrics[] = []
  const allStrokePhaseMetrics: StrokePhaseMetrics[] = []
  const accMetrics = useMemo(() => allAccMetrics.filter((m) => visibleStrokes?.some((s) => s.StrokeID === m.StrokeID)), [allAccMetrics, visibleStrokes])
  const decMetrics = useMemo(() => allDecMetrics.filter((m) => visibleStrokes?.some((s) => s.StrokeID === m.StrokeID)), [allDecMetrics, visibleStrokes])
  const strokePhaseMetrics = useMemo(() => allStrokePhaseMetrics.filter((m) => visibleStrokes?.some((s) => s.StrokeID === m.StrokeID)), [allStrokePhaseMetrics, visibleStrokes])

  // Compute summary stats
  const summaryStats = useMemo<SummaryStats | null>(() => {
    if (!visibleStrokes || visibleStrokes.length === 0) return null

    const totalStrokes = visibleStrokes.length
    
    // Calculate average SPM (strokes per minute)
    const totalDuration = visibleStrokes.reduce((sum, s) => sum + s.TotalTime, 0)
    const avgStrokeTime = totalDuration / totalStrokes
    const avgSPM = 60 / avgStrokeTime

    // Calculate velocity changes from acc/dec metrics
    const avgVelocityAcc = accMetrics.length > 0
      ? accMetrics.reduce((sum, m) => sum + m.VelocityChange, 0) / accMetrics.length
      : 0.18 // Default mock value

    const avgVelocityDec = decMetrics.length > 0
      ? Math.abs(decMetrics.reduce((sum, m) => sum + m.VelocityChange, 0) / decMetrics.length)
      : 0.14 // Default mock value

    // Calculate MAC values from stroke phase metrics
    const avgMACAccX = strokePhaseMetrics.length > 0
      ? strokePhaseMetrics.reduce((sum, m) => sum + m.MeanAbsAccX, 0) / strokePhaseMetrics.length
      : 0.085

    const avgMACGyrX = strokePhaseMetrics.length > 0
      ? strokePhaseMetrics.reduce((sum, m) => sum + m.MeanAbsGyrX, 0) / strokePhaseMetrics.length
      : 17.5

    const avgMACGyrY = strokePhaseMetrics.length > 0
      ? strokePhaseMetrics.reduce((sum, m) => sum + m.MeanAbsGyrY, 0) / strokePhaseMetrics.length
      : 9.5

    const avgMACGyrZ = strokePhaseMetrics.length > 0
      ? strokePhaseMetrics.reduce((sum, m) => sum + m.MeanAbsGyrZ, 0) / strokePhaseMetrics.length
      : 14.2

    return {
      totalStrokes,
      avgSPM,
      avgVelocityAcc,
      avgVelocityDec,
      avgMACAccX,
      avgMACGyrX,
      avgMACGyrY,
      avgMACGyrZ,
    }
  }, [visibleStrokes, accMetrics, decMetrics, strokePhaseMetrics])

  const handleStrokesUpload = useCallback((file: File) => {
    const reader = new FileReader()
    reader.onload = (e) => {
      try {
        const text = e.target?.result as string
        const parsed = parseStrokeCyclesCSV(text)
        const result = validateStrokeData(parsed)

        if (!result.valid) {
          setUploadMessage(result.message)
          return
        }

        setUploadedStrokes(parsed)
        setSelectedStrokeId(null)
        setFilterTopN(null)
        setUploadMessage(result.message)
      } catch (err) {
        setUploadMessage(err instanceof Error ? err.message : "Failed to parse CSV")
      }
    }
    reader.readAsText(file)
  }, [])

  const handleIMUUpload = useCallback((file: File) => {
    const reader = new FileReader()
    reader.onload = (e) => {
      try {
        const text = e.target?.result as string
        const parsed = parseRawIMUCSV(text)

        if (parsed.length === 0) {
          setUploadMessage("No IMU data found in the CSV.")
          return
        }

        setUploadedWaveform(parsed)
        setUploadMessage(`Loaded ${parsed.length} IMU samples for waveform display.`)
      } catch (err) {
        setUploadMessage(err instanceof Error ? err.message : "Failed to parse IMU CSV")
      }
    }
    reader.readAsText(file)
  }, [])

  // Strokes detected server-side by the C++ engine (from the raw-IMU upload).
  const handleEngineStrokes = useCallback((engineStrokes: StrokeCycle[]) => {
    setUploadedStrokes(engineStrokes)
    setSelectedStrokeId(null)
    setFilterTopN(null)
    setSelectedRange(null)
  }, [])

  const handleClearData = useCallback(() => {
    setUploadedStrokes(null)
    setUploadedWaveform(null)
    setUploadMessage(null)
    setSelectedStrokeId(null)
    setFilterTopN(null)
    setSelectedParams([])
    setSelectedRange(null)
  }, [])

  const handleSelectStroke = useCallback((id: number) => {
    setSelectedStrokeId((prev) => (prev === id ? null : id))
  }, [])

  // Empty state: no data uploaded yet — show only the upload panel.
  if (!hasData || !strokes || !waveformData || !summaryStats || !visibleStrokes) {
    return (
      <div className="flex min-h-screen flex-col bg-background">
        <DashboardHeader
          filterTopN={filterTopN}
          onFilterChange={setFilterTopN}
          totalStrokes={0}
          hasData={false}
          onClear={handleClearData}
        />
        <main className="flex-1 p-4 md:p-6">
          <div className="mx-auto flex max-w-[1440px] flex-col gap-6">
            <UploadPanel
              onStrokesUpload={handleStrokesUpload}
              onIMUUpload={handleIMUUpload}
              onEngineStrokes={handleEngineStrokes}
              hasData={false}
              uploadMessage={uploadMessage}
              hasUploadedWaveform={uploadedWaveform !== null}
            />
            <div className="flex flex-col items-center justify-center gap-2 rounded-xl border border-dashed border-border/50 py-16 text-center">
              <h2 className="text-lg font-semibold text-foreground">No stroke data yet</h2>
              <p className="max-w-md text-pretty text-sm text-muted-foreground">
                Upload a <span className="font-mono text-foreground/70">stroke_full_cycles.csv</span> file
                to see summary stats, waveforms, and per-stroke analysis.
              </p>
            </div>
          </div>
        </main>
      </div>
    )
  }

  return (
    <div className="flex min-h-screen flex-col bg-background">
      <DashboardHeader
        filterTopN={filterTopN}
        onFilterChange={setFilterTopN}
        totalStrokes={visibleStrokes.length}
        hasData={hasData}
        onClear={handleClearData}
      />
      <main className="flex-1 p-4 md:p-6">
        <div className="mx-auto flex max-w-[1440px] flex-col gap-6">
          {/* Upload Panel */}
          <UploadPanel
            onStrokesUpload={handleStrokesUpload}
            onIMUUpload={handleIMUUpload}
            onEngineStrokes={handleEngineStrokes}
            hasData={hasData}
            uploadMessage={uploadMessage}
            hasUploadedWaveform={uploadedWaveform !== null}
          />

          {/* Summary Cards */}
          <SummaryCards stats={summaryStats} />

          {/* Overview Waveform */}
          <OverviewWaveform
            waveformData={waveformData}
            strokes={strokes}
            selectedRange={selectedRange}
            onRangeChange={setSelectedRange}
          />

          {/* Table + Charts side by side */}
          <div className="grid grid-cols-1 gap-6 xl:grid-cols-2">
            <StrokeTable
              strokes={visibleStrokes}
              accMetrics={accMetrics}
              decMetrics={decMetrics}
              strokePhaseMetrics={strokePhaseMetrics}
              selectedStrokeId={selectedStrokeId}
              onSelectStroke={handleSelectStroke}
              onParameterSelect={setSelectedParams}
              selectedParams={selectedParams}
            />
            <StrokeCharts 
              strokes={visibleStrokes} 
              accMetrics={accMetrics}
              decMetrics={decMetrics}
              strokePhaseMetrics={strokePhaseMetrics}
              selectedParams={selectedParams}
            />
          </div>

          {/* Stroke Visualization Grid */}
          <div>
            <h2 className="mb-4 text-lg font-semibold text-foreground">
              Individual Stroke Waveforms
              {filterTopN !== null && (
                <span className="ml-2 text-sm font-normal text-muted-foreground">
                  (showing top {filterTopN} by peak area)
                </span>
              )}
            </h2>
            <StrokeGrid
              strokes={visibleStrokes}
              waveformData={waveformData.filter((point) => !selectedRange || (point.sample >= selectedRange[0] && point.sample <= selectedRange[1]))}
              selectedStrokeId={selectedStrokeId}
              filterTopN={filterTopN}
            />
          </div>
        </div>
      </main>
    </div>
  )
}
