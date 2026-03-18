"use client"

import { useState, useCallback } from "react"
import useSWR from "swr"
import { SummaryCards } from "@/components/summary-cards"
import { OverviewWaveform } from "@/components/overview-waveform"
import { StrokeTable } from "@/components/stroke-table"
import { StrokeCharts } from "@/components/stroke-charts"
import { StrokeGrid } from "@/components/stroke-grid"
import { DashboardHeader } from "@/components/dashboard-header"
import { UploadPanel } from "@/components/upload-panel"
import { Skeleton } from "@/components/ui/skeleton"
import type { StrokeCycle } from "@/lib/types"
import {
  parseStrokeCyclesCSV,
  parseRawIMUCSV,
  validateStrokeData,
} from "@/lib/csv-parser"

const fetcher = (url: string) => fetch(url).then((res) => res.json())

type DataSource = "demo" | "uploaded"

export function Dashboard() {
  const [selectedStrokeId, setSelectedStrokeId] = useState<number | null>(null)
  const [filterTopN, setFilterTopN] = useState<number | null>(null)

  // Upload state
  const [dataSource, setDataSource] = useState<DataSource>("demo")
  const [uploadedStrokes, setUploadedStrokes] = useState<StrokeCycle[] | null>(null)
  const [uploadedWaveform, setUploadedWaveform] = useState<{ sample: number; accY: number }[] | null>(null)
  const [uploadMessage, setUploadMessage] = useState<string | null>(null)

  // Demo data from API (fallback)
  const { data: demoStrokes, isLoading: strokesLoading } = useSWR<StrokeCycle[]>(
    "/api/strokes",
    fetcher
  )
  const { data: demoWaveform, isLoading: waveformLoading } = useSWR<
    { sample: number; accY: number }[]
  >("/api/overview", fetcher)

  // Active data based on source
  const strokes = dataSource === "uploaded" && uploadedStrokes ? uploadedStrokes : demoStrokes
  const waveformData = dataSource === "uploaded" && uploadedWaveform ? uploadedWaveform : demoWaveform

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
        setDataSource("uploaded")
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
          setUploadMessage("No Acc_Y data found in the IMU CSV.")
          return
        }

        setUploadedWaveform(parsed)
        if (dataSource !== "uploaded" && uploadedStrokes) {
          setDataSource("uploaded")
        }
        setUploadMessage(`Loaded ${parsed.length} IMU samples for waveform display.`)
      } catch (err) {
        setUploadMessage(err instanceof Error ? err.message : "Failed to parse IMU CSV")
      }
    }
    reader.readAsText(file)
  }, [dataSource, uploadedStrokes])

  const handleResetToDemo = useCallback(() => {
    setDataSource("demo")
    setUploadedStrokes(null)
    setUploadedWaveform(null)
    setUploadMessage(null)
    setSelectedStrokeId(null)
    setFilterTopN(null)
  }, [])

  const handleSelectStroke = useCallback((id: number) => {
    setSelectedStrokeId((prev) => (prev === id ? null : id))
  }, [])

  const isLoading = dataSource === "demo" && (strokesLoading || waveformLoading)

  if (isLoading) {
    return (
      <div className="flex min-h-screen flex-col bg-background">
        <div className="border-b border-border/50 bg-card px-6 py-4">
          <Skeleton className="h-8 w-64" />
        </div>
        <main className="flex-1 p-6">
          <div className="mx-auto flex max-w-[1440px] flex-col gap-6">
            <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-4">
              {Array.from({ length: 4 }).map((_, i) => (
                <Skeleton key={i} className="h-24 rounded-xl" />
              ))}
            </div>
            <Skeleton className="h-[340px] rounded-xl" />
            <div className="grid grid-cols-1 gap-6 lg:grid-cols-2">
              <Skeleton className="h-[500px] rounded-xl" />
              <Skeleton className="h-[500px] rounded-xl" />
            </div>
          </div>
        </main>
      </div>
    )
  }

  if (!strokes || !waveformData) return null

  return (
    <div className="flex min-h-screen flex-col bg-background">
      <DashboardHeader
        filterTopN={filterTopN}
        onFilterChange={setFilterTopN}
        totalStrokes={strokes.length}
        dataSource={dataSource}
        onResetToDemo={handleResetToDemo}
      />
      <main className="flex-1 p-4 md:p-6">
        <div className="mx-auto flex max-w-[1440px] flex-col gap-6">
          {/* Upload Panel */}
          <UploadPanel
            onStrokesUpload={handleStrokesUpload}
            onIMUUpload={handleIMUUpload}
            dataSource={dataSource}
            uploadMessage={uploadMessage}
            hasUploadedWaveform={uploadedWaveform !== null}
          />

          {/* Summary Cards */}
          <SummaryCards strokes={strokes} />

          {/* Overview Waveform */}
          <OverviewWaveform waveformData={waveformData} strokes={strokes} />

          {/* Table + Charts side by side */}
          <div className="grid grid-cols-1 gap-6 xl:grid-cols-2">
            <StrokeTable
              strokes={strokes}
              selectedStrokeId={selectedStrokeId}
              onSelectStroke={handleSelectStroke}
            />
            <StrokeCharts strokes={strokes} />
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
              strokes={strokes}
              waveformData={waveformData}
              selectedStrokeId={selectedStrokeId}
              filterTopN={filterTopN}
            />
          </div>
        </div>
      </main>
    </div>
  )
}
