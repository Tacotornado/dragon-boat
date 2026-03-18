"use client"

import { useState, useCallback } from "react"
import useSWR from "swr"
import { SummaryCards } from "@/components/summary-cards"
import { OverviewWaveform } from "@/components/overview-waveform"
import { StrokeTable } from "@/components/stroke-table"
import { StrokeCharts } from "@/components/stroke-charts"
import { StrokeGrid } from "@/components/stroke-grid"
import { DashboardHeader } from "@/components/dashboard-header"
import { Skeleton } from "@/components/ui/skeleton"
import type { StrokeCycle } from "@/lib/types"

const fetcher = (url: string) => fetch(url).then((res) => res.json())

export function Dashboard() {
  const [selectedStrokeId, setSelectedStrokeId] = useState<number | null>(null)
  const [filterTopN, setFilterTopN] = useState<number | null>(null)

  const { data: strokes, isLoading: strokesLoading, mutate: mutateStrokes } = useSWR<StrokeCycle[]>(
    "/api/strokes",
    fetcher
  )

  const { data: waveformData, isLoading: waveformLoading } = useSWR<
    { sample: number; accY: number }[]
  >("/api/overview", fetcher)

  const handleRunAnalysis = useCallback(async () => {
    const res = await fetch("/api/run", { method: "POST" })
    const result = await res.json()
    if (result.success) {
      mutateStrokes()
    }
    return result
  }, [mutateStrokes])

  const handleSelectStroke = useCallback((id: number) => {
    setSelectedStrokeId((prev) => (prev === id ? null : id))
  }, [])

  const isLoading = strokesLoading || waveformLoading

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
        onRunAnalysis={handleRunAnalysis}
        filterTopN={filterTopN}
        onFilterChange={setFilterTopN}
        totalStrokes={strokes.length}
      />
      <main className="flex-1 p-4 md:p-6">
        <div className="mx-auto flex max-w-[1440px] flex-col gap-6">
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
