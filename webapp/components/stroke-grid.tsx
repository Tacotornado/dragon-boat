"use client"

import { useRef, useEffect, useState } from "react"
import { StrokeVisualization } from "@/components/stroke-visualization"
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"
import type { StrokeCycle, IMUDataPoint, AxisType } from "@/lib/types"

interface StrokeGridProps {
  strokes: StrokeCycle[]
  waveformData: IMUDataPoint[]
  selectedStrokeId: number | null
  filterTopN: number | null
}

const AXIS_OPTIONS: { value: AxisType; label: string }[] = [
  { value: "Acc_Y", label: "Acc Y" },
  { value: "Acc_X", label: "Acc X" },
  { value: "Gyr_X", label: "Gyr X" },
  { value: "Gyr_Y", label: "Gyr Y" },
  { value: "Gyr_Z", label: "Gyr Z" },
]

export function StrokeGrid({
  strokes,
  waveformData,
  selectedStrokeId,
  filterTopN,
}: StrokeGridProps) {
  const refs = useRef<Map<number, HTMLDivElement>>(new Map())
  const [selectedAxis, setSelectedAxis] = useState<AxisType>("Acc_Y")

  // Scroll to selected stroke
  useEffect(() => {
    if (selectedStrokeId !== null) {
      const el = refs.current.get(selectedStrokeId)
      if (el) {
        el.scrollIntoView({ behavior: "smooth", block: "center" })
      }
    }
  }, [selectedStrokeId])

  // Apply top-N filter
  let displayStrokes = strokes
  if (filterTopN !== null) {
    const sorted = [...strokes].sort(
      (a, b) => b.Area_Positive - a.Area_Positive
    )
    const topIds = new Set(sorted.slice(0, filterTopN).map((s) => s.StrokeID))
    displayStrokes = strokes.filter((s) => topIds.has(s.StrokeID))
  }

  return (
    <div className="flex flex-col gap-4">
      {/* Axis selector */}
      <div className="flex items-center justify-between">
        <p className="text-sm text-muted-foreground">
          Showing {displayStrokes.length} of {strokes.length} strokes
        </p>
        <Select value={selectedAxis} onValueChange={(v) => setSelectedAxis(v as AxisType)}>
          <SelectTrigger className="w-[120px]">
            <SelectValue />
          </SelectTrigger>
          <SelectContent>
            {AXIS_OPTIONS.map((axis) => (
              <SelectItem key={axis.value} value={axis.value}>
                {axis.label}
              </SelectItem>
            ))}
          </SelectContent>
        </Select>
      </div>

      {/* Grid */}
      <div className="grid grid-cols-1 gap-3 sm:grid-cols-2 sm:gap-4 lg:grid-cols-3">
        {displayStrokes.map((stroke) => (
          <div
            key={stroke.StrokeID}
            ref={(el) => {
              if (el) refs.current.set(stroke.StrokeID, el)
              else refs.current.delete(stroke.StrokeID)
            }}
            className={`rounded-xl transition-all duration-300 ${
              selectedStrokeId === stroke.StrokeID
                ? "ring-2 ring-primary ring-offset-2 ring-offset-background scale-[1.01]"
                : "ring-0"
            }`}
          >
            <StrokeVisualization
              stroke={stroke}
              waveformData={waveformData}
              selectedAxis={selectedAxis}
            />
          </div>
        ))}
      </div>
    </div>
  )
}
