"use client"

import { useRef, useEffect } from "react"
import { StrokeVisualization } from "@/components/stroke-visualization"
import type { StrokeCycle } from "@/lib/types"

interface StrokeGridProps {
  strokes: StrokeCycle[]
  waveformData: { sample: number; accY: number }[]
  selectedStrokeId: number | null
  filterTopN: number | null
}

export function StrokeGrid({
  strokes,
  waveformData,
  selectedStrokeId,
  filterTopN,
}: StrokeGridProps) {
  const refs = useRef<Map<number, HTMLDivElement>>(new Map())

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
          />
        </div>
      ))}
    </div>
  )
}
