"use client"

import { useMemo, useState } from "react"
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  ReferenceLine,
} from "recharts"
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card"
import { ChartContainer, ChartTooltip, ChartTooltipContent } from "@/components/ui/chart"
import type { StrokeCycle } from "@/lib/types"

interface OverviewWaveformProps {
  waveformData: { sample: number; accY: number }[]
  strokes: StrokeCycle[]
}

const WAVE_COLOR = "oklch(0.72 0.19 190)"
const PEAK_POS_COLOR = "oklch(0.70 0.18 150)"
const GRID_COLOR = "oklch(0.25 0.01 260)"

export function OverviewWaveform({ waveformData, strokes }: OverviewWaveformProps) {
  const [hoveredStroke, setHoveredStroke] = useState<number | null>(null)

  // Downsample waveform for performance (every 3rd sample)
  const downsampled = useMemo(() => {
    if (waveformData.length > 3000) {
      const step = Math.max(1, Math.floor(waveformData.length / 1000))
      return waveformData.filter((_, i) => i % step === 0)
    }
    return waveformData.filter((_, i) => i % 3 === 0)
  }, [waveformData])

  // Check if we have actual waveform coverage for stroke boundaries
  const hasStrokeCoverage = useMemo(() => {
    if (strokes.length === 0 || waveformData.length === 0) return false
    const maxSample = waveformData[waveformData.length - 1]?.sample ?? 0
    return strokes[strokes.length - 1].EndZero <= maxSample
  }, [strokes, waveformData])

  const peakMarkers = useMemo(() => {
    if (!hasStrokeCoverage) return []
    return strokes.map((s) => ({
      strokeId: s.StrokeID,
      startSample: s.StartZero,
    }))
  }, [strokes, hasStrokeCoverage])

  return (
    <Card className="border-border/50">
      <CardHeader>
        <CardTitle className="text-lg">Signal Overview</CardTitle>
        <CardDescription>
          {hasStrokeCoverage
            ? `Filtered Acc_Y waveform with ${strokes.length} detected stroke boundaries`
            : `Acc_Y waveform (${waveformData.length.toLocaleString()} samples)`}
        </CardDescription>
      </CardHeader>
      <CardContent>
        <ChartContainer
          config={{
            accY: {
              label: "Acc Y",
              color: WAVE_COLOR,
            },
          }}
          className="h-[280px] w-full"
        >
          <LineChart
              data={downsampled}
              margin={{ top: 8, right: 12, left: 0, bottom: 4 }}
              onMouseMove={(state) => {
                if (state?.activePayload?.[0] && hasStrokeCoverage) {
                  const sampleIdx = state.activePayload[0].payload.sample
                  const stroke = strokes.find(
                    (s) => sampleIdx >= s.StartZero && sampleIdx <= s.EndZero
                  )
                  setHoveredStroke(stroke?.StrokeID ?? null)
                }
              }}
              onMouseLeave={() => setHoveredStroke(null)}
            >
              <CartesianGrid
                strokeDasharray="3 3"
                stroke={GRID_COLOR}
                vertical={false}
              />
              <XAxis
                dataKey="sample"
                tick={{ fontSize: 10 }}
                tickFormatter={(val) => `${(val / 120).toFixed(1)}s`}
                stroke="oklch(0.65 0.02 260)"
                interval="preserveStartEnd"
              />
              <YAxis
                tick={{ fontSize: 10 }}
                stroke="oklch(0.65 0.02 260)"
                domain={["auto", "auto"]}
              />
              <ReferenceLine y={0} stroke="oklch(0.35 0.01 260)" strokeDasharray="2 4" />
              <ChartTooltip
                content={
                  <ChartTooltipContent
                    formatter={(value) => [
                      `${Number(value).toFixed(4)} m/s\u00B2`,
                      hoveredStroke ? `Stroke ${hoveredStroke}` : "Acc Y",
                    ]}
                  />
                }
              />
              <Line
                type="monotone"
                dataKey="accY"
                stroke={WAVE_COLOR}
                strokeWidth={1.5}
                dot={false}
                activeDot={{
                  r: 3,
                  fill: WAVE_COLOR,
                  stroke: "oklch(0.13 0.005 260)",
                  strokeWidth: 2,
                }}
              />
              {/* Stroke boundary reference lines */}
              {peakMarkers.map((pm) => (
                <ReferenceLine
                  key={`start-${pm.strokeId}`}
                  x={pm.startSample}
                  stroke={
                    hoveredStroke === pm.strokeId
                      ? PEAK_POS_COLOR
                      : "oklch(0.30 0.01 260)"
                  }
                  strokeDasharray="2 3"
                  strokeWidth={hoveredStroke === pm.strokeId ? 1.5 : 0.5}
                />
              ))}
            </LineChart>
        </ChartContainer>
      </CardContent>
    </Card>
  )
}
