"use client"

import { useMemo, useState } from "react"
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  ReferenceLine,
  Brush,
} from "recharts"
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card"
import { ChartContainer, ChartTooltip, ChartTooltipContent } from "@/components/ui/chart"
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"
import type { StrokeCycle, IMUDataPoint, AxisType } from "@/lib/types"

interface OverviewWaveformProps {
  waveformData: IMUDataPoint[]
  strokes: StrokeCycle[]
  selectedRange: [number, number] | null
  onRangeChange: (range: [number, number] | null) => void
}

const AXIS_OPTIONS: { value: AxisType; label: string; unit: string }[] = [
  { value: "Acc_Y", label: "Acc Y", unit: "m/s²" },
  { value: "Acc_X", label: "Acc X", unit: "m/s²" },
  { value: "Gyr_X", label: "Gyr X", unit: "deg/s" },
  { value: "Gyr_Y", label: "Gyr Y", unit: "deg/s" },
  { value: "Gyr_Z", label: "Gyr Z", unit: "deg/s" },
]

const AXIS_COLORS: Record<AxisType, string> = {
  Acc_X: "oklch(0.65 0.20 30)",
  Acc_Y: "oklch(0.72 0.19 190)",
  Gyr_X: "oklch(0.70 0.18 150)",
  Gyr_Y: "oklch(0.78 0.15 85)",
  Gyr_Z: "oklch(0.60 0.22 280)",
}

const PEAK_POS_COLOR = "oklch(0.70 0.18 150)"
const GRID_COLOR = "oklch(0.25 0.01 260)"

export function OverviewWaveform({
  waveformData,
  strokes,
  selectedRange,
  onRangeChange,
}: OverviewWaveformProps) {
  const [hoveredStroke, setHoveredStroke] = useState<number | null>(null)
  const [selectedAxis, setSelectedAxis] = useState<AxisType>("Acc_Y")

  // Downsample waveform for performance
  const downsampled = useMemo(() => {
    if (waveformData.length > 3000) {
      const step = Math.max(1, Math.floor(waveformData.length / 1000))
      return waveformData.filter((_, i) => i % step === 0)
    }
    return waveformData.filter((_, i) => i % 3 === 0)
  }, [waveformData])

  // Transform data for selected axis
  const chartData = useMemo(() => {
    return downsampled.map(d => ({
      sample: d.sample,
      value: d[selectedAxis],
    }))
  }, [downsampled, selectedAxis])

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

  const selectedAxisInfo = AXIS_OPTIONS.find(a => a.value === selectedAxis)
  const axisColor = AXIS_COLORS[selectedAxis]

  return (
    <Card className="border-border/50">
      <CardHeader>
        <div className="flex flex-col gap-3 sm:flex-row sm:items-center sm:justify-between">
          <div>
            <CardTitle className="text-lg">Signal Overview</CardTitle>
            <CardDescription>
              {hasStrokeCoverage
                ? `${selectedAxisInfo?.label} waveform with ${strokes.length} stroke boundaries`
                : `${selectedAxisInfo?.label} waveform (${waveformData.length.toLocaleString()} samples)`}
            </CardDescription>
          </div>
          <Select value={selectedAxis} onValueChange={(v) => setSelectedAxis(v as AxisType)}>
            <SelectTrigger className="w-[140px]">
              <SelectValue />
            </SelectTrigger>
            <SelectContent>
              {AXIS_OPTIONS.map((axis) => (
                <SelectItem key={axis.value} value={axis.value}>
                  {axis.label} ({axis.unit})
                </SelectItem>
              ))}
            </SelectContent>
          </Select>
        </div>
      </CardHeader>
      <CardContent>
        <ChartContainer
          config={{
            value: {
              label: selectedAxisInfo?.label ?? "Value",
              color: axisColor,
            },
          }}
          className="h-[280px] w-full"
        >
          <LineChart
            data={chartData}
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
                    `${Number(value).toFixed(4)} ${selectedAxisInfo?.unit}`,
                    hoveredStroke ? `Stroke ${hoveredStroke}` : selectedAxisInfo?.label,
                  ]}
                />
              }
            />
            <Line
              type="monotone"
              dataKey="value"
              stroke={axisColor}
              strokeWidth={1.5}
              dot={false}
              activeDot={{
                r: 3,
                fill: axisColor,
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
            <Brush
              dataKey="sample"
              height={24}
              travellerWidth={10}
              stroke={axisColor}
              fill="oklch(0.17 0.008 260)"
              startIndex={selectedRange ? Math.max(0, chartData.findIndex((d) => d.sample >= selectedRange[0])) : 0}
              endIndex={selectedRange ? Math.max(0, chartData.findIndex((d) => d.sample >= selectedRange[1])) : Math.max(0, chartData.length - 1)}
              onChange={(range) => {
                if (range.startIndex == null || range.endIndex == null || chartData.length === 0) {
                  onRangeChange(null)
                  return
                }
                const start = chartData[range.startIndex]?.sample
                const end = chartData[range.endIndex]?.sample
                if (typeof start === "number" && typeof end === "number") {
                  onRangeChange([Math.min(start, end), Math.max(start, end)])
                }
              }}
            />
          </LineChart>
        </ChartContainer>
      </CardContent>
    </Card>
  )
}
