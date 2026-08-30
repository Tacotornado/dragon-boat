"use client"

import { useMemo } from "react"
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  ReferenceDot,
} from "recharts"
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card"
import { ChartContainer, ChartTooltip, ChartTooltipContent } from "@/components/ui/chart"
import type { StrokeCycle, IMUDataPoint, AxisType } from "@/lib/types"

interface StrokeVisualizationProps {
  stroke: StrokeCycle
  waveformData: IMUDataPoint[]
  selectedAxis: AxisType
}

const AXIS_COLORS: Record<AxisType, string> = {
  Acc_X: "oklch(0.65 0.20 30)",
  Acc_Y: "oklch(0.72 0.19 190)",
  Gyr_X: "oklch(0.70 0.18 150)",
  Gyr_Y: "oklch(0.78 0.15 85)",
  Gyr_Z: "oklch(0.60 0.22 280)",
}

const AXIS_UNITS: Record<AxisType, string> = {
  Acc_X: "m/s²",
  Acc_Y: "m/s²",
  Gyr_X: "deg/s",
  Gyr_Y: "deg/s",
  Gyr_Z: "deg/s",
}

const POS_PEAK_COLOR = "oklch(0.70 0.18 150)"
const NEG_PEAK_COLOR = "oklch(0.65 0.20 30)"
const ZERO_COLOR = "oklch(0.78 0.15 85)"
const GRID_COLOR = "oklch(0.25 0.01 260)"

export function StrokeVisualization({ stroke, waveformData, selectedAxis }: StrokeVisualizationProps) {
  // Check if waveform data covers this stroke's sample range
  const hasWaveformCoverage = useMemo(() => {
    if (!waveformData || waveformData.length === 0) return false
    const maxSample = waveformData[waveformData.length - 1]?.sample ?? 0
    return stroke.EndZero <= maxSample
  }, [stroke, waveformData])

  const segmentData = useMemo(() => {
    if (!hasWaveformCoverage) return []
    const padSamples = 10
    const start = Math.max(0, stroke.StartZero - padSamples)
    const end = Math.min(waveformData.length - 1, stroke.EndZero + padSamples)
    return waveformData.slice(start, end + 1).map(d => ({
      sample: d.sample,
      value: d[selectedAxis],
    }))
  }, [stroke, waveformData, hasWaveformCoverage, selectedAxis])

  const posPeakY = hasWaveformCoverage ? (waveformData[stroke.PosPeak]?.[selectedAxis] ?? 0) : 0
  const negPeakY = hasWaveformCoverage ? (waveformData[stroke.NegPeak]?.[selectedAxis] ?? 0) : 0

  const axisColor = AXIS_COLORS[selectedAxis]
  const axisUnit = AXIS_UNITS[selectedAxis]

  return (
    <Card className="border-border/50">
      <CardHeader className="pb-2">
        <div className="flex items-center justify-between">
          <CardTitle className="text-base">
            Stroke {stroke.StrokeID}
          </CardTitle>
          <span className="text-xs font-mono text-muted-foreground">
            [{stroke.StartZero} - {stroke.EndZero}]
          </span>
        </div>
        <CardDescription className="text-xs">
          {stroke.TotalTime.toFixed(3)}s | Area+: {stroke.Area_Positive.toFixed(2)} | Area-: {stroke.Area_Negative.toFixed(2)}
        </CardDescription>
      </CardHeader>
      <CardContent>
        {hasWaveformCoverage ? (
          <ChartContainer
            config={{
              value: { label: selectedAxis.replace("_", " "), color: axisColor },
            }}
            className="h-[180px] w-full"
          >
            <LineChart
              data={segmentData}
              margin={{ top: 8, right: 8, left: -10, bottom: 4 }}
            >
              <CartesianGrid
                strokeDasharray="3 3"
                stroke={GRID_COLOR}
                vertical={false}
              />
              <XAxis
                dataKey="sample"
                tick={{ fontSize: 9 }}
                stroke="oklch(0.65 0.02 260)"
                tickFormatter={(val) => `${val}`}
              />
              <YAxis
                tick={{ fontSize: 9 }}
                stroke="oklch(0.65 0.02 260)"
              />
              <ChartTooltip 
                content={
                  <ChartTooltipContent 
                    formatter={(value) => [`${Number(value).toFixed(4)} ${axisUnit}`, selectedAxis.replace("_", " ")]}
                  />
                } 
              />
              <Line
                type="monotone"
                dataKey="value"
                stroke={axisColor}
                strokeWidth={2}
                dot={false}
              />
              {/* Only show peak markers for Acc_Y (where peaks are detected) */}
              {selectedAxis === "Acc_Y" && (
                <>
                  <ReferenceDot
                    x={stroke.PosPeak}
                    y={posPeakY}
                    r={5}
                    fill={POS_PEAK_COLOR}
                    stroke="oklch(0.13 0.005 260)"
                    strokeWidth={2}
                  />
                  <ReferenceDot
                    x={stroke.NegPeak}
                    y={negPeakY}
                    r={5}
                    fill={NEG_PEAK_COLOR}
                    stroke="oklch(0.13 0.005 260)"
                    strokeWidth={2}
                  />
                  <ReferenceDot
                    x={stroke.StartZero}
                    y={waveformData[stroke.StartZero]?.[selectedAxis] ?? 0}
                    r={3}
                    fill={ZERO_COLOR}
                    stroke="oklch(0.13 0.005 260)"
                    strokeWidth={1.5}
                  />
                  <ReferenceDot
                    x={stroke.MidZero}
                    y={waveformData[stroke.MidZero]?.[selectedAxis] ?? 0}
                    r={3}
                    fill={ZERO_COLOR}
                    stroke="oklch(0.13 0.005 260)"
                    strokeWidth={1.5}
                  />
                  <ReferenceDot
                    x={stroke.EndZero}
                    y={waveformData[stroke.EndZero]?.[selectedAxis] ?? 0}
                    r={3}
                    fill={ZERO_COLOR}
                    stroke="oklch(0.13 0.005 260)"
                    strokeWidth={1.5}
                  />
                </>
              )}
            </LineChart>
          </ChartContainer>
        ) : (
          <div className="flex h-[180px] flex-col items-center justify-center gap-2 rounded-lg bg-secondary/30">
            <p className="text-sm text-muted-foreground">No waveform data for this range</p>
            <p className="text-xs text-muted-foreground/70">
              Upload your raw IMU CSV to see stroke waveforms
            </p>
          </div>
        )}
      </CardContent>
    </Card>
  )
}
