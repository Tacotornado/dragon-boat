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
import { ChartContainer, ChartTooltip, ChartTooltipContent } from "@/components/ui/chart"
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select"
import type { IMUDataPoint, StrokeCycle, AxisType } from "@/lib/types"

interface LiveWaveformProps {
  samples: IMUDataPoint[]
  strokes: StrokeCycle[]
}

const AXIS_OPTIONS: { value: AxisType; label: string; color: string }[] = [
  { value: "Acc_Y", label: "Acc Y", color: "oklch(0.72 0.19 190)" },
  { value: "Acc_X", label: "Acc X", color: "oklch(0.65 0.20 30)" },
  { value: "Gyr_X", label: "Gyr X", color: "oklch(0.78 0.15 85)" },
  { value: "Gyr_Y", label: "Gyr Y", color: "oklch(0.60 0.22 280)" },
  { value: "Gyr_Z", label: "Gyr Z", color: "oklch(0.70 0.18 150)" },
]

export function LiveWaveform({ samples, strokes }: LiveWaveformProps) {
  const [selectedAxis, setSelectedAxis] = useState<AxisType>("Acc_Y")
  
  const axisConfig = AXIS_OPTIONS.find(a => a.value === selectedAxis)!
  
  // Prepare chart data
  const chartData = useMemo(() => {
    return samples.map((s, i) => ({
      sample: s.sample || i,
      value: s[selectedAxis],
    }))
  }, [samples, selectedAxis])
  
  // Find stroke boundaries that are within our sample range
  const strokeBoundaries = useMemo(() => {
    if (samples.length === 0) return []
    
    const minSample = samples[0]?.sample || 0
    const maxSample = samples[samples.length - 1]?.sample || samples.length
    
    return strokes
      .filter(s => s.StartZero >= minSample && s.StartZero <= maxSample)
      .map(s => s.StartZero)
  }, [samples, strokes])
  
  if (samples.length === 0) {
    return (
      <div className="flex h-[280px] items-center justify-center rounded-lg border border-dashed border-border/50 bg-secondary/20">
        <p className="text-sm text-muted-foreground">
          Waiting for live data...
        </p>
      </div>
    )
  }
  
  return (
    <div>
      <div className="mb-3 flex items-center justify-between">
        <Select
          value={selectedAxis}
          onValueChange={(v) => setSelectedAxis(v as AxisType)}
        >
          <SelectTrigger className="w-[140px] h-8 text-xs">
            <SelectValue />
          </SelectTrigger>
          <SelectContent>
            {AXIS_OPTIONS.map(axis => (
              <SelectItem key={axis.value} value={axis.value}>
                <span className="flex items-center gap-2">
                  <span 
                    className="size-2 rounded-full" 
                    style={{ backgroundColor: axis.color }}
                  />
                  {axis.label}
                </span>
              </SelectItem>
            ))}
          </SelectContent>
        </Select>
        
        <span className="text-xs text-muted-foreground">
          Last {samples.length} samples (~{(samples.length / 120).toFixed(1)}s)
        </span>
      </div>
      
      <ChartContainer
        config={{
          value: { label: axisConfig.label, color: axisConfig.color },
        }}
        className="h-[250px] w-full"
      >
        <LineChart
          data={chartData}
          margin={{ top: 8, right: 12, left: 0, bottom: 4 }}
        >
          <CartesianGrid
            strokeDasharray="3 3"
            stroke="oklch(0.25 0.01 260)"
            vertical={false}
          />
          <XAxis
            dataKey="sample"
            tick={{ fontSize: 10 }}
            stroke="oklch(0.65 0.02 260)"
            tickFormatter={(v) => `${Math.floor(v / 120)}s`}
          />
          <YAxis
            tick={{ fontSize: 10 }}
            stroke="oklch(0.65 0.02 260)"
            domain={["auto", "auto"]}
          />
          <ChartTooltip content={<ChartTooltipContent />} />
          
          {/* Zero line */}
          <ReferenceLine
            y={0}
            stroke="oklch(0.65 0.02 260)"
            strokeDasharray="3 3"
          />
          
          {/* Stroke boundaries */}
          {strokeBoundaries.map((x, i) => (
            <ReferenceLine
              key={i}
              x={x}
              stroke="oklch(0.72 0.19 190)"
              strokeDasharray="2 2"
              strokeOpacity={0.5}
            />
          ))}
          
          <Line
            type="monotone"
            dataKey="value"
            stroke={axisConfig.color}
            strokeWidth={1.5}
            dot={false}
            isAnimationActive={false} // Disable animation for performance
          />
        </LineChart>
      </ChartContainer>
    </div>
  )
}
