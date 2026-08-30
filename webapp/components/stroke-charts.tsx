"use client"

import { useMemo } from "react"
import {
  LineChart,
  Line,
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Legend,
} from "recharts"
import {
  Card,
  CardContent,
  CardHeader,
  CardTitle,
  CardDescription,
} from "@/components/ui/card"
import { ChartContainer, ChartTooltip, ChartTooltipContent } from "@/components/ui/chart"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import type { StrokeCycle, AccPhaseMetrics, DecPhaseMetrics, StrokePhaseMetrics } from "@/lib/types"

interface StrokeChartsProps {
  strokes: StrokeCycle[]
  accMetrics?: AccPhaseMetrics[]
  decMetrics?: DecPhaseMetrics[]
  strokePhaseMetrics?: StrokePhaseMetrics[]
  selectedParams?: string[]
}

const PARAM_COLORS = [
  "oklch(0.72 0.19 190)",
  "oklch(0.65 0.20 30)",
]

const DURATION_COLOR = "oklch(0.72 0.19 190)"
const POS_AREA_COLOR = "oklch(0.70 0.18 150)"
const NEG_AREA_COLOR = "oklch(0.65 0.20 30)"
const GRID_COLOR = "oklch(0.25 0.01 260)"

// Parameter definitions for custom chart
const PARAM_GETTERS: Record<string, {
  label: string
  unit: string
  getValue: (s: StrokeCycle, acc?: AccPhaseMetrics, dec?: DecPhaseMetrics, sp?: StrokePhaseMetrics) => number
}> = {
  TotalTime: { label: "Duration", unit: "s", getValue: (s) => s.TotalTime },
  TimeToPosPeak: { label: "T+ Peak", unit: "s", getValue: (s) => s.TimeToPosPeak },
  TimeToNegPeak: { label: "T- Peak", unit: "s", getValue: (s) => s.TimeToNegPeak },
  Area_Positive: { label: "Area+", unit: "", getValue: (s) => s.Area_Positive },
  Area_Negative: { label: "Area-", unit: "", getValue: (s) => s.Area_Negative },
  acc_Duration: { label: "Acc Duration", unit: "s", getValue: (_, acc) => acc?.Duration ?? 0 },
  acc_VelocityChange: { label: "ΔVel Acc", unit: "m/s", getValue: (_, acc) => acc?.VelocityChange ?? 0 },
  acc_AccY_PosPeak: { label: "Acc-Y+ Peak", unit: "m/s²", getValue: (_, acc) => acc?.AccY_PosPeak ?? 0 },
  acc_RAD: { label: "RAD", unit: "m/s³", getValue: (_, acc) => acc?.RAD ?? 0 },
  acc_MeanAbsAccX: { label: "|Acc-X|", unit: "m/s²", getValue: (_, acc) => acc?.MeanAbsAccX ?? 0 },
  acc_MeanAbsGyrX: { label: "|Gyr-X|", unit: "deg/s", getValue: (_, acc) => acc?.MeanAbsGyrX ?? 0 },
  dec_Duration: { label: "Dec Duration", unit: "s", getValue: (_, __, dec) => dec?.Duration ?? 0 },
  dec_VelocityChange: { label: "ΔVel Dec", unit: "m/s", getValue: (_, __, dec) => dec?.VelocityChange ?? 0 },
  dec_AccY_NegPeak: { label: "Acc-Y- Peak", unit: "m/s²", getValue: (_, __, dec) => dec?.AccY_NegPeak ?? 0 },
  dec_MeanAbsAccX: { label: "|Acc-X| Dec", unit: "m/s²", getValue: (_, __, dec) => dec?.MeanAbsAccX ?? 0 },
  stroke_GyrX_PosPeak: { label: "Gyr-X+ Peak", unit: "deg/s", getValue: (_, __, ___, sp) => sp?.GyrX_PosPeak ?? 0 },
  stroke_GyrX_NegPeak: { label: "Gyr-X- Peak", unit: "deg/s", getValue: (_, __, ___, sp) => sp?.GyrX_NegPeak ?? 0 },
  stroke_MeanAbsAccX: { label: "MAC Acc-X", unit: "m/s²", getValue: (_, __, ___, sp) => sp?.MeanAbsAccX ?? 0 },
  stroke_MeanAbsGyrX: { label: "MAC Gyr-X", unit: "deg/s", getValue: (_, __, ___, sp) => sp?.MeanAbsGyrX ?? 0 },
  stroke_MeanAbsGyrY: { label: "MAC Gyr-Y", unit: "deg/s", getValue: (_, __, ___, sp) => sp?.MeanAbsGyrY ?? 0 },
  stroke_MeanAbsGyrZ: { label: "MAC Gyr-Z", unit: "deg/s", getValue: (_, __, ___, sp) => sp?.MeanAbsGyrZ ?? 0 },
}

export function StrokeCharts({ 
  strokes, 
  accMetrics, 
  decMetrics, 
  strokePhaseMetrics,
  selectedParams = [] 
}: StrokeChartsProps) {
  const chartData = useMemo(() => {
    return strokes.map((s) => {
      const acc = accMetrics?.find(m => m.StrokeID === s.StrokeID)
      const dec = decMetrics?.find(m => m.StrokeID === s.StrokeID)
      const sp = strokePhaseMetrics?.find(m => m.StrokeID === s.StrokeID)
      
      const data: Record<string, string | number> = {
        stroke: `${s.StrokeID}`,
        duration: s.TotalTime,
        areaPos: s.Area_Positive,
        areaNeg: s.Area_Negative,
        ratio: s.Area_Positive / Math.max(0.01, s.Area_Negative),
      }

      // Add custom parameter values
      selectedParams.forEach((param, idx) => {
        const getter = PARAM_GETTERS[param]
        if (getter) {
          data[`param${idx}`] = getter.getValue(s, acc, dec, sp)
        }
      })

      return data
    })
  }, [strokes, accMetrics, decMetrics, strokePhaseMetrics, selectedParams])

  const hasCustomParams = selectedParams.length > 0
  const param0 = selectedParams[0] ? PARAM_GETTERS[selectedParams[0]] : null
  const param1 = selectedParams[1] ? PARAM_GETTERS[selectedParams[1]] : null

  return (
    <Card className="border-border/50">
      <CardHeader>
        <CardTitle className="text-lg">Performance Charts</CardTitle>
        <CardDescription>
          {hasCustomParams 
            ? `Comparing: ${param0?.label ?? ""}${param1 ? ` vs ${param1.label}` : ""}`
            : "Stroke-by-stroke trends and area comparisons"}
        </CardDescription>
      </CardHeader>
      <CardContent>
        <Tabs defaultValue={hasCustomParams ? "custom" : "duration"}>
          <TabsList className="mb-4 flex-wrap h-auto gap-1">
            {hasCustomParams && <TabsTrigger value="custom">Custom Compare</TabsTrigger>}
            <TabsTrigger value="duration">Duration Trend</TabsTrigger>
            <TabsTrigger value="areas">Area Comparison</TabsTrigger>
            <TabsTrigger value="ratio">Symmetry Ratio</TabsTrigger>
          </TabsList>

          {hasCustomParams && (
            <TabsContent value="custom">
              <ChartContainer
                config={{
                  param0: { label: param0?.label ?? "Param 1", color: PARAM_COLORS[0] },
                  param1: { label: param1?.label ?? "Param 2", color: PARAM_COLORS[1] },
                }}
                className="h-[300px] w-full"
              >
                <LineChart data={chartData} margin={{ top: 8, right: 12, left: 0, bottom: 4 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke={GRID_COLOR} vertical={false} />
                  <XAxis
                    dataKey="stroke"
                    tick={{ fontSize: 10 }}
                    stroke="oklch(0.65 0.02 260)"
                    label={{ value: "Stroke ID", position: "insideBottom", offset: -2, fontSize: 11, fill: "oklch(0.65 0.02 260)" }}
                  />
                  <YAxis
                    yAxisId="left"
                    tick={{ fontSize: 10 }}
                    stroke={PARAM_COLORS[0]}
                  />
                  {param1 && (
                    <YAxis
                      yAxisId="right"
                      orientation="right"
                      tick={{ fontSize: 10 }}
                      stroke={PARAM_COLORS[1]}
                    />
                  )}
                  <ChartTooltip content={<ChartTooltipContent />} />
                  <Legend />
                  <Line
                    yAxisId="left"
                    type="monotone"
                    dataKey="param0"
                    name={`${param0?.label} (${param0?.unit})`}
                    stroke={PARAM_COLORS[0]}
                    strokeWidth={2}
                    dot={{ r: 3, fill: PARAM_COLORS[0], stroke: "oklch(0.13 0.005 260)", strokeWidth: 1.5 }}
                    activeDot={{ r: 5, fill: PARAM_COLORS[0] }}
                  />
                  {param1 && (
                    <Line
                      yAxisId="right"
                      type="monotone"
                      dataKey="param1"
                      name={`${param1.label} (${param1.unit})`}
                      stroke={PARAM_COLORS[1]}
                      strokeWidth={2}
                      dot={{ r: 3, fill: PARAM_COLORS[1], stroke: "oklch(0.13 0.005 260)", strokeWidth: 1.5 }}
                      activeDot={{ r: 5, fill: PARAM_COLORS[1] }}
                    />
                  )}
                </LineChart>
              </ChartContainer>
            </TabsContent>
          )}

          <TabsContent value="duration">
            <ChartContainer
              config={{
                duration: { label: "Duration (s)", color: DURATION_COLOR },
              }}
              className="h-[300px] w-full"
            >
              <LineChart data={chartData} margin={{ top: 8, right: 12, left: 0, bottom: 4 }}>
                <CartesianGrid strokeDasharray="3 3" stroke={GRID_COLOR} vertical={false} />
                <XAxis
                  dataKey="stroke"
                  tick={{ fontSize: 10 }}
                  stroke="oklch(0.65 0.02 260)"
                  label={{ value: "Stroke ID", position: "insideBottom", offset: -2, fontSize: 11, fill: "oklch(0.65 0.02 260)" }}
                />
                <YAxis
                  tick={{ fontSize: 10 }}
                  stroke="oklch(0.65 0.02 260)"
                  domain={["dataMin - 0.01", "dataMax + 0.01"]}
                />
                <ChartTooltip content={<ChartTooltipContent />} />
                <Line
                  type="monotone"
                  dataKey="duration"
                  stroke={DURATION_COLOR}
                  strokeWidth={2}
                  dot={{ r: 3, fill: DURATION_COLOR, stroke: "oklch(0.13 0.005 260)", strokeWidth: 1.5 }}
                  activeDot={{ r: 5, fill: DURATION_COLOR }}
                />
              </LineChart>
            </ChartContainer>
          </TabsContent>

          <TabsContent value="areas">
            <ChartContainer
              config={{
                areaPos: { label: "Area +", color: POS_AREA_COLOR },
                areaNeg: { label: "Area -", color: NEG_AREA_COLOR },
              }}
              className="h-[300px] w-full"
            >
              <BarChart data={chartData} margin={{ top: 8, right: 12, left: 0, bottom: 4 }}>
                <CartesianGrid strokeDasharray="3 3" stroke={GRID_COLOR} vertical={false} />
                <XAxis
                  dataKey="stroke"
                  tick={{ fontSize: 10 }}
                  stroke="oklch(0.65 0.02 260)"
                  label={{ value: "Stroke ID", position: "insideBottom", offset: -2, fontSize: 11, fill: "oklch(0.65 0.02 260)" }}
                />
                <YAxis tick={{ fontSize: 10 }} stroke="oklch(0.65 0.02 260)" />
                <ChartTooltip content={<ChartTooltipContent />} />
                <Bar dataKey="areaPos" fill={POS_AREA_COLOR} radius={[3, 3, 0, 0]} />
                <Bar dataKey="areaNeg" fill={NEG_AREA_COLOR} radius={[3, 3, 0, 0]} />
              </BarChart>
            </ChartContainer>
          </TabsContent>

          <TabsContent value="ratio">
            <ChartContainer
              config={{
                ratio: { label: "Pos/Neg Ratio", color: "oklch(0.78 0.15 85)" },
              }}
              className="h-[300px] w-full"
            >
              <LineChart data={chartData} margin={{ top: 8, right: 12, left: 0, bottom: 4 }}>
                <CartesianGrid strokeDasharray="3 3" stroke={GRID_COLOR} vertical={false} />
                <XAxis
                  dataKey="stroke"
                  tick={{ fontSize: 10 }}
                  stroke="oklch(0.65 0.02 260)"
                  label={{ value: "Stroke ID", position: "insideBottom", offset: -2, fontSize: 11, fill: "oklch(0.65 0.02 260)" }}
                />
                <YAxis tick={{ fontSize: 10 }} stroke="oklch(0.65 0.02 260)" />
                <ChartTooltip content={<ChartTooltipContent />} />
                <Line
                  type="monotone"
                  dataKey="ratio"
                  stroke="oklch(0.78 0.15 85)"
                  strokeWidth={2}
                  dot={{ r: 3, fill: "oklch(0.78 0.15 85)", stroke: "oklch(0.13 0.005 260)", strokeWidth: 1.5 }}
                  activeDot={{ r: 5, fill: "oklch(0.78 0.15 85)" }}
                />
              </LineChart>
            </ChartContainer>
          </TabsContent>
        </Tabs>
      </CardContent>
    </Card>
  )
}
