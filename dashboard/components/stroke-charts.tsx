"use client"

import {
  LineChart,
  Line,
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
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
import type { StrokeCycle } from "@/lib/types"

interface StrokeChartsProps {
  strokes: StrokeCycle[]
}

const DURATION_COLOR = "oklch(0.72 0.19 190)"
const POS_AREA_COLOR = "oklch(0.70 0.18 150)"
const NEG_AREA_COLOR = "oklch(0.65 0.20 30)"
const GRID_COLOR = "oklch(0.25 0.01 260)"

export function StrokeCharts({ strokes }: StrokeChartsProps) {
  const chartData = strokes.map((s) => ({
    stroke: `${s.StrokeID}`,
    duration: s.TotalTime,
    areaPos: s.Area_Positive,
    areaNeg: s.Area_Negative,
    ratio: s.Area_Positive / Math.max(0.01, s.Area_Negative),
  }))

  return (
    <Card className="border-border/50">
      <CardHeader>
        <CardTitle className="text-lg">Performance Charts</CardTitle>
        <CardDescription>
          Stroke-by-stroke trends and area comparisons
        </CardDescription>
      </CardHeader>
      <CardContent>
        <Tabs defaultValue="duration">
          <TabsList className="mb-4">
            <TabsTrigger value="duration">Duration Trend</TabsTrigger>
            <TabsTrigger value="areas">Area Comparison</TabsTrigger>
            <TabsTrigger value="ratio">Symmetry Ratio</TabsTrigger>
          </TabsList>

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
