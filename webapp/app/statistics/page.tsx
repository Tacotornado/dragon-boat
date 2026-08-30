"use client"

import { useCallback, useState } from "react"
import useSWR from "swr"
import Link from "next/link"
import { ArrowLeft, Upload } from "lucide-react"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table"
import { Skeleton } from "@/components/ui/skeleton"
import type { StrokeCycle, AccPhaseMetrics, DecPhaseMetrics, StrokePhaseMetrics } from "@/lib/types"
import {
  parseStrokeCyclesCSV,
  parseAccPhaseCSV,
  parseDecPhaseCSV,
  parseStrokePhaseCSV,
} from "@/lib/csv-parser"

const fetcher = (url: string) => fetch(url).then((res) => res.json())

function StatCard({ label, value, unit }: { label: string; value: string; unit?: string }) {
  return (
    <div className="rounded-lg bg-secondary/50 p-3">
      <p className="text-xs text-muted-foreground">{label}</p>
      <p className="font-mono text-lg font-semibold">
        {value}
        {unit && <span className="ml-1 text-xs text-muted-foreground">{unit}</span>}
      </p>
    </div>
  )
}

function calculateStats(values: number[]) {
  if (values.length === 0) return { mean: 0, std: 0, min: 0, max: 0 }
  const mean = values.reduce((a, b) => a + b, 0) / values.length
  const variance = values.reduce((sum, v) => sum + Math.pow(v - mean, 2), 0) / values.length
  const std = Math.sqrt(variance)
  const min = Math.min(...values)
  const max = Math.max(...values)
  return { mean, std, min, max }
}

export default function StatisticsPage() {
  const [strokes, setStrokes] = useState<StrokeCycle[] | null>(null)
  const [accMetrics, setAccMetrics] = useState<AccPhaseMetrics[] | null>(null)
  const [decMetrics, setDecMetrics] = useState<DecPhaseMetrics[] | null>(null)
  const [strokePhaseMetrics, setStrokePhaseMetrics] = useState<StrokePhaseMetrics[] | null>(null)
  const [uploadMessage, setUploadMessage] = useState<string | null>(null)

  // Demo data
  const { data: demoStrokes, isLoading } = useSWR<StrokeCycle[]>("/api/strokes", fetcher)

  const activeStrokes = strokes ?? demoStrokes
  const hasUploadedData = strokes !== null

  const handleFileUpload = useCallback((file: File, type: "strokes" | "acc" | "dec" | "stroke") => {
    const reader = new FileReader()
    reader.onload = (e) => {
      try {
        const text = e.target?.result as string
        switch (type) {
          case "strokes":
            setStrokes(parseStrokeCyclesCSV(text))
            setUploadMessage(`Loaded stroke cycles from ${file.name}`)
            break
          case "acc":
            setAccMetrics(parseAccPhaseCSV(text))
            setUploadMessage(`Loaded acceleration phase metrics from ${file.name}`)
            break
          case "dec":
            setDecMetrics(parseDecPhaseCSV(text))
            setUploadMessage(`Loaded deceleration phase metrics from ${file.name}`)
            break
          case "stroke":
            setStrokePhaseMetrics(parseStrokePhaseCSV(text))
            setUploadMessage(`Loaded stroke phase metrics from ${file.name}`)
            break
        }
      } catch (err) {
        setUploadMessage(err instanceof Error ? err.message : "Failed to parse CSV")
      }
    }
    reader.readAsText(file)
  }, [])

  if (isLoading) {
    return (
      <div className="flex min-h-screen flex-col bg-background p-6">
        <Skeleton className="h-8 w-64 mb-6" />
        <div className="grid gap-6 md:grid-cols-3">
          <Skeleton className="h-[400px]" />
          <Skeleton className="h-[400px]" />
          <Skeleton className="h-[400px]" />
        </div>
      </div>
    )
  }

  // Calculate statistics for each phase
  const strokeStats = activeStrokes ? {
    totalTime: calculateStats(activeStrokes.map(s => s.TotalTime)),
    timeToPosPeak: calculateStats(activeStrokes.map(s => s.TimeToPosPeak)),
    timeToNegPeak: calculateStats(activeStrokes.map(s => s.TimeToNegPeak)),
    areaPos: calculateStats(activeStrokes.map(s => s.Area_Positive)),
    areaNeg: calculateStats(activeStrokes.map(s => s.Area_Negative)),
  } : null

  const accStats = accMetrics ? {
    duration: calculateStats(accMetrics.map(m => m.Duration)),
    posPeak: calculateStats(accMetrics.map(m => m.AccY_PosPeak)),
    timeToPeak: calculateStats(accMetrics.map(m => m.TimeToPosPeakPercent)),
    rad: calculateStats(accMetrics.map(m => m.RAD)),
    velocity: calculateStats(accMetrics.map(m => m.VelocityChange)),
    macAccX: calculateStats(accMetrics.map(m => m.MeanAbsAccX)),
    macGyrX: calculateStats(accMetrics.map(m => m.MeanAbsGyrX)),
    macGyrY: calculateStats(accMetrics.map(m => m.MeanAbsGyrY)),
    macGyrZ: calculateStats(accMetrics.map(m => m.MeanAbsGyrZ)),
  } : null

  const decStats = decMetrics ? {
    duration: calculateStats(decMetrics.map(m => m.Duration)),
    negPeak: calculateStats(decMetrics.map(m => m.AccY_NegPeak)),
    timeToNegPeak: calculateStats(decMetrics.map(m => m.TimeToNegPeakPercent)),
    velocity: calculateStats(decMetrics.map(m => m.VelocityChange)),
    macAccX: calculateStats(decMetrics.map(m => m.MeanAbsAccX)),
    macGyrX: calculateStats(decMetrics.map(m => m.MeanAbsGyrX)),
    macGyrY: calculateStats(decMetrics.map(m => m.MeanAbsGyrY)),
    macGyrZ: calculateStats(decMetrics.map(m => m.MeanAbsGyrZ)),
  } : null

  const spStats = strokePhaseMetrics ? {
    gyrXPos: calculateStats(strokePhaseMetrics.map(m => m.GyrX_PosPeak)),
    gyrXNeg: calculateStats(strokePhaseMetrics.map(m => m.GyrX_NegPeak)),
    macAccX: calculateStats(strokePhaseMetrics.map(m => m.MeanAbsAccX)),
    macGyrX: calculateStats(strokePhaseMetrics.map(m => m.MeanAbsGyrX)),
    macGyrY: calculateStats(strokePhaseMetrics.map(m => m.MeanAbsGyrY)),
    macGyrZ: calculateStats(strokePhaseMetrics.map(m => m.MeanAbsGyrZ)),
  } : null

  return (
    <div className="flex min-h-screen flex-col bg-background">
      {/* Header */}
      <header className="sticky top-0 z-20 border-b border-border/50 bg-card/95 backdrop-blur supports-[backdrop-filter]:bg-card/60">
        <div className="mx-auto flex max-w-[1440px] items-center gap-4 px-4 py-3 md:px-6">
          <Link href="/">
            <Button variant="ghost" size="sm" className="gap-2">
              <ArrowLeft className="size-4" />
              Back to Dashboard
            </Button>
          </Link>
          <div className="flex-1">
            <h1 className="text-lg font-semibold">Parameter Statistics</h1>
            <p className="text-xs text-muted-foreground">
              {hasUploadedData ? "Showing uploaded data" : "Showing demo data"} 
              {activeStrokes && ` - ${activeStrokes.length} strokes`}
            </p>
          </div>
        </div>
      </header>

      <main className="flex-1 p-4 md:p-6">
        <div className="mx-auto flex max-w-[1440px] flex-col gap-6">
          {/* Upload Section */}
          <Card className="border-border/50">
            <CardHeader>
              <CardTitle className="text-base">Upload Phase Metrics (Optional)</CardTitle>
              <CardDescription>
                Export Excel sheets to CSV and upload for detailed statistics
              </CardDescription>
            </CardHeader>
            <CardContent>
              {uploadMessage && (
                <p className="mb-4 text-sm text-primary">{uploadMessage}</p>
              )}
              <div className="grid grid-cols-2 gap-3 md:grid-cols-4">
                <UploadButton label="Stroke Cycles CSV" onUpload={(f) => handleFileUpload(f, "strokes")} />
                <UploadButton label="Acceleration Phase CSV" onUpload={(f) => handleFileUpload(f, "acc")} />
                <UploadButton label="Deceleration Phase CSV" onUpload={(f) => handleFileUpload(f, "dec")} />
                <UploadButton label="Stroke Phase CSV" onUpload={(f) => handleFileUpload(f, "stroke")} />
              </div>
            </CardContent>
          </Card>

          {/* Statistics Tabs */}
          <Tabs defaultValue="overview">
            <TabsList className="mb-4">
              <TabsTrigger value="overview">Overview</TabsTrigger>
              <TabsTrigger value="acceleration">Acceleration Phase</TabsTrigger>
              <TabsTrigger value="deceleration">Deceleration Phase</TabsTrigger>
              <TabsTrigger value="stroke">Stroke Phase</TabsTrigger>
            </TabsList>

            <TabsContent value="overview">
              <Card className="border-border/50">
                <CardHeader>
                  <CardTitle>Basic Stroke Statistics</CardTitle>
                  <CardDescription>Summary statistics from stroke_full_cycles.csv</CardDescription>
                </CardHeader>
                <CardContent>
                  {strokeStats ? (
                    <div className="overflow-auto">
                      <Table>
                        <TableHeader>
                          <TableRow>
                            <TableHead>Parameter</TableHead>
                            <TableHead className="text-right">Mean</TableHead>
                            <TableHead className="text-right">Std Dev</TableHead>
                            <TableHead className="text-right">Min</TableHead>
                            <TableHead className="text-right">Max</TableHead>
                          </TableRow>
                        </TableHeader>
                        <TableBody>
                          <StatRow label="Total Time (s)" stats={strokeStats.totalTime} />
                          <StatRow label="Time to Pos Peak (s)" stats={strokeStats.timeToPosPeak} />
                          <StatRow label="Time to Neg Peak (s)" stats={strokeStats.timeToNegPeak} />
                          <StatRow label="Area Positive" stats={strokeStats.areaPos} />
                          <StatRow label="Area Negative" stats={strokeStats.areaNeg} />
                        </TableBody>
                      </Table>
                    </div>
                  ) : (
                    <p className="text-muted-foreground">No stroke data available</p>
                  )}
                </CardContent>
              </Card>
            </TabsContent>

            <TabsContent value="acceleration">
              <Card className="border-border/50">
                <CardHeader>
                  <CardTitle>Acceleration Phase Statistics</CardTitle>
                  <CardDescription>Upload acceleration_phase.csv for detailed metrics</CardDescription>
                </CardHeader>
                <CardContent>
                  {accStats ? (
                    <div className="overflow-auto">
                      <Table>
                        <TableHeader>
                          <TableRow>
                            <TableHead>Parameter</TableHead>
                            <TableHead className="text-right">Mean</TableHead>
                            <TableHead className="text-right">Std Dev</TableHead>
                            <TableHead className="text-right">Min</TableHead>
                            <TableHead className="text-right">Max</TableHead>
                          </TableRow>
                        </TableHeader>
                        <TableBody>
                          <StatRow label="Duration (s)" stats={accStats.duration} />
                          <StatRow label="Acc-Y+ Peak (m/s²)" stats={accStats.posPeak} />
                          <StatRow label="Time to Peak (%)" stats={accStats.timeToPeak} />
                          <StatRow label="RAD (m/s³)" stats={accStats.rad} />
                          <StatRow label="ΔVelocity (m/s)" stats={accStats.velocity} />
                          <StatRow label="|Acc-X| (m/s²)" stats={accStats.macAccX} />
                          <StatRow label="|Gyr-X| (deg/s)" stats={accStats.macGyrX} />
                          <StatRow label="|Gyr-Y| (deg/s)" stats={accStats.macGyrY} />
                          <StatRow label="|Gyr-Z| (deg/s)" stats={accStats.macGyrZ} />
                        </TableBody>
                      </Table>
                    </div>
                  ) : (
                    <div className="text-center py-8">
                      <p className="text-muted-foreground mb-2">No acceleration phase data</p>
                      <p className="text-xs text-muted-foreground">Upload the Acceleration Phase sheet exported from Excel</p>
                    </div>
                  )}
                </CardContent>
              </Card>
            </TabsContent>

            <TabsContent value="deceleration">
              <Card className="border-border/50">
                <CardHeader>
                  <CardTitle>Deceleration Phase Statistics</CardTitle>
                  <CardDescription>Upload deceleration_phase.csv for detailed metrics</CardDescription>
                </CardHeader>
                <CardContent>
                  {decStats ? (
                    <div className="overflow-auto">
                      <Table>
                        <TableHeader>
                          <TableRow>
                            <TableHead>Parameter</TableHead>
                            <TableHead className="text-right">Mean</TableHead>
                            <TableHead className="text-right">Std Dev</TableHead>
                            <TableHead className="text-right">Min</TableHead>
                            <TableHead className="text-right">Max</TableHead>
                          </TableRow>
                        </TableHeader>
                        <TableBody>
                          <StatRow label="Duration (s)" stats={decStats.duration} />
                          <StatRow label="Acc-Y- Peak (m/s²)" stats={decStats.negPeak} />
                          <StatRow label="Time to Neg Peak (%)" stats={decStats.timeToNegPeak} />
                          <StatRow label="ΔVelocity (m/s)" stats={decStats.velocity} />
                          <StatRow label="|Acc-X| (m/s²)" stats={decStats.macAccX} />
                          <StatRow label="|Gyr-X| (deg/s)" stats={decStats.macGyrX} />
                          <StatRow label="|Gyr-Y| (deg/s)" stats={decStats.macGyrY} />
                          <StatRow label="|Gyr-Z| (deg/s)" stats={decStats.macGyrZ} />
                        </TableBody>
                      </Table>
                    </div>
                  ) : (
                    <div className="text-center py-8">
                      <p className="text-muted-foreground mb-2">No deceleration phase data</p>
                      <p className="text-xs text-muted-foreground">Upload the Deceleration Phase sheet exported from Excel</p>
                    </div>
                  )}
                </CardContent>
              </Card>
            </TabsContent>

            <TabsContent value="stroke">
              <Card className="border-border/50">
                <CardHeader>
                  <CardTitle>Stroke Phase Statistics</CardTitle>
                  <CardDescription>Upload stroke_phase.csv for detailed metrics</CardDescription>
                </CardHeader>
                <CardContent>
                  {spStats ? (
                    <div className="overflow-auto">
                      <Table>
                        <TableHeader>
                          <TableRow>
                            <TableHead>Parameter</TableHead>
                            <TableHead className="text-right">Mean</TableHead>
                            <TableHead className="text-right">Std Dev</TableHead>
                            <TableHead className="text-right">Min</TableHead>
                            <TableHead className="text-right">Max</TableHead>
                          </TableRow>
                        </TableHeader>
                        <TableBody>
                          <StatRow label="Gyr-X+ Peak (deg/s)" stats={spStats.gyrXPos} />
                          <StatRow label="Gyr-X- Peak (deg/s)" stats={spStats.gyrXNeg} />
                          <StatRow label="|Acc-X| (m/s²)" stats={spStats.macAccX} />
                          <StatRow label="|Gyr-X| (deg/s)" stats={spStats.macGyrX} />
                          <StatRow label="|Gyr-Y| (deg/s)" stats={spStats.macGyrY} />
                          <StatRow label="|Gyr-Z| (deg/s)" stats={spStats.macGyrZ} />
                        </TableBody>
                      </Table>
                    </div>
                  ) : (
                    <div className="text-center py-8">
                      <p className="text-muted-foreground mb-2">No stroke phase data</p>
                      <p className="text-xs text-muted-foreground">Upload the Stroke Phase sheet exported from Excel</p>
                    </div>
                  )}
                </CardContent>
              </Card>
            </TabsContent>
          </Tabs>
        </div>
      </main>
    </div>
  )
}

function UploadButton({ label, onUpload }: { label: string; onUpload: (file: File) => void }) {
  return (
    <label className="flex cursor-pointer flex-col items-center gap-2 rounded-lg border border-dashed border-border/50 p-3 transition-colors hover:border-primary hover:bg-primary/5">
      <Upload className="size-4 text-muted-foreground" />
      <span className="text-xs text-center text-muted-foreground">{label}</span>
      <input
        type="file"
        accept=".csv"
        className="sr-only"
        onChange={(e) => {
          const file = e.target.files?.[0]
          if (file) onUpload(file)
        }}
      />
    </label>
  )
}

function StatRow({ label, stats }: { label: string; stats: { mean: number; std: number; min: number; max: number } }) {
  return (
    <TableRow>
      <TableCell className="font-medium">{label}</TableCell>
      <TableCell className="text-right font-mono">{stats.mean.toFixed(4)}</TableCell>
      <TableCell className="text-right font-mono">{stats.std.toFixed(4)}</TableCell>
      <TableCell className="text-right font-mono">{stats.min.toFixed(4)}</TableCell>
      <TableCell className="text-right font-mono">{stats.max.toFixed(4)}</TableCell>
    </TableRow>
  )
}
