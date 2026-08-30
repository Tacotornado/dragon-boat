"use client"

import { useState, useMemo } from "react"
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table"
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Button } from "@/components/ui/button"
import {
  DropdownMenu,
  DropdownMenuCheckboxItem,
  DropdownMenuContent,
  DropdownMenuLabel,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu"
import type { StrokeCycle, AccPhaseMetrics, DecPhaseMetrics, StrokePhaseMetrics } from "@/lib/types"
import { ArrowUpDown, ChevronDown, ChevronUp, Settings2 } from "lucide-react"

interface StrokeTableProps {
  strokes: StrokeCycle[]
  accMetrics?: AccPhaseMetrics[]
  decMetrics?: DecPhaseMetrics[]
  strokePhaseMetrics?: StrokePhaseMetrics[]
  selectedStrokeId: number | null
  onSelectStroke: (id: number) => void
  onParameterSelect?: (params: string[]) => void
  selectedParams?: string[]
}

type SortDirection = "asc" | "desc"

interface ColumnDef {
  key: string
  label: string
  category: "basic" | "acceleration" | "deceleration" | "stroke"
  format?: (v: number) => string
  getValue: (stroke: StrokeCycle, acc?: AccPhaseMetrics, dec?: DecPhaseMetrics, sp?: StrokePhaseMetrics) => number | string
}

const ALL_COLUMNS: ColumnDef[] = [
  // Basic
  { key: "StrokeID", label: "Stroke", category: "basic", getValue: (s) => s.StrokeID },
  { key: "TotalTime", label: "Duration (s)", category: "basic", format: (v) => v.toFixed(3), getValue: (s) => s.TotalTime },
  { key: "TimeToPosPeak", label: "T+ Peak (s)", category: "basic", format: (v) => v.toFixed(3), getValue: (s) => s.TimeToPosPeak },
  { key: "TimeToNegPeak", label: "T- Peak (s)", category: "basic", format: (v) => v.toFixed(3), getValue: (s) => s.TimeToNegPeak },
  { key: "Area_Positive", label: "Area +", category: "basic", format: (v) => v.toFixed(2), getValue: (s) => s.Area_Positive },
  { key: "Area_Negative", label: "Area -", category: "basic", format: (v) => v.toFixed(2), getValue: (s) => s.Area_Negative },
  // Acceleration phase
  { key: "acc_Duration", label: "Acc Duration", category: "acceleration", format: (v) => v.toFixed(3), getValue: (_, acc) => acc?.Duration ?? 0 },
  { key: "acc_VelocityChange", label: "ΔVel Acc", category: "acceleration", format: (v) => v.toFixed(3), getValue: (_, acc) => acc?.VelocityChange ?? 0 },
  { key: "acc_AccY_PosPeak", label: "Acc-Y+ Peak", category: "acceleration", format: (v) => v.toFixed(3), getValue: (_, acc) => acc?.AccY_PosPeak ?? 0 },
  { key: "acc_RAD", label: "RAD", category: "acceleration", format: (v) => v.toFixed(3), getValue: (_, acc) => acc?.RAD ?? 0 },
  { key: "acc_MeanAbsAccX", label: "|Acc-X|", category: "acceleration", format: (v) => v.toFixed(3), getValue: (_, acc) => acc?.MeanAbsAccX ?? 0 },
  { key: "acc_MeanAbsGyrX", label: "|Gyr-X|", category: "acceleration", format: (v) => v.toFixed(2), getValue: (_, acc) => acc?.MeanAbsGyrX ?? 0 },
  // Deceleration phase
  { key: "dec_Duration", label: "Dec Duration", category: "deceleration", format: (v) => v.toFixed(3), getValue: (_, __, dec) => dec?.Duration ?? 0 },
  { key: "dec_VelocityChange", label: "ΔVel Dec", category: "deceleration", format: (v) => v.toFixed(3), getValue: (_, __, dec) => dec?.VelocityChange ?? 0 },
  { key: "dec_AccY_NegPeak", label: "Acc-Y- Peak", category: "deceleration", format: (v) => v.toFixed(3), getValue: (_, __, dec) => dec?.AccY_NegPeak ?? 0 },
  { key: "dec_MeanAbsAccX", label: "|Acc-X|", category: "deceleration", format: (v) => v.toFixed(3), getValue: (_, __, dec) => dec?.MeanAbsAccX ?? 0 },
  // Stroke phase
  { key: "stroke_GyrX_PosPeak", label: "Gyr-X+ Peak", category: "stroke", format: (v) => v.toFixed(2), getValue: (_, __, ___, sp) => sp?.GyrX_PosPeak ?? 0 },
  { key: "stroke_GyrX_NegPeak", label: "Gyr-X- Peak", category: "stroke", format: (v) => v.toFixed(2), getValue: (_, __, ___, sp) => sp?.GyrX_NegPeak ?? 0 },
  { key: "stroke_MeanAbsAccX", label: "MAC Acc-X", category: "stroke", format: (v) => v.toFixed(3), getValue: (_, __, ___, sp) => sp?.MeanAbsAccX ?? 0 },
  { key: "stroke_MeanAbsGyrX", label: "MAC Gyr-X", category: "stroke", format: (v) => v.toFixed(2), getValue: (_, __, ___, sp) => sp?.MeanAbsGyrX ?? 0 },
  { key: "stroke_MeanAbsGyrY", label: "MAC Gyr-Y", category: "stroke", format: (v) => v.toFixed(2), getValue: (_, __, ___, sp) => sp?.MeanAbsGyrY ?? 0 },
  { key: "stroke_MeanAbsGyrZ", label: "MAC Gyr-Z", category: "stroke", format: (v) => v.toFixed(2), getValue: (_, __, ___, sp) => sp?.MeanAbsGyrZ ?? 0 },
]

const DEFAULT_COLUMNS = ["StrokeID", "TotalTime", "Area_Positive", "Area_Negative", "acc_VelocityChange", "dec_VelocityChange"]

export function StrokeTable({ 
  strokes, 
  accMetrics, 
  decMetrics, 
  strokePhaseMetrics,
  selectedStrokeId, 
  onSelectStroke,
  onParameterSelect,
  selectedParams,
}: StrokeTableProps) {
  const [sortField, setSortField] = useState<string>("StrokeID")
  const [sortDir, setSortDir] = useState<SortDirection>("asc")
  const [visibleColumns, setVisibleColumns] = useState<string[]>(DEFAULT_COLUMNS)

  const toggleSort = (field: string) => {
    if (sortField === field) {
      setSortDir(sortDir === "asc" ? "desc" : "asc")
    } else {
      setSortField(field)
      setSortDir("asc")
    }
  }

  const toggleColumn = (key: string) => {
    setVisibleColumns(prev => {
      if (prev.includes(key)) {
        return prev.filter(k => k !== key)
      }
      return [...prev, key]
    })
  }

  const activeColumns = useMemo(() => {
    return ALL_COLUMNS.filter(col => visibleColumns.includes(col.key))
  }, [visibleColumns])

  const getMetricsForStroke = (strokeId: number) => {
    return {
      acc: accMetrics?.find(m => m.StrokeID === strokeId),
      dec: decMetrics?.find(m => m.StrokeID === strokeId),
      sp: strokePhaseMetrics?.find(m => m.StrokeID === strokeId),
    }
  }

  const sorted = useMemo(() => {
    const col = ALL_COLUMNS.find(c => c.key === sortField)
    if (!col) return strokes

    return [...strokes].sort((a, b) => {
      const metricsA = getMetricsForStroke(a.StrokeID)
      const metricsB = getMetricsForStroke(b.StrokeID)
      const aVal = col.getValue(a, metricsA.acc, metricsA.dec, metricsA.sp)
      const bVal = col.getValue(b, metricsB.acc, metricsB.dec, metricsB.sp)
      if (typeof aVal === "number" && typeof bVal === "number") {
        return sortDir === "asc" ? aVal - bVal : bVal - aVal
      }
      return 0
    })
  }, [strokes, sortField, sortDir, accMetrics, decMetrics, strokePhaseMetrics])

  const SortIcon = ({ field }: { field: string }) => {
    if (sortField !== field) return <ArrowUpDown className="ml-1 size-3 text-muted-foreground" />
    return sortDir === "asc" ? (
      <ChevronUp className="ml-1 size-3 text-primary" />
    ) : (
      <ChevronDown className="ml-1 size-3 text-primary" />
    )
  }

  // Compute percentile thresholds for highlighting
  const posAreas = strokes.map((s) => s.Area_Positive).sort((a, b) => b - a)
  const top20Threshold = posAreas[Math.floor(posAreas.length * 0.2)] ?? 0

  // Handle parameter selection for charts
  const handleParamClick = (key: string) => {
    if (!onParameterSelect || !selectedParams) return
    if (selectedParams.includes(key)) {
      onParameterSelect(selectedParams.filter(p => p !== key))
    } else if (selectedParams.length < 2) {
      onParameterSelect([...selectedParams, key])
    } else {
      onParameterSelect([selectedParams[1], key])
    }
  }

  return (
    <Card className="border-border/50">
      <CardHeader>
        <div className="flex flex-col gap-2 sm:flex-row sm:items-center sm:justify-between">
          <div>
            <CardTitle className="text-lg">Stroke Metrics</CardTitle>
            <CardDescription>
              Click columns to sort. Click rows to highlight stroke. 
              {onParameterSelect && " Click column headers to add to chart."}
            </CardDescription>
          </div>
          <DropdownMenu>
            <DropdownMenuTrigger asChild>
              <Button variant="outline" size="sm" className="gap-2">
                <Settings2 className="size-4" />
                Columns
              </Button>
            </DropdownMenuTrigger>
            <DropdownMenuContent align="end" className="w-56 max-h-[400px] overflow-auto">
              <DropdownMenuLabel>Basic Metrics</DropdownMenuLabel>
              {ALL_COLUMNS.filter(c => c.category === "basic").map(col => (
                <DropdownMenuCheckboxItem
                  key={col.key}
                  checked={visibleColumns.includes(col.key)}
                  onCheckedChange={() => toggleColumn(col.key)}
                >
                  {col.label}
                </DropdownMenuCheckboxItem>
              ))}
              <DropdownMenuSeparator />
              <DropdownMenuLabel>Acceleration Phase</DropdownMenuLabel>
              {ALL_COLUMNS.filter(c => c.category === "acceleration").map(col => (
                <DropdownMenuCheckboxItem
                  key={col.key}
                  checked={visibleColumns.includes(col.key)}
                  onCheckedChange={() => toggleColumn(col.key)}
                >
                  {col.label}
                </DropdownMenuCheckboxItem>
              ))}
              <DropdownMenuSeparator />
              <DropdownMenuLabel>Deceleration Phase</DropdownMenuLabel>
              {ALL_COLUMNS.filter(c => c.category === "deceleration").map(col => (
                <DropdownMenuCheckboxItem
                  key={col.key}
                  checked={visibleColumns.includes(col.key)}
                  onCheckedChange={() => toggleColumn(col.key)}
                >
                  {col.label}
                </DropdownMenuCheckboxItem>
              ))}
              <DropdownMenuSeparator />
              <DropdownMenuLabel>Stroke Phase</DropdownMenuLabel>
              {ALL_COLUMNS.filter(c => c.category === "stroke").map(col => (
                <DropdownMenuCheckboxItem
                  key={col.key}
                  checked={visibleColumns.includes(col.key)}
                  onCheckedChange={() => toggleColumn(col.key)}
                >
                  {col.label}
                </DropdownMenuCheckboxItem>
              ))}
            </DropdownMenuContent>
          </DropdownMenu>
        </div>
      </CardHeader>
      <CardContent className="px-2 sm:px-6">
        <div className="max-h-[420px] overflow-auto rounded-lg border border-border/50 -mx-2 sm:mx-0">
          <Table>
            <TableHeader className="sticky top-0 z-10 bg-card">
              <TableRow className="border-border/50 hover:bg-transparent">
                {activeColumns.map((col) => (
                  <TableHead
                    key={col.key}
                    className={`cursor-pointer select-none whitespace-nowrap px-2 text-xs font-medium transition-colors hover:text-foreground sm:px-4 ${
                      selectedParams?.includes(col.key) 
                        ? "text-primary bg-primary/10" 
                        : "text-muted-foreground"
                    }`}
                    onClick={(e) => {
                      if (e.shiftKey && onParameterSelect) {
                        handleParamClick(col.key)
                      } else {
                        toggleSort(col.key)
                      }
                    }}
                  >
                    <span className="flex items-center">
                      {col.label}
                      <SortIcon field={col.key} />
                    </span>
                  </TableHead>
                ))}
                <TableHead className="whitespace-nowrap px-2 text-xs font-medium text-muted-foreground sm:px-4">
                  Rating
                </TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {sorted.map((stroke) => {
                const metrics = getMetricsForStroke(stroke.StrokeID)
                const isSelected = selectedStrokeId === stroke.StrokeID
                const isTopStroke = stroke.Area_Positive >= top20Threshold
                return (
                  <TableRow
                    key={stroke.StrokeID}
                    className={`cursor-pointer border-border/30 transition-colors ${
                      isSelected
                        ? "bg-primary/10 hover:bg-primary/15"
                        : "hover:bg-secondary/50"
                    }`}
                    onClick={() => onSelectStroke(stroke.StrokeID)}
                  >
                    {activeColumns.map((col) => {
                      const value = col.getValue(stroke, metrics.acc, metrics.dec, metrics.sp)
                      return (
                        <TableCell
                          key={col.key}
                          className={`whitespace-nowrap px-2 font-mono text-xs sm:px-4 sm:text-sm ${
                            col.key === "StrokeID" ? "font-semibold" : ""
                          }`}
                        >
                          {typeof value === "number" && col.format
                            ? col.format(value)
                            : value}
                        </TableCell>
                      )
                    })}
                    <TableCell className="px-2 sm:px-4">
                      {isTopStroke ? (
                        <Badge className="bg-chart-5/15 text-chart-5 border-chart-5/30 text-xs">
                          Top 20%
                        </Badge>
                      ) : (
                        <span className="text-xs text-muted-foreground">--</span>
                      )}
                    </TableCell>
                  </TableRow>
                )
              })}
            </TableBody>
          </Table>
        </div>
      </CardContent>
    </Card>
  )
}
