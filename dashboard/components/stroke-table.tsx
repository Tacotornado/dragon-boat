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
import type { StrokeCycle } from "@/lib/types"
import { ArrowUpDown, ChevronDown, ChevronUp } from "lucide-react"

interface StrokeTableProps {
  strokes: StrokeCycle[]
  selectedStrokeId: number | null
  onSelectStroke: (id: number) => void
}

type SortField = keyof StrokeCycle
type SortDirection = "asc" | "desc"

export function StrokeTable({ strokes, selectedStrokeId, onSelectStroke }: StrokeTableProps) {
  const [sortField, setSortField] = useState<SortField>("StrokeID")
  const [sortDir, setSortDir] = useState<SortDirection>("asc")

  const toggleSort = (field: SortField) => {
    if (sortField === field) {
      setSortDir(sortDir === "asc" ? "desc" : "asc")
    } else {
      setSortField(field)
      setSortDir("asc")
    }
  }

  const sorted = useMemo(() => {
    return [...strokes].sort((a, b) => {
      const aVal = a[sortField]
      const bVal = b[sortField]
      if (typeof aVal === "number" && typeof bVal === "number") {
        return sortDir === "asc" ? aVal - bVal : bVal - aVal
      }
      return 0
    })
  }, [strokes, sortField, sortDir])

  const SortIcon = ({ field }: { field: SortField }) => {
    if (sortField !== field) return <ArrowUpDown className="ml-1 size-3 text-muted-foreground" />
    return sortDir === "asc" ? (
      <ChevronUp className="ml-1 size-3 text-primary" />
    ) : (
      <ChevronDown className="ml-1 size-3 text-primary" />
    )
  }

  const columns: { key: SortField; label: string; format?: (v: number) => string }[] = [
    { key: "StrokeID", label: "Stroke" },
    { key: "TotalTime", label: "Duration (s)", format: (v) => v.toFixed(3) },
    { key: "TimeToPosPeak", label: "T+ Peak (s)", format: (v) => v.toFixed(3) },
    { key: "TimeToNegPeak", label: "T- Peak (s)", format: (v) => v.toFixed(3) },
    { key: "Area_Positive", label: "Area +", format: (v) => v.toFixed(2) },
    { key: "Area_Negative", label: "Area -", format: (v) => v.toFixed(2) },
  ]

  // Compute percentile thresholds for highlighting
  const posAreas = strokes.map((s) => s.Area_Positive).sort((a, b) => b - a)
  const top20Threshold = posAreas[Math.floor(posAreas.length * 0.2)] ?? 0

  return (
    <Card className="border-border/50">
      <CardHeader>
        <CardTitle className="text-lg">Stroke Metrics</CardTitle>
        <CardDescription>
          Click any row to highlight the corresponding stroke visualization. Sortable columns.
        </CardDescription>
      </CardHeader>
      <CardContent className="px-2 sm:px-6">
        <div className="max-h-[420px] overflow-auto rounded-lg border border-border/50 -mx-2 sm:mx-0">
          <Table>
            <TableHeader className="sticky top-0 z-10 bg-card">
              <TableRow className="border-border/50 hover:bg-transparent">
                {columns.map((col) => (
                  <TableHead
                    key={col.key}
                    className="cursor-pointer select-none whitespace-nowrap px-2 text-xs font-medium text-muted-foreground transition-colors hover:text-foreground sm:px-4"
                    onClick={() => toggleSort(col.key)}
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
                    {columns.map((col) => (
                      <TableCell
                        key={col.key}
                        className={`whitespace-nowrap px-2 font-mono text-xs sm:px-4 sm:text-sm ${
                          col.key === "StrokeID" ? "font-semibold" : ""
                        }`}
                      >
                        {col.format
                          ? col.format(stroke[col.key] as number)
                          : stroke[col.key]}
                      </TableCell>
                    ))}
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
