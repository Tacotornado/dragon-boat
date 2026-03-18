"use client"

import { Button } from "@/components/ui/button"
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select"
import { Badge } from "@/components/ui/badge"
import { Activity, RotateCcw } from "lucide-react"

interface DashboardHeaderProps {
  filterTopN: number | null
  onFilterChange: (v: number | null) => void
  totalStrokes: number
  dataSource: "demo" | "uploaded"
  onResetToDemo: () => void
}

export function DashboardHeader({
  filterTopN,
  onFilterChange,
  totalStrokes,
  dataSource,
  onResetToDemo,
}: DashboardHeaderProps) {
  return (
    <header className="sticky top-0 z-30 border-b border-border/50 bg-card/80 backdrop-blur-md">
      <div className="mx-auto flex max-w-[1440px] flex-col gap-3 px-4 py-3 sm:flex-row sm:items-center sm:justify-between md:px-6">
        <div className="flex items-center gap-3">
          <div className="flex size-9 shrink-0 items-center justify-center rounded-lg bg-primary text-primary-foreground">
            <Activity className="size-5" />
          </div>
          <div className="min-w-0">
            <h1 className="truncate text-base font-semibold leading-tight text-foreground sm:text-lg">
              Stroke Analysis
            </h1>
            <p className="hidden text-xs text-muted-foreground sm:block">
              IMU Performance Dashboard
            </p>
          </div>
          <Badge
            variant={dataSource === "uploaded" ? "default" : "secondary"}
            className="shrink-0 text-xs"
          >
            {dataSource === "uploaded" ? "Live" : "Demo"}
          </Badge>
        </div>

        <div className="flex items-center gap-2 sm:gap-3">
          <Select
            value={filterTopN?.toString() ?? "all"}
            onValueChange={(v) =>
              onFilterChange(v === "all" ? null : parseInt(v))
            }
          >
            <SelectTrigger className="h-9 w-full min-w-[120px] flex-1 border-border/50 bg-secondary text-xs sm:w-[160px] sm:flex-none">
              <SelectValue placeholder="Filter strokes" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="all">All {totalStrokes} strokes</SelectItem>
              <SelectItem value="5">Top 5 strokes</SelectItem>
              <SelectItem value="10">Top 10 strokes</SelectItem>
              <SelectItem value="15">Top 15 strokes</SelectItem>
            </SelectContent>
          </Select>

          {dataSource === "uploaded" && (
            <Button
              variant="outline"
              size="sm"
              onClick={onResetToDemo}
              className="shrink-0 gap-1.5 text-xs sm:gap-2"
            >
              <RotateCcw className="size-3.5" />
              <span className="hidden sm:inline">Reset to Demo</span>
              <span className="sm:hidden">Reset</span>
            </Button>
          )}
        </div>
      </div>
    </header>
  )
}
