"use client"

import { useState } from "react"
import { Button } from "@/components/ui/button"
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select"
import { Activity, Play, Loader2 } from "lucide-react"

interface DashboardHeaderProps {
  onRunAnalysis: () => Promise<{ success: boolean; message: string }>
  filterTopN: number | null
  onFilterChange: (v: number | null) => void
  totalStrokes: number
}

export function DashboardHeader({
  onRunAnalysis,
  filterTopN,
  onFilterChange,
  totalStrokes,
}: DashboardHeaderProps) {
  const [isRunning, setIsRunning] = useState(false)
  const [status, setStatus] = useState<string | null>(null)

  const handleRun = async () => {
    setIsRunning(true)
    setStatus(null)
    try {
      const result = await onRunAnalysis()
      setStatus(result.message)
      setTimeout(() => setStatus(null), 3000)
    } catch {
      setStatus("Analysis failed")
      setTimeout(() => setStatus(null), 3000)
    } finally {
      setIsRunning(false)
    }
  }

  return (
    <header className="sticky top-0 z-30 border-b border-border/50 bg-card/80 backdrop-blur-md">
      <div className="mx-auto flex max-w-[1440px] items-center justify-between px-4 py-3 md:px-6">
        <div className="flex items-center gap-3">
          <div className="flex size-9 items-center justify-center rounded-lg bg-primary text-primary-foreground">
            <Activity className="size-5" />
          </div>
          <div>
            <h1 className="text-lg font-semibold leading-tight text-foreground">
              Stroke Analysis
            </h1>
            <p className="text-xs text-muted-foreground">
              IMU Performance Dashboard
            </p>
          </div>
        </div>

        <div className="flex items-center gap-3">
          {status && (
            <span className="text-xs text-primary animate-in fade-in">
              {status}
            </span>
          )}

          <Select
            value={filterTopN?.toString() ?? "all"}
            onValueChange={(v) =>
              onFilterChange(v === "all" ? null : parseInt(v))
            }
          >
            <SelectTrigger className="w-[160px] h-9 text-xs bg-secondary border-border/50">
              <SelectValue placeholder="Filter strokes" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="all">All {totalStrokes} strokes</SelectItem>
              <SelectItem value="5">Top 5 strokes</SelectItem>
              <SelectItem value="10">Top 10 strokes</SelectItem>
              <SelectItem value="15">Top 15 strokes</SelectItem>
            </SelectContent>
          </Select>

          <Button
            size="sm"
            onClick={handleRun}
            disabled={isRunning}
            className="gap-2"
          >
            {isRunning ? (
              <Loader2 className="size-4 animate-spin" />
            ) : (
              <Play className="size-4" />
            )}
            {isRunning ? "Running..." : "Run Analysis"}
          </Button>
        </div>
      </div>
    </header>
  )
}
