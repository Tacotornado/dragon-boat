"use client"

import { Card, CardContent } from "@/components/ui/card"
import type { StrokeCycle } from "@/lib/types"
import { Activity, Timer, TrendingUp, TrendingDown } from "lucide-react"

interface SummaryCardsProps {
  strokes: StrokeCycle[]
}

export function SummaryCards({ strokes }: SummaryCardsProps) {
  const totalStrokes = strokes.length
  const avgDuration =
    strokes.reduce((sum, s) => sum + s.TotalTime, 0) / totalStrokes
  const maxPosPeak = Math.max(...strokes.map((s) => s.Area_Positive))
  const maxNegPeak = Math.max(...strokes.map((s) => s.Area_Negative))

  const cards = [
    {
      label: "Total Strokes",
      value: totalStrokes.toString(),
      subtitle: "Detected cycles",
      icon: Activity,
      accentClass: "text-primary",
    },
    {
      label: "Avg Duration",
      value: `${avgDuration.toFixed(3)}s`,
      subtitle: "Per stroke cycle",
      icon: Timer,
      accentClass: "text-chart-1",
    },
    {
      label: "Max Positive Area",
      value: maxPosPeak.toFixed(2),
      subtitle: "Acceleration impulse",
      icon: TrendingUp,
      accentClass: "text-chart-5",
    },
    {
      label: "Max Negative Area",
      value: maxNegPeak.toFixed(2),
      subtitle: "Deceleration impulse",
      icon: TrendingDown,
      accentClass: "text-chart-2",
    },
  ]

  return (
    <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-4">
      {cards.map((card) => (
        <Card key={card.label} className="border-border/50">
          <CardContent className="flex items-center gap-4 pt-0">
            <div className={`rounded-lg bg-secondary p-3 ${card.accentClass}`}>
              <card.icon className="size-5" />
            </div>
            <div className="flex flex-col gap-0.5">
              <span className="text-xs text-muted-foreground">{card.label}</span>
              <span className="text-2xl font-semibold tracking-tight font-mono">
                {card.value}
              </span>
              <span className="text-xs text-muted-foreground">{card.subtitle}</span>
            </div>
          </CardContent>
        </Card>
      ))}
    </div>
  )
}
