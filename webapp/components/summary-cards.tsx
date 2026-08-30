"use client"

import { Card, CardContent } from "@/components/ui/card"
import type { SummaryStats } from "@/lib/types"
import { 
  Activity, 
  Gauge, 
  TrendingUp, 
  TrendingDown, 
  Zap,
  RotateCw
} from "lucide-react"

interface SummaryCardsProps {
  stats: SummaryStats
}

export function SummaryCards({ stats }: SummaryCardsProps) {
  // Primary metrics row
  const primaryCards = [
    {
      label: "Total Strokes",
      value: stats.totalStrokes.toString(),
      subtitle: "Detected cycles",
      icon: Activity,
      accentClass: "text-primary",
    },
    {
      label: "Avg SPM",
      value: stats.avgSPM.toFixed(1),
      subtitle: "Strokes per minute",
      icon: Gauge,
      accentClass: "text-chart-1",
    },
    {
      label: "Avg ΔVel (Acc)",
      value: `${stats.avgVelocityAcc.toFixed(3)}`,
      subtitle: "m/s in acc phase",
      icon: TrendingUp,
      accentClass: "text-chart-5",
    },
    {
      label: "Avg ΔVel (Dec)",
      value: `${stats.avgVelocityDec.toFixed(3)}`,
      subtitle: "m/s in dec phase",
      icon: TrendingDown,
      accentClass: "text-chart-2",
    },
  ]

  // Secondary metrics row (MAC - Mean Absolute Change)
  const secondaryCards = [
    {
      label: "MAC Acc-X",
      value: stats.avgMACAccX.toFixed(3),
      subtitle: "m/s²",
      icon: Zap,
      accentClass: "text-primary",
    },
    {
      label: "MAC Gyr-X",
      value: stats.avgMACGyrX.toFixed(2),
      subtitle: "deg/s",
      icon: RotateCw,
      accentClass: "text-chart-1",
    },
    {
      label: "MAC Gyr-Y",
      value: stats.avgMACGyrY.toFixed(2),
      subtitle: "deg/s",
      icon: RotateCw,
      accentClass: "text-chart-5",
    },
    {
      label: "MAC Gyr-Z",
      value: stats.avgMACGyrZ.toFixed(2),
      subtitle: "deg/s",
      icon: RotateCw,
      accentClass: "text-chart-2",
    },
  ]

  return (
    <div className="flex flex-col gap-4">
      {/* Primary metrics */}
      <div className="grid grid-cols-2 gap-3 sm:gap-4 lg:grid-cols-4">
        {primaryCards.map((card) => (
          <Card key={card.label} className="border-border/50">
            <CardContent className="flex items-center gap-3 p-3 sm:gap-4 sm:p-4">
              <div className={`shrink-0 rounded-lg bg-secondary p-2 sm:p-3 ${card.accentClass}`}>
                <card.icon className="size-4 sm:size-5" />
              </div>
              <div className="flex min-w-0 flex-col gap-0.5">
                <span className="truncate text-xs text-muted-foreground">{card.label}</span>
                <span className="truncate font-mono text-lg font-semibold tracking-tight sm:text-2xl">
                  {card.value}
                </span>
                <span className="hidden text-xs text-muted-foreground sm:block">{card.subtitle}</span>
              </div>
            </CardContent>
          </Card>
        ))}
      </div>

      {/* Secondary metrics - MAC values */}
      <div className="grid grid-cols-2 gap-3 sm:gap-4 lg:grid-cols-4">
        {secondaryCards.map((card) => (
          <Card key={card.label} className="border-border/50 bg-secondary/30">
            <CardContent className="flex items-center gap-3 p-3 sm:gap-4 sm:p-4">
              <div className={`shrink-0 rounded-lg bg-secondary p-2 ${card.accentClass}`}>
                <card.icon className="size-4" />
              </div>
              <div className="flex min-w-0 flex-col gap-0.5">
                <span className="truncate text-xs text-muted-foreground">{card.label}</span>
                <span className="truncate font-mono text-base font-semibold tracking-tight sm:text-xl">
                  {card.value}
                </span>
                <span className="text-xs text-muted-foreground">{card.subtitle}</span>
              </div>
            </CardContent>
          </Card>
        ))}
      </div>
    </div>
  )
}
