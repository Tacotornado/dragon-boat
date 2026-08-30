"use client"

import { useEffect, useRef } from "react"
import { Badge } from "@/components/ui/badge"
import type { StrokeCycle } from "@/lib/types"

interface LiveStrokeListProps {
  strokes: StrokeCycle[]
}

export function LiveStrokeList({ strokes }: LiveStrokeListProps) {
  const scrollRef = useRef<HTMLDivElement>(null)
  
  // Auto-scroll to bottom when new strokes arrive
  useEffect(() => {
    if (scrollRef.current) {
      scrollRef.current.scrollTop = scrollRef.current.scrollHeight
    }
  }, [strokes.length])
  
  if (strokes.length === 0) {
    return (
      <div className="flex h-[200px] items-center justify-center rounded-lg border border-dashed border-border/50 bg-secondary/20">
        <p className="text-sm text-muted-foreground">
          No strokes detected yet...
        </p>
      </div>
    )
  }
  
  // Show last 20 strokes
  const recentStrokes = strokes.slice(-20).reverse()
  
  return (
    <div 
      ref={scrollRef}
      className="max-h-[300px] overflow-auto rounded-lg border border-border/50"
    >
      <table className="w-full text-sm">
        <thead className="sticky top-0 bg-card">
          <tr className="border-b border-border/50 text-left text-xs text-muted-foreground">
            <th className="px-3 py-2 font-medium">ID</th>
            <th className="px-3 py-2 font-medium">Duration</th>
            <th className="px-3 py-2 font-medium">Area+</th>
            <th className="px-3 py-2 font-medium">Area-</th>
            <th className="px-3 py-2 font-medium">Balance</th>
          </tr>
        </thead>
        <tbody>
          {recentStrokes.map((stroke, i) => {
            const balance = stroke.Area_Positive / (stroke.Area_Positive + stroke.Area_Negative)
            const isRecent = i === 0
            
            return (
              <tr 
                key={stroke.StrokeID}
                className={`border-b border-border/30 transition-colors ${
                  isRecent ? "bg-primary/5" : "hover:bg-secondary/50"
                }`}
              >
                <td className="px-3 py-2">
                  <span className="font-mono font-medium">
                    #{stroke.StrokeID}
                  </span>
                  {isRecent && (
                    <Badge 
                      variant="secondary" 
                      className="ml-2 bg-green-500/20 text-green-400 text-[10px]"
                    >
                      NEW
                    </Badge>
                  )}
                </td>
                <td className="px-3 py-2 font-mono">
                  {stroke.TotalTime.toFixed(3)}s
                </td>
                <td className="px-3 py-2 font-mono text-green-400">
                  {stroke.Area_Positive.toFixed(2)}
                </td>
                <td className="px-3 py-2 font-mono text-red-400">
                  {stroke.Area_Negative.toFixed(2)}
                </td>
                <td className="px-3 py-2">
                  <div className="flex items-center gap-2">
                    <div className="h-1.5 w-16 overflow-hidden rounded-full bg-secondary">
                      <div 
                        className="h-full bg-primary transition-all"
                        style={{ width: `${balance * 100}%` }}
                      />
                    </div>
                    <span className="font-mono text-xs text-muted-foreground">
                      {(balance * 100).toFixed(0)}%
                    </span>
                  </div>
                </td>
              </tr>
            )
          })}
        </tbody>
      </table>
    </div>
  )
}
