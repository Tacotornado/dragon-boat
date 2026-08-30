"use client"

import { use } from "react"
import Link from "next/link"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Skeleton } from "@/components/ui/skeleton"
import { 
  ArrowLeft, 
  Radio,
  Wifi,
  WifiOff,
  Activity,
  Clock,
  Waves,
  Gauge,
  Zap
} from "lucide-react"
import { useLiveSession } from "@/hooks/use-live-session"
import { LiveWaveform } from "@/components/live-waveform"
import { LiveStrokeList } from "@/components/live-stroke-list"

export default function LiveSessionPage({ 
  params 
}: { 
  params: Promise<{ id: string }> 
}) {
  const { id } = use(params)
  const {
    session,
    isConnected,
    isLoading,
    error,
    recentSamples,
    strokes,
    stats,
    sampleCount,
  } = useLiveSession(id)
  
  if (isLoading) {
    return (
      <div className="min-h-screen bg-background p-4">
        <Skeleton className="mb-4 h-12 w-full" />
        <div className="grid gap-4 md:grid-cols-2">
          <Skeleton className="h-[300px]" />
          <Skeleton className="h-[300px]" />
        </div>
      </div>
    )
  }
  
  if (error || !session) {
    return (
      <div className="flex min-h-screen items-center justify-center bg-background">
        <Card className="max-w-md">
          <CardContent className="p-6 text-center">
            <p className="text-destructive">{error || "Session not found"}</p>
            <Link href="/sessions">
              <Button variant="outline" className="mt-4">
                Back to Sessions
              </Button>
            </Link>
          </CardContent>
        </Card>
      </div>
    )
  }
  
  const duration = sampleCount / 120 // seconds at 120Hz
  const spm = strokes.length > 0 && duration > 0 
    ? (strokes.length / duration) * 60 
    : 0
  
  return (
    <div className="min-h-screen bg-background">
      {/* Header */}
      <header className="sticky top-0 z-50 border-b border-border/50 bg-background/95 backdrop-blur">
        <div className="mx-auto flex max-w-[1440px] items-center justify-between px-4 py-3 md:px-6">
          <div className="flex items-center gap-3">
            <Link href="/sessions" className="text-muted-foreground hover:text-foreground">
              <ArrowLeft className="size-5" />
            </Link>
            <div className="flex size-9 items-center justify-center rounded-lg bg-primary text-primary-foreground">
              <Radio className="size-5" />
            </div>
            <div>
              <h1 className="text-lg font-semibold leading-tight text-foreground">
                {session.name}
              </h1>
              <p className="text-xs text-muted-foreground">
                {session.athleteName || "Live Session"}
              </p>
            </div>
          </div>
          
          <div className="flex items-center gap-3">
            <Badge 
              variant={isConnected ? "default" : "secondary"}
              className={isConnected ? "bg-green-500/20 text-green-400" : ""}
            >
              {isConnected ? (
                <><Wifi className="mr-1.5 size-3" /> Connected</>
              ) : (
                <><WifiOff className="mr-1.5 size-3" /> Disconnected</>
              )}
            </Badge>
            <Badge 
              className={
                session.status === "recording" 
                  ? "bg-red-500/20 text-red-400 animate-pulse" 
                  : "bg-muted text-muted-foreground"
              }
            >
              {session.status === "recording" && (
                <span className="mr-1.5 inline-block size-2 rounded-full bg-red-500 animate-pulse" />
              )}
              {session.status}
            </Badge>
          </div>
        </div>
      </header>
      
      <main className="mx-auto max-w-[1440px] p-4 md:p-6">
        {/* Stats Row */}
        <div className="mb-6 grid grid-cols-2 gap-3 sm:grid-cols-4 lg:grid-cols-6">
          <StatCard
            icon={<Clock className="size-4" />}
            label="Duration"
            value={formatDuration(duration)}
          />
          <StatCard
            icon={<Waves className="size-4" />}
            label="Strokes"
            value={strokes.length.toString()}
          />
          <StatCard
            icon={<Gauge className="size-4" />}
            label="SPM"
            value={spm.toFixed(1)}
          />
          <StatCard
            icon={<Activity className="size-4" />}
            label="Samples"
            value={`${(sampleCount / 1000).toFixed(1)}k`}
          />
          {stats && (
            <>
              <StatCard
                icon={<Zap className="size-4" />}
                label="Avg ΔV Acc"
                value={`${stats.avgVelocityAcc.toFixed(2)} m/s`}
              />
              <StatCard
                icon={<Zap className="size-4" />}
                label="Avg ΔV Dec"
                value={`${stats.avgVelocityDec.toFixed(2)} m/s`}
              />
            </>
          )}
        </div>
        
        {/* Main Content */}
        <div className="grid gap-6 lg:grid-cols-2">
          {/* Live Waveform */}
          <Card className="border-border/50 bg-card lg:col-span-2">
            <CardHeader className="pb-2">
              <CardTitle className="flex items-center gap-2 text-sm font-medium">
                <Activity className="size-4 text-primary" />
                Live Signal
                {recentSamples.length > 0 && (
                  <Badge variant="secondary" className="ml-auto text-xs">
                    {recentSamples.length} samples
                  </Badge>
                )}
              </CardTitle>
            </CardHeader>
            <CardContent>
              <LiveWaveform 
                samples={recentSamples} 
                strokes={strokes}
              />
            </CardContent>
          </Card>
          
          {/* Recent Strokes */}
          <Card className="border-border/50 bg-card lg:col-span-2">
            <CardHeader className="pb-2">
              <CardTitle className="flex items-center gap-2 text-sm font-medium">
                <Waves className="size-4 text-primary" />
                Detected Strokes
                <Badge variant="secondary" className="ml-auto text-xs">
                  {strokes.length} total
                </Badge>
              </CardTitle>
            </CardHeader>
            <CardContent>
              <LiveStrokeList strokes={strokes} />
            </CardContent>
          </Card>
        </div>
        
        {/* Instructions */}
        {session.status === "idle" && strokes.length === 0 && (
          <Card className="mt-6 border-border/50 bg-secondary/30">
            <CardContent className="p-6">
              <h3 className="mb-2 font-semibold text-foreground">
                Waiting for IMU Data
              </h3>
              <p className="text-sm text-muted-foreground">
                This session is ready to receive data. Send IMU samples to:
              </p>
              <code className="mt-2 block rounded bg-secondary px-3 py-2 font-mono text-xs text-foreground">
                POST /api/sessions/{id}/ingest
              </code>
              <p className="mt-3 text-xs text-muted-foreground">
                Body: {"{"} samples: [{"{"} sample, Acc_X, Acc_Y, Gyr_X, Gyr_Y, Gyr_Z {"}"}] {"}"}
              </p>
            </CardContent>
          </Card>
        )}
      </main>
    </div>
  )
}

function StatCard({ 
  icon, 
  label, 
  value 
}: { 
  icon: React.ReactNode
  label: string
  value: string 
}) {
  return (
    <div className="rounded-lg border border-border/50 bg-card p-3">
      <div className="mb-1 flex items-center gap-2 text-muted-foreground">
        {icon}
        <span className="text-xs">{label}</span>
      </div>
      <p className="font-mono text-lg font-semibold text-foreground">{value}</p>
    </div>
  )
}

function formatDuration(seconds: number): string {
  const mins = Math.floor(seconds / 60)
  const secs = Math.floor(seconds % 60)
  return `${mins}:${secs.toString().padStart(2, "0")}`
}
