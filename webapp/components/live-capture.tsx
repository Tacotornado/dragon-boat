"use client"

import Link from "next/link"
import { useLiveCapture } from "@/hooks/use-live-capture"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Activity, ArrowLeft, Bluetooth, CircleStop, Radio, RotateCcw, Waves } from "lucide-react"

function formatDuration(seconds: number) {
  return `${Math.floor(seconds / 60)}:${String(Math.floor(seconds % 60)).padStart(2, "0")}`
}

export function LiveCapture({ sessionId }: { sessionId?: string }) {
  const capture = useLiveCapture(sessionId)
  const isRecording = capture.status === "recording"
  const isConnected = capture.status === "connected" || isRecording

  return (
    <main className="min-h-screen bg-background text-foreground">
      <header className="border-b border-border/50 bg-background/95 px-4 py-4 backdrop-blur">
        <div className="mx-auto flex max-w-2xl items-center justify-between">
          <div className="flex items-center gap-3">
            <Link href="/" className="rounded-md p-2 text-muted-foreground hover:bg-secondary hover:text-foreground" aria-label="Back to dashboard">
              <ArrowLeft className="size-5" />
            </Link>
            <div className="flex size-10 items-center justify-center rounded-xl bg-primary text-primary-foreground">
              <Radio className="size-5" />
            </div>
            <div>
              <h1 className="font-semibold tracking-tight">Live capture</h1>
              <p className="text-xs text-muted-foreground">Xsens DOT field monitor</p>
            </div>
          </div>
          <Badge variant={isRecording ? "default" : "secondary"} className={isRecording ? "bg-destructive text-destructive-foreground" : ""}>
            {isRecording ? "Recording" : capture.status === "error" ? "Error" : "Ready"}
          </Badge>
        </div>
      </header>

      <div className="mx-auto flex max-w-2xl flex-col gap-4 p-4 pb-10 sm:p-6">
        <Card className="border-border/50 bg-card">
          <CardContent className="flex items-center justify-between gap-4 p-4">
            <div className="flex min-w-0 items-center gap-3">
              <span className={`flex size-11 shrink-0 items-center justify-center rounded-full ${isConnected ? "bg-primary/15 text-primary" : "bg-secondary text-muted-foreground"}`}>
                <Bluetooth className="size-5" />
              </span>
              <div className="min-w-0">
                <p className="truncate font-medium">{capture.deviceName || "No sensor connected"}</p>
                <p className="text-sm text-muted-foreground">{isConnected ? "Bluetooth link active" : "Connect an Xsens DOT to begin"}</p>
              </div>
            </div>
  <div className="flex flex-wrap gap-2">
  {isConnected && (
  <Button variant="outline" onClick={() => void capture.identify()}>
  Identify sensor
  </Button>
  )}
  <Button variant={isConnected ? "outline" : "default"} onClick={isConnected ? capture.disconnect : capture.connect} disabled={capture.status === "connecting"}>
  {capture.status === "connecting" ? "Connecting..." : isConnected ? "Disconnect" : "Connect Xsens DOT"}
  </Button>
  </div>
          </CardContent>
        </Card>

        {capture.error && <div role="alert" className="rounded-lg border border-destructive/40 bg-destructive/10 px-4 py-3 text-sm text-destructive">{capture.error}</div>}

        <section className="grid grid-cols-2 gap-3 sm:grid-cols-4">
          {[
            ["Samples", capture.sampleCount.toLocaleString()],
            ["Rate", `${capture.sampleRate} Hz`],
            ["Buffered", capture.samples.length.toLocaleString()],
            ["Duration", formatDuration((capture.samples.at(-1)?.sample || 0) / 120)],
          ].map(([label, value]) => (
            <Card key={label} className="border-border/50 bg-card">
              <CardContent className="p-4">
                <p className="text-xs uppercase tracking-wide text-muted-foreground">{label}</p>
                <p className="mt-1 text-xl font-semibold tabular-nums">{value}</p>
              </CardContent>
            </Card>
          ))}
        </section>

        <Card className="border-border/50 bg-card">
          <CardHeader className="pb-3"><CardTitle className="flex items-center gap-2 text-base"><Activity className="size-4 text-primary" />Signal status</CardTitle></CardHeader>
          <CardContent>
            <div className="flex h-28 items-end gap-px overflow-hidden rounded-lg bg-secondary/60 p-3" aria-label="Recent sensor samples">
              {capture.samples.slice(-100).map((sample, index) => {
                const height = Math.max(5, Math.min(100, 50 + sample.Acc_X * 8))
                return <span key={`${sample.timestamp}-${index}`} className="flex-1 rounded-t-sm bg-primary/70" style={{ height: `${height}%` }} />
              })}
              {capture.samples.length === 0 && <p className="m-auto text-sm text-muted-foreground">Waiting for sensor packets</p>}
            </div>
          </CardContent>
        </Card>

        <div className="grid grid-cols-1 gap-3 sm:grid-cols-2">
          <Button size="lg" className="h-14 gap-2" disabled={!isConnected || isRecording} onClick={() => void capture.start()}>
            <Waves className="size-5" /> Start recording
          </Button>
          <Button size="lg" variant="destructive" className="h-14 gap-2" disabled={!isRecording} onClick={() => void capture.stop()}>
            <CircleStop className="size-5" /> Stop and save
          </Button>
        </div>

        {capture.engineResult && (
          <div className="rounded-lg border border-primary/30 bg-primary/10 px-4 py-3 text-sm">
            {capture.engineResult.ran
              ? `C++ engine completed: ${capture.engineResult.strokeCount} strokes detected.`
              : capture.engineResult.error || "C++ engine did not run."}
          </div>
        )}

        <p className="text-center text-xs leading-relaxed text-muted-foreground">
          Web Bluetooth requires Chrome or Edge, HTTPS or localhost, and a nearby powered-on Xsens DOT.
        </p>
      </div>
    </main>
  )
}
