"use client"

import useSWR from "swr"
import Link from "next/link"
import { useState } from "react"
import { ArrowLeft, CalendarDays, Download, FileText, Trash2, Upload, UserRound } from "lucide-react"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { RunArtifacts } from "@/components/run-artifacts"

const fetcher = (url: string) => fetch(url).then((response) => response.json())

export default function RecordsPage() {
  const { data, mutate } = useSWR("/api/runs", fetcher)
  const [deleting, setDeleting] = useState<string | null>(null)
  const runs = Array.isArray(data) ? data : []

  async function removeRun(id: string) {
    if (!window.confirm("Delete this stored run?")) return
    setDeleting(id)
    await fetch(`/api/runs/${id}`, { method: "DELETE" })
    await mutate()
    setDeleting(null)
  }

  return (
    <main className="min-h-screen bg-background px-4 py-6 text-foreground sm:px-6 lg:px-10">
      <div className="mx-auto max-w-6xl space-y-6">
        <header className="flex flex-col gap-4 sm:flex-row sm:items-center sm:justify-between">
          <div className="flex items-center gap-3">
            <Link href="/">
              <Button variant="ghost" size="icon" aria-label="Back to dashboard"><ArrowLeft className="size-4" /></Button>
            </Link>
            <div>
              <p className="text-xs font-medium uppercase tracking-[0.18em] text-primary">Stroke archive</p>
              <h1 className="text-2xl font-semibold tracking-tight sm:text-3xl">Past records</h1>
            </div>
          </div>
          <Link href="/">
            <Button className="gap-2"><Upload className="size-4" />Upload another run</Button>
          </Link>
        </header>

        {runs.length === 0 ? (
          <Card className="border-dashed">
            <CardContent className="flex flex-col items-center gap-3 py-16 text-center">
              <FileText className="size-10 text-muted-foreground" />
              <h2 className="text-lg font-semibold">No saved runs yet</h2>
              <p className="max-w-md text-sm text-muted-foreground">Upload a raw Xsens IMU CSV from the dashboard and complete the racer metadata form to create a historical record.</p>
            </CardContent>
          </Card>
        ) : (
          <div className="grid gap-4 md:grid-cols-2">
            {runs.map((run: { id: string; name: string; racerName: string; uploaderName?: string | null; note?: string | null; runDate: string; raceDate?: string | null; uploadedAt: string; sampleCount: number; strokeCount: number }) => (
              <Card key={run.id} className="border-border/70 transition-colors hover:border-primary/50">
                <CardHeader className="flex flex-row items-start justify-between gap-4 space-y-0">
                  <div className="min-w-0">
                    <CardTitle className="truncate text-lg">{run.name}</CardTitle>
                    <p className="mt-1 flex items-center gap-1.5 text-sm text-muted-foreground"><UserRound className="size-3.5" />{run.racerName}</p>
                  </div>
                  <Badge variant="secondary">{run.strokeCount || 0} strokes</Badge>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div className="grid grid-cols-2 gap-3 text-sm">
                    <div><p className="text-xs text-muted-foreground">Run date</p><p className="mt-1 flex items-center gap-1.5"><CalendarDays className="size-3.5 text-primary" />{run.runDate}</p></div>
                    <div><p className="text-xs text-muted-foreground">Samples</p><p className="mt-1 font-mono">{run.sampleCount.toLocaleString()}</p></div>
                  </div>
                  {run.note && <p className="rounded-md bg-secondary/60 px-3 py-2 text-sm text-muted-foreground">{run.note}</p>}
                  <RunArtifacts runId={run.id} runName={run.name} />
                  <div className="flex items-center justify-between gap-3 border-t border-border/60 pt-3">
                    <p className="text-xs text-muted-foreground">Uploaded {new Date(run.uploadedAt).toLocaleString()}</p>
                    <div className="flex gap-2"><a href={`/api/runs/${run.id}/export`} download><Button size="sm" variant="outline" aria-label={`Export ${run.name}`}><Download className="size-4" /></Button></a><Link href={`/?run=${run.id}`}><Button size="sm">Open</Button></Link><Button size="sm" variant="ghost" onClick={() => removeRun(run.id)} disabled={deleting === run.id} aria-label={`Delete ${run.name}`}><Trash2 className="size-4 text-destructive" /></Button></div>
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>
        )}
      </div>
    </main>
  )
}
