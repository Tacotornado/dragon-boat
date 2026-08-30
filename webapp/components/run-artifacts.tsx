"use client"

import useSWR from "swr"
import { useState } from "react"
import Image from "next/image"
import { ChevronDown, FileSpreadsheet, FileText, ImageIcon } from "lucide-react"
import { Button } from "@/components/ui/button"

interface Artifact {
  id: number
  kind: "cycles_csv" | "phase_metrics_xlsx" | "overview_plot" | "stroke_plot"
  filename: string
  byteSize: number
  strokeId: number | null
}

const fetcher = (url: string) => fetch(url).then((response) => response.json())

function formatSize(bytes: number) {
  if (bytes < 1024) return `${bytes} B`
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(0)} KB`
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`
}

export function RunArtifacts({ runId, runName }: { runId: string; runName: string }) {
  const [open, setOpen] = useState(false)
  const { data, isLoading } = useSWR<Artifact[]>(open ? `/api/runs/${runId}/artifacts` : null, fetcher)
  const artifacts = Array.isArray(data) ? data : []

  const overview = artifacts.find((artifact) => artifact.kind === "overview_plot")
  const documents = artifacts.filter(
    (artifact) => artifact.kind === "cycles_csv" || artifact.kind === "phase_metrics_xlsx",
  )
  const strokePlots = artifacts.filter((artifact) => artifact.kind === "stroke_plot")

  const href = (artifact: Artifact, inline = false) =>
    `/api/runs/${runId}/artifacts/${artifact.id}${inline ? "?inline=1" : ""}`

  return (
    <div className="border-t border-border/60 pt-3">
      <Button
        variant="ghost"
        size="sm"
        className="h-auto w-full justify-between px-0 text-xs font-medium uppercase tracking-[0.14em] text-muted-foreground hover:bg-transparent hover:text-foreground"
        onClick={() => setOpen((value) => !value)}
        aria-expanded={open}
      >
        Engine outputs
        <ChevronDown className={`size-4 transition-transform ${open ? "rotate-180" : ""}`} />
      </Button>

      {open && (
        <div className="mt-3 flex flex-col gap-3">
          {isLoading && <p className="text-xs text-muted-foreground">Loading artifacts…</p>}

          {!isLoading && artifacts.length === 0 && (
            <p className="text-xs text-muted-foreground">
              No engine outputs stored. Enable the C++ engine and re-upload this run.
            </p>
          )}

          {overview && (
            <a
              href={href(overview, true)}
              target="_blank"
              rel="noreferrer"
              className="group relative block overflow-hidden rounded-lg border border-border/60 bg-secondary/40"
            >
              <Image
                src={href(overview, true)}
                alt={`Detected stroke cycles overview for ${runName}`}
                width={1600}
                height={600}
                unoptimized
                className="h-auto w-full transition-opacity group-hover:opacity-90"
              />
            </a>
          )}

          {documents.length > 0 && (
            <div className="flex flex-wrap gap-2">
              {documents.map((artifact) => (
                <a key={artifact.id} href={href(artifact)} download>
                  <Button variant="outline" size="sm" className="gap-2">
                    {artifact.kind === "phase_metrics_xlsx" ? (
                      <FileSpreadsheet className="size-3.5 text-primary" />
                    ) : (
                      <FileText className="size-3.5 text-primary" />
                    )}
                    <span className="font-mono text-xs">{artifact.filename}</span>
                    <span className="text-xs text-muted-foreground">{formatSize(artifact.byteSize)}</span>
                  </Button>
                </a>
              ))}
            </div>
          )}

          {strokePlots.length > 0 && (
            <div className="flex flex-col gap-2">
              <p className="flex items-center gap-1.5 text-xs text-muted-foreground">
                <ImageIcon className="size-3.5" />
                {strokePlots.length} per-stroke plots
              </p>
              <div className="flex flex-wrap gap-1.5">
                {strokePlots.map((artifact) => (
                  <a key={artifact.id} href={href(artifact, true)} target="_blank" rel="noreferrer">
                    <Button variant="secondary" size="sm" className="h-7 px-2 font-mono text-xs">
                      {artifact.strokeId ?? "?"}
                    </Button>
                  </a>
                ))}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  )
}
