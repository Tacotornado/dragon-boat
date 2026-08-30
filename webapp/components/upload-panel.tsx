"use client"

import { useCallback, useRef, useState } from "react"
import { Card, CardContent } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Upload, FileText, Waves, CheckCircle2, AlertCircle } from "lucide-react"

interface UploadPanelProps {
  onStrokesUpload: (file: File) => void
  onIMUUpload: (file: File) => void
  dataSource: "demo" | "uploaded"
  uploadMessage: string | null
  hasUploadedWaveform: boolean
}

export function UploadPanel({
  onStrokesUpload,
  onIMUUpload,
  dataSource,
  uploadMessage,
  hasUploadedWaveform,
}: UploadPanelProps) {
  const strokesInputRef = useRef<HTMLInputElement>(null)
  const imuInputRef = useRef<HTMLInputElement>(null)
  const [isDraggingStrokes, setIsDraggingStrokes] = useState(false)
  const [isDraggingIMU, setIsDraggingIMU] = useState(false)
  const [strokesFileName, setStrokesFileName] = useState<string | null>(null)
  const [imuFileName, setIMUFileName] = useState<string | null>(null)

  const handleStrokesFile = useCallback(
    (file: File) => {
      setStrokesFileName(file.name)
      onStrokesUpload(file)
    },
    [onStrokesUpload]
  )

  const handleIMUFile = useCallback(
    async (file: File) => {
      setIMUFileName(file.name)
      onIMUUpload(file)

      const racerName = window.prompt("Enter the racer name for this run:")?.trim()
      if (!racerName) return
      const uploaderName = window.prompt("Uploader name (optional):")?.trim() ?? ""
      const note = window.prompt("Run note (optional):")?.trim() ?? ""
      const raceDate = window.prompt("Race date (YYYY-MM-DD, optional):")?.trim() ?? ""
      const form = new FormData()
      form.append("file", file)
      form.append("racerName", racerName)
      if (uploaderName) form.append("uploaderName", uploaderName)
      if (note) form.append("note", note)
      if (raceDate) form.append("raceDate", raceDate)
      try {
        const response = await fetch("/api/runs", { method: "POST", body: form })
        const responseText = await response.text()
        let result: {
          error?: string
          name?: string
          sampleCount?: number
          strokeCount?: number
          engine?: { enabled: boolean; ran: boolean; error: string | null }
        } = {}
        try {
          result = responseText ? JSON.parse(responseText) : {}
        } catch {
          result = { error: responseText || `Upload failed (${response.status})` }
        }
        if (!response.ok) throw new Error(result.error || `Upload failed (${response.status})`)

        const summary = [`Saved as ${result.name}`, `${result.sampleCount ?? 0} samples`]
        if (result.strokeCount) summary.push(`${result.strokeCount} strokes detected by the C++ engine`)
        else if (result.engine && !result.engine.enabled) summary.push("C++ engine disabled (USE_CPP_ENGINE=false)")
        else if (result.engine?.error) summary.push(`Engine error: ${result.engine.error}`)
        window.alert(summary.join("\n"))
      } catch (error) {
        window.alert(error instanceof Error ? error.message : "Could not save run")
      }
    },
    [onIMUUpload]
  )

  const handleDrop = useCallback(
    (e: React.DragEvent, handler: (file: File) => void) => {
      e.preventDefault()
      setIsDraggingStrokes(false)
      setIsDraggingIMU(false)
      const file = e.dataTransfer.files[0]
      if (file && file.name.endsWith(".csv")) {
        handler(file)
      }
    },
    []
  )

  return (
    <Card className="border-border/50">
      <CardContent className="pt-0">
        <div className="flex flex-col gap-4">
          {/* Status message */}
          {uploadMessage && (
            <div
              className={`flex items-center gap-2 rounded-lg px-3 py-2 text-sm ${
                uploadMessage.toLowerCase().includes("error") ||
                uploadMessage.toLowerCase().includes("missing") ||
                uploadMessage.toLowerCase().includes("no ")
                  ? "bg-destructive/10 text-destructive"
                  : "bg-primary/10 text-primary"
              }`}
            >
              {uploadMessage.toLowerCase().includes("error") ||
              uploadMessage.toLowerCase().includes("missing") ||
              uploadMessage.toLowerCase().includes("no ") ? (
                <AlertCircle className="size-4 shrink-0" />
              ) : (
                <CheckCircle2 className="size-4 shrink-0" />
              )}
              {uploadMessage}
            </div>
          )}

          <div className="grid grid-cols-1 gap-4 md:grid-cols-2">
            {/* Stroke results upload */}
            <div
              className={`relative flex flex-col items-center gap-2 rounded-xl border-2 border-dashed p-4 transition-colors sm:gap-3 sm:p-6 ${
                isDraggingStrokes
                  ? "border-primary bg-primary/5"
                  : strokesFileName
                    ? "border-primary/40 bg-primary/5"
                    : "border-border/50 hover:border-border"
              }`}
              onDragOver={(e) => {
                e.preventDefault()
                setIsDraggingStrokes(true)
              }}
              onDragLeave={() => setIsDraggingStrokes(false)}
              onDrop={(e) => handleDrop(e, handleStrokesFile)}
            >
              <div className="flex size-8 items-center justify-center rounded-lg bg-secondary text-primary sm:size-10">
                {strokesFileName ? <CheckCircle2 className="size-5" /> : <FileText className="size-5" />}
              </div>
              <div className="text-center">
                <p className="text-sm font-medium text-foreground">
                  {strokesFileName ? strokesFileName : "Stroke Results CSV"}
                </p>
                <p className="mt-1 hidden text-xs text-muted-foreground sm:block">
                  Upload <span className="font-mono text-foreground/70">stroke_full_cycles.csv</span>
                </p>
              </div>
              <Button
                variant="outline"
                size="sm"
                className="gap-2"
                onClick={() => strokesInputRef.current?.click()}
              >
                <Upload className="size-3.5" />
                {strokesFileName ? "Replace" : "Choose File"}
              </Button>
              <input
                ref={strokesInputRef}
                type="file"
                accept=".csv"
                className="sr-only"
                onChange={(e) => {
                  const file = e.target.files?.[0]
                  if (file) handleStrokesFile(file)
                }}
              />
            </div>

            {/* Raw IMU data upload (optional) */}
            <div
              className={`relative flex flex-col items-center gap-2 rounded-xl border-2 border-dashed p-4 transition-colors sm:gap-3 sm:p-6 ${
                isDraggingIMU
                  ? "border-primary bg-primary/5"
                  : imuFileName
                    ? "border-primary/40 bg-primary/5"
                    : "border-border/50 hover:border-border"
              }`}
              onDragOver={(e) => {
                e.preventDefault()
                setIsDraggingIMU(true)
              }}
              onDragLeave={() => setIsDraggingIMU(false)}
              onDrop={(e) => handleDrop(e, handleIMUFile)}
            >
              <div className="flex size-8 items-center justify-center rounded-lg bg-secondary text-muted-foreground sm:size-10">
                {imuFileName ? <CheckCircle2 className="size-5 text-primary" /> : <Waves className="size-5" />}
              </div>
              <div className="text-center">
                <p className="text-sm font-medium text-foreground">
                  {imuFileName ? imuFileName : "Raw IMU Data (Optional)"}
                </p>
                <p className="mt-1 hidden text-xs text-muted-foreground sm:block">
                  Upload <span className="font-mono text-foreground/70">Canoe_xsens_dot.csv</span>
                </p>
              </div>
              <Button
                variant="outline"
                size="sm"
                className="gap-2"
                onClick={() => imuInputRef.current?.click()}
              >
                <Upload className="size-3.5" />
                {imuFileName ? "Replace" : "Choose File"}
              </Button>
              <input
                ref={imuInputRef}
                type="file"
                accept=".csv"
                className="sr-only"
                onChange={(e) => {
                  const file = e.target.files?.[0]
                  if (file) handleIMUFile(file)
                }}
              />
            </div>
          </div>

          {/* How-to hint - hidden on mobile */}
          {dataSource === "demo" && (
            <div className="hidden rounded-lg bg-secondary/50 px-4 py-3 sm:block">
              <p className="text-xs leading-relaxed text-muted-foreground">
                <span className="font-medium text-foreground">How it works:</span>{" "}
                Run your C++ <span className="font-mono">stroke_analysis</span> binary locally, then upload the resulting{" "}
                <span className="font-mono">stroke_full_cycles.csv</span> here. Currently showing demo data.
              </p>
            </div>
          )}
        </div>
      </CardContent>
    </Card>
  )
}
