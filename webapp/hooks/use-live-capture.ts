"use client"

import { useCallback, useEffect, useRef, useState } from "react"
import type { IMUDataPoint } from "@/lib/types"
import {
  connectXsensDot,
  disconnectXsensDot,
  identifyXsensDot,
  startXsensDotMeasurement,
  type XsensDotDevice,
} from "@/lib/xsens-dot-ble"

const BATCH_SIZE = 24
const FLUSH_INTERVAL_MS = 250

export type CaptureStatus = "idle" | "connecting" | "connected" | "recording" | "stopped" | "error"

export function useLiveCapture(existingSessionId?: string) {
  const [status, setStatus] = useState<CaptureStatus>("idle")
  const [error, setError] = useState<string | null>(null)
  const [samples, setSamples] = useState<IMUDataPoint[]>([])
  const [sampleRate, setSampleRate] = useState(0)
  const [deviceName, setDeviceName] = useState<string | null>(null)
  const [sessionId, setSessionId] = useState<string | null>(null)
  const [engineResult, setEngineResult] = useState<{ ran: boolean; strokeCount: number; error?: string } | null>(null)
  const connectionRef = useRef<XsensDotDevice | null>(null)
  const batchRef = useRef<IMUDataPoint[]>([])
  const totalRef = useRef(0)
  const firstTimestampRef = useRef<number | null>(null)
  const lastFlushRef = useRef(0)
  const rateWindowRef = useRef<number[]>([])

  const flush = useCallback(async (force = false) => {
    if (!sessionId || (!force && batchRef.current.length === 0)) return
    const batch = batchRef.current.splice(0)
    if (batch.length === 0) return
    try {
      await fetch(`/api/sessions/${sessionId}/ingest`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ samples: batch }),
      })
    } catch (cause) {
      setError(cause instanceof Error ? cause.message : "Could not send sensor data")
    }
  }, [sessionId])

  const handleSamples = useCallback((incoming: IMUDataPoint[]) => {
    const normalized = incoming.map((sample) => ({ ...sample }))
    totalRef.current += normalized.length
    batchRef.current.push(...normalized)
    setSamples((previous) => [...previous, ...normalized].slice(-1200))
    const now = performance.now()
    rateWindowRef.current.push(now)
    rateWindowRef.current = rateWindowRef.current.filter((time) => now - time < 1000)
    setSampleRate(rateWindowRef.current.length)
    if (batchRef.current.length >= BATCH_SIZE || now - lastFlushRef.current > FLUSH_INTERVAL_MS) {
      lastFlushRef.current = now
      void flush()
    }
  }, [flush])

  const connect = useCallback(async () => {
    setStatus("connecting")
    setError(null)
    try {
      const connection = await connectXsensDot((packet) => handleSamples(packet.samples))
      connectionRef.current = connection
      setDeviceName(connection.device.name || "Xsens DOT")
      setStatus("connected")
    } catch (cause) {
      setStatus("error")
      setError(cause instanceof Error ? cause.message : "Could not connect to sensor")
    }
  }, [handleSamples])

  const start = useCallback(async () => {
    if (status !== "connected") return
    try {
    let recordingSessionId = existingSessionId
    if (!recordingSessionId) {
      const response = await fetch("/api/sessions", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ name: `Live capture ${new Date().toLocaleString()}` }),
      })
      if (!response.ok) throw new Error("Could not create live session")
      const session = await response.json()
      recordingSessionId = session.id
    }
    if (!recordingSessionId) throw new Error("Could not determine recording session")
    if (!connectionRef.current) throw new Error("Sensor disconnected before recording started")
    await startXsensDotMeasurement(connectionRef.current)
    setSessionId(recordingSessionId)
    setSamples([])
    totalRef.current = 0
    batchRef.current = []
    firstTimestampRef.current = null
    setStatus("recording")
    } catch (cause) {
      setError(cause instanceof Error ? cause.message : "Could not start sensor recording")
      setStatus("error")
    }
  }, [status])

  const stop = useCallback(async () => {
    await flush(true)
    if (sessionId) {
      const response = await fetch(`/api/sessions/${sessionId}`, {
        method: "PATCH",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ status: "completed" }),
      })
      const result = await response.json()
      if (!response.ok) throw new Error(result.error || "Could not finalize session")
      setEngineResult({
        ran: Boolean(result.engine?.ran),
        strokeCount: result.engine?.cycles?.length ?? 0,
        error: result.engine?.error,
      })
      if (result.engine?.error) setError(`Final engine run failed: ${result.engine.error}`)
    }
    setStatus("stopped")
  }, [flush, sessionId])

  const identify = useCallback(async () => {
    if (!connectionRef.current) return
    try {
      setError(null)
      await identifyXsensDot(connectionRef.current)
    } catch (cause) {
      setError(cause instanceof Error ? cause.message : "Could not identify sensor")
    }
  }, [])

  const disconnect = useCallback(async () => {
    await disconnectXsensDot(connectionRef.current)
    connectionRef.current = null
    setStatus("idle")
  }, [])

  useEffect(() => () => { void disconnectXsensDot(connectionRef.current) }, [])

  return {
    status,
    error,
    samples,
    sampleRate,
    sampleCount: totalRef.current,
    deviceName,
    sessionId,
    engineResult,
    connect,
    identify,
    start,
    stop,
    disconnect,
  }
}
