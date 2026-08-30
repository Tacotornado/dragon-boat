"use client"

import { useState, useEffect, useCallback, useRef } from "react"
import type {
  Session,
  SessionWithData,
  StrokeCycle,
  IMUDataPoint,
  SummaryStats,
} from "@/lib/types"

interface LiveSessionState {
  session: Session | null
  isConnected: boolean
  isLoading: boolean
  error: string | null
  // Live data
  recentSamples: IMUDataPoint[]
  strokes: StrokeCycle[]
  stats: SummaryStats | null
  sampleCount: number
}

interface LiveSessionActions {
  connect: () => void
  disconnect: () => void
  refresh: () => Promise<void>
}

const MAX_RECENT_SAMPLES = 1200 // 10 seconds at 120Hz

export function useLiveSession(sessionId: string | null): LiveSessionState & LiveSessionActions {
  const [state, setState] = useState<LiveSessionState>({
    session: null,
    isConnected: false,
    isLoading: true,
    error: null,
    recentSamples: [],
    strokes: [],
    stats: null,
    sampleCount: 0,
  })
  
  const eventSourceRef = useRef<EventSource | null>(null)
  const reconnectTimeoutRef = useRef<NodeJS.Timeout | null>(null)
  
  // Fetch initial session data
  const refresh = useCallback(async () => {
    if (!sessionId) return
    
    setState(prev => ({ ...prev, isLoading: true, error: null }))
    
    try {
      const res = await fetch(`/api/sessions/${sessionId}`)
      if (!res.ok) throw new Error("Failed to fetch session")
      
      const session: SessionWithData = await res.json()
      
      setState(prev => ({
        ...prev,
        session,
        strokes: session.strokes || [],
        recentSamples: (session.rawData || []).slice(-MAX_RECENT_SAMPLES),
        sampleCount: session.sampleCount || 0,
        isLoading: false,
      }))
    } catch (error) {
      setState(prev => ({
        ...prev,
        error: error instanceof Error ? error.message : "Unknown error",
        isLoading: false,
      }))
    }
  }, [sessionId])
  
  // Connect to SSE stream
  const connect = useCallback(() => {
    if (!sessionId || eventSourceRef.current) return
    
    const eventSource = new EventSource(`/api/sessions/${sessionId}/stream`)
    eventSourceRef.current = eventSource
    
    eventSource.onopen = () => {
      setState(prev => ({ ...prev, isConnected: true, error: null }))
    }
    
    eventSource.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data)
        
        switch (data.type) {
          case "connected":
            setState(prev => ({
              ...prev,
              sampleCount: data.sampleCount || prev.sampleCount,
            }))
            break
            
          case "imu_data":
            setState(prev => {
              const newSamples = [...prev.recentSamples, ...data.samples]
              return {
                ...prev,
                recentSamples: newSamples.slice(-MAX_RECENT_SAMPLES),
                sampleCount: prev.sampleCount + data.sampleCount,
              }
            })
            break
            
          case "stroke_detected":
            setState(prev => ({
              ...prev,
              strokes: [...prev.strokes, ...data.newStrokes],
              stats: data.stats || prev.stats,
            }))
            break
        }
      } catch (error) {
        console.error("SSE parse error:", error)
      }
    }
    
    eventSource.onerror = () => {
      setState(prev => ({ ...prev, isConnected: false }))
      eventSource.close()
      eventSourceRef.current = null
      
      // Attempt to reconnect after 3 seconds
      reconnectTimeoutRef.current = setTimeout(() => {
        connect()
      }, 3000)
    }
  }, [sessionId])
  
  // Disconnect from SSE stream
  const disconnect = useCallback(() => {
    if (reconnectTimeoutRef.current) {
      clearTimeout(reconnectTimeoutRef.current)
      reconnectTimeoutRef.current = null
    }
    
    if (eventSourceRef.current) {
      eventSourceRef.current.close()
      eventSourceRef.current = null
    }
    
    setState(prev => ({ ...prev, isConnected: false }))
  }, [])
  
  // Auto-connect and fetch on mount
  useEffect(() => {
    if (sessionId) {
      refresh()
      connect()
    }
    
    return () => {
      disconnect()
    }
  }, [sessionId, refresh, connect, disconnect])
  
  return {
    ...state,
    connect,
    disconnect,
    refresh,
  }
}
