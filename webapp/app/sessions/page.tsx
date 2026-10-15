"use client"

import { useState, useEffect } from "react"
import Link from "next/link"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Input } from "@/components/ui/input"
import { 
  Activity, 
  Plus, 
  Play, 
  Clock, 
  Waves,
  ArrowLeft,
  Trash2,
  Radio
} from "lucide-react"
import type { Session } from "@/lib/types"

export default function SessionsPage() {
  const [sessions, setSessions] = useState<Session[]>([])
  const [isLoading, setIsLoading] = useState(true)
  const [newSessionName, setNewSessionName] = useState("")
  const [athleteName, setAthleteName] = useState("")
  const [isCreating, setIsCreating] = useState(false)
  
  useEffect(() => {
    fetchSessions()
  }, [])
  
  async function fetchSessions() {
    try {
      const res = await fetch("/api/sessions")
      const data = await res.json()
      setSessions(data)
    } catch (error) {
      console.error("Failed to fetch sessions:", error)
    } finally {
      setIsLoading(false)
    }
  }
  
  async function createSession() {
    if (!newSessionName.trim()) return
    
    setIsCreating(true)
    try {
      const res = await fetch("/api/sessions", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ 
          name: newSessionName,
          athleteName: athleteName || undefined
        }),
      })
      
      if (res.ok) {
        setNewSessionName("")
        setAthleteName("")
        fetchSessions()
      }
    } catch (error) {
      console.error("Failed to create session:", error)
    } finally {
      setIsCreating(false)
    }
  }
  
  async function deleteSession(id: string) {
    try {
      await fetch(`/api/sessions/${id}`, { method: "DELETE" })
      fetchSessions()
    } catch (error) {
      console.error("Failed to delete session:", error)
    }
  }
  
  function formatDuration(seconds: number): string {
    const mins = Math.floor(seconds / 60)
    const secs = Math.floor(seconds % 60)
    return `${mins}:${secs.toString().padStart(2, "0")}`
  }
  
  function formatDate(iso: string): string {
    return new Date(iso).toLocaleDateString("en-US", {
      month: "short",
      day: "numeric",
      hour: "2-digit",
      minute: "2-digit",
    })
  }
  
  const statusColors: Record<string, string> = {
    idle: "bg-muted text-muted-foreground",
    recording: "bg-red-500/20 text-red-400 animate-pulse",
    processing: "bg-yellow-500/20 text-yellow-400",
    completed: "bg-green-500/20 text-green-400",
    error: "bg-destructive/20 text-destructive",
  }
  
  return (
    <div className="min-h-screen bg-background">
      <header className="sticky top-0 z-50 border-b border-border/50 bg-background/95 backdrop-blur">
        <div className="mx-auto flex max-w-[1440px] items-center justify-between px-4 py-3 md:px-6">
          <div className="flex items-center gap-3">
            <Link href="/" className="text-muted-foreground hover:text-foreground">
              <ArrowLeft className="size-5" />
            </Link>
            <div className="flex size-9 items-center justify-center rounded-lg bg-primary text-primary-foreground">
              <Radio className="size-5" />
            </div>
            <div>
              <h1 className="text-lg font-semibold leading-tight text-foreground">
                Live Sessions
              </h1>
              <p className="text-xs text-muted-foreground">
                Real-time stroke analysis
              </p>
            </div>
          </div>
          <Link href="/live">
            <Button variant="outline" className="gap-2">
              <Radio className="size-4" />
              Connect sensor
            </Button>
          </Link>
        </div>
      </header>
      
      <main className="mx-auto max-w-[1440px] p-4 md:p-6">
        {/* Create New Session */}
        <Card className="mb-6 border-border/50 bg-card">
          <CardHeader className="pb-3">
            <CardTitle className="text-base">New Session</CardTitle>
            <CardDescription>
              Create a recording session for live IMU data
            </CardDescription>
          </CardHeader>
          <CardContent>
            <div className="flex flex-col gap-3 sm:flex-row">
              <Input
                placeholder="Session name (e.g., Morning Practice)"
                value={newSessionName}
                onChange={(e) => setNewSessionName(e.target.value)}
                className="flex-1"
              />
              <Input
                placeholder="Athlete name (optional)"
                value={athleteName}
                onChange={(e) => setAthleteName(e.target.value)}
                className="flex-1"
              />
              <Button 
                onClick={createSession} 
                disabled={!newSessionName.trim() || isCreating}
                className="gap-2"
              >
                <Plus className="size-4" />
                Create
              </Button>
            </div>
          </CardContent>
        </Card>
        
        {/* Sessions List */}
        <div className="grid gap-4 sm:grid-cols-2 lg:grid-cols-3">
          {isLoading ? (
            Array.from({ length: 3 }).map((_, i) => (
              <Card key={i} className="animate-pulse border-border/50 bg-card">
                <CardContent className="p-6">
                  <div className="h-4 w-3/4 rounded bg-secondary" />
                  <div className="mt-2 h-3 w-1/2 rounded bg-secondary" />
                </CardContent>
              </Card>
            ))
          ) : sessions.length === 0 ? (
            <Card className="col-span-full border-border/50 bg-card">
              <CardContent className="flex flex-col items-center justify-center py-12 text-center">
                <Activity className="mb-3 size-12 text-muted-foreground/50" />
                <p className="text-sm text-muted-foreground">
                  No sessions yet. Create one to start recording.
                </p>
              </CardContent>
            </Card>
          ) : (
            sessions.map((session) => (
              <Card 
                key={session.id} 
                className="group border-border/50 bg-card transition-colors hover:border-primary/30"
              >
                <CardContent className="p-4">
                  <div className="mb-3 flex items-start justify-between">
                    <div className="flex-1 pr-2">
                      <h3 className="font-semibold text-foreground line-clamp-1">
                        {session.name}
                      </h3>
                      {session.athleteName && (
                        <p className="text-xs text-muted-foreground">
                          {session.athleteName}
                        </p>
                      )}
                    </div>
                    <Badge className={statusColors[session.status]}>
                      {session.status}
                    </Badge>
                  </div>
                  
                  <div className="mb-4 grid grid-cols-3 gap-2 text-center">
                    <div className="rounded-lg bg-secondary/50 p-2">
                      <Clock className="mx-auto mb-1 size-4 text-muted-foreground" />
                      <p className="font-mono text-sm font-medium">
                        {formatDuration(session.duration)}
                      </p>
                      <p className="text-[10px] text-muted-foreground">Duration</p>
                    </div>
                    <div className="rounded-lg bg-secondary/50 p-2">
                      <Waves className="mx-auto mb-1 size-4 text-muted-foreground" />
                      <p className="font-mono text-sm font-medium">
                        {session.strokeCount}
                      </p>
                      <p className="text-[10px] text-muted-foreground">Strokes</p>
                    </div>
                    <div className="rounded-lg bg-secondary/50 p-2">
                      <Activity className="mx-auto mb-1 size-4 text-muted-foreground" />
                      <p className="font-mono text-sm font-medium">
                        {(session.sampleCount / 1000).toFixed(1)}k
                      </p>
                      <p className="text-[10px] text-muted-foreground">Samples</p>
                    </div>
                  </div>
                  
                  <p className="mb-3 text-xs text-muted-foreground">
                    Created {formatDate(session.createdAt)}
                  </p>
                  
                  <div className="flex gap-2">
                    <Link href={`/live?sessionId=${encodeURIComponent(session.id)}`} className="flex-1">
                      <Button variant="default" size="sm" className="w-full gap-2">
                        <Radio className="size-3.5" />
                        Record
                      </Button>
                    </Link>
                    <Link href={`/sessions/${session.id}`}>
                      <Button variant="outline" size="sm" aria-label={`Open ${session.name}`}>
                        <Play className="size-3.5" />
                      </Button>
                    </Link>
                    <Button 
                      variant="ghost" 
                      size="sm"
                      onClick={() => deleteSession(session.id)}
                      className="text-muted-foreground hover:text-destructive"
                    >
                      <Trash2 className="size-4" />
                    </Button>
                  </div>
                </CardContent>
              </Card>
            ))
          )}
        </div>
      </main>
    </div>
  )
}
