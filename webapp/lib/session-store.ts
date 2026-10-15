"use server"

import { 
  Session, 
  SessionWithData, 
  StrokeCycle, 
  IMUDataPoint,
  AccPhaseMetrics,
  DecPhaseMetrics,
  StrokePhaseMetrics,
  SummaryStats
} from "./types"

// In-memory store for development
// In production, replace with Redis or PostgreSQL
const sessions = new Map<string, SessionWithData>()

function generateId(): string {
  return `sess_${Date.now()}_${Math.random().toString(36).slice(2, 9)}`
}

export async function createSession(name: string, athleteName?: string): Promise<Session> {
  const id = generateId()
  const now = new Date().toISOString()
  
  const session: SessionWithData = {
    id,
    name,
    athleteName,
    createdAt: now,
    updatedAt: now,
    status: "idle",
    duration: 0,
    sampleCount: 0,
    strokeCount: 0,
    rawData: [],
    strokes: [],
    accMetrics: [],
    decMetrics: [],
    strokePhaseMetrics: [],
  }
  
  sessions.set(id, session)
  return sessionToPublic(session)
}

export async function getSession(id: string): Promise<SessionWithData | null> {
  return sessions.get(id) || null
}

export async function listSessions(): Promise<Session[]> {
  return Array.from(sessions.values())
    .map(sessionToPublic)
    .sort((a, b) => new Date(b.createdAt).getTime() - new Date(a.createdAt).getTime())
}

export async function updateSessionStatus(
  id: string, 
  status: Session["status"]
): Promise<Session | null> {
  const session = sessions.get(id)
  if (!session) return null
  
  session.status = status
  session.updatedAt = new Date().toISOString()
  return sessionToPublic(session)
}

export async function appendIMUData(
  id: string, 
  samples: IMUDataPoint[]
): Promise<boolean> {
  const session = sessions.get(id)
  if (!session) return false
  
  session.rawData.push(...samples)
  session.sampleCount = session.rawData.length
  session.updatedAt = new Date().toISOString()
  
  // Update duration based on sample count (assuming 120Hz)
  session.duration = session.sampleCount / 120
  
  return true
}

export async function addStroke(
  id: string,
  stroke: StrokeCycle,
  accPhase?: AccPhaseMetrics,
  decPhase?: DecPhaseMetrics,
  strokePhase?: StrokePhaseMetrics
): Promise<boolean> {
  const session = sessions.get(id)
  if (!session) return false
  
  session.strokes.push(stroke)
  if (accPhase) session.accMetrics.push(accPhase)
  if (decPhase) session.decMetrics.push(decPhase)
  if (strokePhase) session.strokePhaseMetrics.push(strokePhase)
  
  session.strokeCount = session.strokes.length
  session.updatedAt = new Date().toISOString()
  
  return true
}

export async function replaceSessionAnalysis(
  id: string,
  strokes: StrokeCycle[],
): Promise<boolean> {
  const session = sessions.get(id)
  if (!session) return false
  session.strokes = strokes
  session.strokeCount = strokes.length
  session.updatedAt = new Date().toISOString()
  return true
}

export async function deleteSession(id: string): Promise<boolean> {
  return sessions.delete(id)
}

export async function getSessionStats(id: string): Promise<SummaryStats | null> {
  const session = sessions.get(id)
  if (!session || session.strokes.length === 0) return null
  
  const strokes = session.strokes
  const accMetrics = session.accMetrics
  const strokePhaseMetrics = session.strokePhaseMetrics
  
  // Calculate SPM (strokes per minute)
  const totalTime = strokes.reduce((sum, s) => sum + s.TotalTime, 0)
  const avgSPM = totalTime > 0 ? (strokes.length / totalTime) * 60 : 0
  
  // Calculate average velocity changes
  const avgVelocityAcc = accMetrics.length > 0
    ? accMetrics.reduce((sum, m) => sum + (m.VelocityChange || 0), 0) / accMetrics.length
    : 0
  
  const avgVelocityDec = session.decMetrics.length > 0
    ? session.decMetrics.reduce((sum, m) => sum + Math.abs(m.VelocityChange || 0), 0) / session.decMetrics.length
    : 0
  
  // Calculate MAC values from stroke phase metrics
  const avgMACAccX = strokePhaseMetrics.length > 0
    ? strokePhaseMetrics.reduce((sum, m) => sum + (m.MeanAbsAccX || 0), 0) / strokePhaseMetrics.length
    : 0
  
  const avgMACGyrX = strokePhaseMetrics.length > 0
    ? strokePhaseMetrics.reduce((sum, m) => sum + (m.MeanAbsGyrX || 0), 0) / strokePhaseMetrics.length
    : 0
  
  const avgMACGyrY = strokePhaseMetrics.length > 0
    ? strokePhaseMetrics.reduce((sum, m) => sum + (m.MeanAbsGyrY || 0), 0) / strokePhaseMetrics.length
    : 0
  
  const avgMACGyrZ = strokePhaseMetrics.length > 0
    ? strokePhaseMetrics.reduce((sum, m) => sum + (m.MeanAbsGyrZ || 0), 0) / strokePhaseMetrics.length
    : 0
  
  return {
    totalStrokes: strokes.length,
    avgSPM,
    avgVelocityAcc,
    avgVelocityDec,
    avgMACAccX,
    avgMACGyrX,
    avgMACGyrY,
    avgMACGyrZ,
  }
}

function sessionToPublic(session: SessionWithData): Session {
  const { rawData, strokes, accMetrics, decMetrics, strokePhaseMetrics, ...publicData } = session
  return publicData
}
