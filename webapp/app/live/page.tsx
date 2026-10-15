import { LiveCapture } from "@/components/live-capture"

export const metadata = {
  title: "Live Capture | Stroke Analysis",
  description: "Connect an Xsens DOT and record IMU data in real time.",
}

export default async function LivePage({
  searchParams,
}: {
  searchParams: Promise<{ sessionId?: string }>
}) {
  const { sessionId } = await searchParams
  return <LiveCapture sessionId={sessionId} />
}
