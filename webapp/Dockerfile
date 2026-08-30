# All stages use Debian (glibc) so the C++ binary built here runs in the runner.

# Stage 1: Node dependencies
FROM node:22-bookworm-slim AS deps
WORKDIR /app
RUN corepack enable && corepack prepare pnpm@latest --activate
COPY package.json pnpm-lock.yaml* ./
RUN pnpm install --frozen-lockfile --ignore-scripts

# Stage 2: Next.js build
FROM node:22-bookworm-slim AS builder
WORKDIR /app
RUN corepack enable && corepack prepare pnpm@latest --activate
COPY --from=deps /app/node_modules ./node_modules
COPY . .
ENV NEXT_TELEMETRY_DISABLED=1
ENV NODE_ENV=production
RUN pnpm build

# Stage 3: C++ engine build (stroke_analysis + utilities)
FROM debian:bookworm-slim AS engine
WORKDIR /engine
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential \
      python3 \
      python3-dev \
      python3-numpy \
      python3-matplotlib \
      libxlnt-dev \
      libcxxopts-dev \
      ca-certificates \
    && rm -rf /var/lib/apt/lists/*
COPY engine/ ./
# Tolerate a missing utilities/ or vendor/ drop-in: the app falls back to the
# TypeScript path when /app/engine/stroke_analysis is absent.
RUN mkdir -p /engine-out \
    && (make all && cp stroke_analysis /engine-out/ && echo "engine built") \
       || echo "WARNING: C++ engine not built (missing sources or headers)"

# Stage 4: Runner
FROM node:22-bookworm-slim AS runner
WORKDIR /app

ENV NODE_ENV=production
ENV NEXT_TELEMETRY_DISABLED=1
# matplotlib-cpp must render headlessly inside the container.
ENV MPLBACKEND=Agg

# Runtime libraries the binary links against (python + matplotlib + xlnt).
RUN apt-get update && apt-get install -y --no-install-recommends \
      python3 \
      python3-numpy \
      python3-matplotlib \
      wget \
    && (apt-get install -y --no-install-recommends libxlnt1.5 \
        || apt-get install -y --no-install-recommends libxlnt-dev) \
    && rm -rf /var/lib/apt/lists/*

RUN groupadd --system --gid 1001 nodejs && useradd --system --uid 1001 --gid nodejs nextjs

COPY --from=builder /app/public ./public
COPY --from=builder /app/.next/standalone ./
COPY --from=builder /app/.next/static ./.next/static
COPY --from=engine /engine-out/ ./engine/

# data/ holds per-run plots and workbooks produced by the engine.
RUN mkdir -p /app/data/runs && chown -R nextjs:nodejs /app/data /app/engine

USER nextjs

EXPOSE 3000
ENV PORT=3000
ENV HOSTNAME="0.0.0.0"

CMD ["node", "server.js"]
