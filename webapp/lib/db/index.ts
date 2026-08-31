import { Pool } from "pg"
import { drizzle } from "drizzle-orm/node-postgres"
import * as schema from "./schema"

const globalForDb = globalThis as unknown as { strokePool?: Pool }

export const pool = globalForDb.strokePool ?? new Pool({ connectionString: process.env.DATABASE_URL })
if (process.env.NODE_ENV !== "production") globalForDb.strokePool = pool

export const db = drizzle(pool, { schema })
