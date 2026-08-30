-- Artifacts produced by the C++ stroke_analysis binary for each uploaded run.
CREATE TABLE IF NOT EXISTS stroke_run_artifacts (
  id BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
  run_id TEXT NOT NULL,
  kind TEXT NOT NULL,
  filename TEXT NOT NULL,
  relative_path TEXT NOT NULL,
  byte_size INTEGER NOT NULL DEFAULT 0,
  stroke_id INTEGER,
  created_at TIMESTAMPTZ NOT NULL DEFAULT now()
);

CREATE INDEX IF NOT EXISTS stroke_run_artifacts_run_id_idx ON stroke_run_artifacts (run_id);
CREATE UNIQUE INDEX IF NOT EXISTS stroke_run_artifacts_path_idx ON stroke_run_artifacts (relative_path);
