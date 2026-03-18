import polars as pl
from pathlib import Path

# ---------------------------
# Dataset path
# ---------------------------
DATA_DIR = Path("/mnt/MPLL/dragonboat/NAS/data/waterscenes/imu/imu")

if not DATA_DIR.exists():
    raise FileNotFoundError(f"{DATA_DIR} not found!")

pattern = str(DATA_DIR / "*.csv")

print("Scanning all CSV files...")

# ---------------------------
# Define schema
# ---------------------------
SCHEMA = {
    "timestamp": pl.Float64,
    "pitch": pl.Float64,
    "roll": pl.Float64,
    "yaw": pl.Float64,
    "angular_velocity_x": pl.Float64,
    "angular_velocity_y": pl.Float64,
    "angular_velocity_z": pl.Float64,
    "linear_acceleration_x": pl.Float64,
    "linear_acceleration_y": pl.Float64,
    "linear_acceleration_z": pl.Float64,
    "magnetic_field_x": pl.Float64,
    "magnetic_field_y": pl.Float64,
    "magnetic_field_z": pl.Float64,
}

# ---------------------------
# Scan ALL files at once
# ---------------------------
lf = (
    pl.scan_csv(
        pattern,
        schema=SCHEMA,
        null_values=["", "null", "NaN"],
        ignore_errors=True,
        include_file_paths="file_path",  # automatically add source file path
    )
)

# ---------------------------
# Extract file_id from filename
# ---------------------------
lf = lf.with_columns(
    pl.col("file_path")
    .str.extract(r"(\d+)\.csv$", 1)
    .cast(pl.Int32)
    .alias("file_id")
).drop("file_path")

# ---------------------------
# Reorder columns (same as dataset + file_id)
# ---------------------------
lf = lf.select([
    "file_id",
    "timestamp",
    "pitch",
    "roll",
    "yaw",
    "angular_velocity_x",
    "angular_velocity_y",
    "angular_velocity_z",
    "linear_acceleration_x",
    "linear_acceleration_y",
    "linear_acceleration_z",
    "magnetic_field_x",
    "magnetic_field_y",
    "magnetic_field_z",
])

# ---------------------------
# Save directly using streaming
# ---------------------------
lf.sink_parquet(
    "all_imu_data.parquet",
    compression="zstd"
)

print("Saved all_imu_data.parquet")

# ---------------------------
# Convert to CSV for C++
# ---------------------------
df = pl.read_parquet("all_imu_data.parquet")

df.write_csv("all_imu_data.csv")

print("Saved all_imu_data.csv")