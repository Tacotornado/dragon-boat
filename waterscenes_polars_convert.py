import polars as pl

# Load Parquet
df = pl.read_parquet("all_imu_data.parquet")

# Save to CSV for the C++ program
df.write_csv("all_imu_data.csv")
print("Saved all_imu_data.csv")