# Deploying Stroke Analysis Dashboard on Ubuntu Server (No Sudo)

This guide walks you through deploying the dashboard on an Ubuntu server where you **do not have sudo access**. It assumes Docker is already installed and your user has permission to run it.

## Prerequisites

- SSH access to your server
- Docker available to your user (check with `docker --version`)
- At least 1GB RAM and 10GB disk space

---

## Step 1: Check Docker Access

```bash
ssh your-username@your-server-ip

# Check if Docker is available
docker --version
docker compose version
```

If Docker isn't available or you get permission errors, contact your system administrator to:
1. Install Docker
2. Add your user to the `docker` group: `sudo usermod -aG docker your-username`

---

## Step 2: Transfer the Project

**Option A: From v0 - Download and transfer**

1. In v0, click the three dots (top right) and select "Download ZIP"
2. On your local machine:
```bash
scp stroke-analysis.zip your-username@your-server-ip:~/
```

3. On the server:
```bash
mkdir -p ~/stroke-analysis
cd ~/stroke-analysis
unzip ~/stroke-analysis.zip
rm ~/stroke-analysis.zip
```

**Option B: Clone from GitHub**
```bash
git clone https://github.com/your-username/stroke-analysis.git ~/stroke-analysis
cd ~/stroke-analysis
```

---

## Step 3: Configure Environment

```bash
cd ~/stroke-analysis

# Create environment file
cp .env.example .env

# Edit if needed (optional)
nano .env
```

Default settings work out of the box. The TypeScript stroke processor runs by default (no C++ needed).

---

## Step 4: Deploy

```bash
cd ~/stroke-analysis

# Build and start (this takes 2-3 minutes first time)
docker compose up -d --build

# Watch the build progress
docker compose logs -f
```

Press `Ctrl+C` to exit logs once you see "Ready" or the server starts.

---

## Step 5: Access the Dashboard

Get your server's IP or use the hostname:
```bash
hostname -I
```

Open in your browser:
```
http://your-server-ip:3000
```

**If port 3000 is blocked**, ask your admin to open it, or change the port in `docker-compose.yml`:
```yaml
ports:
  - "8080:3000"  # Use port 8080 instead
```

Then rebuild: `docker compose up -d --build`

---

## Connecting Your Xsens IMU

### 1. Create a Session

Go to `http://your-server-ip:3000/sessions` and click "New Session"

Or via curl:
```bash
curl -X POST http://your-server-ip:3000/api/sessions \
  -H "Content-Type: application/json" \
  -d '{"name": "Morning Practice", "athleteName": "John Doe"}'
```

Response:
```json
{"id": "abc123", "name": "Morning Practice", ...}
```

### 2. Send IMU Data

From your Xsens integration, POST data to the ingest endpoint:

```bash
curl -X POST http://your-server-ip:3000/api/sessions/abc123/ingest \
  -H "Content-Type: application/json" \
  -d '{
    "samples": [
      {"sample": 0, "Acc_X": 0.1, "Acc_Y": 2.3, "Gyr_X": 10.5, "Gyr_Y": -5.2, "Gyr_Z": 3.1},
      {"sample": 1, "Acc_X": 0.2, "Acc_Y": 2.5, "Gyr_X": 11.0, "Gyr_Y": -5.0, "Gyr_Z": 3.0}
    ]
  }'
```

**Tips:**
- Send batches of 10-50 samples at 10-20 Hz for smooth real-time updates
- Sample rate should match your IMU (typically 120 Hz for Xsens)

### 3. View Live Dashboard

Open `http://your-server-ip:3000/sessions/abc123` to see real-time waveforms and stroke detection.

---

## Example: Python Script for Xsens Integration

```python
import requests
import time

SERVER = "http://your-server-ip:3000"
SESSION_ID = "abc123"  # From step 1

def send_imu_batch(samples):
    """Send a batch of IMU samples to the dashboard."""
    resp = requests.post(
        f"{SERVER}/api/sessions/{SESSION_ID}/ingest",
        json={"samples": samples}
    )
    return resp.json()

# Example: sending data from your Xsens callback
def on_xsens_data(acc_x, acc_y, gyr_x, gyr_y, gyr_z):
    global sample_buffer, sample_count
    
    sample_buffer.append({
        "sample": sample_count,
        "Acc_X": acc_x,
        "Acc_Y": acc_y,
        "Gyr_X": gyr_x,
        "Gyr_Y": gyr_y,
        "Gyr_Z": gyr_z
    })
    sample_count += 1
    
    # Send every 12 samples (~100ms at 120Hz)
    if len(sample_buffer) >= 12:
        send_imu_batch(sample_buffer)
        sample_buffer = []

sample_buffer = []
sample_count = 0
```

---

## Common Commands

```bash
cd ~/stroke-analysis

# Check status
docker compose ps

# View logs
docker compose logs -f

# Stop
docker compose down

# Restart
docker compose restart

# Rebuild after code changes
docker compose up -d --build

# Clean up disk space
docker system prune -f
```

---

## Troubleshooting

**"permission denied" when running docker:**
```bash
# Check if you're in the docker group
groups
# Should include "docker". If not, ask admin to run:
# sudo usermod -aG docker your-username
```

**Container won't start:**
```bash
docker compose logs stroke-analysis
```

**Port already in use:**
```bash
# Find what's using port 3000
ss -tlnp | grep 3000

# Change port in docker-compose.yml to something else (e.g., 3001)
```

**Can't connect from browser:**
- Check firewall rules with your admin
- Try a different port
- Make sure the container is running: `docker compose ps`

**Out of disk space:**
```bash
# Check disk usage
df -h

# Clean Docker
docker system prune -f
docker image prune -a -f
```

---

## URL Reference

| URL | Purpose |
|-----|---------|
| `/` | Main dashboard - upload CSV results |
| `/sessions` | Manage live recording sessions |
| `/sessions/{id}` | Real-time session view |
| `/statistics` | Detailed statistics overview |

## API Reference

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/api/sessions` | GET | List all sessions |
| `/api/sessions` | POST | Create new session |
| `/api/sessions/{id}` | GET | Get session details |
| `/api/sessions/{id}` | PATCH | Update session (start/stop) |
| `/api/sessions/{id}/ingest` | POST | Send IMU data |
| `/api/sessions/{id}/stream` | GET | SSE stream for live updates |
