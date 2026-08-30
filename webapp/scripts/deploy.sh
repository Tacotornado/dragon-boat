#!/bin/bash
set -e

# ============================================
# Stroke Analysis Dashboard - Deployment Script
# No sudo required - assumes Docker is available
# ============================================

echo "=== Stroke Analysis Dashboard Deployment ==="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if Docker is installed and accessible
if ! command -v docker &> /dev/null; then
    echo -e "${RED}Error: Docker is not installed or not in PATH.${NC}"
    echo "Ask your system administrator to install Docker and add you to the docker group."
    exit 1
fi

# Check if user can run Docker without sudo
if ! docker ps &> /dev/null 2>&1; then
    echo -e "${RED}Error: Cannot connect to Docker daemon.${NC}"
    echo "Ask your admin to add you to the docker group:"
    echo "  sudo usermod -aG docker $USER"
    echo "Then log out and back in."
    exit 1
fi

# Check if Docker Compose is available
if ! docker compose version &> /dev/null; then
    echo -e "${RED}Error: Docker Compose is not available.${NC}"
    echo "Ask your admin to update Docker (Compose is included in recent versions)."
    exit 1
fi

echo -e "${GREEN}Docker is available and accessible.${NC}"
echo ""

# Create directories
mkdir -p engine
mkdir -p data

# Check for .env file
if [ ! -f ".env" ]; then
    echo -e "${YELLOW}No .env file found. Creating from .env.example...${NC}"
    if [ -f ".env.example" ]; then
        cp .env.example .env
        echo -e "${GREEN}Created .env file.${NC}"
    else
        echo "NODE_ENV=production" > .env
        echo "PORT=3000" >> .env
        echo "USE_CPP_ENGINE=false" >> .env
        echo -e "${GREEN}Created default .env file.${NC}"
    fi
fi

echo ""
echo "=== Building and Starting Services ==="
echo "This may take 2-3 minutes on first run..."
echo ""

# Build and start
docker compose up -d --build

# Wait for container to be ready
echo ""
echo "Waiting for container to start..."
sleep 5

# Check if running
if docker compose ps | grep -q "running"; then
    echo ""
    echo -e "${GREEN}=== Deployment Complete ===${NC}"
    echo ""
    
    # Try to get IP address
    IP=$(hostname -I 2>/dev/null | awk '{print $1}')
    if [ -z "$IP" ]; then
        IP="your-server-ip"
    fi
    
    echo "Dashboard is available at: http://${IP}:3000"
    echo ""
    echo "Quick commands:"
    echo "  View logs:     docker compose logs -f"
    echo "  Stop:          docker compose down"
    echo "  Restart:       docker compose restart"
    echo ""
else
    echo -e "${RED}Container may not have started correctly.${NC}"
    echo "Check logs with: docker compose logs"
fi
