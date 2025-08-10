# Docker Deployment Guide - Lipid Docking Platform

## 🐳 Single-Server Docker Architecture

This guide covers the Docker deployment for the **unified single-server architecture** where Django serves both the API and the React frontend.

## 📋 Quick Start

### Prerequisites
- Docker installed and running
- Docker Compose (optional, for easier management)

### 1. Build and Run (Production)
```bash
# Build and run production container
./scripts/docker_build.sh

# Or using PowerShell on Windows
./scripts/docker_build.ps1
```

### 2. Access the Application
- **Frontend**: http://localhost:8000/
- **API Documentation**: http://localhost:8000/api/info
- **Django Admin**: http://localhost:8000/admin/
- **Health Check**: http://localhost:8000/api/healthz

## 🏗️ Architecture Overview

### Multi-Stage Build Process
1. **Frontend Stage**: Builds React application using Node.js
2. **Backend Stage**: Sets up Django with the built frontend assets

### File Structure
```
├── server/
│   ├── Dockerfile              # Production container
│   ├── Dockerfile.dev           # Development container
│   ├── .dockerignore           # Docker ignore rules
│   └── requirements.txt        # Python dependencies
├── docker-compose.yml          # Container orchestration
├── scripts/
│   ├── docker_build.sh         # Linux/Mac build script
│   ├── docker_build.ps1        # Windows PowerShell script
│   └── validate_docker_setup.py # Setup validation
└── DOCKER_README.md            # This file
```

## 🚀 Deployment Options

### Option 1: Direct Docker Commands

#### Production Build
```bash
# Build the image
docker build -t lipid-docking:latest -f server/Dockerfile .

# Run the container
docker run -p 8000:8000 --name lipid-docking-app lipid-docking:latest
```

#### Development Build
```bash
# Build development image
docker build -t lipid-docking:dev -f server/Dockerfile.dev .

# Run with live reload
docker run -p 8000:8000 \
  -v $(pwd)/server:/app \
  -v $(pwd)/lipid_viewer/dist:/app/staticfiles/frontend \
  -e DEBUG=True \
  --name lipid-docking-dev \
  lipid-docking:dev
```

### Option 2: Docker Compose

#### Production
```bash
docker-compose up app
```

#### Development
```bash
docker-compose --profile dev up app-dev
```

### Option 3: Automated Scripts

#### Linux/Mac
```bash
# Production
./scripts/docker_build.sh

# Development
./scripts/docker_build.sh --environment development

# Build only (no run)
./scripts/docker_build.sh --no-build

# Verbose output
./scripts/docker_build.sh --verbose
```

#### Windows PowerShell
```powershell
# Production
./scripts/docker_build.ps1

# Development
./scripts/docker_build.ps1 -Environment development

# Build only
./scripts/docker_build.ps1 -NoBuild

# Verbose output
./scripts/docker_build.ps1 -Verbose
```

## 🔧 Configuration

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `DEBUG` | `False` | Django debug mode |
| `DJANGO_SETTINGS_MODULE` | `core.settings` | Django settings module |
| `PYTHONDONTWRITEBYTECODE` | `1` | Prevent .pyc files |
| `PYTHONUNBUFFERED` | `1` | Unbuffered output |

### Production Configuration
- **Server**: Gunicorn with 3 workers
- **Static Files**: Served by Whitenoise
- **Security**: Non-root user (`appuser`)
- **Health Check**: Automatic health monitoring

### Development Configuration
- **Server**: Django development server
- **Live Reload**: Volume mounts for code changes
- **Debug Mode**: Enabled by default

## 🧪 Testing Docker Setup

### Validate Configuration
```bash
python scripts/validate_docker_setup.py
```

This script checks:
- ✅ Dockerfile configuration
- ✅ Docker Compose setup
- ✅ Django settings
- ✅ Frontend build
- ✅ Requirements and dependencies

### Run Container Tests
```bash
# Build and test the container
docker build -t lipid-docking-test -f server/Dockerfile .
docker run --rm lipid-docking-test python -m pytest api/test_docker_deployment.py
```

## 🐛 Troubleshooting

### Common Issues

#### 1. Container Fails to Start
```bash
# Check logs
docker logs lipid-docking-app

# Common causes:
# - Missing frontend build: run `npm run build` in lipid_viewer/
# - Port already in use: check if another service uses port 8000
# - Permission issues: ensure Docker has proper permissions
```

#### 2. Frontend Not Loading
```bash
# Ensure frontend is built
cd lipid_viewer
npm run build

# Check static files in container
docker exec -it lipid-docking-app ls -la /app/staticfiles/
```

#### 3. API Endpoints Not Working
```bash
# Test health endpoint
curl http://localhost:8000/api/healthz

# Check Django settings
docker exec -it lipid-docking-app python manage.py check
```

#### 4. Docker Build Fails
```bash
# Clear Docker cache
docker system prune -a

# Build with no cache
docker build --no-cache -t lipid-docking:latest -f server/Dockerfile .
```

### Performance Optimization

#### Build Time
- Use `.dockerignore` to exclude unnecessary files
- Multi-stage build reduces final image size
- Layer caching optimizes rebuild times

#### Runtime
- Gunicorn with multiple workers for production
- Static files served efficiently by Whitenoise
- Health checks ensure container reliability

## 📊 Monitoring

### Health Checks
```bash
# Container health status
docker ps --format "table {{.Names}}\t{{.Status}}"

# Detailed health info
docker inspect lipid-docking-app | grep -A 10 Health
```

### Logs
```bash
# Follow logs
docker logs -f lipid-docking-app

# Last 100 lines
docker logs --tail 100 lipid-docking-app
```

### Resource Usage
```bash
# Container stats
docker stats lipid-docking-app

# Detailed resource info
docker exec -it lipid-docking-app ps aux
```

## 🔄 Development Workflow

### 1. Code Changes
```bash
# For development container with volume mounts:
# - Edit code in ./server/ or ./lipid_viewer/
# - Changes automatically reflected in container
# - Django auto-reloads on Python file changes
```

### 2. Frontend Changes
```bash
# Rebuild frontend
cd lipid_viewer
npm run build

# Restart container if using production image
docker restart lipid-docking-app
```

### 3. Dependency Changes
```bash
# Rebuild image after changing requirements.txt
docker build -t lipid-docking:latest -f server/Dockerfile .
```

## 🚢 Production Deployment

### Security Considerations
- ✅ Non-root user in container
- ✅ Minimal base image (python:slim)
- ✅ No sensitive data in image
- ✅ Health checks for reliability

### Scaling
```bash
# Multiple instances with different ports
docker run -p 8001:8000 --name lipid-docking-app-1 lipid-docking:latest
docker run -p 8002:8000 --name lipid-docking-app-2 lipid-docking:latest

# Or use Docker Compose scaling
docker-compose up --scale app=3
```

### Persistence
```bash
# Mount volume for database persistence
docker run -p 8000:8000 \
  -v lipid-docking-data:/app/data \
  --name lipid-docking-app \
  lipid-docking:latest
```

## 📚 Additional Resources

- [Django Deployment Checklist](https://docs.djangoproject.com/en/stable/howto/deployment/checklist/)
- [Docker Best Practices](https://docs.docker.com/develop/best-practices/)
- [Gunicorn Configuration](https://docs.gunicorn.org/en/stable/configure.html)
- [Whitenoise Documentation](http://whitenoise.evans.io/en/stable/)

## 🆘 Support

If you encounter issues:

1. **Check the validation script**: `python scripts/validate_docker_setup.py`
2. **Review container logs**: `docker logs lipid-docking-app`
3. **Test health endpoint**: `curl http://localhost:8000/api/healthz`
4. **Verify frontend build**: Check `lipid_viewer/dist/` directory exists

For additional help, ensure all tests pass:
```bash
cd server
python -m pytest api/test_docker_deployment.py -v
```
