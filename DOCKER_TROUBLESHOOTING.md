# Docker Deployment Troubleshooting Guide

## Issue: White Blank Page with Asset Loading Errors

### Symptoms
- Browser shows white blank page
- Console errors: `Failed to load resource: 404 (Not Found)`
- MIME type errors: `'text/html' is not a supported stylesheet MIME type`
- Browser trying to load outdated asset filenames

### Root Cause Analysis

1. **Asset Hash Mismatch**: Vite generates unique hashes for each build (e.g., `index-VADf9BTh.css`), but browser cache may have older HTML with different hashes (e.g., `index-D9tVEuBa.css`)

2. **Static File Serving**: Django static file configuration may not be properly serving frontend assets

3. **Docker Volume Mounting**: Frontend assets may not be properly copied/mounted in container

## Solution Steps

### Step 1: Verify Container Status

```bash
# Check if container is running
docker ps

# Check container logs
docker-compose -f docker-compose.test.yml logs app-test

# If container stopped, restart it
docker-compose -f docker-compose.test.yml up -d app-test
```

### Step 2: Clear Browser Cache

The most common cause is browser cache serving old HTML with outdated asset references.

**Chrome/Edge:**
- Press `Ctrl + Shift + R` (hard refresh)
- Or `F12` → Network tab → Disable cache checkbox → `Ctrl + R`

**Firefox:**
- Press `Ctrl + F5`
- Or `Ctrl + Shift + R`

**Alternative: Incognito/Private Mode**
- Open `http://localhost:8000` in incognito/private browsing mode

### Step 3: Verify Asset Files in Container

```bash
# Check if assets exist in container
docker exec lipidrendering-app-test-1 ls -la /app/staticfiles/frontend/assets/

# Check HTML content
docker exec lipidrendering-app-test-1 cat /app/staticfiles/frontend/index.html | grep assets

# Verify static files are served
curl -I http://localhost:8000/assets/index-VADf9BTh.css
```

### Step 4: Test API Endpoints

```bash
# Test health endpoint
curl http://localhost:8000/api/healthz

# Test capabilities
curl http://localhost:8000/api/dock/capabilities

# Test main page
curl -I http://localhost:8000/
```

### Step 5: Run Automated Test Script

```bash
# Run the comprehensive test script
python test_docker_deployment.py
```

## Common Issues and Fixes

### Issue 1: Assets Return 404

**Cause**: Static file URL routing not working

**Fix**: 
1. Check Django `STATIC_URL` and `STATIC_ROOT` settings
2. Verify `collectstatic` ran successfully
3. Check URL patterns in `core/urls.py`

```python
# In core/urls.py
urlpatterns += static('/assets/', document_root=os.path.join(str(settings.STATIC_ROOT), 'frontend', 'assets'))
```

### Issue 2: Wrong MIME Types

**Cause**: Django serving HTML error pages instead of static files

**Fix**:
1. Ensure static files exist at expected paths
2. Check Django static file configuration
3. Verify container file permissions

```bash
# Fix permissions if needed
docker exec lipidrendering-app-test-1 chmod -R 755 /app/staticfiles/
```

### Issue 3: Container Won't Start

**Cause**: Port conflicts, build errors, or configuration issues

**Fix**:
```bash
# Stop conflicting containers
docker stop $(docker ps -q --filter "expose=8000")

# Remove old containers
docker-compose -f docker-compose.test.yml down

# Rebuild and start fresh
docker-compose -f docker-compose.test.yml build app-test
docker-compose -f docker-compose.test.yml up -d app-test
```

### Issue 4: Frontend Build Outdated

**Cause**: Frontend assets not rebuilt before Docker build

**Fix**:
```bash
# Rebuild frontend
cd lipid_viewer
npm run build

# Rebuild Docker container
cd ..
docker-compose -f docker-compose.test.yml build app-test
docker-compose -f docker-compose.test.yml up -d app-test
```

## Manual Testing Checklist

### ✅ Container Health
- [ ] Container is running (`docker ps`)
- [ ] No error logs (`docker-compose logs`)
- [ ] Health endpoint responds (`/api/healthz`)

### ✅ Static Files
- [ ] Assets exist in container (`/app/staticfiles/frontend/assets/`)
- [ ] HTML references correct asset hashes
- [ ] Static URL routing works

### ✅ Frontend Loading
- [ ] Main page loads without 404s
- [ ] CSS files have correct MIME type (`text/css`)
- [ ] JS files have correct MIME type (`application/javascript`)
- [ ] Browser cache cleared

### ✅ API Functionality
- [ ] Health check works
- [ ] Capabilities endpoint works
- [ ] Mock docking enabled in CI mode

## Production Deployment Notes

For production deployment, consider these improvements:

1. **Use Production WSGI Server**:
   ```dockerfile
   CMD ["gunicorn", "--bind", "0.0.0.0:8000", "--workers", "3", "core.wsgi:application"]
   ```

2. **Serve Static Files with Nginx** (for high traffic):
   ```nginx
   location /static/ {
       alias /app/staticfiles/;
       expires 1y;
       add_header Cache-Control "public, immutable";
   }
   ```

3. **Use Production Settings**:
   ```python
   DJANGO_SETTINGS_MODULE=core.settings  # Not settings_ci
   DOCKING_ALLOW_MOCK=False             # For real docking
   ```

4. **Health Check Configuration**:
   ```yaml
   healthcheck:
     test: ["CMD", "python", "-c", "import requests; requests.get('http://localhost:8000/api/healthz', timeout=10)"]
     interval: 30s
     timeout: 10s
     retries: 3
   ```

## Quick Fix Commands

```bash
# Complete reset and restart
docker-compose -f docker-compose.test.yml down
docker-compose -f docker-compose.test.yml build app-test
docker-compose -f docker-compose.test.yml up -d app-test

# Wait for startup
sleep 10

# Test the deployment
python test_docker_deployment.py

# If successful, open browser to http://localhost:8000 and hard refresh
```

## Expected Working State

When everything is working correctly:

1. **Container Status**: `docker ps` shows healthy container
2. **API Response**: `curl http://localhost:8000/api/healthz` returns `{"status": "healthy"}`
3. **Main Page**: `curl http://localhost:8000` returns HTML with correct asset references
4. **Assets Load**: CSS/JS files return with correct MIME types
5. **Browser**: Shows React application without console errors

## Support

If issues persist after following this guide:

1. **Collect Debug Information**:
   ```bash
   # Container logs
   docker-compose -f docker-compose.test.yml logs app-test > container.log
   
   # Container file structure
   docker exec lipidrendering-app-test-1 find /app/staticfiles -type f > files.log
   
   # Network test
   python test_docker_deployment.py > test_results.log
   ```

2. **Check Browser Developer Tools**:
   - Network tab for failed requests
   - Console for JavaScript errors
   - Sources tab to verify asset loading

3. **Verify Docker Environment**:
   - Docker version: `docker --version`
   - Available memory: `docker system info`
   - Port conflicts: `netstat -tulpn | grep 8000`
