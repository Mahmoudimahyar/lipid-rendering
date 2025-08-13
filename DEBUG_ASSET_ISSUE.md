# 🚨 CRITICAL ASSET ISSUE DEBUG

## 🎯 Current Problem
Browser requesting OLD asset hashes that don't exist:
- ❌ `index-tqUhe6Oj.js` (404 error)
- ❌ `index-D9tVEuBa.css` (404 error)

## 🔍 Root Cause Analysis

The issue is that **Django is serving an old cached HTML file** with outdated asset references. This is happening because:

1. **Docker Volume Mount**: The container might be using old cached files
2. **Django Template Cache**: Django might be caching the old HTML
3. **Collectstatic Issue**: Old assets might not have been properly replaced

## 🚀 STEP-BY-STEP MANUAL FIX

### Step 1: Complete Container Cleanup
```cmd
docker-compose -f docker-compose.test.yml down --remove-orphans
docker system prune -f
```

### Step 2: Fresh Frontend Build
```cmd
cd lipid_viewer
rmdir /s dist
npm run build
dir dist\assets
cd ..
```

### Step 3: Verify Local Assets
Check that `lipid_viewer/dist/index.html` contains the correct asset hashes:
- Should reference `index-VADf9BTh.css`
- Should reference `index-aL4CIk45.js`

### Step 4: Clean Docker Build
```cmd
docker-compose -f docker-compose.test.yml build --no-cache app-test
```

### Step 5: Start and Verify Container
```cmd
docker-compose -f docker-compose.test.yml up -d app-test
timeout /t 15
docker exec lipidrendering-app-test-1 ls /app/staticfiles/frontend/assets/
docker exec lipidrendering-app-test-1 cat /app/staticfiles/frontend/index.html | findstr "assets"
```

### Step 6: Test Browser Access
1. Open incognito/private mode
2. Go to `http://localhost:8000`
3. Check Network tab for 404 errors

## 🔧 ALTERNATIVE DIAGNOSIS

If the manual steps don't work, run this diagnosis:

### Check What's Actually Being Served
```cmd
curl -s http://localhost:8000 | findstr "assets"
```

This will show what asset hashes Django is actually serving vs. what exists.

### Check Container Contents
```cmd
docker exec lipidrendering-app-test-1 find /app -name "index.html" -exec head -50 {} \;
```

This shows all index.html files in the container and their contents.

## 🎯 LIKELY SOLUTIONS

### Solution A: Template Cache Issue
The Django template cache might be serving old HTML. Add this to your Docker rebuild:

```dockerfile
# In Dockerfile, after copying files
RUN python manage.py collectstatic --noinput --clear
```

### Solution B: Volume Mount Issue
The `docker-compose.test.yml` might be mounting old volumes. Check for:

```yaml
volumes:
  - ./lipid_viewer/dist:/app/staticfiles/frontend
```

If this exists, it's overriding the container's built-in files with your local files.

### Solution C: Multiple HTML Files
There might be multiple `index.html` files in the container. Check:

```cmd
docker exec lipidrendering-app-test-1 find /app -name "index.html"
```

## 🚨 EMERGENCY WORKAROUND

If nothing else works, you can manually fix the asset references:

1. Get the correct asset hashes from your local build:
   ```cmd
   type lipid_viewer\dist\index.html | findstr "assets"
   ```

2. Update the container's HTML file:
   ```cmd
   docker exec -it lipidrendering-app-test-1 /bin/bash
   cd /app/staticfiles/frontend
   # Edit index.html to have correct asset references
   ```

## 🔍 DEBUGGING CHECKLIST

- [ ] Local `dist/` directory has correct asset hashes
- [ ] Container `/app/staticfiles/frontend/assets/` has same files as local
- [ ] Container HTML references match actual asset files
- [ ] Django serves the same HTML as in container
- [ ] Browser cache completely cleared (try different browser)
- [ ] No volume mounts overriding container files

## ⚡ QUICK TEST

Run this single command to test everything:

```cmd
echo "=== LOCAL ASSETS ===" && dir lipid_viewer\dist\assets && echo "=== CONTAINER ASSETS ===" && docker exec lipidrendering-app-test-1 ls /app/staticfiles/frontend/assets/ && echo "=== DJANGO SERVES ===" && curl -s http://localhost:8000 | findstr "assets"
```

This will show you exactly where the mismatch is occurring.
