# 🔧 Docker Asset Loading Issue - Complete Solution

## 🎯 Problem Identified

**Browser Error Messages:**
```
Refused to apply style from 'http://localhost:8000/assets/index-D9tVEuBa.css' 
because its MIME type ('text/html') is not a supported stylesheet MIME type

GET http://localhost:8000/assets/index-tqUhe6Oj.js net::ERR_ABORTED 404 (Not Found)
```

**Root Cause:** Browser cache serving old HTML with outdated asset hashes.

## 📊 Analysis

### What Browser Requests (OLD):
- ❌ `index-D9tVEuBa.css` (doesn't exist)
- ❌ `index-tqUhe6Oj.js` (doesn't exist)

### What Actually Exists (NEW):
- ✅ `index-VADf9BTh.css` (current build)
- ✅ `index-aL4CIk45.js` (current build)

### Why This Happens:
1. **Vite Build Process**: Generates unique hashes for each build
2. **Browser Caching**: Cached old HTML with old asset references
3. **Container Mismatch**: Docker has new build, browser has old cache

## 🚀 IMMEDIATE SOLUTION

### Option 1: Quick Fix (Browser Cache Clear)
```bash
# Open browser to http://localhost:8000
# Press Ctrl + Shift + R (hard refresh)
# OR open in incognito/private mode
```

### Option 2: Complete Rebuild (If Option 1 fails)
```bash
# Run the automated fix script
fix_asset_issue.bat

# OR manually:
docker-compose -f docker-compose.test.yml down
cd lipid_viewer && npm run build && cd ..
docker-compose -f docker-compose.test.yml build --no-cache app-test
docker-compose -f docker-compose.test.yml up -d app-test
```

## 🛠️ SYSTEMATIC DEBUGGING STEPS

### Step 1: Check Container Status
```bash
docker-compose -f docker-compose.test.yml ps
docker-compose -f docker-compose.test.yml logs app-test
```

### Step 2: Verify Assets in Container
```bash
docker exec lipidrendering-app-test-1 ls -la /app/staticfiles/frontend/assets/
docker exec lipidrendering-app-test-1 cat /app/staticfiles/frontend/index.html | findstr "assets"
```

### Step 3: Test API Endpoints
```bash
curl http://localhost:8000/api/healthz
curl -I http://localhost:8000
```

### Step 4: Verify Browser Behavior
1. Open Browser Developer Tools (F12)
2. Go to Network tab
3. Reload page and check for 404 errors
4. Look at what assets are being requested vs. what exists

## 🔧 PERMANENT FIXES IMPLEMENTED

### 1. Cache-Busting Headers
Added to Django `ReactAppView`:
```python
@method_decorator(never_cache, name='dispatch')
class ReactAppView(TemplateView):
    def get(self, request, *args, **kwargs):
        response = super().get(request, *args, **kwargs)
        response['Cache-Control'] = 'no-cache, no-store, must-revalidate'
        response['Pragma'] = 'no-cache'
        response['Expires'] = '0'
        return response
```

### 2. Proper Static File Routing
```python
urlpatterns += static('/assets/', 
    document_root=os.path.join(str(settings.STATIC_ROOT), 'frontend', 'assets'))
```

### 3. Automated Fix Script
Created `fix_asset_issue.bat` for complete rebuild automation.

## 🧪 TESTING CHECKLIST

### ✅ Container Health
- [ ] Container is running (`docker ps`)
- [ ] No errors in logs (`docker logs`)
- [ ] Health endpoint responds (`/api/healthz`)

### ✅ Asset Verification
- [ ] Assets exist in container (`/app/staticfiles/frontend/assets/`)
- [ ] HTML references match actual files
- [ ] Static files served with correct MIME types

### ✅ Browser Testing
- [ ] Hard refresh clears cache (`Ctrl+Shift+R`)
- [ ] Incognito mode works
- [ ] No 404 errors in Network tab
- [ ] CSS/JS load correctly

## 🐛 COMMON ISSUES & SOLUTIONS

### Issue: Container Not Running
```bash
# Solution: Start container
docker-compose -f docker-compose.test.yml up -d app-test
```

### Issue: Port 8000 Occupied
```bash
# Solution: Find and kill process
netstat -ano | findstr :8000
taskkill /PID <PID> /F
```

### Issue: Build Assets Don't Match
```bash
# Solution: Rebuild everything
cd lipid_viewer
npm run build
cd ..
docker-compose -f docker-compose.test.yml build --no-cache app-test
```

### Issue: Static Files Not Served
```bash
# Solution: Check Django collectstatic
docker exec lipidrendering-app-test-1 python manage.py collectstatic --noinput
```

## 📈 VERIFICATION COMMANDS

### Test Main Page
```bash
curl -I http://localhost:8000
# Should return 200 OK with HTML content
```

### Test Actual Assets
```bash
# Get actual CSS filename from container
CSS_FILE=$(docker exec lipidrendering-app-test-1 ls /app/staticfiles/frontend/assets/ | grep '\.css')
curl -I http://localhost:8000/assets/$CSS_FILE
# Should return 200 OK with text/css content-type
```

### Test API
```bash
curl http://localhost:8000/api/healthz
# Should return {"status": "healthy"}
```

## 🎯 SUCCESS INDICATORS

When everything works correctly:

1. **Container**: Shows as "healthy" or "running"
2. **API**: `/api/healthz` returns JSON status
3. **Main Page**: Returns HTML without 404s
4. **Assets**: CSS/JS files load with correct MIME types
5. **Browser**: React app loads without console errors

## 📞 NEXT STEPS

1. **Run the fix script**: `fix_asset_issue.bat`
2. **Hard refresh browser**: `Ctrl + Shift + R`
3. **Test in incognito mode**: Verify without cache interference
4. **Check developer tools**: Confirm no 404 errors

If issues persist, the cache-busting headers will prevent this problem in future deployments.

## 🔄 PREVENTION

To prevent this issue in future:

1. **Always hard refresh** after Docker rebuilds
2. **Use incognito mode** for testing new deployments
3. **Clear Docker cache** when making significant changes
4. **Rebuild frontend** before rebuilding Docker container

The implemented cache-busting headers should eliminate this issue going forward!
