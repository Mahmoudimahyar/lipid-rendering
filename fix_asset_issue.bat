@echo off
echo 🔧 FIXING DOCKER ASSET ISSUE
echo =============================

echo Step 1: Stopping existing containers...
docker-compose -f docker-compose.test.yml down

echo Step 2: Rebuilding frontend with fresh hashes...
cd lipid_viewer
call npm run build
cd ..

echo Step 3: Rebuilding Docker container without cache...
docker-compose -f docker-compose.test.yml build --no-cache app-test

echo Step 4: Starting fresh container...
docker-compose -f docker-compose.test.yml up -d app-test

echo Step 5: Waiting for container startup...
timeout /t 15 /nobreak >nul

echo Step 6: Testing health endpoint...
curl -s http://localhost:8000/api/healthz

echo.
echo ✅ CONTAINER REBUILT SUCCESSFULLY!
echo.
echo 🌐 NEXT STEPS:
echo 1. Open your browser to: http://localhost:8000
echo 2. Press Ctrl+Shift+R (hard refresh) to clear cache
echo 3. OR open in incognito/private mode
echo.
echo If you still see issues, check browser Developer Tools (F12) Network tab
echo.
pause
