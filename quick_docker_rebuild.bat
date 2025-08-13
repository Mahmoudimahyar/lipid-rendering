@echo off
echo 🔧 QUICK DOCKER REBUILD TO FIX ASSETS
echo =====================================

echo Step 1: Stopping container...
docker-compose -f docker-compose.test.yml down

echo Step 2: Ensuring frontend is built...
cd lipid_viewer
call npm run build
cd ..

echo Step 3: Building only essential layers (fast)...
docker-compose -f docker-compose.test.yml build --no-cache --parallel app-test

echo Step 4: Starting container...
docker-compose -f docker-compose.test.yml up -d app-test

echo Step 5: Waiting for startup...
timeout /t 15 /nobreak >nul

echo Step 6: Running migrations...
docker exec lipidrendering-app-test-1 python manage.py migrate --settings=core.settings_ci

echo Step 7: Copying correct template...
docker exec lipidrendering-app-test-1 cp /app/staticfiles/frontend/index.html /app/templates/index.html

echo.
echo ✅ REBUILD COMPLETE!
echo Open browser in INCOGNITO mode to: http://localhost:8000
echo.
pause
