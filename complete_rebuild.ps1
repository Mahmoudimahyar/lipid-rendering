#!/usr/bin/env powershell
# Complete Docker Rebuild Script
# Fixes asset hash mismatch issue

Write-Host "🔧 COMPLETE DOCKER REBUILD - FIXING ASSET ISSUE" -ForegroundColor Cyan
Write-Host "=================================================" -ForegroundColor Cyan

# Step 1: Clean up everything
Write-Host "`n1. Cleaning up existing containers and images..." -ForegroundColor Yellow
docker-compose -f docker-compose.test.yml down --remove-orphans
docker system prune -f

# Step 2: Rebuild frontend with fresh hashes
Write-Host "`n2. Rebuilding frontend with fresh asset hashes..." -ForegroundColor Yellow
Set-Location lipid_viewer
if (Test-Path "dist") {
    Remove-Item -Recurse -Force dist
    Write-Host "   Removed old dist directory" -ForegroundColor Green
}
npm run build
Write-Host "   Frontend build complete" -ForegroundColor Green

# Check what assets were generated
Write-Host "`n   📁 New assets generated:" -ForegroundColor Cyan
Get-ChildItem dist/assets/ | ForEach-Object { Write-Host "      - $($_.Name)" -ForegroundColor White }

# Check HTML asset references
$htmlContent = Get-Content dist/index.html -Raw
$cssMatches = [regex]::Matches($htmlContent, '/assets/(index-[a-zA-Z0-9]+\.css)')
$jsMatches = [regex]::Matches($htmlContent, '/assets/(index-[a-zA-Z0-9]+\.js)')

Write-Host "`n   📄 HTML references:" -ForegroundColor Cyan
foreach ($match in $cssMatches) {
    Write-Host "      CSS: $($match.Groups[1].Value)" -ForegroundColor White
}
foreach ($match in $jsMatches) {
    Write-Host "      JS: $($match.Groups[1].Value)" -ForegroundColor White
}

Set-Location ..

# Step 3: Completely rebuild Docker container
Write-Host "`n3. Rebuilding Docker container without cache..." -ForegroundColor Yellow
docker-compose -f docker-compose.test.yml build --no-cache app-test
Write-Host "   Docker build complete" -ForegroundColor Green

# Step 4: Start fresh container
Write-Host "`n4. Starting fresh container..." -ForegroundColor Yellow
docker-compose -f docker-compose.test.yml up -d app-test

# Step 5: Wait for startup and verify
Write-Host "`n5. Waiting for container startup..." -ForegroundColor Yellow
Start-Sleep -Seconds 20

# Check container status
$containerStatus = docker ps --filter name=lipidrendering-app-test-1 --format "{{.Status}}"
if ($containerStatus) {
    Write-Host "   ✅ Container status: $containerStatus" -ForegroundColor Green
} else {
    Write-Host "   ❌ Container not running!" -ForegroundColor Red
    Write-Host "`n   Container logs:" -ForegroundColor Yellow
    docker-compose -f docker-compose.test.yml logs app-test
    exit 1
}

# Step 6: Verify assets in container
Write-Host "`n6. Verifying assets in container..." -ForegroundColor Yellow
$containerAssets = docker exec lipidrendering-app-test-1 ls /app/staticfiles/frontend/assets/ 2>$null
if ($containerAssets) {
    Write-Host "   📁 Container assets:" -ForegroundColor Cyan
    $containerAssets.Split("`n") | ForEach-Object { 
        if ($_.Trim()) { Write-Host "      - $($_.Trim())" -ForegroundColor White }
    }
} else {
    Write-Host "   ❌ Cannot access container assets!" -ForegroundColor Red
}

# Step 7: Check container HTML
Write-Host "`n7. Checking container HTML references..." -ForegroundColor Yellow
$containerHtml = docker exec lipidrendering-app-test-1 cat /app/staticfiles/frontend/index.html 2>$null
if ($containerHtml) {
    $containerCssMatches = [regex]::Matches($containerHtml, '/assets/(index-[a-zA-Z0-9]+\.css)')
    $containerJsMatches = [regex]::Matches($containerHtml, '/assets/(index-[a-zA-Z0-9]+\.js)')
    
    Write-Host "   📄 Container HTML references:" -ForegroundColor Cyan
    foreach ($match in $containerCssMatches) {
        Write-Host "      CSS: $($match.Groups[1].Value)" -ForegroundColor White
    }
    foreach ($match in $containerJsMatches) {
        Write-Host "      JS: $($match.Groups[1].Value)" -ForegroundColor White
    }
} else {
    Write-Host "   ❌ Cannot read container HTML!" -ForegroundColor Red
}

# Step 8: Test API
Write-Host "`n8. Testing API health..." -ForegroundColor Yellow
try {
    $response = Invoke-WebRequest -Uri "http://localhost:8000/api/healthz" -UseBasicParsing -TimeoutSec 10
    if ($response.StatusCode -eq 200) {
        Write-Host "   ✅ API health: $($response.StatusCode)" -ForegroundColor Green
    } else {
        Write-Host "   ⚠️  API health: $($response.StatusCode)" -ForegroundColor Yellow
    }
} catch {
    Write-Host "   ❌ API health check failed: $($_.Exception.Message)" -ForegroundColor Red
}

# Step 9: Test asset serving
Write-Host "`n9. Testing asset serving..." -ForegroundColor Yellow
if ($containerCssMatches.Count -gt 0) {
    $cssFile = $containerCssMatches[0].Groups[1].Value
    try {
        $cssResponse = Invoke-WebRequest -Uri "http://localhost:8000/assets/$cssFile" -UseBasicParsing -TimeoutSec 5
        $contentType = $cssResponse.Headers.'Content-Type'
        Write-Host "   ✅ CSS file ($cssFile): $($cssResponse.StatusCode), MIME: $contentType" -ForegroundColor Green
    } catch {
        Write-Host "   ❌ CSS file ($cssFile): Failed - $($_.Exception.Message)" -ForegroundColor Red
    }
}

# Final instructions
Write-Host "`n🎉 REBUILD COMPLETE!" -ForegroundColor Green
Write-Host "=" * 50 -ForegroundColor Green
Write-Host "`n📋 NEXT STEPS:" -ForegroundColor Cyan
Write-Host "1. Open browser to: http://localhost:8000" -ForegroundColor White
Write-Host "2. Press Ctrl+Shift+Del to open Clear browsing data" -ForegroundColor White
Write-Host "3. Clear 'Cached images and files' for 'All time'" -ForegroundColor White
Write-Host "4. OR use incognito/private browsing mode" -ForegroundColor White
Write-Host "5. Check browser Developer Tools (F12) Network tab" -ForegroundColor White
Write-Host "`n⚠️  If you still see old asset names, your browser has aggressive caching." -ForegroundColor Yellow
Write-Host "   Try different browser or clear all browser data." -ForegroundColor Yellow

Write-Host "`nPress any key to continue..." -ForegroundColor Gray
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
