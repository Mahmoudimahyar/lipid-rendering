#!/usr/bin/env powershell
# Verify Docker deployment and debug asset issues

Write-Host "🔍 DOCKER DEPLOYMENT VERIFICATION" -ForegroundColor Cyan
Write-Host "=================================" -ForegroundColor Cyan

# Check container status
Write-Host "`n1. Container Status:" -ForegroundColor Yellow
$containers = docker ps --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}"
Write-Host $containers

# Check specific container
$ourContainer = docker ps --filter name=lipidrendering-app-test-1 --format "{{.Names}}"
if ($ourContainer) {
    Write-Host "   ✅ Our container is running: $ourContainer" -ForegroundColor Green
    
    # Get container logs
    Write-Host "`n2. Container Logs (last 10 lines):" -ForegroundColor Yellow
    docker logs $ourContainer --tail 10
    
    # Check assets in container
    Write-Host "`n3. Assets in Container:" -ForegroundColor Yellow
    $assets = docker exec $ourContainer ls /app/staticfiles/frontend/assets/ 2>$null
    if ($assets) {
        $assets.Split("`n") | ForEach-Object { 
            if ($_.Trim()) { Write-Host "   - $($_.Trim())" -ForegroundColor White }
        }
    } else {
        Write-Host "   ❌ Cannot access container assets" -ForegroundColor Red
    }
    
    # Check HTML in container
    Write-Host "`n4. HTML Asset References in Container:" -ForegroundColor Yellow
    $html = docker exec $ourContainer cat /app/staticfiles/frontend/index.html 2>$null
    if ($html) {
        $cssMatches = [regex]::Matches($html, '/assets/(index-[a-zA-Z0-9]+\.css)')
        $jsMatches = [regex]::Matches($html, '/assets/(index-[a-zA-Z0-9]+\.js)')
        
        foreach ($match in $cssMatches) {
            Write-Host "   CSS: $($match.Groups[1].Value)" -ForegroundColor Cyan
        }
        foreach ($match in $jsMatches) {
            Write-Host "   JS: $($match.Groups[1].Value)" -ForegroundColor Cyan
        }
    } else {
        Write-Host "   ❌ Cannot read container HTML" -ForegroundColor Red
    }
    
} else {
    Write-Host "   ❌ Container not running!" -ForegroundColor Red
    Write-Host "`n   All containers:" -ForegroundColor Yellow
    docker ps -a --format "table {{.Names}}\t{{.Status}}"
    exit
}

# Test API endpoints
Write-Host "`n5. API Testing:" -ForegroundColor Yellow
try {
    $health = Invoke-WebRequest -Uri "http://localhost:8000/api/healthz" -UseBasicParsing -TimeoutSec 5
    Write-Host "   ✅ Health: $($health.StatusCode)" -ForegroundColor Green
} catch {
    Write-Host "   ❌ Health: Failed - $($_.Exception.Message)" -ForegroundColor Red
}

try {
    $main = Invoke-WebRequest -Uri "http://localhost:8000" -UseBasicParsing -TimeoutSec 5
    Write-Host "   ✅ Main page: $($main.StatusCode)" -ForegroundColor Green
    
    # Check what HTML is being served
    $htmlContent = $main.Content
    $cssMatches = [regex]::Matches($htmlContent, '/assets/(index-[a-zA-Z0-9]+\.css)')
    $jsMatches = [regex]::Matches($htmlContent, '/assets/(index-[a-zA-Z0-9]+\.js)')
    
    Write-Host "`n6. HTML Being Served by Django:" -ForegroundColor Yellow
    foreach ($match in $cssMatches) {
        Write-Host "   CSS: $($match.Groups[1].Value)" -ForegroundColor Cyan
    }
    foreach ($match in $jsMatches) {
        Write-Host "   JS: $($match.Groups[1].Value)" -ForegroundColor Cyan
    }
    
} catch {
    Write-Host "   ❌ Main page: Failed - $($_.Exception.Message)" -ForegroundColor Red
}

# Test actual asset files
if ($cssMatches.Count -gt 0) {
    Write-Host "`n7. Testing Asset Files:" -ForegroundColor Yellow
    $cssFile = $cssMatches[0].Groups[1].Value
    try {
        $cssTest = Invoke-WebRequest -Uri "http://localhost:8000/assets/$cssFile" -UseBasicParsing -TimeoutSec 5
        Write-Host "   ✅ CSS ($cssFile): $($cssTest.StatusCode)" -ForegroundColor Green
    } catch {
        Write-Host "   ❌ CSS ($cssFile): $($_.Exception.Response.StatusCode)" -ForegroundColor Red
    }
}

if ($jsMatches.Count -gt 0) {
    $jsFile = $jsMatches[0].Groups[1].Value
    try {
        $jsTest = Invoke-WebRequest -Uri "http://localhost:8000/assets/$jsFile" -UseBasicParsing -TimeoutSec 5
        Write-Host "   ✅ JS ($jsFile): $($jsTest.StatusCode)" -ForegroundColor Green
    } catch {
        Write-Host "   ❌ JS ($jsFile): $($_.Exception.Response.StatusCode)" -ForegroundColor Red
    }
}

Write-Host "`n📋 SUMMARY:" -ForegroundColor Cyan
Write-Host "If you see ❌ for asset files, the issue is in the container build." -ForegroundColor White
Write-Host "If you see ✅ for asset files but browser still shows 404, it's browser cache." -ForegroundColor White
Write-Host "`nFor browser cache issues:" -ForegroundColor Yellow
Write-Host "1. Try incognito/private mode" -ForegroundColor White
Write-Host "2. Clear all browser data (Ctrl+Shift+Del)" -ForegroundColor White
Write-Host "3. Try different browser" -ForegroundColor White
