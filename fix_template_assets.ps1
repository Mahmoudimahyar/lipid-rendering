# Permanent fix for Docker template asset issue

Write-Host "🔧 FIXING TEMPLATE ASSET REFERENCES" -ForegroundColor Cyan
Write-Host "===================================" -ForegroundColor Cyan

# Step 1: Copy correct HTML to templates
Write-Host "`n1. Copying correct HTML to templates directory..." -ForegroundColor Yellow
docker exec lipidrendering-app-test-1 cp /app/staticfiles/frontend/index.html /app/templates/index.html

# Step 2: Verify the fix
Write-Host "`n2. Verifying asset references..." -ForegroundColor Yellow
$templateAssets = docker exec lipidrendering-app-test-1 grep "assets/" /app/templates/index.html
$staticAssets = docker exec lipidrendering-app-test-1 grep "assets/" /app/staticfiles/frontend/index.html

Write-Host "   📄 Template HTML assets:" -ForegroundColor Cyan
$templateAssets | ForEach-Object { Write-Host "      $_" -ForegroundColor White }

Write-Host "   📄 Static HTML assets:" -ForegroundColor Cyan  
$staticAssets | ForEach-Object { Write-Host "      $_" -ForegroundColor White }

# Step 3: Test Django serves correct HTML
Write-Host "`n3. Testing what Django serves..." -ForegroundColor Yellow
try {
    $response = Invoke-WebRequest -Uri "http://localhost:8000" -UseBasicParsing
    $servedAssets = $response.Content | Select-String "assets/"
    Write-Host "   📄 Django serves:" -ForegroundColor Cyan
    $servedAssets | ForEach-Object { Write-Host "      $($_.Line.Trim())" -ForegroundColor White }
    
    if ($response.Content -match "index-VADf9BTh\.css") {
        Write-Host "   ✅ SUCCESS: Django serving correct CSS!" -ForegroundColor Green
    } else {
        Write-Host "   ❌ PROBLEM: Django still serving old CSS!" -ForegroundColor Red
    }
    
    if ($response.Content -match "index-aL4CIk45\.js") {
        Write-Host "   ✅ SUCCESS: Django serving correct JS!" -ForegroundColor Green
    } else {
        Write-Host "   ❌ PROBLEM: Django still serving old JS!" -ForegroundColor Red
    }
    
} catch {
    Write-Host "   ❌ Error testing Django: $($_.Exception.Message)" -ForegroundColor Red
}

Write-Host "`n🎯 INSTRUCTIONS:" -ForegroundColor Cyan
Write-Host "1. Open browser in incognito mode: http://localhost:8000" -ForegroundColor White
Write-Host "2. Hard refresh (Ctrl+Shift+R) to clear cache" -ForegroundColor White
Write-Host "3. Check Network tab - should see index-VADf9BTh.css and index-aL4CIk45.js" -ForegroundColor White
Write-Host "4. If still broken, try different browser" -ForegroundColor White
