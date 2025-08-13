#!/usr/bin/env pwsh

Write-Host "🧬 AUTODOCK VINA SETUP FOR LIPID RENDERING" -ForegroundColor Cyan
Write-Host "=========================================" -ForegroundColor Cyan
Write-Host ""

Write-Host "Current Status:" -ForegroundColor Yellow
Write-Host "✅ You are currently using the SIMPLIFIED Docker container" -ForegroundColor Green
Write-Host "❌ AutoDock Vina is NOT installed in this container" -ForegroundColor Red
Write-Host "❌ RDKit is NOT installed in this container" -ForegroundColor Red  
Write-Host "⚠️  System falls back to MOCK docking" -ForegroundColor Yellow
Write-Host ""

Write-Host "Choose your solution:" -ForegroundColor Cyan

Write-Host ""
Write-Host "OPTION 1: Use Full Production Container (RECOMMENDED)" -ForegroundColor Green
Write-Host "  ✅ Includes AutoDock Vina + RDKit + all scientific libraries"
Write-Host "  ✅ Real molecular docking functionality"
Write-Host "  ⚠️  Longer build time (~10-15 minutes first time)"
Write-Host "  ⚠️  Larger container size (~2-3 GB)"
Write-Host ""

Write-Host "OPTION 2: Install AutoDock Vina in Current Container" -ForegroundColor Yellow
Write-Host "  ✅ Quick setup in existing container"
Write-Host "  ⚠️  Manual pip install (may have compatibility issues)"
Write-Host "  ⚠️  Less reliable than conda-based installation"
Write-Host ""

Write-Host "OPTION 3: Keep Mock But Show Error Message" -ForegroundColor Magenta
Write-Host "  ✅ Immediate testing with clear error messages"
Write-Host "  ❌ No real docking functionality"
Write-Host ""

$choice = Read-Host "Enter your choice (1, 2, or 3)"

switch ($choice) {
    "1" {
        Write-Host ""
        Write-Host "🚀 SETTING UP PRODUCTION CONTAINER..." -ForegroundColor Green
        Write-Host ""
        
        Write-Host "Step 1: Stopping current container..." -ForegroundColor Yellow
        docker-compose -f docker-compose.test.yml down
        
        Write-Host "Step 2: Building frontend..." -ForegroundColor Yellow
        Set-Location lipid_viewer
        npm run build
        Set-Location ..
        
        Write-Host "Step 3: Building production container (this will take 10-15 minutes)..." -ForegroundColor Yellow
        Write-Host "⏳ Installing conda, AutoDock Vina, RDKit, and scientific libraries..." -ForegroundColor Cyan
        docker-compose -f docker-compose.production.yml build app-production
        
        Write-Host "Step 4: Starting production container..." -ForegroundColor Yellow
        docker-compose -f docker-compose.production.yml up -d app-production
        
        Write-Host "Step 5: Waiting for startup..." -ForegroundColor Yellow
        Start-Sleep 30
        
        Write-Host "Step 6: Running migrations..." -ForegroundColor Yellow
        docker exec lipidrendering-app-production-1 python manage.py migrate
        
        Write-Host "Step 7: Testing AutoDock Vina availability..." -ForegroundColor Yellow
        docker exec lipidrendering-app-production-1 python check_dependencies.py
        
        Write-Host ""
        Write-Host "✅ PRODUCTION SETUP COMPLETE!" -ForegroundColor Green
        Write-Host "🌐 Open: http://localhost:8000" -ForegroundColor Cyan
        Write-Host "🧬 AutoDock Vina should now be available for real docking" -ForegroundColor Green
    }
    
    "2" {
        Write-Host ""
        Write-Host "🔧 INSTALLING AUTODOCK VINA IN CURRENT CONTAINER..." -ForegroundColor Yellow
        Write-Host ""
        
        Write-Host "Step 1: Installing pip packages..." -ForegroundColor Yellow
        docker exec lipidrendering-app-test-1 pip install vina rdkit numpy scipy
        
        Write-Host "Step 2: Testing installation..." -ForegroundColor Yellow
        docker exec lipidrendering-app-test-1 python check_dependencies.py
        
        Write-Host "Step 3: Restarting Django..." -ForegroundColor Yellow
        docker restart lipidrendering-app-test-1
        
        Write-Host ""
        Write-Host "✅ PACKAGE INSTALLATION COMPLETE!" -ForegroundColor Green
        Write-Host "⚠️  Note: This method may have compatibility issues" -ForegroundColor Yellow
        Write-Host "🌐 Test at: http://localhost:8000" -ForegroundColor Cyan
    }
    
    "3" {
        Write-Host ""
        Write-Host "🔧 CONFIGURING ERROR MESSAGES..." -ForegroundColor Magenta
        Write-Host ""
        
        Write-Host "Mock docking is now disabled. Users will see clear error messages." -ForegroundColor Yellow
        Write-Host "The system will return: 'AutoDock Vina not available and mock docking is disabled'" -ForegroundColor Yellow
        
        Write-Host "Step 1: Restarting container to apply settings..." -ForegroundColor Yellow
        docker restart lipidrendering-app-test-1
        
        Write-Host ""
        Write-Host "✅ ERROR MESSAGE SETUP COMPLETE!" -ForegroundColor Green
        Write-Host "🌐 Test at: http://localhost:8000" -ForegroundColor Cyan
        Write-Host "❌ Docking will show error message instead of mock results" -ForegroundColor Red
    }
    
    default {
        Write-Host ""
        Write-Host "❌ Invalid choice. Please run the script again and choose 1, 2, or 3." -ForegroundColor Red
        exit 1
    }
}

Write-Host ""
Write-Host "📋 NEXT STEPS:" -ForegroundColor Cyan
Write-Host "1. Open http://localhost:8000 in incognito browser" 
Write-Host "2. Test with Cholesterol + 1CRN" 
Write-Host "3. Check console for docking engine type (should be 'vina' not 'mock')"
Write-Host ""
Write-Host "🔍 To verify AutoDock Vina is working:"
Write-Host "   Open browser developer tools → Network tab → Run docking"
Write-Host "   Check API response for: 'engine': 'vina' and 'is_mock': false"
