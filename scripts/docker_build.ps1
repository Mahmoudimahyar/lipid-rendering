# Docker build script for Lipid Docking single-server deployment
param(
    [string]$Environment = "production",
    [switch]$NoBuild,
    [switch]$Verbose
)

$ErrorActionPreference = "Stop"

Write-Host "🐳 Docker Build Script for Lipid Docking Platform" -ForegroundColor Green
Write-Host "Environment: $Environment" -ForegroundColor Yellow

# Change to project root
$ProjectRoot = Split-Path -Parent $PSScriptRoot
Set-Location $ProjectRoot

# Build the React frontend first
Write-Host "📦 Building React frontend..." -ForegroundColor Blue
Set-Location "lipid_viewer"

if (Test-Path "node_modules") {
    if ($Verbose) { Write-Host "Running npm run build..." }
    npm run build
    if ($LASTEXITCODE -ne 0) {
        Write-Error "❌ Frontend build failed"
        exit 1
    }
} else {
    Write-Host "⚠️  Node modules not found. Installing dependencies..." -ForegroundColor Yellow
    npm install
    npm run build
    if ($LASTEXITCODE -ne 0) {
        Write-Error "❌ Frontend build failed"
        exit 1
    }
}

Write-Host "✅ Frontend build completed" -ForegroundColor Green

# Return to project root
Set-Location $ProjectRoot

if (-not $NoBuild) {
    # Build Docker image
    Write-Host "🔨 Building Docker image..." -ForegroundColor Blue
    
    $ImageTag = "lipid-docking:latest"
    if ($Environment -eq "development") {
        $DockerFile = "server/Dockerfile.dev"
        $ImageTag = "lipid-docking:dev"
    } else {
        $DockerFile = "server/Dockerfile"
    }
    
    $BuildArgs = @(
        "docker", "build",
        "-t", $ImageTag,
        "-f", $DockerFile,
        "."
    )
    
    if ($Verbose) {
        $BuildArgs += "--progress=plain"
        Write-Host "Running: $($BuildArgs -join ' ')" -ForegroundColor Gray
    }
    
    & $BuildArgs[0] $BuildArgs[1..($BuildArgs.Length-1)]
    
    if ($LASTEXITCODE -ne 0) {
        Write-Error "❌ Docker build failed"
        exit 1
    }
    
    Write-Host "✅ Docker image built successfully: $ImageTag" -ForegroundColor Green
}

# Run the container
Write-Host "🚀 Starting Docker container..." -ForegroundColor Blue

$RunArgs = @(
    "docker", "run",
    "--rm",
    "-p", "8000:8000",
    "--name", "lipid-docking-app"
)

if ($Environment -eq "development") {
    $RunArgs += @(
        "-v", "${ProjectRoot}/server:/app",
        "-v", "${ProjectRoot}/lipid_viewer/dist:/app/staticfiles/frontend",
        "-e", "DEBUG=True"
    )
    $ImageTag = "lipid-docking:dev"
} else {
    $RunArgs += @(
        "-e", "DEBUG=False"
    )
    $ImageTag = "lipid-docking:latest"
}

$RunArgs += $ImageTag

if ($Verbose) {
    Write-Host "Running: $($RunArgs -join ' ')" -ForegroundColor Gray
}

Write-Host "🌐 Application will be available at: http://localhost:8000" -ForegroundColor Cyan
Write-Host "🔍 API documentation at: http://localhost:8000/api/info" -ForegroundColor Cyan
Write-Host "Press Ctrl+C to stop the container" -ForegroundColor Yellow

& $RunArgs[0] $RunArgs[1..($RunArgs.Length-1)]
