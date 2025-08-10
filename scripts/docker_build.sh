#!/bin/bash
# Docker build script for Lipid Docking single-server deployment

set -e

ENVIRONMENT="production"
NO_BUILD=false
VERBOSE=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --environment)
            ENVIRONMENT="$2"
            shift 2
            ;;
        --no-build)
            NO_BUILD=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        *)
            echo "Unknown option $1"
            exit 1
            ;;
    esac
done

echo "🐳 Docker Build Script for Lipid Docking Platform"
echo "Environment: $ENVIRONMENT"

# Change to project root
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

# Build the React frontend first
echo "📦 Building React frontend..."
cd "lipid_viewer"

if [ -d "node_modules" ]; then
    if [ "$VERBOSE" = true ]; then
        echo "Running npm run build..."
    fi
    npm run build
else
    echo "⚠️  Node modules not found. Installing dependencies..."
    npm install
    npm run build
fi

echo "✅ Frontend build completed"

# Return to project root
cd "$PROJECT_ROOT"

if [ "$NO_BUILD" = false ]; then
    # Build Docker image
    echo "🔨 Building Docker image..."
    
    IMAGE_TAG="lipid-docking:latest"
    if [ "$ENVIRONMENT" = "development" ]; then
        DOCKER_FILE="server/Dockerfile.dev"
        IMAGE_TAG="lipid-docking:dev"
    else
        DOCKER_FILE="server/Dockerfile"
    fi
    
    BUILD_ARGS=("docker" "build" "-t" "$IMAGE_TAG" "-f" "$DOCKER_FILE" ".")
    
    if [ "$VERBOSE" = true ]; then
        BUILD_ARGS+=("--progress=plain")
        echo "Running: ${BUILD_ARGS[*]}"
    fi
    
    "${BUILD_ARGS[@]}"
    
    echo "✅ Docker image built successfully: $IMAGE_TAG"
fi

# Run the container
echo "🚀 Starting Docker container..."

RUN_ARGS=("docker" "run" "--rm" "-p" "8000:8000" "--name" "lipid-docking-app")

if [ "$ENVIRONMENT" = "development" ]; then
    RUN_ARGS+=("-v" "$PROJECT_ROOT/server:/app")
    RUN_ARGS+=("-v" "$PROJECT_ROOT/lipid_viewer/dist:/app/staticfiles/frontend")
    RUN_ARGS+=("-e" "DEBUG=True")
    IMAGE_TAG="lipid-docking:dev"
else
    RUN_ARGS+=("-e" "DEBUG=False")
    IMAGE_TAG="lipid-docking:latest"
fi

RUN_ARGS+=("$IMAGE_TAG")

if [ "$VERBOSE" = true ]; then
    echo "Running: ${RUN_ARGS[*]}"
fi

echo "🌐 Application will be available at: http://localhost:8000"
echo "🔍 API documentation at: http://localhost:8000/api/info"
echo "Press Ctrl+C to stop the container"

"${RUN_ARGS[@]}"
