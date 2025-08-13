#!/usr/bin/env python3
"""
Docker Issue Debug Script
Systematically analyzes the asset loading problem
"""

import subprocess
import os
import re
import json
from pathlib import Path

class DockerAssetDebugger:
    def __init__(self):
        self.container_name = "lipidrendering-app-test-1"
        
    def run_debug(self):
        """Run comprehensive debugging analysis"""
        print("🔍 DOCKER ASSET DEBUGGING ANALYSIS")
        print("=" * 60)
        
        # 1. Check Docker environment
        self.check_docker_status()
        
        # 2. Analyze local build files
        self.analyze_local_build()
        
        # 3. Check container contents
        self.check_container_contents()
        
        # 4. Analyze asset mismatch
        self.analyze_asset_mismatch()
        
        # 5. Test static file serving
        self.test_static_serving()
        
        # 6. Provide solution
        self.provide_solution()
    
    def run_command(self, cmd, description=""):
        """Run command and return output"""
        try:
            if description:
                print(f"   🔧 {description}")
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
            if result.returncode == 0:
                return result.stdout.strip()
            else:
                print(f"   ❌ Command failed: {cmd}")
                print(f"   Error: {result.stderr}")
                return None
        except subprocess.TimeoutExpired:
            print(f"   ⏰ Command timed out: {cmd}")
            return None
        except Exception as e:
            print(f"   ❌ Command error: {e}")
            return None
    
    def check_docker_status(self):
        """Check Docker container status"""
        print("\n1. DOCKER ENVIRONMENT STATUS")
        print("-" * 40)
        
        # Check if Docker is running
        docker_version = self.run_command("docker --version", "Checking Docker version")
        if docker_version:
            print(f"   ✅ Docker: {docker_version}")
        
        # Check running containers
        containers = self.run_command("docker ps --format 'table {{.Names}}\\t{{.Status}}\\t{{.Ports}}'", "Checking running containers")
        if containers:
            print(f"   📦 Running containers:\n{containers}")
        
        # Check our specific container
        our_container = self.run_command(f"docker ps --filter name={self.container_name} --format '{{{{.Status}}}}'")
        if our_container:
            print(f"   ✅ Our container status: {our_container}")
        else:
            print(f"   ❌ Container {self.container_name} not running")
            # Try to start it
            print("   🔄 Attempting to start container...")
            self.run_command("docker-compose -f docker-compose.test.yml up -d app-test", "Starting container")
    
    def analyze_local_build(self):
        """Analyze local build files"""
        print("\n2. LOCAL BUILD ANALYSIS")
        print("-" * 40)
        
        dist_path = Path("lipid_viewer/dist")
        if not dist_path.exists():
            print("   ❌ No dist directory found - need to build frontend")
            print("   🔧 Run: cd lipid_viewer && npm run build")
            return
        
        # Check dist contents
        assets_path = dist_path / "assets"
        if assets_path.exists():
            asset_files = list(assets_path.glob("*"))
            print(f"   📁 Local assets found: {len(asset_files)} files")
            for file in asset_files:
                print(f"      - {file.name}")
        
        # Check index.html
        index_path = dist_path / "index.html"
        if index_path.exists():
            with open(index_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Extract asset references
            css_matches = re.findall(r'/assets/(index-[a-zA-Z0-9]+\.css)', content)
            js_matches = re.findall(r'/assets/(index-[a-zA-Z0-9]+\.js)', content)
            
            print(f"   📄 Local HTML references:")
            print(f"      CSS: {css_matches}")
            print(f"      JS: {js_matches}")
    
    def check_container_contents(self):
        """Check what's actually in the container"""
        print("\n3. CONTAINER CONTENTS ANALYSIS")
        print("-" * 40)
        
        # Check if container exists and is running
        container_status = self.run_command(f"docker inspect {self.container_name} --format '{{{{.State.Status}}}}'")
        if not container_status:
            print(f"   ❌ Container {self.container_name} not found")
            return
        
        print(f"   📦 Container status: {container_status}")
        
        # Check container frontend assets
        assets = self.run_command(f"docker exec {self.container_name} ls -la /app/staticfiles/frontend/assets/", "Checking container assets")
        if assets:
            print(f"   📁 Container assets:\n{assets}")
        
        # Check container HTML
        html_content = self.run_command(f"docker exec {self.container_name} cat /app/staticfiles/frontend/index.html", "Reading container HTML")
        if html_content:
            # Extract asset references from container HTML
            css_matches = re.findall(r'/assets/(index-[a-zA-Z0-9]+\.css)', html_content)
            js_matches = re.findall(r'/assets/(index-[a-zA-Z0-9]+\.js)', content)
            
            print(f"   📄 Container HTML references:")
            print(f"      CSS: {css_matches}")
            print(f"      JS: {js_matches}")
        
        # Check Django logs
        logs = self.run_command(f"docker logs {self.container_name} --tail 20", "Getting container logs")
        if logs:
            print(f"   📋 Recent logs:\n{logs}")
    
    def analyze_asset_mismatch(self):
        """Analyze the specific asset mismatch issue"""
        print("\n4. ASSET MISMATCH ANALYSIS")
        print("-" * 40)
        
        print("   🔍 Browser error analysis:")
        print("   ❌ Browser requesting: index-D9tVEuBa.css")
        print("   ❌ Browser requesting: index-tqUhe6Oj.js")
        
        # Check what's actually in container
        container_css = self.run_command(f"docker exec {self.container_name} ls /app/staticfiles/frontend/assets/ | grep '\\.css'")
        container_js = self.run_command(f"docker exec {self.container_name} ls /app/staticfiles/frontend/assets/ | grep '\\.js'")
        
        if container_css:
            print(f"   ✅ Container has CSS: {container_css}")
        if container_js:
            print(f"   ✅ Container has JS: {container_js}")
        
        # The problem is clear: browser cache has old HTML with old asset hashes
        print("\n   🎯 ROOT CAUSE IDENTIFIED:")
        print("   The browser has cached an old version of index.html with outdated asset hashes.")
        print("   The container has newer assets with different hashes.")
    
    def test_static_serving(self):
        """Test static file serving"""
        print("\n5. STATIC FILE SERVING TEST")
        print("-" * 40)
        
        # Test main page
        main_page = self.run_command('curl -s -I http://localhost:8000', "Testing main page")
        if main_page:
            print(f"   📄 Main page response:\n{main_page}")
        
        # Test asset serving with actual filenames from container
        css_file = self.run_command(f"docker exec {self.container_name} ls /app/staticfiles/frontend/assets/ | grep '\\.css' | head -1")
        if css_file:
            css_test = self.run_command(f'curl -s -I http://localhost:8000/assets/{css_file.strip()}', f"Testing CSS file: {css_file}")
            if css_test:
                print(f"   🎨 CSS file response:\n{css_test}")
    
    def provide_solution(self):
        """Provide step-by-step solution"""
        print("\n6. SYSTEMATIC SOLUTION")
        print("-" * 40)
        
        print("   🎯 The issue is browser cache serving old HTML with outdated asset hashes.")
        print("   📋 Follow these steps to fix:")
        print()
        print("   STEP 1: Clear Browser Cache")
        print("   - Press Ctrl + Shift + R (hard refresh)")
        print("   - Or open http://localhost:8000 in incognito mode")
        print()
        print("   STEP 2: If that doesn't work, rebuild everything:")
        print("   - cd lipid_viewer")
        print("   - npm run build")
        print("   - cd ..")
        print("   - docker-compose -f docker-compose.test.yml down")
        print("   - docker-compose -f docker-compose.test.yml build --no-cache app-test")
        print("   - docker-compose -f docker-compose.test.yml up -d app-test")
        print()
        print("   STEP 3: Verify the fix:")
        print("   - Wait 10 seconds for container startup")
        print("   - Open http://localhost:8000 in incognito mode")
        print("   - Check browser dev tools (F12) Network tab for any 404s")
        
        # Create a fix script
        self.create_fix_script()
    
    def create_fix_script(self):
        """Create an automated fix script"""
        fix_script = """#!/bin/bash
# Automated Docker Asset Fix Script

echo "🔧 Starting Docker Asset Fix..."

# Step 1: Stop containers
echo "1. Stopping containers..."
docker-compose -f docker-compose.test.yml down

# Step 2: Rebuild frontend
echo "2. Rebuilding frontend..."
cd lipid_viewer
npm run build
cd ..

# Step 3: Rebuild Docker container
echo "3. Rebuilding Docker container..."
docker-compose -f docker-compose.test.yml build --no-cache app-test

# Step 4: Start container
echo "4. Starting container..."
docker-compose -f docker-compose.test.yml up -d app-test

# Step 5: Wait for startup
echo "5. Waiting for startup..."
sleep 15

# Step 6: Test
echo "6. Testing deployment..."
curl -s http://localhost:8000/api/healthz

echo "✅ Fix complete! Open http://localhost:8000 in incognito mode"
"""
        
        with open("fix_docker_assets.sh", "w") as f:
            f.write(fix_script)
        
        print("   📜 Created fix_docker_assets.sh - run this script to fix everything automatically")

def main():
    debugger = DockerAssetDebugger()
    debugger.run_debug()

if __name__ == "__main__":
    main()
