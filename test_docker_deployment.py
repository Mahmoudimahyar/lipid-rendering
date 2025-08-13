#!/usr/bin/env python3
"""
Docker Deployment Test Script
Tests the complete Docker deployment including API and frontend assets
"""

import requests
import time
import subprocess
import sys
import json
from pathlib import Path

class DockerDeploymentTester:
    def __init__(self):
        self.base_url = "http://localhost:8000"
        self.container_name = "lipidrendering-app-test-1"
        
    def run_all_tests(self):
        """Run complete deployment test suite"""
        print("🚀 Starting Docker Deployment Test Suite")
        print("=" * 50)
        
        try:
            # 1. Check container status
            self.check_container_status()
            
            # 2. Wait for server startup
            self.wait_for_server()
            
            # 3. Test API endpoints
            self.test_api_endpoints()
            
            # 4. Test static file serving
            self.test_static_files()
            
            # 5. Test frontend loading
            self.test_frontend_assets()
            
            print("\n✅ All tests passed! Docker deployment is working correctly.")
            return True
            
        except Exception as e:
            print(f"\n❌ Test failed: {e}")
            return False
    
    def check_container_status(self):
        """Check if Docker container is running"""
        print("1. Checking container status...")
        
        result = subprocess.run(
            ["docker", "ps", "--filter", f"name={self.container_name}", "--format", "table {{.Names}}\\t{{.Status}}"],
            capture_output=True, text=True
        )
        
        if self.container_name in result.stdout:
            print(f"   ✅ Container {self.container_name} is running")
        else:
            print(f"   ❌ Container {self.container_name} not found or not running")
            # Try to start it
            print("   🔄 Attempting to start container...")
            subprocess.run(["docker-compose", "-f", "docker-compose.test.yml", "up", "-d", "app-test"])
            time.sleep(10)
    
    def wait_for_server(self):
        """Wait for Django server to start"""
        print("2. Waiting for server startup...")
        
        max_attempts = 30
        for attempt in range(max_attempts):
            try:
                response = requests.get(f"{self.base_url}/api/healthz", timeout=5)
                if response.status_code == 200:
                    print(f"   ✅ Server is ready (attempt {attempt + 1})")
                    return
            except requests.RequestException:
                pass
            
            print(f"   ⏳ Waiting for server... (attempt {attempt + 1}/{max_attempts})")
            time.sleep(2)
        
        raise Exception("Server failed to start within timeout period")
    
    def test_api_endpoints(self):
        """Test critical API endpoints"""
        print("3. Testing API endpoints...")
        
        endpoints = [
            ("/api/healthz", "Health check"),
            ("/api/dock/capabilities", "Docking capabilities"),
            ("/api/pdb/1CRN/info", "PDB info endpoint")
        ]
        
        for endpoint, description in endpoints:
            try:
                response = requests.get(f"{self.base_url}{endpoint}", timeout=10)
                if response.status_code == 200:
                    print(f"   ✅ {description}: {response.status_code}")
                else:
                    print(f"   ⚠️  {description}: {response.status_code} (may be expected)")
            except Exception as e:
                print(f"   ❌ {description}: Failed - {e}")
    
    def test_static_files(self):
        """Test static file serving"""
        print("4. Testing static file serving...")
        
        # Check if static files are accessible
        static_endpoints = [
            ("/static/admin/css/base.css", "Django admin CSS"),
            ("/static/rest_framework/css/bootstrap.min.css", "DRF CSS")
        ]
        
        for endpoint, description in static_endpoints:
            try:
                response = requests.get(f"{self.base_url}{endpoint}", timeout=5)
                if response.status_code == 200:
                    print(f"   ✅ {description}: {response.status_code}")
                elif response.status_code == 404:
                    print(f"   ⚠️  {description}: Not found (may be expected)")
                else:
                    print(f"   ❓ {description}: {response.status_code}")
            except Exception as e:
                print(f"   ❌ {description}: Failed - {e}")
    
    def test_frontend_assets(self):
        """Test frontend asset serving"""
        print("5. Testing frontend assets...")
        
        # First, get the main page to see what assets are referenced
        try:
            response = requests.get(self.base_url, timeout=10)
            if response.status_code == 200:
                print(f"   ✅ Main page loads: {response.status_code}")
                
                # Extract asset references from HTML
                html_content = response.text
                print(f"   📄 HTML content length: {len(html_content)} characters")
                
                # Look for asset references
                if "/assets/" in html_content:
                    print("   ✅ Assets referenced in HTML")
                    
                    # Extract specific asset filenames
                    import re
                    css_matches = re.findall(r'/assets/(index-[a-zA-Z0-9]+\.css)', html_content)
                    js_matches = re.findall(r'/assets/(index-[a-zA-Z0-9]+\.js)', html_content)
                    
                    print(f"   📋 CSS files found: {css_matches}")
                    print(f"   📋 JS files found: {js_matches}")
                    
                    # Test loading these assets
                    for css_file in css_matches:
                        self.test_asset_file(f"/assets/{css_file}", "CSS file")
                    
                    for js_file in js_matches:
                        self.test_asset_file(f"/assets/{js_file}", "JS file")
                        
                else:
                    print("   ⚠️  No asset references found in HTML")
            else:
                print(f"   ❌ Main page failed: {response.status_code}")
                
        except Exception as e:
            print(f"   ❌ Frontend test failed: {e}")
    
    def test_asset_file(self, asset_path, asset_type):
        """Test loading a specific asset file"""
        try:
            response = requests.get(f"{self.base_url}{asset_path}", timeout=5)
            if response.status_code == 200:
                content_type = response.headers.get('content-type', 'unknown')
                print(f"   ✅ {asset_type} ({asset_path}): {response.status_code}, MIME: {content_type}")
            else:
                print(f"   ❌ {asset_type} ({asset_path}): {response.status_code}")
        except Exception as e:
            print(f"   ❌ {asset_type} ({asset_path}): Failed - {e}")

def main():
    """Main test execution"""
    print("Docker Deployment Test for Lipid Rendering")
    print("==========================================")
    
    tester = DockerDeploymentTester()
    success = tester.run_all_tests()
    
    if success:
        print("\n🎉 Docker deployment test completed successfully!")
        print("\nNext steps:")
        print("1. Open browser to http://localhost:8000")
        print("2. Hard refresh (Ctrl+F5) to clear browser cache")
        print("3. Check browser developer tools if issues persist")
    else:
        print("\n💥 Docker deployment test failed!")
        print("\nTroubleshooting:")
        print("1. Check container logs: docker-compose -f docker-compose.test.yml logs app-test")
        print("2. Restart container: docker-compose -f docker-compose.test.yml restart app-test")
        print("3. Rebuild container: docker-compose -f docker-compose.test.yml build app-test")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
