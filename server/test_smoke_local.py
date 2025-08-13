#!/usr/bin/env python3
"""
Local Smoke Test Validation

This script runs the smoke test locally to validate it works
before committing to CI. It can run with or without scientific libraries.
"""

import os
import sys
import subprocess
import time
import signal
from pathlib import Path

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Setup Django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'core.settings')
import django
django.setup()

from api.test_smoke_seeded import SmokeTestRunner


class LocalSmokeTestValidator:
    """Runs and validates smoke test in local environment"""
    
    def __init__(self):
        self.django_process = None
        self.server_port = 8001  # Use different port to avoid conflicts
        
    def start_django_server(self):
        """Start Django development server"""
        print("🚀 Starting Django development server...")
        
        # Run migrations first
        print("   Running migrations...")
        subprocess.run([
            sys.executable, 'manage.py', 'migrate', '--run-syncdb'
        ], check=True)
        
        # Start server
        self.django_process = subprocess.Popen([
            sys.executable, 'manage.py', 'runserver', f'127.0.0.1:{self.server_port}'
        ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Wait for server to start
        print("   Waiting for server startup...")
        time.sleep(8)
        
        # Check if server is running
        try:
            import requests
            response = requests.get(f'http://127.0.0.1:{self.server_port}/api/healthz', timeout=5)
            if response.status_code == 200:
                print("✅ Django server started successfully")
                return True
            else:
                print(f"❌ Server health check failed: {response.status_code}")
                return False
        except Exception as e:
            print(f"❌ Failed to connect to server: {e}")
            return False
    
    def stop_django_server(self):
        """Stop Django development server"""
        if self.django_process:
            print("🛑 Stopping Django server...")
            self.django_process.terminate()
            try:
                self.django_process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.django_process.kill()
                self.django_process.wait()
            print("✅ Django server stopped")
    
    def check_scientific_libraries(self):
        """Check which scientific libraries are available"""
        libraries = {}
        
        try:
            import rdkit
            libraries['rdkit'] = rdkit.__version__
        except ImportError:
            libraries['rdkit'] = None
        
        try:
            import vina
            libraries['vina'] = getattr(vina, '__version__', 'unknown')
        except ImportError:
            libraries['vina'] = None
        
        try:
            import openbabel
            libraries['openbabel'] = 'available'
        except ImportError:
            libraries['openbabel'] = None
        
        try:
            import meeko
            libraries['meeko'] = meeko.__version__
        except ImportError:
            libraries['meeko'] = None
        
        print("📚 Scientific Libraries Status:")
        for lib, version in libraries.items():
            status = f"✅ {version}" if version else "❌ Not available"
            print(f"   {lib}: {status}")
        
        return libraries
    
    def run_local_smoke_test(self):
        """Run the smoke test against local server"""
        print("\n🧪 Running Local Smoke Test...")
        
        runner = SmokeTestRunner(base_url=f"http://127.0.0.1:{self.server_port}")
        result = runner.run_full_smoke_test()
        
        print(f"\n📊 Smoke Test Result: {result['status']}")
        
        if result['status'] == 'PASSED':
            print("🎉 Local smoke test PASSED!")
            return True
        else:
            print("💥 Local smoke test FAILED!")
            if 'validation' in result:
                validation = result['validation']
                for error in validation.get('errors', []):
                    print(f"   ❌ {error}")
            if 'error' in result:
                print(f"   💥 Error: {result['error']}")
            return False
    
    def run_validation(self):
        """Run complete local validation"""
        success = True
        
        try:
            print("🔍 Local Smoke Test Validation")
            print("=" * 40)
            
            # Check scientific libraries
            libraries = self.check_scientific_libraries()
            
            # Start Django server
            if not self.start_django_server():
                print("❌ Failed to start Django server")
                return False
            
            # Run smoke test
            success = self.run_local_smoke_test()
            
        except KeyboardInterrupt:
            print("\n⚠️  Test interrupted by user")
            success = False
        except Exception as e:
            print(f"\n💥 Validation failed with exception: {e}")
            success = False
        finally:
            # Always stop server
            self.stop_django_server()
        
        print("\n" + "=" * 40)
        if success:
            print("✅ Local validation PASSED - Ready for CI!")
        else:
            print("❌ Local validation FAILED - Fix issues before CI")
        
        return success


def main():
    """Main entry point"""
    validator = LocalSmokeTestValidator()
    success = validator.run_validation()
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
