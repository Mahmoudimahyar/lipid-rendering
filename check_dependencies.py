#!/usr/bin/env python3

import sys
import pkg_resources

print("=== DEPENDENCY CHECK ===")
print(f"Python version: {sys.version}")
print(f"Python path: {sys.path}")

# Check for scientific libraries
packages_to_check = [
    'vina', 'rdkit', 'numpy', 'scipy', 'pandas', 
    'django', 'djangorestframework', 'openbabel'
]

print("\n=== INSTALLED PACKAGES ===")
installed_packages = {pkg.project_name.lower(): pkg.version for pkg in pkg_resources.working_set}

for package in packages_to_check:
    if package.lower() in installed_packages:
        print(f"✅ {package}: {installed_packages[package.lower()]}")
    else:
        print(f"❌ {package}: NOT INSTALLED")

print("\n=== IMPORT TEST ===")
# Test imports
try:
    import vina
    print("✅ Vina import: SUCCESS")
    print(f"   Version: {getattr(vina, '__version__', 'unknown')}")
except ImportError as e:
    print(f"❌ Vina import: FAILED - {e}")

try:
    from rdkit import Chem
    print("✅ RDKit import: SUCCESS")
except ImportError as e:
    print(f"❌ RDKit import: FAILED - {e}")

try:
    import numpy
    print("✅ NumPy import: SUCCESS")
except ImportError as e:
    print(f"❌ NumPy import: FAILED - {e}")

print("\n=== DJANGO SETTINGS ===")
try:
    import django
    from django.conf import settings
    print(f"Django settings module: {settings.SETTINGS_MODULE}")
    print(f"DOCKING_ALLOW_MOCK: {getattr(settings, 'DOCKING_ALLOW_MOCK', 'NOT_SET')}")
except Exception as e:
    print(f"Django check failed: {e}")

print("\n=== AUTODOCK VINA AVAILABILITY ===")
try:
    from api.real_docking_engine import RealDockingEngine
    print(f"RealDockingEngine.is_available(): {RealDockingEngine.is_available()}")
except Exception as e:
    print(f"RealDockingEngine check failed: {e}")

try:
    from api.docking_utils import DockingEngine
    print(f"DockingEngine.is_real_docking_available(): {DockingEngine.is_real_docking_available()}")
except Exception as e:
    print(f"DockingEngine check failed: {e}")
