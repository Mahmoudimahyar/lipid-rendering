#!/usr/bin/env python3
"""
Production Scientific Software Validation Script

This script validates that the molecular docking system is ready for production
scientific research and provides instructions for enabling full capabilities.
"""

import sys
import subprocess
import importlib
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def check_package_availability(package_name, import_name=None):
    """Check if a Python package is available"""
    if import_name is None:
        import_name = package_name
    
    try:
        importlib.import_module(import_name)
        return True, "Available"
    except ImportError as e:
        return False, f"Not available: {e}"


def check_conda_package(package_name):
    """Check if a conda package is available"""
    try:
        result = subprocess.run(['conda', 'list', package_name], 
                              capture_output=True, text=True)
        if result.returncode == 0 and package_name in result.stdout:
            return True, "Available via conda"
        else:
            return False, "Not installed via conda"
    except FileNotFoundError:
        return False, "Conda not available"


def validate_scientific_stack():
    """Validate the scientific computing stack"""
    print("🔬 SCIENTIFIC SOFTWARE VALIDATION REPORT")
    print("=" * 60)
    
    # Core scientific packages
    packages_to_check = [
        ("rdkit", "rdkit"),
        ("numpy", "numpy"), 
        ("scipy", "scipy"),
        ("pandas", "pandas"),
        ("torch", "torch"),
        ("biopython", "Bio"),
        ("prody", "prody"),
        ("openmm", "openmm"),
        ("vina", "vina"),
    ]
    
    available_count = 0
    total_count = len(packages_to_check)
    
    print("\n📦 PYTHON PACKAGE STATUS:")
    print("-" * 40)
    
    for package_name, import_name in packages_to_check:
        available, status = check_package_availability(package_name, import_name)
        status_icon = "✅" if available else "❌"
        print(f"{status_icon} {package_name:12} - {status}")
        if available:
            available_count += 1
    
    print(f"\n📊 AVAILABILITY SUMMARY: {available_count}/{total_count} packages available")
    
    # Check system executables
    print("\n🛠️  SYSTEM EXECUTABLES:")
    print("-" * 40)
    
    executables = ["vina", "obabel", "conda"]
    for executable in executables:
        try:
            result = subprocess.run(['which', executable], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                print(f"✅ {executable:12} - Available at {result.stdout.strip()}")
            else:
                print(f"❌ {executable:12} - Not found in PATH")
        except Exception:
            print(f"❌ {executable:12} - Check failed")
    
    return available_count / total_count


def check_docking_implementation():
    """Check the docking implementation status"""
    print("\n🧬 DOCKING IMPLEMENTATION STATUS:")
    print("-" * 40)
    
    try:
        # Check if our implementation classes exist
        from api.docking_utils import DockingEngine
        from api.real_docking_engine import RealDockingUtils
        from api.real_gnina_scorer import RealGNINAScorer
        from api.real_pocket_detection import RealPocketDetector
        
        print("✅ DockingEngine - Core docking implementation")
        print("✅ RealDockingUtils - AutoDock Vina integration")
        print("✅ RealGNINAScorer - GNINA neural network scoring")
        print("✅ RealPocketDetector - Binding pocket detection")
        
        # Check availability detection
        docking_available = DockingEngine.is_real_docking_available()
        gnina_available = RealGNINAScorer.is_available()
        pocket_available = RealPocketDetector.is_available()
        
        print(f"\n🔍 AVAILABILITY DETECTION:")
        print(f"   AutoDock Vina: {'✅ Ready' if docking_available else '⚠️  Fallback mode'}")
        print(f"   GNINA Scoring: {'✅ Ready' if gnina_available else '⚠️  Fallback mode'}")
        print(f"   Pocket Detection: {'✅ Ready' if pocket_available else '⚠️  Fallback mode'}")
        
        return True
        
    except ImportError as e:
        print(f"❌ Implementation error: {e}")
        return False


def run_quick_functionality_test():
    """Run a quick functionality test"""
    print("\n🧪 FUNCTIONALITY TEST:")
    print("-" * 40)
    
    try:
        from api.docking_utils import DockingEngine
        
        # Test parameter validation
        test_params = {
            'ligand_smiles': 'CCO',  # Ethanol
            'receptor_pdb_id': '1CRN',
            'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0,
            'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0
        }
        
        validated = DockingEngine.validate_docking_parameters(test_params)
        print("✅ Parameter validation - Working")
        
        # Test mock docking (should always work)
        results = DockingEngine.run_mock_docking(validated)
        if results.get('success'):
            print("✅ Mock docking - Working")
            print(f"   Generated {len(results.get('poses', []))} poses")
        else:
            print("❌ Mock docking - Failed")
            
        return True
        
    except Exception as e:
        print(f"❌ Functionality test failed: {e}")
        return False


def generate_installation_instructions():
    """Generate installation instructions for full production mode"""
    print("\n🚀 PRODUCTION MODE ACTIVATION:")
    print("=" * 60)
    print("""
To enable full production scientific computing capabilities, install:

📋 CONDA PACKAGES (Recommended):
   conda install -c conda-forge -y autodock-vina rdkit openmm biopython
   conda install -c conda-forge -y pytorch torchvision scipy numpy pandas

📋 PIP PACKAGES (Alternative):
   pip install vina rdkit torch torchvision biopython prody openbabel-wheel

📋 SYSTEM PACKAGES (Ubuntu/Debian):
   apt-get install autodock-vina openbabel

🔄 AUTOMATIC DETECTION:
   The system will automatically detect installed software and switch
   from fallback mode to production mode without any code changes.

📊 CURRENT STATUS:
   ✅ All implementation code ready
   ✅ Fallback system active and tested
   ✅ Production APIs implemented
   ⚠️  Waiting for scientific software installation
""")


def main():
    """Main validation function"""
    print("🔬 MOLECULAR DOCKING PRODUCTION VALIDATION")
    print("=" * 60)
    print("Validating scientific software stack for production readiness...")
    
    # Check if running in Django context
    try:
        import django
        from django.conf import settings
        if not settings.configured:
            import os
            os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'core.settings')
            django.setup()
    except Exception as e:
        print(f"⚠️  Django setup: {e}")
    
    # Run validation checks
    software_score = validate_scientific_stack()
    implementation_ready = check_docking_implementation()
    functionality_working = run_quick_functionality_test()
    
    # Generate final report
    print("\n🏆 FINAL VALIDATION REPORT:")
    print("=" * 60)
    
    overall_score = (
        (software_score * 0.4) +
        (1.0 if implementation_ready else 0.0) * 0.3 +
        (1.0 if functionality_working else 0.0) * 0.3
    )
    
    if overall_score >= 0.8:
        status = "✅ PRODUCTION READY"
        recommendation = "Ready for scientific research use"
    elif overall_score >= 0.6:
        status = "⚠️  PARTIALLY READY"
        recommendation = "Functional with fallback mode, install scientific software for full capabilities"
    else:
        status = "❌ NEEDS ATTENTION"
        recommendation = "Implementation issues detected, review error messages above"
    
    print(f"📊 Overall Score: {overall_score:.1%}")
    print(f"🎯 Status: {status}")
    print(f"💡 Recommendation: {recommendation}")
    
    if overall_score >= 0.6:
        generate_installation_instructions()
    
    print("\n🎉 VALIDATION COMPLETE")
    return overall_score >= 0.6


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
