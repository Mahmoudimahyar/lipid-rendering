#!/usr/bin/env python3
"""
Validate Docker setup for single-server deployment.
This script checks that all Docker-related files are properly configured.
"""

import os
import sys
import yaml
import json
from pathlib import Path


def validate_docker_files():
    """Validate that all required Docker files exist and are properly configured"""
    project_root = Path(__file__).parent.parent
    server_dir = project_root / "server"
    
    results = []
    
    # Check Dockerfile
    dockerfile = server_dir / "Dockerfile"
    if dockerfile.exists():
        content = dockerfile.read_text()
        checks = [
            ("Multi-stage build", "FROM node:" in content and "FROM python:" in content),
            ("Gunicorn server", "gunicorn" in content),
            ("Health check", "HEALTHCHECK" in content),
            ("Non-root user", "USER appuser" in content),
            ("Static files", "collectstatic" in content),
            ("Port exposure", "EXPOSE 8000" in content)
        ]
        
        for check_name, passed in checks:
            results.append(("Dockerfile", check_name, passed))
    else:
        results.append(("Dockerfile", "File exists", False))
    
    # Check Dockerfile.dev
    dockerfile_dev = server_dir / "Dockerfile.dev"
    if dockerfile_dev.exists():
        content = dockerfile_dev.read_text()
        results.append(("Dockerfile.dev", "Development config", "runserver" in content))
    else:
        results.append(("Dockerfile.dev", "File exists", False))
    
    # Check docker-compose.yml
    compose_file = project_root / "docker-compose.yml"
    if compose_file.exists():
        try:
            with open(compose_file, 'r') as f:
                compose_data = yaml.safe_load(f)
            
            checks = [
                ("Has services", "services" in compose_data),
                ("Has app service", "app" in compose_data.get("services", {})),
                ("Has dev profile", any("dev" in service.get("profiles", []) 
                                      for service in compose_data.get("services", {}).values())),
                ("Port mapping", any("8000:8000" in service.get("ports", []) 
                                   for service in compose_data.get("services", {}).values()))
            ]
            
            for check_name, passed in checks:
                results.append(("docker-compose.yml", check_name, passed))
                
        except yaml.YAMLError as e:
            results.append(("docker-compose.yml", "Valid YAML", False))
    else:
        results.append(("docker-compose.yml", "File exists", False))
    
    # Check requirements.txt
    requirements_file = server_dir / "requirements.txt"
    if requirements_file.exists():
        content = requirements_file.read_text()
        checks = [
            ("Has Django", "Django==" in content),
            ("Has Gunicorn", "gunicorn==" in content),
            ("Has Whitenoise", "whitenoise==" in content),
            ("Has test dependencies", "pytest==" in content)
        ]
        
        for check_name, passed in checks:
            results.append(("requirements.txt", check_name, passed))
    else:
        results.append(("requirements.txt", "File exists", False))
    
    # Check .dockerignore
    dockerignore = server_dir / ".dockerignore"
    if dockerignore.exists():
        content = dockerignore.read_text()
        checks = [
            ("Ignores Python cache", "__pycache__" in content),
            ("Ignores virtual env", "venv/" in content),
            ("Ignores node modules", "node_modules/" in content)
        ]
        
        for check_name, passed in checks:
            results.append((".dockerignore", check_name, passed))
    else:
        results.append((".dockerignore", "File exists", False))
    
    return results


def validate_django_settings():
    """Validate Django settings for Docker deployment"""
    sys.path.append(str(Path(__file__).parent.parent / "server"))
    
    try:
        os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'core.settings')
        import django
        from django.conf import settings
        django.setup()
        
        results = []
        
        # Check static files configuration
        results.append(("Django Settings", "Static URL configured", 
                       hasattr(settings, 'STATIC_URL') and settings.STATIC_URL))
        results.append(("Django Settings", "Static root configured", 
                       hasattr(settings, 'STATIC_ROOT') and settings.STATIC_ROOT))
        results.append(("Django Settings", "Whitenoise in middleware", 
                       'whitenoise' in str(settings.MIDDLEWARE).lower()))
        
        # Check template configuration
        results.append(("Django Settings", "Templates configured", 
                       len(settings.TEMPLATES) > 0))
        
        # Check for React integration
        if hasattr(settings, 'STATICFILES_DIRS'):
            results.append(("Django Settings", "Frontend integration", 
                           any('dist' in str(d) for d in settings.STATICFILES_DIRS)))
        else:
            results.append(("Django Settings", "Frontend integration", False))
        
        return results
        
    except Exception as e:
        return [("Django Settings", f"Configuration error: {e}", False)]


def check_frontend_build():
    """Check that frontend is built and ready"""
    project_root = Path(__file__).parent.parent
    frontend_dist = project_root / "lipid_viewer" / "dist"
    
    results = []
    
    if frontend_dist.exists():
        # Check for essential files
        index_html = frontend_dist / "index.html"
        assets_dir = frontend_dist / "assets"
        
        results.append(("Frontend Build", "index.html exists", index_html.exists()))
        results.append(("Frontend Build", "Assets directory exists", assets_dir.exists()))
        
        if assets_dir.exists():
            js_files = list(assets_dir.glob("*.js"))
            results.append(("Frontend Build", "JavaScript assets exist", len(js_files) > 0))
    else:
        results.append(("Frontend Build", "Build directory exists", False))
    
    return results


def print_results(results):
    """Print validation results in a formatted way"""
    print("🐳 Docker Setup Validation Results")
    print("=" * 50)
    
    categories = {}
    for category, check, passed in results:
        if category not in categories:
            categories[category] = []
        categories[category].append((check, passed))
    
    total_checks = len(results)
    passed_checks = sum(1 for _, _, passed in results if passed)
    
    for category, checks in categories.items():
        print(f"\n📁 {category}:")
        for check, passed in checks:
            status = "✅" if passed else "❌"
            print(f"  {status} {check}")
    
    print(f"\n📊 Summary: {passed_checks}/{total_checks} checks passed")
    
    if passed_checks == total_checks:
        print("🎉 All checks passed! Docker setup is ready.")
        return True
    else:
        print("⚠️  Some checks failed. Please review the configuration.")
        return False


def main():
    """Main validation function"""
    print("Starting Docker setup validation...")
    
    all_results = []
    
    # Validate Docker files
    all_results.extend(validate_docker_files())
    
    # Validate Django settings
    all_results.extend(validate_django_settings())
    
    # Check frontend build
    all_results.extend(check_frontend_build())
    
    # Print results
    success = print_results(all_results)
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
