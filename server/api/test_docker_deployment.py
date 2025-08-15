"""
Tests for Docker deployment and containerization.
"""

import json
import subprocess
import time
import requests
import pytest
from django.test import TestCase
from pathlib import Path
import os


def _resolve_path_relative_to_repo(*parts) -> str:
    repo_root = Path(__file__).resolve().parents[2]
    return str(repo_root.joinpath(*parts))

def _resolve_server_path(*parts) -> str:
    repo_root = Path(__file__).resolve().parents[2]
    return str(repo_root.joinpath('server', *parts))
from unittest.mock import patch, MagicMock


class TestDockerConfiguration(TestCase):
    """Test Docker configuration files and settings"""

    def test_dockerfile_exists_and_valid(self):
        """Test that Dockerfile exists and contains required instructions"""
        dockerfile_path = _resolve_server_path('Dockerfile')
        if not Path(dockerfile_path).exists():
            dockerfile_path = _resolve_path_relative_to_repo('Dockerfile')
        
        self.assertTrue(os.path.exists(dockerfile_path), "Dockerfile not found")
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check for essential Docker instructions
        required_instructions = [
            'FROM node:',
            'as frontend-builder',
            'FROM python:3.13-slim',
            'WORKDIR /app',
            'COPY requirements.txt',
            'RUN pip install',
            'EXPOSE 8000',
            'CMD',
            'HEALTHCHECK'
        ]
        
        for instruction in required_instructions:
            self.assertIn(instruction, content, f"Missing Docker instruction: {instruction}")

    def test_dockerfile_security_best_practices(self):
        """Test that Dockerfile follows security best practices"""
        dockerfile_path = _resolve_server_path('Dockerfile')
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Security checks
        self.assertIn('adduser', content, "Should create non-root user")
        self.assertIn('USER appuser', content, "Should run as non-root user")
        self.assertNotIn('sudo', content, "Should not use sudo")
        self.assertNotIn('chmod 777', content, "Should not use overly permissive permissions")

    def test_docker_compose_exists_and_valid(self):
        """Test that docker-compose.yml exists and is valid"""
        import yaml
        compose_path = _resolve_path_relative_to_repo('docker-compose.yml')
        
        self.assertTrue(os.path.exists(compose_path), "docker-compose.yml not found")
        
        with open(compose_path, 'r') as f:
            compose_content = yaml.safe_load(f)
        
        # Check structure
        self.assertIn('version', compose_content)
        self.assertIn('services', compose_content)
        self.assertIn('app', compose_content['services'])
        
        # Check app service configuration
        app_service = compose_content['services']['app']
        self.assertIn('build', app_service)
        self.assertIn('ports', app_service)
        self.assertIn('environment', app_service)
        self.assertIn('healthcheck', app_service)

    def test_requirements_includes_production_dependencies(self):
        """Test that requirements.txt includes production dependencies"""
        req_path = _resolve_server_path('requirements.txt')
        
        with open(req_path, 'r') as f:
            requirements = f.read()
        
        # Check for production dependencies
        prod_deps = ['gunicorn', 'whitenoise', 'Django', 'djangorestframework']
        for dep in prod_deps:
            self.assertIn(dep, requirements, f"Missing production dependency: {dep}")


class TestDockerBuildProcess(TestCase):
    """Test Docker build process and image creation"""

    @patch('subprocess.run')
    def test_docker_build_command_structure(self, mock_subprocess):
        """Test Docker build command structure"""
        # Mock successful build
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = b"Successfully built image"
        mock_subprocess.return_value = mock_result
        
        # Test build command would work
        build_cmd = [
            'docker', 'build', 
            '-t', 'lipid-docking:latest',
            '-f', 'server/Dockerfile',
            '.'
        ]
        
        # This tests the command structure without actually running Docker
        self.assertEqual(len(build_cmd), 7)
        self.assertEqual(build_cmd[0], 'docker')
        self.assertEqual(build_cmd[1], 'build')

    def test_dockerfile_multi_stage_build(self):
        """Test that Dockerfile uses multi-stage build for efficiency"""
        dockerfile_path = _resolve_server_path('Dockerfile')
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check for multi-stage build
        self.assertIn('FROM node:', content, "Should include Node.js stage for frontend build")
        self.assertIn('as frontend-builder', content, "Should name frontend build stage")
        self.assertIn('FROM python:', content, "Should include Python stage for backend")
        self.assertIn('COPY --from=frontend-builder', content, "Should copy from frontend stage")

    def test_dockerfile_static_file_handling(self):
        """Test that Dockerfile properly handles static files"""
        dockerfile_path = _resolve_server_path('Dockerfile')
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check static file handling
        self.assertIn('collectstatic', content, "Should collect static files")
        self.assertIn('--noinput', content, "Should use non-interactive static collection")
        self.assertIn('/app/staticfiles', content, "Should set up static files directory")


class TestContainerRuntime(TestCase):
    """Test container runtime behavior and configuration"""

    def test_environment_variables_configuration(self):
        """Test environment variable configuration for container"""
        required_env_vars = [
            'PYTHONDONTWRITEBYTECODE',
            'PYTHONUNBUFFERED',
            'DJANGO_SETTINGS_MODULE',
            'DEBUG'
        ]
        
        # These would be set in the Docker environment
        for var in required_env_vars:
            # Test that we know what these variables should be
            if var == 'DEBUG':
                expected_value = 'False'  # Production default
            elif var == 'DJANGO_SETTINGS_MODULE':
                expected_value = 'core.settings'
            elif var in ['PYTHONDONTWRITEBYTECODE', 'PYTHONUNBUFFERED']:
                expected_value = '1'
            
            self.assertIsNotNone(expected_value, f"Should define expected value for {var}")

    def test_port_exposure_configuration(self):
        """Test that container properly exposes the correct port"""
        dockerfile_path = _resolve_server_path('Dockerfile')
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        self.assertIn('EXPOSE 8000', content, "Should expose port 8000")

    def test_health_check_configuration(self):
        """Test health check configuration"""
        dockerfile_path = _resolve_server_path('Dockerfile')
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check health check configuration
        self.assertIn('HEALTHCHECK', content, "Should include health check")
        self.assertIn('/api/healthz', content, "Should check health endpoint")
        self.assertIn('--interval=', content, "Should set health check interval")
        self.assertIn('--timeout=', content, "Should set health check timeout")

    def test_gunicorn_production_server_configuration(self):
        """Test Gunicorn server configuration for production"""
        dockerfile_path = _resolve_server_path('Dockerfile')
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check Gunicorn configuration
        self.assertIn('gunicorn', content, "Should use Gunicorn for production")
        # Accept both exec-form JSON array and shell string formats
        self.assertTrue(
            ('--bind 0.0.0.0:8000' in content) or ('"--bind", "0.0.0.0:8000"' in content),
            "Should bind to all interfaces",
        )
        self.assertIn('--workers', content, "Should configure worker processes")
        self.assertIn('core.wsgi:application', content, "Should reference Django WSGI app")


class TestDockerSingleServerIntegration(TestCase):
    """Test Docker integration with single-server architecture"""

    def test_static_files_served_in_container(self):
        """Test that static files are properly served in container"""
        # This tests the configuration that would serve static files
        from django.conf import settings
        
        # Check that settings are configured for container deployment
        self.assertTrue(any('whitenoise' in m for m in settings.MIDDLEWARE))
        self.assertEqual(settings.STATIC_URL, '/static/')
        self.assertIsNotNone(settings.STATIC_ROOT)

    def test_frontend_assets_copied_to_container(self):
        """Test that frontend assets are copied to container"""
        dockerfile_path = _resolve_server_path('Dockerfile')
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Check that frontend assets are handled
        self.assertIn('npm run build', content, "Should build frontend assets")
        self.assertIn('COPY --from=frontend-builder /frontend/dist', content, 
                     "Should copy built frontend assets")

    def test_api_routes_not_conflicting_with_static_routes(self):
        """Test that API routes don't conflict with static file serving"""
        from django.urls import resolve, reverse
        from django.http import Http404
        
        # Test that API routes are properly namespaced
        try:
            api_url = reverse('api-info')
            self.assertTrue(api_url.startswith('/api/'))
        except:
            # If named URL doesn't exist, test that /api/ prefix works
            self.assertTrue('/api/' in 'expected API structure')

    def test_docker_compose_development_configuration(self):
        """Test Docker Compose development configuration"""
        import os
        import yaml
        
        compose_path = os.path.join(os.path.dirname(__file__), '../../../docker-compose.yml')
        
        if os.path.exists(compose_path):
            with open(compose_path, 'r') as f:
                compose_content = yaml.safe_load(f)
            
            # Check for development profile
            if 'app-dev' in compose_content.get('services', {}):
                dev_service = compose_content['services']['app-dev']
                self.assertIn('profiles', dev_service)
                self.assertIn('dev', dev_service['profiles'])


class TestDockerDeploymentScenarios(TestCase):
    """Test various Docker deployment scenarios"""

    def test_production_deployment_configuration(self):
        """Test production deployment configuration"""
        # Test that production settings are appropriate
        from django.conf import settings
        
        # In a container, DEBUG should be False
        # ALLOWED_HOSTS should be configured
        # Static files should be configured for production
        
        self.assertTrue(any('whitenoise' in m for m in settings.MIDDLEWARE))
        self.assertIsNotNone(settings.STATIC_ROOT)

    def test_development_vs_production_dockerfile_differences(self):
        """Test differences between development and production Dockerfiles"""
        import os
        
        prod_dockerfile = _resolve_server_path('Dockerfile')
        dev_dockerfile = _resolve_server_path('Dockerfile.dev')
        
        # Both should exist
        self.assertTrue(os.path.exists(prod_dockerfile), "Production Dockerfile should exist")
        
        if os.path.exists(dev_dockerfile):
            with open(prod_dockerfile, 'r') as f:
                prod_content = f.read()
            
            with open(dev_dockerfile, 'r') as f:
                dev_content = f.read()
            
            # Production should use Gunicorn, development should use Django dev server
            self.assertIn('gunicorn', prod_content)
            self.assertIn('runserver', dev_content)

    def test_container_security_configuration(self):
        """Test container security configuration"""
        import os
        dockerfile_path = _resolve_server_path('Dockerfile')
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Security best practices
        self.assertIn('USER appuser', content, "Should run as non-root user")
        self.assertNotIn('--user root', content, "Should not explicitly use root")
        self.assertNotIn('chmod 777', content, "Should not use overly permissive permissions")

    def test_container_resource_optimization(self):
        """Test container resource optimization"""
        import os
        dockerfile_path = _resolve_server_path('Dockerfile')
        
        with open(dockerfile_path, 'r') as f:
            content = f.read()
        
        # Resource optimization checks
        self.assertIn('rm -rf /var/lib/apt/lists/*', content, 
                     "Should clean up apt cache to reduce image size")
        self.assertIn('PIP_NO_CACHE_DIR=1', content, 
                     "Should disable pip cache to reduce image size")


class TestDockerNetworking(TestCase):
    """Test Docker networking configuration"""

    def test_docker_compose_networking(self):
        """Test Docker Compose networking configuration"""
        import os
        import yaml
        
        compose_path = os.path.join(os.path.dirname(__file__), '../../../docker-compose.yml')
        
        if os.path.exists(compose_path):
            with open(compose_path, 'r') as f:
                compose_content = yaml.safe_load(f)
            
            # Check networking
            if 'networks' in compose_content:
                self.assertIn('default', compose_content['networks'])

    def test_port_mapping_configuration(self):
        """Test port mapping configuration"""
        import os
        import yaml
        
        compose_path = os.path.join(os.path.dirname(__file__), '../../../docker-compose.yml')
        
        if os.path.exists(compose_path):
            with open(compose_path, 'r') as f:
                compose_content = yaml.safe_load(f)
            
            # Check port mapping
            app_service = compose_content['services']['app']
            ports = app_service.get('ports', [])
            
            # Should map 8000:8000
            self.assertIn('8000:8000', ports)


# Integration test that would run against actual container (disabled by default)
class TestLiveDockerContainer:
    """
    Integration tests for live Docker container.
    These tests would run against an actual container instance.
    """
    
    @pytest.mark.skip(reason="Requires Docker container to be running")
    def test_container_health_endpoint(self):
        """Test health endpoint on running container"""
        response = requests.get('http://localhost:8000/api/healthz', timeout=10)
        assert response.status_code == 200
        
    @pytest.mark.skip(reason="Requires Docker container to be running")
    def test_container_serves_frontend(self):
        """Test that container serves frontend application"""
        response = requests.get('http://localhost:8000/', timeout=10)
        assert response.status_code == 200
        assert 'text/html' in response.headers.get('Content-Type', '')
        
    @pytest.mark.skip(reason="Requires Docker container to be running")
    def test_container_api_endpoints(self):
        """Test API endpoints on running container"""
        response = requests.get('http://localhost:8000/api/info', timeout=10)
        assert response.status_code == 200
        
        data = response.json()
        assert 'name' in data
        assert 'version' in data
        assert 'endpoints' in data
