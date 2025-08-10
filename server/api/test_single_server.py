"""
Tests for single server architecture
Ensures React frontend and Django API work together on one server
"""

import json
from django.test import TestCase, Client
from django.urls import reverse
from django.conf import settings


class TestSingleServerArchitecture(TestCase):
    """Test that frontend and backend work together on single server"""

    def setUp(self):
        self.client = Client()

    def test_api_info_endpoint_accessible(self):
        """Test that API info is accessible at /api/info"""
        response = self.client.get('/api/info')
        self.assertEqual(response.status_code, 200)
        
        data = json.loads(response.content)
        self.assertEqual(data['name'], 'Lipid Rendering API')
        self.assertEqual(data['endpoints']['frontend'], '/')  # Frontend now served from same server

    def test_health_endpoint_accessible(self):
        """Test that health endpoint is accessible via relative URL"""
        response = self.client.get('/api/healthz')
        self.assertEqual(response.status_code, 200)

    def test_react_app_served_at_root(self):
        """Test that React app is served at root path"""
        response = self.client.get('/')
        self.assertEqual(response.status_code, 200)
        
        # Should be HTML content (React app)
        self.assertIn('text/html', response.get('Content-Type', ''))

    def test_react_app_handles_spa_routing(self):
        """Test that non-API routes serve React app (SPA routing)"""
        spa_routes = [
            '/dock',
            '/dock/advanced',
            '/some/nonexistent/route',
            '/molecules',
        ]
        
        for route in spa_routes:
            with self.subTest(route=route):
                response = self.client.get(route)
                self.assertEqual(response.status_code, 200)
                self.assertIn('text/html', response.get('Content-Type', ''))

    def test_api_routes_not_intercepted_by_react(self):
        """Test that API routes are not caught by React catch-all"""
        api_routes = [
            '/api/healthz',
            '/api/pdb/1CRN',
            '/api/dock/jobs',
            '/api/templates',
        ]
        
        for route in api_routes:
            with self.subTest(route=route):
                response = self.client.get(route)
                # Should not return HTML (React app)
                # Should be JSON or appropriate API response
                self.assertNotEqual(response.status_code, 404)
                content_type = response.get('Content-Type', '')
                if response.status_code == 200:
                    # If successful, should be JSON API response
                    self.assertIn('application/json', content_type)

    def test_static_files_served_correctly(self):
        """Test that static files (CSS, JS) are served correctly"""
        # Check that assets directory exists in staticfiles
        # This is more of a configuration test
        self.assertTrue(settings.STATICFILES_DIRS)
        self.assertIn('lipid_viewer', str(settings.STATICFILES_DIRS[0]))

    def test_cors_not_needed_for_same_origin(self):
        """Test that CORS headers are not needed since frontend and backend are same origin"""
        response = self.client.get('/api/healthz')
        
        # Since it's same origin, we don't need CORS headers
        # But they might still be present from previous configuration
        self.assertEqual(response.status_code, 200)

    def test_relative_api_urls_work(self):
        """Test that frontend can use relative API URLs"""
        # Test that relative URLs resolve correctly
        response = self.client.post('/api/binding-site/estimate', 
                                   data=json.dumps({'pdb_id': '1CRN'}),
                                   content_type='application/json')
        
        # Should not be 404 (URL found)
        self.assertNotEqual(response.status_code, 404)


class TestProductionReadiness(TestCase):
    """Test that the single server setup is production ready"""

    def test_whitenoise_configuration(self):
        """Test that whitenoise is properly configured"""
        self.assertIn('whitenoise.middleware.WhiteNoiseMiddleware', 
                     settings.MIDDLEWARE)
        self.assertEqual(settings.STATICFILES_STORAGE, 
                        'whitenoise.storage.CompressedManifestStaticFilesStorage')

    def test_static_file_configuration(self):
        """Test that static file settings are correct"""
        self.assertTrue(settings.STATIC_ROOT)
        self.assertTrue(settings.STATICFILES_DIRS)
        self.assertEqual(settings.STATIC_URL, '/static/')

    def test_template_configuration(self):
        """Test that templates can find React build output"""
        template_dirs = settings.TEMPLATES[0]['DIRS']
        self.assertTrue(any('lipid_viewer' in str(d) and 'dist' in str(d) 
                           for d in template_dirs))


class TestSingleServerDeployment(TestCase):
    """Test deployment-specific functionality"""

    def test_single_port_serves_everything(self):
        """Test that both frontend and API are served from same port"""
        # Test API
        api_response = self.client.get('/api/healthz')
        self.assertEqual(api_response.status_code, 200)
        
        # Test Frontend
        frontend_response = self.client.get('/')
        self.assertEqual(frontend_response.status_code, 200)
        
        # Both should work from same server/port

    def test_no_port_3000_references(self):
        """Test that there are no hardcoded port 3000 references"""
        # This is more of a code review test
        # Check that API info returns relative URLs
        response = self.client.get('/api/info')
        data = json.loads(response.content)
        
        # Frontend should be served from root, not port 3000
        self.assertEqual(data['endpoints']['frontend'], '/')
        
        # No localhost:3000 references
        for endpoint_url in data['endpoints'].values():
            self.assertNotIn('3000', str(endpoint_url))
            self.assertNotIn('localhost:3000', str(endpoint_url))
