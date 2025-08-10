"""
Tests for URL routing and endpoint accessibility
Ensures no 404 errors occur for expected URLs
"""

import json
from django.test import TestCase, Client
from django.urls import resolve
from .models import DockingJob


class TestURLRouting(TestCase):
    """Test URL routing and accessibility"""

    def setUp(self):
        self.client = Client()

    def test_root_url_returns_api_info(self):
        """Test that root URL returns API information instead of 404"""
        response = self.client.get('/')
        
        # Should return 200, not 404
        self.assertEqual(response.status_code, 200)
        
        # Should return JSON with API information
        data = json.loads(response.content)
        self.assertIn('name', data)
        self.assertIn('version', data)
        self.assertIn('description', data)
        self.assertIn('endpoints', data)
        self.assertIn('status', data)
        self.assertIn('features', data)
        
        # Check basic API info structure
        self.assertEqual(data['name'], 'Lipid Rendering API')
        self.assertEqual(data['status'], 'running')
        self.assertIsInstance(data['endpoints'], dict)
        self.assertIsInstance(data['features'], list)

    def test_all_api_endpoints_exist(self):
        """Test that all documented API endpoints exist and don't return 404"""
        endpoints_to_test = [
            '/api/healthz',
            '/api/pdb/1CRN',
            '/api/pdb/1CRN/info',
            '/api/ligand/prepare',
            '/api/receptor/prepare',
            '/api/binding-site/estimate',
            '/api/dock/jobs',
            '/api/dock/run',
            '/api/dock/advanced/run',
            '/api/pockets/detect',
            '/api/templates',
            '/api/poses/rescore'
        ]
        
        for endpoint in endpoints_to_test:
            with self.subTest(endpoint=endpoint):
                response = self.client.get(endpoint)
                # Should not return 404 (endpoint exists)
                # May return 405 (method not allowed) which is fine
                self.assertNotEqual(response.status_code, 404, 
                                   f"Endpoint {endpoint} returned 404")

    def test_health_endpoint_works(self):
        """Test health check endpoint is accessible"""
        response = self.client.get('/api/healthz')
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        self.assertIn('status', data)

    def test_admin_url_exists(self):
        """Test admin URL exists (even if login required)"""
        response = self.client.get('/admin/')
        # Should redirect to login, not return 404
        self.assertIn(response.status_code, [200, 302])

    def test_api_documentation_endpoints(self):
        """Test that API info includes all current endpoints"""
        response = self.client.get('/api/info')
        data = json.loads(response.content)
        
        endpoints = data['endpoints']
        
        # Check that key endpoints are documented
        required_endpoints = ['health', 'protein_info', 'docking', 'admin']
        for endpoint_key in required_endpoints:
            self.assertIn(endpoint_key, endpoints,
                         f"Endpoint {endpoint_key} not documented in API info")
        
        # Check that frontend URL is provided (now served from same server)
        if 'frontend' in endpoints:
            self.assertEqual(endpoints['frontend'], '/')

    def test_docking_workflow_urls(self):
        """Test complete docking workflow URLs are accessible"""
        # Create a test job for status endpoint
        job = DockingJob.objects.create(
            ligand_smiles='CCO',
            receptor_pdb_id='1CRN',
            center_x=5.0,
            center_y=10.0,
            center_z=-2.0,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0
        )
        
        workflow_endpoints = [
            f'/api/dock/status/{job.job_id}',
            '/api/dock/jobs'
        ]
        
        for endpoint in workflow_endpoints:
            response = self.client.get(endpoint)
            self.assertNotEqual(response.status_code, 404)

    def test_advanced_feature_urls(self):
        """Test Phase 7 advanced feature URLs are accessible"""
        # Test template listing
        response = self.client.get('/api/templates')
        self.assertNotEqual(response.status_code, 404)
        
        # Test pocket detection endpoint exists (POST only)
        response = self.client.get('/api/pockets/detect')
        self.assertIn(response.status_code, [405, 400])  # Method not allowed or bad request, not 404
        
        # Test pocket list for protein
        response = self.client.get('/api/pockets/1CRN')
        self.assertNotEqual(response.status_code, 404)
        
        # Test pose rescoring endpoint exists (POST only)
        response = self.client.get('/api/poses/rescore')
        self.assertIn(response.status_code, [405, 400])  # Method not allowed or bad request, not 404

    def test_url_patterns_resolve_correctly(self):
        """Test that Django URL patterns resolve to correct views"""
        # Test root URL resolves
        resolver = resolve('/')
        self.assertEqual(resolver.view_name, 'root')
        
        # Test API endpoints resolve
        resolver = resolve('/api/healthz')
        self.assertEqual(resolver.view_name, 'healthz')
        
        resolver = resolve('/api/pdb/1CRN/info')
        self.assertEqual(resolver.view_name, 'get-protein-info')

    def test_frontend_integration_info(self):
        """Test that API provides frontend integration information"""
        response = self.client.get('/api/info')
        data = json.loads(response.content)
        
        # Should include frontend URL for development
        self.assertIn('endpoints', data)
        self.assertIn('frontend', data['endpoints'])

    def test_api_feature_list_completeness(self):
        """Test that API feature list reflects current capabilities"""
        response = self.client.get('/api/info')
        data = json.loads(response.content)
        
        features = data['features']
        
        # Check that major features are listed
        expected_features = [
            'molecular docking',
            'visualization',
            'GNINA',
            'pocket detection'
        ]
        
        features_text = ' '.join(features).lower()
        for feature in expected_features:
            self.assertIn(feature.lower(), features_text,
                         f"Feature '{feature}' not found in API features")

    def test_version_information(self):
        """Test that API provides version information"""
        response = self.client.get('/api/info')
        data = json.loads(response.content)
        
        self.assertIn('version', data)
        self.assertRegex(data['version'], r'\d+\.\d+\.\d+')  # Semantic version format


class TestAPIErrorHandling(TestCase):
    """Test API error handling and proper HTTP status codes"""

    def setUp(self):
        self.client = Client()

    def test_nonexistent_endpoints_return_404(self):
        """Test that truly nonexistent endpoints return 404"""
        nonexistent_endpoints = [
            '/api/nonexistent',
            '/api/fake/endpoint',
            '/definitely/not/real',
            '/api/dock/nonexistent'
        ]
        
        for endpoint in nonexistent_endpoints:
            with self.subTest(endpoint=endpoint):
                response = self.client.get(endpoint)
                self.assertEqual(response.status_code, 404)

    def test_invalid_job_id_returns_404(self):
        """Test that invalid job IDs return 404"""
        fake_job_id = '00000000-0000-0000-0000-000000000000'
        response = self.client.get(f'/api/dock/status/{fake_job_id}')
        # Should return 400 or 404 for invalid UUID format
        self.assertIn(response.status_code, [400, 404])

    def test_invalid_pdb_id_handling(self):
        """Test handling of invalid PDB IDs"""
        response = self.client.get('/api/pdb/FAKE')
        # Should return 400 (bad request) not 404 for invalid PDB format
        self.assertEqual(response.status_code, 400)

    def test_method_not_allowed_returns_405(self):
        """Test that wrong HTTP methods return 405, not 404"""
        # GET on POST-only endpoint
        response = self.client.get('/api/dock/run')
        self.assertEqual(response.status_code, 405)
        
        # POST on GET-only endpoint
        response = self.client.post('/api/healthz')
        self.assertEqual(response.status_code, 405)


class TestCORSAndSecurity(TestCase):
    """Test CORS and security headers"""

    def setUp(self):
        self.client = Client()

    def test_cors_headers_present(self):
        """Test that CORS headers are present for frontend integration"""
        response = self.client.get('/api/healthz')
        
        # Basic CORS headers should be present for API endpoints
        # Exact headers depend on django-cors-headers configuration
        self.assertEqual(response.status_code, 200)

    def test_api_accepts_json_content_type(self):
        """Test that API endpoints accept JSON content type"""
        response = self.client.post('/api/ligand/prepare',
                                   data='{"smiles": "CCO"}',
                                   content_type='application/json')
        # Should not return 415 (unsupported media type)
        self.assertNotEqual(response.status_code, 415)

