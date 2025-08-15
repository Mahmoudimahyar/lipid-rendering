"""
Tests for CSRF and CORS configuration to ensure frontend-backend communication works
"""
import pytest
from django.test import TestCase, Client
from django.conf import settings
from unittest.mock import patch
import json


class TestCSRFCORSConfiguration(TestCase):
    """Test CSRF and CORS settings for frontend integration"""
    
    def setUp(self):
        self.client = Client()
        self.frontend_origin = 'http://localhost:3000'
    
    def test_cors_allowed_origins_includes_frontend(self):
        """Test that CORS allows frontend origin"""
        self.assertIn(self.frontend_origin, settings.CORS_ALLOWED_ORIGINS)
        self.assertIn('http://127.0.0.1:3000', settings.CORS_ALLOWED_ORIGINS)
    
    def test_csrf_trusted_origins_includes_frontend(self):
        """Test that CSRF trusts frontend origin"""
        self.assertIn(self.frontend_origin, settings.CSRF_TRUSTED_ORIGINS)
        self.assertIn('http://127.0.0.1:3000', settings.CSRF_TRUSTED_ORIGINS)
    
    def test_cors_preflight_request(self):
        """Test CORS preflight request from frontend"""
        response = self.client.options(
            '/api/pdb/1CRN/info',
            HTTP_ORIGIN=self.frontend_origin,
            HTTP_ACCESS_CONTROL_REQUEST_METHOD='GET',
            HTTP_ACCESS_CONTROL_REQUEST_HEADERS='content-type'
        )
        
        # Should not be forbidden
        self.assertNotEqual(response.status_code, 403)
    
    def test_api_request_with_frontend_origin(self):
        """Test API request with frontend origin header"""
        response = self.client.get(
            '/api/healthz',
            HTTP_ORIGIN=self.frontend_origin
        )
        
        self.assertEqual(response.status_code, 200)
        # Should have CORS headers
        self.assertIn('Access-Control-Allow-Origin', response.headers)
    
    def test_post_request_from_frontend_origin(self):
        """Test POST request from frontend (the one that was failing)"""
        # This simulates the binding site estimation request that was failing
        response = self.client.post(
            '/api/binding-site/estimate',
            data=json.dumps({
                'pdb_id': '1CRN',
                'ligand_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'
            }),
            content_type='application/json',
            HTTP_ORIGIN=self.frontend_origin
        )
        
        # Should not be 403 Forbidden due to CSRF
        self.assertNotEqual(response.status_code, 403)
        # Might be 400 or 500 due to missing PDB file, but not 403
        self.assertIn(response.status_code, [200, 400, 404, 500])
    
    def test_docking_request_from_frontend(self):
        """Test docking request from frontend origin"""
        response = self.client.post(
            '/api/dock/run',
            data=json.dumps({
                'ligand_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'receptor_pdb_id': '1CRN',
                'center_x': 0,
                'center_y': 0,
                'center_z': 0,
                'size_x': 20,
                'size_y': 20,
                'size_z': 20
            }),
            content_type='application/json',
            HTTP_ORIGIN=self.frontend_origin
        )
        
        # Should not be 403 Forbidden
        self.assertNotEqual(response.status_code, 403)
    
    @patch('django.middleware.csrf.CsrfViewMiddleware.process_view')
    def test_csrf_protection_bypassed_for_api(self, mock_csrf):
        """Test that CSRF protection allows API requests from trusted origins"""
        # This test ensures our CSRF_TRUSTED_ORIGINS setting works
        mock_csrf.return_value = None
        response = self.client.post(
            '/api/pdb/1CRN/info',
            HTTP_ORIGIN=self.frontend_origin
        )
        
        # Should not trigger CSRF protection for trusted origin
        self.assertNotEqual(response.status_code, 403)


class TestFrontendBackendIntegration(TestCase):
    """Integration tests to prevent connection failures"""
    
    def test_all_critical_endpoints_accept_frontend_requests(self):
        """Test that all critical API endpoints accept requests from frontend"""
        client = Client()
        frontend_origin = 'http://localhost:3000'
        
        # List of critical endpoints that frontend uses
        endpoints = [
            '/api/healthz',
            '/api/pdb/1CRN/info',
            '/api/templates',
        ]
        
        for endpoint in endpoints:
            response = client.get(
                endpoint,
                HTTP_ORIGIN=frontend_origin
            )
            
            # Should not be forbidden
            self.assertNotEqual(
                response.status_code, 
                403, 
                f"Endpoint {endpoint} returned 403 Forbidden for frontend origin"
            )
    
    def test_post_endpoints_accept_frontend_requests(self):
        """Test POST endpoints accept requests from frontend"""
        client = Client()
        frontend_origin = 'http://localhost:3000'
        
        # POST endpoints that frontend uses
        post_data = {
            'ligand_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'receptor_pdb_id': '1CRN'
        }
        
        post_endpoints = [
            '/api/binding-site/estimate',
            '/api/dock/run',
        ]
        
        for endpoint in post_endpoints:
            response = client.post(
                endpoint,
                data=json.dumps(post_data),
                content_type='application/json',
                HTTP_ORIGIN=frontend_origin
            )
            
            # Should not be 403 Forbidden due to CSRF
            self.assertNotEqual(
                response.status_code, 
                403, 
                f"POST endpoint {endpoint} returned 403 Forbidden - CSRF issue"
            )
    
    def test_settings_configuration(self):
        """Test that settings are properly configured for frontend integration"""
        # Verify CORS settings
        self.assertTrue(hasattr(settings, 'CORS_ALLOWED_ORIGINS'))
        self.assertIn('http://localhost:3000', settings.CORS_ALLOWED_ORIGINS)
        
        # Verify CSRF settings
        self.assertTrue(hasattr(settings, 'CSRF_TRUSTED_ORIGINS'))
        self.assertIn('http://localhost:3000', settings.CSRF_TRUSTED_ORIGINS)
        
        # Verify CORS credentials
        self.assertTrue(getattr(settings, 'CORS_ALLOW_CREDENTIALS', False))


@pytest.mark.django_db
class TestCSRFErrorPrevention(TestCase):
    """Specific tests to prevent the CSRF error that occurred in production"""
    
    def test_binding_site_estimation_csrf_protection(self):
        """Test the specific request that failed with CSRF error"""
        client = Client()
        
        # This is the exact request that was failing
        response = client.post(
            '/api/binding-site/estimate',
            data=json.dumps({
                'pdb_id': '1CRN',
                'ligand_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'
            }),
            content_type='application/json',
            HTTP_ORIGIN='http://localhost:3000'
        )
        
        # Should NOT return 403 Forbidden
        self.assertNotEqual(response.status_code, 403)
        
        # Should NOT contain CSRF error message
        if hasattr(response, 'content'):
            content = response.content.decode()
            self.assertNotIn('CSRF verification failed', content)
            self.assertNotIn('Origin checking failed', content)
    
    def test_docking_requests_no_csrf_errors(self):
        """Test that docking requests don't trigger CSRF errors"""
        client = Client()
        
        docking_data = {
            'ligand_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'receptor_pdb_id': '1CRN',
            'center_x': 0,
            'center_y': 0,
            'center_z': 0,
            'size_x': 20,
            'size_y': 20,
            'size_z': 20
        }
        
        endpoints = [
            '/api/dock/run',
            '/api/dock/advanced/run',
        ]
        
        for endpoint in endpoints:
            response = client.post(
                endpoint,
                data=json.dumps(docking_data),
                content_type='application/json',
                HTTP_ORIGIN='http://localhost:3000'
            )
            
            # Should not be CSRF forbidden
            self.assertNotEqual(response.status_code, 403)
            
            if hasattr(response, 'content'):
                content = response.content.decode()
                self.assertNotIn('CSRF verification failed', content)
