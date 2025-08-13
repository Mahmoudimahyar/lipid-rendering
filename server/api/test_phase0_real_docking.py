"""
Tests for Phase 0: Enforce REAL_DOCKING & Deprecate Mock by Default

This module tests the DOCKING_ALLOW_MOCK configuration and capabilities endpoint.
"""

import pytest
import json
import sys
from django.test import TestCase, Client
from django.urls import reverse
from django.conf import settings
from unittest.mock import patch, MagicMock

from api.models import DockingJob
from api.docking_utils import DockingEngine


class TestDockingCapabilities(TestCase):
    """Test the /api/dock/capabilities endpoint"""
    
    def setUp(self):
        self.client = Client()
        self.capabilities_url = reverse('dock-capabilities')
    
    def test_capabilities_endpoint_exists(self):
        """Test that the capabilities endpoint is accessible"""
        response = self.client.get(self.capabilities_url)
        self.assertEqual(response.status_code, 200)
        
        data = json.loads(response.content)
        self.assertIn('vina_available', data)
        self.assertIn('vina_version', data)
        self.assertIn('engine_default', data)
        self.assertIn('mock_allowed', data)
        self.assertIn('advanced_features', data)
        self.assertIn('status', data)
    
    def test_capabilities_with_vina_unavailable_mock_disabled(self):
        """Test capabilities when Vina unavailable and mock disabled"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': False}):
                response = self.client.get(self.capabilities_url)
                self.assertEqual(response.status_code, 200)
                
                data = json.loads(response.content)
                self.assertFalse(data['vina_available'])
                self.assertIsNone(data['vina_version'])
                self.assertEqual(data['engine_default'], 'unavailable')
                self.assertFalse(data['mock_allowed'])
                self.assertEqual(data['status'], 'limited')
    
    def test_capabilities_with_vina_unavailable_mock_enabled(self):
        """Test capabilities when Vina unavailable but mock enabled"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                response = self.client.get(self.capabilities_url)
                self.assertEqual(response.status_code, 200)
                
                data = json.loads(response.content)
                self.assertFalse(data['vina_available'])
                self.assertIsNone(data['vina_version'])
                self.assertEqual(data['engine_default'], 'mock')
                self.assertTrue(data['mock_allowed'])
                self.assertEqual(data['status'], 'operational')
    
    def test_capabilities_with_vina_available(self):
        """Test capabilities when Vina is available"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=True):
            # Mock the import of vina within the view function
            mock_vina = MagicMock()
            mock_vina.__version__ = '1.2.5'
            with patch.dict('sys.modules', {'vina': mock_vina}):
                response = self.client.get(self.capabilities_url)
                self.assertEqual(response.status_code, 200)
                
                data = json.loads(response.content)
                self.assertTrue(data['vina_available'])
                self.assertEqual(data['vina_version'], '1.2.5')
                self.assertEqual(data['engine_default'], 'vina')
                self.assertEqual(data['status'], 'operational')


class TestDockingAllowMockConfiguration(TestCase):
    """Test DOCKING_ALLOW_MOCK configuration behavior"""
    
    def setUp(self):
        self.client = Client()
        self.dock_url = reverse('run-docking')
        self.test_data = {
            'ligand_smiles': 'CCO',  # Simple ethanol
            'receptor_pdb_id': '1CRN',
            'center_x': 0.0,
            'center_y': 0.0,
            'center_z': 0.0,
            'size_x': 20.0,
            'size_y': 20.0,
            'size_z': 20.0,
            'exhaustiveness': 8,
            'num_modes': 9
        }
    
    def test_dock_with_vina_unavailable_mock_disabled_returns_503(self):
        """Test that docking returns 503 when Vina unavailable and mock disabled"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': False}):
                response = self.client.post(
                    self.dock_url,
                    data=json.dumps(self.test_data),
                    content_type='application/json'
                )
                
                self.assertEqual(response.status_code, 503)
                data = json.loads(response.content)
                self.assertIn('AutoDock Vina not available', data['error'])
                self.assertEqual(data['engine'], 'unavailable')
    
    def test_dock_with_vina_unavailable_mock_enabled_succeeds(self):
        """Test that docking succeeds with mock when Vina unavailable but mock enabled"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                response = self.client.post(
                    self.dock_url,
                    data=json.dumps(self.test_data),
                    content_type='application/json'
                )
                
                self.assertEqual(response.status_code, 200)
                data = json.loads(response.content)
                self.assertIn('job_id', data)
                self.assertEqual(data['status'], 'pending')
                
                # Check that job was created
                job = DockingJob.objects.get(job_id=data['job_id'])
                self.assertEqual(job.status, 'pending')
    
    def test_dock_with_vina_available_ignores_mock_setting(self):
        """Test that docking uses Vina when available, regardless of mock setting"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=True):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': False}):
                response = self.client.post(
                    self.dock_url,
                    data=json.dumps(self.test_data),
                    content_type='application/json'
                )
                
                self.assertEqual(response.status_code, 200)
                data = json.loads(response.content)
                self.assertIn('job_id', data)
                self.assertEqual(data['status'], 'pending')


class TestEngineMetadataInResults(TestCase):
    """Test that engine and is_mock metadata is included in results"""
    
    def setUp(self):
        self.test_params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 0.0,
            'center_y': 0.0,
            'center_z': 0.0,
            'size_x': 20.0,
            'size_y': 20.0,
            'size_z': 20.0,
            'exhaustiveness': 8,
            'num_modes': 9
        }
    
    def test_mock_docking_includes_engine_metadata(self):
        """Test that mock docking results include engine='mock' and is_mock=True"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result = DockingEngine.run_production_docking(self.test_params)
                
                self.assertTrue(result.get('success', False))
                self.assertEqual(result.get('engine'), 'mock')
                self.assertTrue(result.get('is_mock'))
    
    def test_real_docking_includes_engine_metadata(self):
        """Test that real docking results include engine='vina' and is_mock=False"""
        mock_real_result = {
            'success': True,
            'poses': [{'mode': 1, 'affinity': -8.5}],
            'calculation_time': 30.5,
            'method': 'AutoDock Vina',
            'version': '1.2.5'
        }
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=True):
            with patch('api.docking_utils.RealDockingUtils.run_production_docking', return_value=mock_real_result):
                result = DockingEngine.run_production_docking(self.test_params)
                
                self.assertTrue(result.get('success', False))
                self.assertEqual(result.get('engine'), 'vina')
                self.assertFalse(result.get('is_mock'))
    
    def test_unavailable_docking_includes_engine_metadata(self):
        """Test that unavailable docking includes proper engine metadata"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': False}):
                result = DockingEngine.run_production_docking(self.test_params)
                
                self.assertFalse(result.get('success', True))
                self.assertEqual(result.get('engine'), 'unavailable')
                self.assertFalse(result.get('is_mock'))
                self.assertIn('AutoDock Vina not available', result.get('error', ''))


class TestJobDatabaseMetadata(TestCase):
    """Test that engine and is_mock fields are stored in database"""
    
    def setUp(self):
        self.job_data = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 0.0,
            'center_y': 0.0,
            'center_z': 0.0,
            'size_x': 20.0,
            'size_y': 20.0,
            'size_z': 20.0,
            'exhaustiveness': 8,
            'num_modes': 9,
            'status': 'pending'
        }
    
    def test_job_model_has_engine_fields(self):
        """Test that DockingJob model has engine and is_mock fields"""
        job = DockingJob.objects.create(**self.job_data)
        
        # Test that fields exist and have proper defaults
        self.assertTrue(hasattr(job, 'engine'))
        self.assertTrue(hasattr(job, 'is_mock'))
        self.assertIsNone(job.engine)  # Default is None
        self.assertFalse(job.is_mock)  # Default is False
    
    def test_job_stores_engine_metadata_on_completion(self):
        """Test that job stores engine metadata when completed"""
        job = DockingJob.objects.create(**self.job_data)
        
        # Simulate completion with engine metadata
        job.status = 'completed'
        job.engine = 'mock'
        job.is_mock = True
        job.results_json = {
            'success': True,
            'engine': 'mock',
            'is_mock': True,
            'poses': []
        }
        job.save()
        
        # Verify storage
        saved_job = DockingJob.objects.get(job_id=job.job_id)
        self.assertEqual(saved_job.engine, 'mock')
        self.assertTrue(saved_job.is_mock)
        self.assertEqual(saved_job.results_json['engine'], 'mock')
        self.assertTrue(saved_job.results_json['is_mock'])
    
    def test_job_status_api_includes_engine_metadata(self):
        """Test that job status API includes engine metadata"""
        # Create job with engine metadata
        job_data = self.job_data.copy()
        job_data['status'] = 'completed'
        job_data['engine'] = 'vina'
        job_data['is_mock'] = False
        job_data['results_json'] = {
            'success': True,
            'engine': 'vina',
            'is_mock': False,
            'poses': [{'mode': 1, 'affinity': -8.5}]
        }
        job = DockingJob.objects.create(**job_data)
        
        client = Client()
        status_url = reverse('get-docking-status', args=[str(job.job_id)])
        response = client.get(status_url)
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        
        # Check that results include engine metadata
        self.assertIn('results', data)
        self.assertEqual(data['results']['engine'], 'vina')
        self.assertFalse(data['results']['is_mock'])


class TestDefaultSettingsBehavior(TestCase):
    """Test default settings behavior for DOCKING_ALLOW_MOCK"""
    
    def test_default_docking_allow_mock_is_false(self):
        """Test that DOCKING_ALLOW_MOCK defaults to False"""
        # Check the actual settings
        self.assertFalse(getattr(settings, 'DOCKING_ALLOW_MOCK', True))  # True as fallback to test it's explicitly False
    
    def test_capabilities_reflects_default_setting(self):
        """Test that capabilities endpoint reflects the default setting"""
        client = Client()
        capabilities_url = reverse('dock-capabilities')
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            response = client.get(capabilities_url)
            self.assertEqual(response.status_code, 200)
            
            data = json.loads(response.content)
            self.assertFalse(data['mock_allowed'])
            self.assertEqual(data['status'], 'limited')


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
