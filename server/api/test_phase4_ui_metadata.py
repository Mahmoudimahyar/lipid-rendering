"""
Tests for Phase 4: UI Visibility & Metadata

This module tests API metadata completeness, export functionality,
and comprehensive job information visibility.
"""

import pytest
import json
import uuid
from datetime import datetime, timezone, timedelta
from django.test import TestCase, Client
from django.urls import reverse
from unittest.mock import patch, MagicMock

from api.models import DockingJob


class TestMetadataAPI(TestCase):
    """Test API returns all metadata correctly"""
    
    def setUp(self):
        self.client = Client()
        self.job_id = str(uuid.uuid4())
        
        # Create a test job with comprehensive data
        self.test_job = DockingJob.objects.create(
            job_id=self.job_id,
            ligand_smiles='CCO',
            receptor_pdb_id='1CRN',
            center_x=0.0,
            center_y=0.0,
            center_z=0.0,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0,
            exhaustiveness=8,
            num_modes=5,
            seed=42,
            engine='vina',
            is_mock=False,
            status='completed',
            started_at=datetime.now(timezone.utc) - timedelta(minutes=5),
            completed_at=datetime.now(timezone.utc),
            results_json={
                'success': True,
                'method': 'AutoDock Vina',
                'version': '1.2.5',
                'calculation_time': 120.5,
                'poses': [
                    {
                        'mode': 1,
                        'affinity': -8.5,
                        'rmsd_lb': 0.0,
                        'rmsd_ub': 0.0,
                        'sdf': 'Mock SDF content for pose 1',
                        'center_x': 0.1,
                        'center_y': -0.2,
                        'center_z': 0.3,
                    },
                    {
                        'mode': 2,
                        'affinity': -7.8,
                        'rmsd_lb': 1.2,
                        'rmsd_ub': 1.5,
                        'sdf': 'Mock SDF content for pose 2',
                        'center_x': 0.5,
                        'center_y': -0.1,
                        'center_z': 0.8,
                    }
                ],
                'summary': {
                    'best_affinity': -8.5,
                    'num_poses': 2,
                    'calculation_time': 120.5
                }
            }
        )
    
    def test_get_docking_status_includes_comprehensive_metadata(self):
        """Test that docking status API returns comprehensive metadata"""
        url = reverse('get-docking-status', kwargs={'job_id': self.job_id})
        response = self.client.get(url)
        
        self.assertEqual(response.status_code, 200)
        data = response.json()
        
        # Check basic job information
        self.assertEqual(data['job_id'], self.job_id)
        self.assertEqual(data['status'], 'completed')
        self.assertEqual(data['ligand_smiles'], 'CCO')
        self.assertEqual(data['receptor_pdb_id'], '1CRN')
        
        # Check docking parameters are included
        self.assertIn('docking_parameters', data)
        params = data['docking_parameters']
        self.assertEqual(params['center_x'], 0.0)
        self.assertEqual(params['center_y'], 0.0)
        self.assertEqual(params['center_z'], 0.0)
        self.assertEqual(params['size_x'], 20.0)
        self.assertEqual(params['size_y'], 20.0)
        self.assertEqual(params['size_z'], 20.0)
        self.assertEqual(params['exhaustiveness'], 8)
        self.assertEqual(params['num_modes'], 5)
        self.assertEqual(params['seed'], 42)
        
        # Check engine metadata is included
        self.assertIn('engine_metadata', data)
        engine = data['engine_metadata']
        self.assertEqual(engine['engine'], 'vina')
        self.assertEqual(engine['is_mock'], False)
        self.assertEqual(engine['method'], 'AutoDock Vina')
        self.assertEqual(engine['version'], '1.2.5')
        self.assertEqual(engine['calculation_time'], 120.5)
        
        # Check performance metrics are included
        self.assertIn('performance_metrics', data)
        metrics = data['performance_metrics']
        self.assertEqual(metrics['best_affinity'], -8.5)
        self.assertEqual(metrics['num_poses_generated'], 2)
        self.assertEqual(metrics['calculation_time'], 120.5)
        
        # Check timing information
        self.assertIn('started_at', data)
        self.assertIn('completed_at', data)
        self.assertIn('duration', data)
    
    def test_get_docking_status_mock_engine_indicators(self):
        """Test that mock engine is properly indicated in metadata"""
        # Update job to mock engine
        self.test_job.engine = 'mock'
        self.test_job.is_mock = True
        self.test_job.results_json['method'] = 'Mock Docking Engine'
        self.test_job.results_json['version'] = 'mock-1.0'
        self.test_job.save()
        
        url = reverse('get-docking-status', kwargs={'job_id': self.job_id})
        response = self.client.get(url)
        
        self.assertEqual(response.status_code, 200)
        data = response.json()
        
        # Check engine metadata shows mock
        engine = data['engine_metadata']
        self.assertEqual(engine['engine'], 'mock')
        self.assertEqual(engine['is_mock'], True)
        self.assertEqual(engine['method'], 'Mock Docking Engine')
        self.assertEqual(engine['version'], 'mock-1.0')
    
    def test_get_docking_status_handles_incomplete_jobs(self):
        """Test metadata for incomplete jobs (running, failed)"""
        # Test running job
        self.test_job.status = 'running'
        self.test_job.completed_at = None
        self.test_job.results_json = None
        self.test_job.save()
        
        url = reverse('get-docking-status', kwargs={'job_id': self.job_id})
        response = self.client.get(url)
        
        self.assertEqual(response.status_code, 200)
        data = response.json()
        
        self.assertEqual(data['status'], 'running')
        self.assertIn('docking_parameters', data)
        self.assertIn('engine_metadata', data)
        self.assertNotIn('completed_at', data)
        self.assertNotIn('performance_metrics', data)
    
    def test_get_docking_status_nonexistent_job(self):
        """Test API response for nonexistent job"""
        fake_job_id = str(uuid.uuid4())
        url = reverse('get-docking-status', kwargs={'job_id': fake_job_id})
        response = self.client.get(url)
        
        self.assertEqual(response.status_code, 400)


class TestExportFunctionality(TestCase):
    """Test export includes metadata JSON alongside pose SDF"""
    
    def setUp(self):
        self.client = Client()
        self.job_id = str(uuid.uuid4())
        
        # Create a completed test job with results
        self.test_job = DockingJob.objects.create(
            job_id=self.job_id,
            ligand_smiles='CCO',
            receptor_pdb_id='1CRN',
            center_x=1.0,
            center_y=2.0,
            center_z=3.0,
            size_x=15.0,
            size_y=16.0,
            size_z=17.0,
            exhaustiveness=12,
            num_modes=3,
            seed=123,
            engine='vina',
            is_mock=False,
            status='completed',
            started_at=datetime.now(timezone.utc) - timedelta(minutes=10),
            completed_at=datetime.now(timezone.utc),
            results_json={
                'success': True,
                'method': 'AutoDock Vina',
                'version': '1.2.5',
                'calculation_time': 250.0,
                'poses': [
                    {
                        'mode': 1,
                        'affinity': -9.2,
                        'rmsd_lb': 0.0,
                        'rmsd_ub': 0.0,
                        'sdf': 'SDF content for best pose',
                        'center_x': 1.1,
                        'center_y': 2.1,
                        'center_z': 3.1,
                    },
                    {
                        'mode': 2,
                        'affinity': -8.7,
                        'rmsd_lb': 1.8,
                        'rmsd_ub': 2.1,
                        'sdf': 'SDF content for second pose',
                        'center_x': 1.5,
                        'center_y': 2.5,
                        'center_z': 3.5,
                    }
                ]
            }
        )
    
    def test_export_docking_results_success(self):
        """Test successful export with comprehensive metadata"""
        url = reverse('export-docking-results', kwargs={'job_id': self.job_id})
        response = self.client.get(url)
        
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response['Content-Type'], 'application/json')
        self.assertIn('attachment', response['Content-Disposition'])
        self.assertIn(f'docking_export_{self.job_id}.json', response['Content-Disposition'])
        
        data = response.json()
        
        # Check top-level structure
        self.assertIn('metadata', data)
        self.assertIn('poses_sdf', data)
        
        # Check metadata structure
        metadata = data['metadata']
        
        # Job information
        self.assertIn('job_information', metadata)
        job_info = metadata['job_information']
        self.assertEqual(job_info['job_id'], self.job_id)
        self.assertEqual(job_info['status'], 'completed')
        self.assertIn('created_at', job_info)
        self.assertIn('started_at', job_info)
        self.assertIn('completed_at', job_info)
        self.assertIn('duration_seconds', job_info)
        
        # Input parameters
        self.assertIn('input_parameters', metadata)
        input_params = metadata['input_parameters']
        self.assertEqual(input_params['ligand_smiles'], 'CCO')
        self.assertEqual(input_params['receptor_pdb_id'], '1CRN')
        self.assertEqual(input_params['binding_site']['center'], [1.0, 2.0, 3.0])
        self.assertEqual(input_params['binding_site']['size'], [15.0, 16.0, 17.0])
        self.assertEqual(input_params['docking_parameters']['exhaustiveness'], 12)
        self.assertEqual(input_params['docking_parameters']['num_modes'], 3)
        self.assertEqual(input_params['docking_parameters']['seed'], 123)
        
        # Engine information
        self.assertIn('engine_information', metadata)
        engine_info = metadata['engine_information']
        self.assertEqual(engine_info['engine'], 'vina')
        self.assertEqual(engine_info['is_mock'], False)
        self.assertEqual(engine_info['method'], 'AutoDock Vina')
        self.assertEqual(engine_info['version'], '1.2.5')
        self.assertEqual(engine_info['calculation_time'], 250.0)
        
        # Results summary
        self.assertIn('results_summary', metadata)
        results = metadata['results_summary']
        self.assertEqual(results['success'], True)
        self.assertEqual(results['num_poses'], 2)
        self.assertEqual(results['best_affinity'], -9.2)
        self.assertIsNotNone(results['mean_affinity'])
        
        # Export information
        self.assertIn('export_information', metadata)
        export_info = metadata['export_information']
        self.assertIn('exported_at', export_info)
        self.assertEqual(export_info['export_format'], 'SDF with metadata')
        self.assertEqual(export_info['software_version'], 'Lipid Rendering v1.0')
        
        # Check SDF data
        poses_sdf = data['poses_sdf']
        self.assertIn('pose_1', poses_sdf)
        self.assertIn('pose_2', poses_sdf)
        
        pose_1 = poses_sdf['pose_1']
        self.assertEqual(pose_1['sdf_content'], 'SDF content for best pose')
        self.assertEqual(pose_1['affinity'], -9.2)
        self.assertEqual(pose_1['rmsd_lb'], 0.0)
        self.assertEqual(pose_1['coordinates']['center_x'], 1.1)
        self.assertEqual(pose_1['coordinates']['center_y'], 2.1)
        self.assertEqual(pose_1['coordinates']['center_z'], 3.1)
    
    def test_export_incomplete_job_fails(self):
        """Test export fails for incomplete jobs"""
        # Set job to running status
        self.test_job.status = 'running'
        self.test_job.save()
        
        url = reverse('export-docking-results', kwargs={'job_id': self.job_id})
        response = self.client.get(url)
        
        self.assertEqual(response.status_code, 400)
    
    def test_export_job_without_results_fails(self):
        """Test export fails for jobs without results"""
        # Remove results
        self.test_job.results_json = None
        self.test_job.save()
        
        url = reverse('export-docking-results', kwargs={'job_id': self.job_id})
        response = self.client.get(url)
        
        self.assertEqual(response.status_code, 400)
    
    def test_export_nonexistent_job_fails(self):
        """Test export fails for nonexistent job"""
        fake_job_id = str(uuid.uuid4())
        url = reverse('export-docking-results', kwargs={'job_id': fake_job_id})
        response = self.client.get(url)
        
        self.assertEqual(response.status_code, 400)


class TestJobCapabilities(TestCase):
    """Test dock/capabilities endpoint includes engine status"""
    
    def test_capabilities_includes_engine_information(self):
        """Test capabilities endpoint returns engine information"""
        url = reverse('dock-capabilities')
        
        # Test the endpoint as-is without complex mocking
        # The endpoint should work even without real vina available
        response = self.client.get(url)
        
        self.assertEqual(response.status_code, 200)
        data = response.json()
        
        # Check basic structure exists
        self.assertIn('vina_available', data)
        self.assertIn('engine_default', data)
        self.assertIn('mock_allowed', data)
        
        # The response should be valid regardless of whether vina is actually available
        # In the test environment, vina_available will likely be False
        self.assertIsInstance(data['vina_available'], bool)
        self.assertIn(data['engine_default'], ['vina', 'mock', 'unavailable'])
    
    def test_capabilities_with_mock_only(self):
        """Test capabilities when only mock engine is available"""
        url = reverse('dock-capabilities')
        
        with patch('api.docking_utils.DockingEngine.is_real_docking_available') as mock_available:
            mock_available.return_value = False
            
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                response = self.client.get(url)
        
        self.assertEqual(response.status_code, 200)
        data = response.json()
        
        self.assertFalse(data['vina_available'])
        self.assertEqual(data['engine_default'], 'mock')


class TestAPIIntegration(TestCase):
    """Test complete API metadata integration"""
    
    def test_metadata_consistency_across_endpoints(self):
        """Test that metadata is consistent across different API endpoints"""
        job_id = str(uuid.uuid4())
        
        # Create test job
        test_job = DockingJob.objects.create(
            job_id=job_id,
            ligand_smiles='CCO',
            receptor_pdb_id='1CRN',
            center_x=0.0,
            center_y=0.0,
            center_z=0.0,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0,
            exhaustiveness=8,
            num_modes=5,
            seed=42,
            engine='vina',
            is_mock=False,
            status='completed',
            results_json={
                'success': True,
                'method': 'AutoDock Vina',
                'version': '1.2.5',
                'calculation_time': 120.5,
                'poses': [{
                    'mode': 1,
                    'affinity': -8.5,
                    'sdf': 'Mock SDF content'
                }]
            }
        )
        
        # Get status metadata
        status_url = reverse('get-docking-status', kwargs={'job_id': job_id})
        status_response = self.client.get(status_url)
        status_data = status_response.json()
        
        # Get export metadata
        export_url = reverse('export-docking-results', kwargs={'job_id': job_id})
        export_response = self.client.get(export_url)
        export_data = export_response.json()
        
        # Check consistency between endpoints
        self.assertEqual(status_data['job_id'], export_data['metadata']['job_information']['job_id'])
        self.assertEqual(status_data['ligand_smiles'], export_data['metadata']['input_parameters']['ligand_smiles'])
        self.assertEqual(status_data['engine_metadata']['engine'], export_data['metadata']['engine_information']['engine'])
        self.assertEqual(status_data['engine_metadata']['is_mock'], export_data['metadata']['engine_information']['is_mock'])
        self.assertEqual(status_data['docking_parameters']['seed'], export_data['metadata']['input_parameters']['docking_parameters']['seed'])
    
    def test_metadata_includes_timing_information(self):
        """Test that metadata includes comprehensive timing information"""
        job_id = str(uuid.uuid4())
        start_time = datetime.now(timezone.utc) - timedelta(minutes=5)
        end_time = datetime.now(timezone.utc)
        
        test_job = DockingJob.objects.create(
            job_id=job_id,
            ligand_smiles='CCO',
            receptor_pdb_id='1CRN',
            center_x=0.0, center_y=0.0, center_z=0.0,
            size_x=20.0, size_y=20.0, size_z=20.0,
            exhaustiveness=8, num_modes=5, seed=42,
            engine='vina', is_mock=False, status='completed',
            started_at=start_time,
            completed_at=end_time,
            results_json={'success': True, 'calculation_time': 280.0}
        )
        
        url = reverse('get-docking-status', kwargs={'job_id': job_id})
        response = self.client.get(url)
        data = response.json()
        
        # Check timing fields are present and reasonable
        self.assertIn('started_at', data)
        self.assertIn('completed_at', data)
        self.assertIn('duration', data)
        
        # Duration should be approximately 5 minutes (300 seconds)
        self.assertGreater(data['duration'], 250)
        self.assertLess(data['duration'], 350)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
