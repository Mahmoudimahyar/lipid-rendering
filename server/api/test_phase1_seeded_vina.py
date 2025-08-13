"""
Tests for Phase 1: Seeded & Deterministic Vina

This module tests the seed parameter functionality and deterministic behavior.
"""

import pytest
import json
from django.test import TestCase, Client
from django.urls import reverse
from unittest.mock import patch, MagicMock

from api.models import DockingJob
from api.docking_utils import DockingEngine


class TestSeedParameterValidation(TestCase):
    """Test seed parameter validation"""
    
    def setUp(self):
        self.test_params_base = {
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
    
    def test_seed_parameter_validation_valid_integer(self):
        """Test that valid integer seeds are accepted"""
        test_cases = [0, 1, 42, 12345, 999999]
        
        for seed in test_cases:
            with self.subTest(seed=seed):
                params = self.test_params_base.copy()
                params['seed'] = seed
                
                result = DockingEngine.validate_docking_parameters(params)
                self.assertEqual(result['seed'], seed)
    
    def test_seed_parameter_validation_none_allowed(self):
        """Test that None seed is allowed (auto-random)"""
        params = self.test_params_base.copy()
        # Don't include seed parameter
        
        result = DockingEngine.validate_docking_parameters(params)
        self.assertIsNone(result['seed'])
    
    def test_seed_parameter_validation_negative_rejected(self):
        """Test that negative seeds are rejected"""
        params = self.test_params_base.copy()
        params['seed'] = -1
        
        with self.assertRaises(ValueError) as cm:
            DockingEngine.validate_docking_parameters(params)
        
        self.assertIn("non-negative integer", str(cm.exception))
    
    def test_seed_parameter_validation_non_integer_rejected(self):
        """Test that non-integer seeds are rejected"""
        # Note: Python's int() accepts some values we might consider non-integer
        # True -> 1, False -> 0, 3.14 -> 3 (truncated)
        # So we test truly invalid values
        invalid_seeds = ['abc', 'not_a_number', [], {}, object()]
        
        for seed in invalid_seeds:
            with self.subTest(seed=seed):
                params = self.test_params_base.copy()
                params['seed'] = seed
                
                with self.assertRaises(ValueError) as cm:
                    DockingEngine.validate_docking_parameters(params)
                
                self.assertIn("integer", str(cm.exception))
    
    def test_seed_parameter_validation_string_integer_accepted(self):
        """Test that string representations of integers are accepted"""
        params = self.test_params_base.copy()
        params['seed'] = "42"
        
        result = DockingEngine.validate_docking_parameters(params)
        self.assertEqual(result['seed'], 42)


class TestSeedInAPIEndpoints(TestCase):
    """Test seed parameter in API endpoints"""
    
    def setUp(self):
        self.client = Client()
        self.dock_url = reverse('run-docking')
        self.test_data = {
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
    
    def test_api_accepts_seed_parameter(self):
        """Test that API accepts seed parameter"""
        data = self.test_data.copy()
        data['seed'] = 42
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                response = self.client.post(
                    self.dock_url,
                    data=json.dumps(data),
                    content_type='application/json'
                )
                
                self.assertEqual(response.status_code, 200)
                result = json.loads(response.content)
                self.assertIn('job_id', result)
                
                # Check that job was created with seed
                job = DockingJob.objects.get(job_id=result['job_id'])
                self.assertEqual(job.seed, 42)
    
    def test_api_accepts_no_seed_parameter(self):
        """Test that API works without seed parameter"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                response = self.client.post(
                    self.dock_url,
                    data=json.dumps(self.test_data),
                    content_type='application/json'
                )
                
                self.assertEqual(response.status_code, 200)
                result = json.loads(response.content)
                self.assertIn('job_id', result)
                
                # Check that job was created with None seed
                job = DockingJob.objects.get(job_id=result['job_id'])
                self.assertIsNone(job.seed)
    
    def test_api_rejects_invalid_seed(self):
        """Test that API rejects invalid seeds"""
        data = self.test_data.copy()
        data['seed'] = -1
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                response = self.client.post(
                    self.dock_url,
                    data=json.dumps(data),
                    content_type='application/json'
                )
                
                self.assertEqual(response.status_code, 400)
                self.assertIn('non-negative integer', response.content.decode())


class TestDeterministicBehavior(TestCase):
    """Test deterministic behavior with seeds"""
    
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
    
    def test_deterministic_mock_docking_same_seed(self):
        """Test that same seed produces identical results in mock docking"""
        seed = 42
        params1 = self.test_params.copy()
        params1['seed'] = seed
        params2 = self.test_params.copy()
        params2['seed'] = seed
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result1 = DockingEngine.run_production_docking(params1)
                result2 = DockingEngine.run_production_docking(params2)
                
                self.assertTrue(result1.get('success'))
                self.assertTrue(result2.get('success'))
                
                # Results should be identical for same seed
                poses1 = result1.get('poses', [])
                poses2 = result2.get('poses', [])
                
                self.assertEqual(len(poses1), len(poses2))
                
                # Compare pose affinities and coordinates
                for i, (pose1, pose2) in enumerate(zip(poses1, poses2)):
                    with self.subTest(pose=i):
                        self.assertEqual(pose1['affinity'], pose2['affinity'])
                        self.assertEqual(pose1['center_x'], pose2['center_x'])
                        self.assertEqual(pose1['center_y'], pose2['center_y'])
                        self.assertEqual(pose1['center_z'], pose2['center_z'])
                        self.assertEqual(pose1['rmsd_lb'], pose2['rmsd_lb'])
                        self.assertEqual(pose1['rmsd_ub'], pose2['rmsd_ub'])
    
    def test_nondeterministic_mock_docking_different_seeds(self):
        """Test that different seeds produce different results"""
        params1 = self.test_params.copy()
        params1['seed'] = 42
        params2 = self.test_params.copy()
        params2['seed'] = 123
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result1 = DockingEngine.run_production_docking(params1)
                result2 = DockingEngine.run_production_docking(params2)
                
                self.assertTrue(result1.get('success'))
                self.assertTrue(result2.get('success'))
                
                poses1 = result1.get('poses', [])
                poses2 = result2.get('poses', [])
                
                self.assertEqual(len(poses1), len(poses2))
                
                # At least one pose should have different properties
                differences_found = False
                for pose1, pose2 in zip(poses1, poses2):
                    if (pose1['affinity'] != pose2['affinity'] or
                        pose1['center_x'] != pose2['center_x'] or
                        pose1['center_y'] != pose2['center_y'] or
                        pose1['center_z'] != pose2['center_z']):
                        differences_found = True
                        break
                
                self.assertTrue(differences_found, "Different seeds should produce different results")
    
    def test_nondeterministic_mock_docking_no_seed(self):
        """Test that no seed produces different results across runs"""
        params1 = self.test_params.copy()
        params1['seed'] = None
        params2 = self.test_params.copy()
        params2['seed'] = None
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result1 = DockingEngine.run_production_docking(params1)
                result2 = DockingEngine.run_production_docking(params2)
                
                self.assertTrue(result1.get('success'))
                self.assertTrue(result2.get('success'))
                
                poses1 = result1.get('poses', [])
                poses2 = result2.get('poses', [])
                
                # Results should be different for random seeds
                # Note: There's a small chance they could be the same, but very unlikely
                differences_found = False
                for pose1, pose2 in zip(poses1, poses2):
                    if (pose1['affinity'] != pose2['affinity'] or
                        pose1['center_x'] != pose2['center_x']):
                        differences_found = True
                        break
                
                # This test might occasionally fail due to randomness, but very unlikely
                self.assertTrue(differences_found, "Random seeds should produce different results")


class TestSeedInResults(TestCase):
    """Test that seed is reflected in results"""
    
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
    
    def test_seed_reflected_in_mock_results(self):
        """Test that seed is included in mock docking results"""
        seed = 42
        params = self.test_params.copy()
        params['seed'] = seed
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result = DockingEngine.run_production_docking(params)
                
                self.assertTrue(result.get('success'))
                self.assertEqual(result.get('seed_used'), seed)
                self.assertEqual(result.get('parameters', {}).get('seed'), seed)
    
    def test_none_seed_reflected_in_mock_results(self):
        """Test that None seed is included in mock docking results"""
        params = self.test_params.copy()
        params['seed'] = None
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result = DockingEngine.run_production_docking(params)
                
                self.assertTrue(result.get('success'))
                self.assertIsNone(result.get('seed_used'))
                self.assertIsNone(result.get('parameters', {}).get('seed'))
    
    def test_seed_reflected_in_real_vina_results(self):
        """Test that seed is included in real Vina docking results (mocked)"""
        seed = 42
        params = self.test_params.copy()
        params['seed'] = seed
        
        mock_vina_result = {
            'success': True,
            'poses': [{'mode': 1, 'affinity': -8.5}],
            'calculation_time': 30.5,
            'method': 'AutoDock Vina',
            'version': '1.2.5',
            'parameters': {
                'center': (0.0, 0.0, 0.0),
                'size': (20.0, 20.0, 20.0),
                'exhaustiveness': 8,
                'num_modes': 9,
                'seed': seed
            }
        }
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=True):
            with patch('api.docking_utils.RealDockingUtils.run_production_docking', return_value=mock_vina_result):
                result = DockingEngine.run_production_docking(params)
                
                self.assertTrue(result.get('success'))
                self.assertEqual(result.get('parameters', {}).get('seed'), seed)


class TestJobStatusAPISeedReporting(TestCase):
    """Test that job status API reports seed information"""
    
    def setUp(self):
        self.client = Client()
    
    def test_job_status_includes_seed(self):
        """Test that job status API includes seed information"""
        # Create job with seed
        job = DockingJob.objects.create(
            ligand_smiles='CCO',
            receptor_pdb_id='1CRN',
            center_x=0.0,
            center_y=0.0,
            center_z=0.0,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0,
            exhaustiveness=8,
            num_modes=9,
            seed=42,
            status='completed',
            engine='mock',
            is_mock=True,
            results_json={
                'success': True,
                'seed_used': 42,
                'parameters': {'seed': 42},
                'poses': []
            }
        )
        
        status_url = reverse('get-docking-status', args=[str(job.job_id)])
        response = self.client.get(status_url)
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        
        # Check that results include seed information
        self.assertIn('results', data)
        self.assertEqual(data['results']['seed_used'], 42)
        self.assertEqual(data['results']['parameters']['seed'], 42)
    
    def test_job_status_includes_none_seed(self):
        """Test that job status API handles None seed correctly"""
        # Create job without seed
        job = DockingJob.objects.create(
            ligand_smiles='CCO',
            receptor_pdb_id='1CRN',
            center_x=0.0,
            center_y=0.0,
            center_z=0.0,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0,
            exhaustiveness=8,
            num_modes=9,
            seed=None,
            status='completed',
            engine='mock',
            is_mock=True,
            results_json={
                'success': True,
                'seed_used': None,
                'parameters': {'seed': None},
                'poses': []
            }
        )
        
        status_url = reverse('get-docking-status', args=[str(job.job_id)])
        response = self.client.get(status_url)
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        
        # Check that results include None seed information
        self.assertIn('results', data)
        self.assertIsNone(data['results']['seed_used'])
        self.assertIsNone(data['results']['parameters']['seed'])


class TestStructuredLogging(TestCase):
    """Test structured logging functionality"""
    
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
            'num_modes': 9,
            'seed': 42
        }
        self.client = Client()
        self.dock_url = reverse('run-docking')
    
    @patch('api.views.logger')
    def test_structured_logging_start_docking(self, mock_logger):
        """Test that structured logging captures docking start"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                # Mock the background task execution to be synchronous
                with patch('threading.Thread') as mock_thread:
                    def run_task_sync(target, args, **kwargs):
                        target(*args)
                        mock_obj = MagicMock()
                        mock_obj.start = MagicMock()
                        return mock_obj
                    mock_thread.side_effect = run_task_sync
                    
                    response = self.client.post(
                        self.dock_url,
                        data=json.dumps(self.test_params),
                        content_type='application/json'
                    )
                    
                    
                    self.assertEqual(response.status_code, 200)
                    
                    # Check that start logging was called
                    start_calls = [call for call in mock_logger.info.call_args_list 
                                  if 'Starting docking job' in call[0][0]]
                    self.assertTrue(len(start_calls) > 0, "Should log docking start")
                    
                    # Check structured data in start logging
                    start_call = start_calls[0]
                    extra_data = start_call[1]['extra']
                    self.assertEqual(extra_data['seed'], 42)
                    self.assertEqual(extra_data['ligand_smiles'], 'CCO')
                    self.assertEqual(extra_data['receptor_pdb_id'], '1CRN')
                    self.assertEqual(extra_data['exhaustiveness'], 8)
    
    @patch('api.views.logger')
    def test_structured_logging_completed_docking(self, mock_logger):
        """Test that structured logging captures docking completion"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                # Mock the background task execution to be synchronous
                with patch('threading.Thread') as mock_thread:
                    def run_task_sync(target, args, **kwargs):
                        target(*args)
                        mock_obj = MagicMock()
                        mock_obj.start = MagicMock()
                        return mock_obj
                    mock_thread.side_effect = run_task_sync
                    
                    response = self.client.post(
                        self.dock_url,
                        data=json.dumps(self.test_params),
                        content_type='application/json'
                    )
                    
                    self.assertEqual(response.status_code, 200)
                    
                    # Check that completion logging was called
                    complete_calls = [call for call in mock_logger.info.call_args_list 
                                     if 'completed successfully' in call[0][0]]
                    self.assertTrue(len(complete_calls) > 0, "Should log docking completion")
                    
                    # Check structured data in completion logging
                    complete_call = complete_calls[0]
                    extra_data = complete_call[1]['extra']
                    self.assertEqual(extra_data['status'], 'completed')
                    self.assertEqual(extra_data['seed_used'], 42)
                    self.assertIn('elapsed_time_seconds', extra_data)
                    self.assertIn('best_affinity', extra_data)
                    self.assertIn('num_poses', extra_data)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
