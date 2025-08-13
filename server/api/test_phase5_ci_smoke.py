"""
Tests for Phase 5: CI with Seeded Smoke Test

This module tests the CI infrastructure and smoke test functionality.
"""

import pytest
import json
import tempfile
import os
from unittest.mock import patch, MagicMock
from django.test import TestCase

from api.test_smoke_seeded import SmokeTestRunner


class TestSmokeTestRunner(TestCase):
    """Test the seeded smoke test runner functionality"""
    
    def setUp(self):
        self.runner = SmokeTestRunner(base_url="http://test-server:8000")
    
    def test_runner_initialization(self):
        """Test smoke test runner initializes correctly"""
        self.assertEqual(self.runner.base_url, "http://test-server:8000")
        self.assertEqual(self.runner.api_url, "http://test-server:8000/api")
        self.assertIsNotNone(self.runner.test_timestamp)
        
        # Check test parameters are deterministic
        params = self.runner.test_params
        self.assertEqual(params['seed'], 42)
        self.assertEqual(params['ligand_smiles'], 'CCO')
        self.assertEqual(params['receptor_pdb_id'], '1CRN')
        self.assertEqual(params['exhaustiveness'], 4)
        self.assertEqual(params['num_modes'], 3)
    
    def test_test_parameters_are_fast(self):
        """Test that smoke test parameters are optimized for speed"""
        params = self.runner.test_params
        
        # Should use small search space for speed
        self.assertLessEqual(params['size_x'], 20.0)
        self.assertLessEqual(params['size_y'], 20.0)
        self.assertLessEqual(params['size_z'], 20.0)
        
        # Should use low exhaustiveness for speed
        self.assertLessEqual(params['exhaustiveness'], 8)
        
        # Should use minimal modes for speed
        self.assertLessEqual(params['num_modes'], 5)
        
        # Should use simple molecule (ethanol)
        self.assertEqual(params['ligand_smiles'], 'CCO')
        
        # Should use small protein (crambin)
        self.assertEqual(params['receptor_pdb_id'], '1CRN')
    
    @patch('requests.get')
    def test_api_health_check_success(self, mock_get):
        """Test successful API health check"""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_get.return_value = mock_response
        
        result = self.runner.test_api_health()
        
        self.assertTrue(result)
        mock_get.assert_called_once_with(
            "http://test-server:8000/api/healthz",
            timeout=10
        )
    
    @patch('requests.get')
    def test_api_health_check_failure(self, mock_get):
        """Test failed API health check"""
        mock_response = MagicMock()
        mock_response.status_code = 500
        mock_get.return_value = mock_response
        
        result = self.runner.test_api_health()
        
        self.assertFalse(result)
    
    @patch('requests.get')
    def test_capabilities_check_success(self, mock_get):
        """Test successful capabilities check"""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            'vina_available': True,
            'engine_default': 'vina',
            'vina_version': '1.2.5'
        }
        mock_get.return_value = mock_response
        
        result = self.runner.test_capabilities()
        
        self.assertEqual(result['vina_available'], True)
        self.assertEqual(result['engine_default'], 'vina')
        mock_get.assert_called_once_with(
            "http://test-server:8000/api/dock/capabilities",
            timeout=10
        )
    
    def test_validate_results_success(self):
        """Test results validation with valid results"""
        mock_results = {
            'results': {
                'poses': [
                    {
                        'mode': 1,
                        'affinity': -8.5,
                        'sdf': 'Mock SDF Content\n$$$$'
                    },
                    {
                        'mode': 2,
                        'affinity': -7.8,
                        'sdf': 'Mock SDF Content 2\n$$$$'
                    },
                    {
                        'mode': 3,
                        'affinity': -7.2,
                        'sdf': 'Mock SDF Content 3\n$$$$'
                    }
                ]
            },
            'engine_metadata': {
                'engine': 'vina',
                'is_mock': False
            },
            'docking_parameters': {
                'seed': 42
            }
        }
        
        validation = self.runner.validate_results(mock_results)
        
        self.assertTrue(validation['passed'])
        self.assertEqual(len(validation['errors']), 0)
        
        # Check specific validations
        check_names = [check['name'] for check in validation['checks']]
        self.assertIn('Has results', check_names)
        self.assertIn('Has poses', check_names)
        self.assertIn('Correct number of poses', check_names)
        self.assertIn('SDF not empty', check_names)
        self.assertIn('Parameters preserved', check_names)
    
    def test_validate_results_missing_poses(self):
        """Test results validation with missing poses"""
        mock_results = {
            'results': {
                'poses': []
            },
            'engine_metadata': {
                'engine': 'vina',
                'is_mock': False
            },
            'docking_parameters': {
                'seed': 42
            }
        }
        
        validation = self.runner.validate_results(mock_results)
        
        self.assertFalse(validation['passed'])
        self.assertGreater(len(validation['errors']), 0)
        
        # Should fail pose count check
        error_messages = ' '.join(validation['errors'])
        self.assertIn('poses', error_messages.lower())
    
    def test_validate_results_wrong_seed(self):
        """Test results validation with wrong seed"""
        mock_results = {
            'results': {
                'poses': [
                    {
                        'mode': 1,
                        'affinity': -8.5,
                        'sdf': 'Mock SDF Content\n$$$$'
                    }
                ]
            },
            'engine_metadata': {
                'engine': 'vina',
                'is_mock': False
            },
            'docking_parameters': {
                'seed': 123  # Wrong seed
            }
        }
        
        validation = self.runner.validate_results(mock_results)
        
        self.assertFalse(validation['passed'])
        
        # Should fail seed check
        error_messages = ' '.join(validation['errors'])
        self.assertIn('seed', error_messages.lower())
    
    def test_calculate_result_hash_deterministic(self):
        """Test that result hash calculation is deterministic"""
        mock_results = {
            'results': {
                'poses': [
                    {'affinity': -8.5},
                    {'affinity': -7.8}
                ]
            },
            'docking_parameters': {
                'seed': 42,
                'exhaustiveness': 4,
                'num_modes': 3
            },
            'engine_metadata': {
                'engine': 'vina',
                'is_mock': False
            }
        }
        
        # Calculate hash twice - should be identical
        hash1 = self.runner.calculate_result_hash(mock_results)
        hash2 = self.runner.calculate_result_hash(mock_results)
        
        self.assertEqual(hash1, hash2)
        self.assertIsInstance(hash1, str)
        self.assertEqual(len(hash1), 16)  # SHA256 truncated to 16 chars
    
    def test_calculate_result_hash_different_results(self):
        """Test that different results produce different hashes"""
        mock_results1 = {
            'results': {
                'poses': [{'affinity': -8.5}]
            },
            'docking_parameters': {'seed': 42},
            'engine_metadata': {'engine': 'vina', 'is_mock': False}
        }
        
        mock_results2 = {
            'results': {
                'poses': [{'affinity': -7.5}]  # Different affinity
            },
            'docking_parameters': {'seed': 42},
            'engine_metadata': {'engine': 'vina', 'is_mock': False}
        }
        
        hash1 = self.runner.calculate_result_hash(mock_results1)
        hash2 = self.runner.calculate_result_hash(mock_results2)
        
        self.assertNotEqual(hash1, hash2)
    
    def test_save_artifacts_creates_files(self):
        """Test that save_artifacts creates expected files"""
        mock_results = {
            'results': {
                'poses': [
                    {
                        'mode': 1,
                        'affinity': -8.5,
                        'sdf': 'Mock SDF Content for pose 1\n$$$$'
                    },
                    {
                        'mode': 2,
                        'affinity': -7.8,
                        'sdf': 'Mock SDF Content for pose 2\n$$$$'
                    }
                ]
            }
        }
        
        mock_validation = {
            'passed': True,
            'checks': [],
            'errors': []
        }
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Change to temp directory
            original_cwd = os.getcwd()
            os.chdir(temp_dir)
            
            try:
                output_file = self.runner.save_artifacts(mock_results, mock_validation)
                
                # Check main results file exists
                self.assertTrue(os.path.exists(output_file))
                
                # Check content of main results file
                with open(output_file, 'r') as f:
                    saved_data = json.load(f)
                
                self.assertIn('smoke_test_metadata', saved_data)
                self.assertIn('docking_results', saved_data)
                self.assertIn('validation_report', saved_data)
                self.assertIn('reproducibility_hash', saved_data)
                
                # Check SDF files were created
                sdf_files = [f for f in os.listdir('.') if f.endswith('.sdf')]
                self.assertEqual(len(sdf_files), 2)  # Two poses
                
            finally:
                os.chdir(original_cwd)


class TestCIIntegration(TestCase):
    """Test CI integration aspects"""
    
    def test_environment_variables_handling(self):
        """Test that runner handles CI environment variables"""
        with patch.dict(os.environ, {
            'GITHUB_SHA': 'abc123',
            'GITHUB_REF_NAME': 'main', 
            'GITHUB_RUN_NUMBER': '42'
        }):
            runner = SmokeTestRunner()
            
            # These should be captured in metadata when saving artifacts
            mock_results = {'results': {'poses': []}}
            mock_validation = {'passed': True, 'checks': [], 'errors': []}
            
            with tempfile.TemporaryDirectory() as temp_dir:
                original_cwd = os.getcwd()
                os.chdir(temp_dir)
                
                try:
                    output_file = runner.save_artifacts(mock_results, mock_validation)
                    
                    with open(output_file, 'r') as f:
                        saved_data = json.load(f)
                    
                    metadata = saved_data['smoke_test_metadata']
                    self.assertEqual(metadata['git_commit'], 'abc123')
                    self.assertEqual(metadata['git_ref'], 'main')
                    self.assertEqual(metadata['build_number'], '42')
                    
                finally:
                    os.chdir(original_cwd)
    
    def test_logging_setup(self):
        """Test that logging is properly configured"""
        runner = SmokeTestRunner()
        
        self.assertIsNotNone(runner.logger)
        self.assertIsNotNone(runner.log_file)
        self.assertTrue(runner.log_file.startswith('smoke_test_'))
        self.assertTrue(runner.log_file.endswith('.log'))
    
    def test_fast_test_parameters_for_ci(self):
        """Test that parameters are optimized for CI speed"""
        runner = SmokeTestRunner()
        params = runner.test_params
        
        # CI-optimized parameters
        self.assertEqual(params['exhaustiveness'], 4)  # Fast
        self.assertEqual(params['num_modes'], 3)  # Minimal
        self.assertEqual(params['ligand_smiles'], 'CCO')  # Simple molecule
        self.assertEqual(params['receptor_pdb_id'], '1CRN')  # Small protein
        
        # Search space should be reasonable for CI
        total_search_volume = params['size_x'] * params['size_y'] * params['size_z']
        self.assertLessEqual(total_search_volume, 8000)  # 20^3 max


class TestSmokeTestReproducibility(TestCase):
    """Test reproducibility aspects of smoke test"""
    
    def test_fixed_parameters_ensure_reproducibility(self):
        """Test that all parameters needed for reproducibility are fixed"""
        runner = SmokeTestRunner()
        params = runner.test_params
        
        # Key parameters for reproducibility must be fixed
        self.assertIsNotNone(params['seed'])
        self.assertIsInstance(params['seed'], int)
        self.assertEqual(params['seed'], 42)  # Specific fixed value
        
        # Input must be fixed
        self.assertEqual(params['ligand_smiles'], 'CCO')
        self.assertEqual(params['receptor_pdb_id'], '1CRN')
        
        # Docking parameters must be fixed
        self.assertIsInstance(params['exhaustiveness'], int)
        self.assertIsInstance(params['num_modes'], int)
        
        # Binding site must be fixed
        for coord in ['center_x', 'center_y', 'center_z']:
            self.assertIsInstance(params[coord], (int, float))
        for size in ['size_x', 'size_y', 'size_z']:
            self.assertIsInstance(params[size], (int, float))
    
    def test_reproducibility_hash_includes_key_components(self):
        """Test that reproducibility hash includes all key components"""
        runner = SmokeTestRunner()
        
        mock_results = {
            'results': {
                'poses': [
                    {'affinity': -8.5},
                    {'affinity': -7.8}
                ]
            },
            'docking_parameters': {
                'seed': 42,
                'exhaustiveness': 4,
                'num_modes': 3
            },
            'engine_metadata': {
                'engine': 'vina',
                'is_mock': False
            }
        }
        
        # Mock the hash calculation to capture what's being hashed
        with patch('api.test_smoke_seeded.json.dumps') as mock_dumps:
            mock_dumps.return_value = '{"mocked": "json"}'  # Return valid JSON string
            
            runner.calculate_result_hash(mock_results)
            
            # Check that dumps was called (hash calculation)
            mock_dumps.assert_called_once()
            hash_components = mock_dumps.call_args[0][0]
            
            # Verify key components are included
            self.assertIn('affinities', hash_components)
            self.assertIn('parameters', hash_components)
            self.assertIn('engine', hash_components)
            
            # Verify specific parameter values
            self.assertEqual(hash_components['parameters']['seed'], 42)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
