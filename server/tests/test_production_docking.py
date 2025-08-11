"""
Tests for production docking implementation with real scientific software.
Validates that the system can fall back gracefully when real software is not available.
"""

import unittest
from unittest.mock import patch, MagicMock
from django.test import TestCase, Client
from django.urls import reverse
import json

from api.docking_utils import DockingEngine
from api.real_docking_engine import RealDockingEngine, RealDockingUtils
from api.real_gnina_scorer import RealGNINAScorer, GNINAScorer
from api.real_pocket_detection import RealPocketDetector, PocketDetector


class ProductionDockingAvailabilityTest(TestCase):
    """Test availability and fallback behavior of production docking software"""
    
    def test_real_docking_availability_check(self):
        """Test availability checking for real docking software"""
        # Test availability check method
        availability = DockingEngine.is_real_docking_available()
        self.assertIsInstance(availability, bool)
        
        # Test that it returns consistent results
        availability2 = DockingEngine.is_real_docking_available()
        self.assertEqual(availability, availability2)
    
    def test_real_gnina_availability_check(self):
        """Test availability checking for real GNINA scoring"""
        availability = RealGNINAScorer.is_available()
        self.assertIsInstance(availability, bool)
    
    def test_real_pocket_detection_availability_check(self):
        """Test availability checking for real pocket detection"""
        availability = RealPocketDetector.is_available()
        self.assertIsInstance(availability, bool)


class ProductionDockingFallbackTest(TestCase):
    """Test fallback behavior when production software is not available"""
    
    @patch('api.docking_utils.REAL_DOCKING_AVAILABLE', False)
    def test_docking_fallback_to_mock(self):
        """Test that docking falls back to mock when real software unavailable"""
        params = {
            'ligand_smiles': 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C',
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
        
        validated_params = DockingEngine.validate_docking_parameters(params)
        results = DockingEngine.run_production_docking(validated_params)
        
        # Should succeed with mock implementation
        self.assertTrue(results.get('success', False))
        self.assertIn('poses', results)
        self.assertIn('summary', results)
        
        # Should indicate mock usage
        summary = results.get('summary', {})
        software = summary.get('software', '')
        self.assertIn('Mock', software)
    
    @patch('api.advanced_docking.REAL_ADVANCED_AVAILABLE', False)
    def test_gnina_fallback_to_mock(self):
        """Test that GNINA scoring falls back to mock when real software unavailable"""
        poses = [
            {'mode': 1, 'affinity': -8.5, 'center_x': 1.0, 'center_y': 2.0, 'center_z': 3.0},
            {'mode': 2, 'affinity': -7.2, 'center_x': 1.5, 'center_y': 2.5, 'center_z': 3.5}
        ]
        
        rescored_poses = GNINAScorer.rescore_poses(poses, "mock_sdf", "mock_pdbqt")
        
        # Should succeed with mock implementation
        self.assertEqual(len(rescored_poses), len(poses))
        for pose in rescored_poses:
            self.assertIn('gnina_score', pose)
            self.assertIn('rescoring_method', pose)
    
    @patch('api.advanced_docking.REAL_ADVANCED_AVAILABLE', False)
    def test_pocket_detection_fallback_to_mock(self):
        """Test that pocket detection falls back to mock when real software unavailable"""
        pockets = PocketDetector.detect_pockets('1CRN', 'mock_pdb_content')
        
        # Should succeed with mock implementation
        self.assertIsInstance(pockets, list)
        self.assertGreater(len(pockets), 0)
        
        for pocket in pockets:
            self.assertIn('pocket_number', pocket)
            self.assertIn('center_x', pocket)
            self.assertIn('center_y', pocket)
            self.assertIn('center_z', pocket)
            self.assertIn('druggability_score', pocket)


class ProductionDockingValidationTest(TestCase):
    """Test parameter validation with production implementations"""
    
    def test_enhanced_smiles_validation(self):
        """Test enhanced SMILES validation with RDKit if available"""
        # Valid SMILES
        valid_params = {
            'ligand_smiles': 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C',  # Cholesterol
            'receptor_pdb_id': '1CRN',
            'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0,
            'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0
        }
        
        try:
            validated = DockingEngine.validate_docking_parameters(valid_params)
            self.assertIn('ligand_smiles', validated)
            self.assertEqual(validated['ligand_smiles'], valid_params['ligand_smiles'])
        except ValueError:
            self.fail("Valid SMILES should not raise ValueError")
        
        # Invalid SMILES
        invalid_params = valid_params.copy()
        invalid_params['ligand_smiles'] = 'INVALID_SMILES_123XYZ'
        
        with self.assertRaises(ValueError):
            DockingEngine.validate_docking_parameters(invalid_params)
    
    def test_enhanced_parameter_validation(self):
        """Test enhanced parameter validation with production requirements"""
        # Test exhaustiveness limits
        params = {
            'ligand_smiles': 'CCO',  # Ethanol
            'receptor_pdb_id': '1CRN',
            'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0,
            'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0,
            'exhaustiveness': 50,  # Too high
            'num_modes': 25  # Too high
        }
        
        # Should either validate with clamping or raise appropriate errors
        try:
            validated = DockingEngine.validate_docking_parameters(params)
            # If validation succeeds, values should be clamped
            self.assertLessEqual(validated['exhaustiveness'], 32)
            self.assertLessEqual(validated['num_modes'], 20)
        except ValueError as e:
            # If validation fails, error message should be informative
            self.assertIn('exhaustiveness', str(e).lower())


class ProductionIntegrationTest(TestCase):
    """Integration tests for production docking through API"""
    
    def setUp(self):
        self.client = Client()
    
    def test_production_docking_api_call(self):
        """Test that API calls work with production docking implementation"""
        data = {
            'smiles': 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C',  # Cholesterol
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
        
        response = self.client.post(
            '/api/dock/run',
            data=json.dumps(data),
            content_type='application/json'
        )
        
        self.assertEqual(response.status_code, 200)
        result = json.loads(response.content)
        
        # Should contain job ID for background processing
        self.assertIn('job_id', result)
        self.assertIn('message', result)
    
    def test_advanced_docking_api_call(self):
        """Test advanced docking API with production features"""
        data = {
            'smiles': 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C',
            'receptor_pdb_id': '1CRN',
            'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0,
            'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0,
            'use_pocket_detection': True,
            'use_gnina_rescoring': True
        }
        
        response = self.client.post(
            '/api/dock/advanced/run',
            data=json.dumps(data),
            content_type='application/json'
        )
        
        self.assertEqual(response.status_code, 200)
        result = json.loads(response.content)
        
        # Should handle advanced features
        self.assertIn('job_id', result)
        self.assertIn('features', result)
        
        features = result.get('features', {})
        self.assertIn('pocket_detection', features)
        self.assertIn('gnina_rescoring', features)


class ProductionPerformanceTest(TestCase):
    """Test performance characteristics of production implementations"""
    
    def test_docking_performance_logging(self):
        """Test that performance metrics are logged appropriately"""
        params = {
            'ligand_smiles': 'CCO',  # Simple molecule for faster testing
            'receptor_pdb_id': '1CRN',
            'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0,
            'size_x': 15.0, 'size_y': 15.0, 'size_z': 15.0,
            'exhaustiveness': 4,  # Lower for faster testing
            'num_modes': 5
        }
        
        validated_params = DockingEngine.validate_docking_parameters(params)
        results = DockingEngine.run_production_docking(validated_params)
        
        # Should include timing information
        if results.get('success'):
            summary = results.get('summary', {})
            self.assertIn('calculation_time', summary)
            
            calc_time = summary.get('calculation_time', 0)
            self.assertGreater(calc_time, 0)  # Should take some time
            self.assertLess(calc_time, 600)   # But not too long for testing
    
    def test_memory_efficient_processing(self):
        """Test that large molecule processing doesn't consume excessive memory"""
        # Test with a moderately complex molecule
        complex_smiles = 'NC(C1=CC=C[N+]([C@H]2[C@H](OC(=O)CCCCCCCCCCCCCCCCC)[C@H](OC(=O)CCCCCCCCCCCCCCCCC)[C@@H](CO)C2)=C1)=O'
        
        params = {
            'ligand_smiles': complex_smiles,
            'receptor_pdb_id': '1CRN',
            'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0,
            'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0,
            'exhaustiveness': 4,
            'num_modes': 5
        }
        
        # Should complete without memory errors
        try:
            validated_params = DockingEngine.validate_docking_parameters(params)
            results = DockingEngine.run_production_docking(validated_params)
            
            # Basic validation that it completed
            self.assertIn('success', results)
            
        except MemoryError:
            self.fail("Processing should not cause memory errors")
        except Exception as e:
            # Other exceptions might be expected if software not available
            # but should be handled gracefully
            self.assertIsInstance(e, (ValueError, RuntimeError))


class ProductionAccuracyTest(TestCase):
    """Test scientific accuracy of production implementations"""
    
    def test_pose_coordinate_accuracy(self):
        """Test that pose coordinates are within reasonable ranges"""
        params = {
            'ligand_smiles': 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C',
            'receptor_pdb_id': '1CRN',
            'center_x': 5.0, 'center_y': 5.0, 'center_z': 5.0,
            'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0
        }
        
        validated_params = DockingEngine.validate_docking_parameters(params)
        results = DockingEngine.run_production_docking(validated_params)
        
        if results.get('success'):
            poses = results.get('poses', [])
            self.assertGreater(len(poses), 0)
            
            for pose in poses:
                # Check pose coordinates are near the binding site center
                pose_x = pose.get('center_x', 0)
                pose_y = pose.get('center_y', 0)
                pose_z = pose.get('center_z', 0)
                
                # Should be within reasonable distance of binding site
                distance_from_center = (
                    (pose_x - 5.0)**2 + (pose_y - 5.0)**2 + (pose_z - 5.0)**2
                )**0.5
                
                self.assertLess(distance_from_center, 15.0, 
                               f"Pose too far from binding site: {distance_from_center:.2f}Å")
                
                # Check affinity values are reasonable
                affinity = pose.get('affinity', 0)
                self.assertLess(affinity, 0)  # Should be negative (favorable)
                self.assertGreater(affinity, -20)  # But not unrealistically strong
    
    def test_gnina_scoring_accuracy(self):
        """Test that GNINA scoring produces reasonable results"""
        poses = [
            {
                'mode': 1, 'affinity': -8.5,
                'center_x': 1.0, 'center_y': 2.0, 'center_z': 3.0
            },
            {
                'mode': 2, 'affinity': -7.2,
                'center_x': 1.5, 'center_y': 2.5, 'center_z': 3.5
            }
        ]
        
        rescored_poses = GNINAScorer.rescore_poses(poses, "mock_sdf", "mock_pdbqt")
        
        for rescored_pose in rescored_poses:
            # GNINA score should be reasonable
            gnina_score = rescored_pose.get('gnina_score', 0)
            self.assertLess(gnina_score, 0)  # Should be negative
            self.assertGreater(gnina_score, -20)  # But not unrealistic
            
            # Should preserve original data
            self.assertIn('mode', rescored_pose)
            self.assertIn('affinity', rescored_pose)
    
    def test_pocket_detection_accuracy(self):
        """Test that pocket detection produces reasonable results"""
        pockets = PocketDetector.detect_pockets('1CRN', 'mock_pdb_content')
        
        self.assertIsInstance(pockets, list)
        self.assertGreater(len(pockets), 0)
        
        for pocket in pockets:
            # Check druggability scores are in valid range
            druggability = pocket.get('druggability_score', 0)
            self.assertGreaterEqual(druggability, 0.0)
            self.assertLessEqual(druggability, 1.0)
            
            # Check volume is reasonable
            volume = pocket.get('volume', 0)
            self.assertGreater(volume, 0)
            self.assertLess(volume, 10000)  # Not unrealistically large


if __name__ == '__main__':
    unittest.main()
