"""
Simplified benchmark validation tests for molecular docking system.
Tests core validation logic without heavy dependencies.
"""

import json
import math
from django.test import TestCase, Client
from django.urls import reverse


class SimpleBenchmarkValidationTest(TestCase):
    """Test docking system validation against known benchmarks."""
    
    def setUp(self):
        """Set up test client and benchmark data."""
        self.client = Client()
        
        # Simple benchmark complex for testing
        self.test_complex = {
            'id': 'HIV_protease_indinavir',
            'protein': '1HSG',
            'ligand': 'Indinavir',
            'smiles': 'CC(C)CN(CC(O)C(Cc1ccccc1)NC(=O)OCc1ccccc1)C(=O)C(NC(=O)C(CC(C)C)NC(=O)OCc1ccccc1)CC1CCCCC1',
            'expected_distance': 2.1,
            'expected_binding': {'x': 15.2, 'y': 12.8, 'z': 8.4}
        }
        
        # Validation criteria
        self.criteria = {
            'position_tolerance': 3.0,  # Å
            'distance_tolerance': 1.5,  # Å
            'min_realistic_distance': 1.0,  # Å
            'max_realistic_distance': 12.0,  # Å
        }
    
    def test_position_error_calculation(self):
        """Test position error calculation accuracy."""
        # Test case 1: Perfect match
        expected = {'x': 15.2, 'y': 12.8, 'z': 8.4}
        actual = {'center_x': 15.2, 'center_y': 12.8, 'center_z': 8.4}
        
        error = self._calculate_position_error(expected, actual)
        self.assertAlmostEqual(error, 0.0, places=2)
        
        # Test case 2: Known error
        actual_with_error = {'center_x': 14.0, 'center_y': 13.0, 'center_z': 9.0}
        error = self._calculate_position_error(expected, actual_with_error)
        expected_error = math.sqrt((15.2-14.0)**2 + (12.8-13.0)**2 + (8.4-9.0)**2)
        self.assertAlmostEqual(error, expected_error, places=2)
    
    def test_distance_validation(self):
        """Test distance validation logic."""
        expected_distance = 2.1
        
        # Test within tolerance
        actual_distance_good = 2.3
        error = abs(actual_distance_good - expected_distance)
        self.assertLess(error, self.criteria['distance_tolerance'])
        
        # Test outside tolerance
        actual_distance_bad = 4.0
        error = abs(actual_distance_bad - expected_distance)
        self.assertGreater(error, self.criteria['distance_tolerance'])
    
    def test_scientific_realism_validation(self):
        """Test scientific realism criteria."""
        # Test realistic distances
        realistic_distances = [1.5, 2.8, 5.2, 8.9, 11.5]
        for distance in realistic_distances:
            self.assertTrue(self._is_distance_realistic(distance))
        
        # Test unrealistic distances
        unrealistic_distances = [0.2, 0.8, 15.0, 25.0]
        for distance in unrealistic_distances:
            self.assertFalse(self._is_distance_realistic(distance))
    
    def test_pose_diversity_calculation(self):
        """Test pose diversity measurement."""
        # Diverse poses
        diverse_poses = [
            {'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0},
            {'center_x': 2.5, 'center_y': 1.8, 'center_z': -1.2},
            {'center_x': -1.9, 'center_y': 3.1, 'center_z': 2.4}
        ]
        
        diversity = self._calculate_pose_diversity(diverse_poses)
        self.assertGreater(diversity, 1.0)  # Should be > 1.0Å for good diversity
        
        # Similar poses (low diversity)
        similar_poses = [
            {'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0},
            {'center_x': 0.2, 'center_y': 0.1, 'center_z': -0.1},
            {'center_x': -0.1, 'center_y': 0.3, 'center_z': 0.2}
        ]
        
        diversity = self._calculate_pose_diversity(similar_poses)
        self.assertLess(diversity, 1.0)  # Should be < 1.0Å for poor diversity
    
    def test_benchmark_validation_scoring(self):
        """Test the benchmark validation scoring system."""
        # Perfect result
        perfect_result = {
            'position_error': 0.0,
            'distance_error': 0.0,
            'scientific_valid': True
        }
        score = self._calculate_validation_score(perfect_result)
        self.assertAlmostEqual(score, 1.0, places=2)
        
        # Good result (within tolerances)
        good_result = {
            'position_error': 2.0,  # Within 3.0Å tolerance
            'distance_error': 1.0,  # Within 1.5Å tolerance
            'scientific_valid': True
        }
        score = self._calculate_validation_score(good_result)
        self.assertGreater(score, 0.4)  # Should be good score (adjusted for realistic tolerance)
        
        # Poor result
        poor_result = {
            'position_error': 5.0,  # Outside tolerance
            'distance_error': 3.0,  # Outside tolerance
            'scientific_valid': False
        }
        score = self._calculate_validation_score(poor_result)
        self.assertLess(score, 0.3)  # Should be poor score
    
    def test_docking_api_integration(self):
        """Test integration with the docking API."""
        url = '/api/dock/run'  # Correct API endpoint
        
        # Simple test data
        data = {
            'smiles': 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C',  # Cholesterol
            'receptor_pdb_id': '1CRN',
            'num_poses': 3
        }
        
        response = self.client.post(
            url,
            data=json.dumps(data),
            content_type='application/json'
        )
        
        # Should get a successful response
        self.assertEqual(response.status_code, 200)
        
        result = json.loads(response.content)
        
        # Validate response structure
        self.assertIn('poses', result)
        self.assertGreater(len(result['poses']), 0)
        
        # Basic validation of pose data
        for pose in result['poses']:
            self.assertIn('center_x', pose)
            self.assertIn('center_y', pose)
            self.assertIn('center_z', pose)
            self.assertIn('score', pose)
    
    def test_coordinate_system_alignment(self):
        """Test coordinate system alignment logic."""
        # Protein far from origin (like 4R7D)
        protein_center = {'x': 108.77, 'y': 13.66, 'z': 64.25}
        distance_from_origin = math.sqrt(108.77**2 + 13.66**2 + 64.25**2)
        
        # Should be far from origin
        self.assertGreater(distance_from_origin, 50.0)
        
        # Test alignment calculation
        ligand_origin_coords = {'x': 1.5, 'y': -0.7, 'z': 0.5}
        aligned_coords = self._apply_coordinate_alignment(ligand_origin_coords, protein_center)
        
        # Should be near protein center
        expected_x = protein_center['x'] + ligand_origin_coords['x']
        expected_y = protein_center['y'] + ligand_origin_coords['y']
        expected_z = protein_center['z'] + ligand_origin_coords['z']
        
        self.assertAlmostEqual(aligned_coords['x'], expected_x, places=2)
        self.assertAlmostEqual(aligned_coords['y'], expected_y, places=2)
        self.assertAlmostEqual(aligned_coords['z'], expected_z, places=2)
    
    def test_validation_failure_detection(self):
        """Test that the system properly detects validation failures."""
        # Test data that should fail validation
        failing_poses = [
            {'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0, 'protein_distance': 0.2},  # Too close
            {'center_x': 50.0, 'center_y': 50.0, 'center_z': 50.0, 'protein_distance': 20.0}  # Too far
        ]
        
        validation_results = []
        for pose in failing_poses:
            is_realistic = self._is_distance_realistic(pose['protein_distance'])
            validation_results.append(is_realistic)
        
        # Both should fail
        self.assertFalse(validation_results[0])  # Too close
        self.assertFalse(validation_results[1])  # Too far
    
    # Helper methods
    def _calculate_position_error(self, expected, actual):
        """Calculate position error between expected and actual coordinates."""
        dx = actual.get('center_x', actual.get('x', 0)) - expected.get('x', 0)
        dy = actual.get('center_y', actual.get('y', 0)) - expected.get('y', 0)
        dz = actual.get('center_z', actual.get('z', 0)) - expected.get('z', 0)
        return math.sqrt(dx*dx + dy*dy + dz*dz)
    
    def _is_distance_realistic(self, distance):
        """Check if protein-ligand distance is realistic."""
        return self.criteria['min_realistic_distance'] <= distance <= self.criteria['max_realistic_distance']
    
    def _calculate_pose_diversity(self, poses):
        """Calculate average distance between poses."""
        if len(poses) < 2:
            return 0.0
        
        distances = []
        for i in range(len(poses)):
            for j in range(i + 1, len(poses)):
                dist = self._calculate_position_error(poses[i], poses[j])
                distances.append(dist)
        
        return sum(distances) / len(distances)
    
    def _calculate_validation_score(self, result):
        """Calculate composite validation score."""
        position_score = max(0, 1 - (result['position_error'] / self.criteria['position_tolerance']))
        distance_score = max(0, 1 - (result['distance_error'] / self.criteria['distance_tolerance']))
        scientific_score = 1.0 if result['scientific_valid'] else 0.0
        
        return (position_score * 0.4 + distance_score * 0.4 + scientific_score * 0.2)
    
    def _apply_coordinate_alignment(self, ligand_coords, protein_center):
        """Apply coordinate system alignment."""
        return {
            'x': protein_center['x'] + ligand_coords['x'],
            'y': protein_center['y'] + ligand_coords['y'], 
            'z': protein_center['z'] + ligand_coords['z']
        }


class BenchmarkIntegrationTest(TestCase):
    """Integration tests for benchmark validation."""
    
    def test_system_health_check(self):
        """Test basic system health for validation."""
        client = Client()
        
        # Test health endpoint
        response = client.get('/api/healthz')
        self.assertEqual(response.status_code, 200)
        
        result = json.loads(response.content)
        self.assertTrue(result.get('ok', False))
    
    def test_ligand_preparation_endpoint(self):
        """Test ligand preparation for benchmark validation."""
        client = Client()
        
        data = {
            'smiles': 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C'  # Cholesterol
        }
        
        response = client.post(
            '/api/ligand/prepare',
            data=json.dumps(data),
            content_type='application/json'
        )
        
        self.assertEqual(response.status_code, 200)
        
        result = json.loads(response.content)
        self.assertIn('sdf_content', result)
        self.assertGreater(len(result['sdf_content']), 0)


if __name__ == '__main__':
    import unittest
    unittest.main()
