"""
Test suite for Phase 2: Real GNINA and Pocket Detection Implementation
"""

import pytest
import json
from django.test import TestCase, Client
from django.urls import reverse
from unittest.mock import patch, MagicMock

from .models import DockingJob, BindingPocket
from .real_gnina_scorer import RealGNINAScorer, GNINAScorer
from .real_pocket_detection import RealPocketDetector, PocketDetector, GeometricPocketDetector


class TestRealGNINAScorer(TestCase):
    """Test real GNINA scorer implementation"""
    
    def setUp(self):
        self.test_ligand_sdf = """
  Mrv2014 01012021

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"""
        self.test_receptor_pdb = "ATOM      1  CA  ALA A   1      20.154  16.967  11.897  1.00 20.00           C"
        self.test_pose = {
            'center_x': 20.0,
            'center_y': 17.0,
            'center_z': 12.0,
            'affinity': -8.5
        }
    
    def test_gnina_scorer_requires_real_implementation(self):
        """Test that GNINAScorer raises RuntimeError when real implementation not available"""
        with patch('api.real_gnina_scorer.RealGNINAScorer.is_available', return_value=False):
            with self.assertRaises(RuntimeError) as context:
                GNINAScorer.rescore_poses([self.test_pose], self.test_ligand_sdf, self.test_receptor_pdb)
            
            self.assertIn("GNINA libraries are not available", str(context.exception))
    
    @patch('api.real_gnina_scorer.GNINA_LIBS_AVAILABLE', True)
    def test_real_gnina_scorer_initialization(self):
        """Test RealGNINAScorer initialization"""
        scorer = RealGNINAScorer()
        self.assertIsNotNone(scorer)
        self.assertIsNotNone(scorer.device)
        self.assertIsNotNone(scorer.temp_dir)
    
    @patch('api.real_gnina_scorer.GNINA_LIBS_AVAILABLE', True)
    def test_real_gnina_scorer_availability(self):
        """Test GNINA availability check"""
        self.assertTrue(RealGNINAScorer.is_available())
    
    @patch('api.real_gnina_scorer.GNINA_LIBS_AVAILABLE', False)
    def test_real_gnina_scorer_unavailability(self):
        """Test GNINA unavailability check"""
        self.assertFalse(RealGNINAScorer.is_available())
    
    @patch('api.real_gnina_scorer.GNINA_LIBS_AVAILABLE', True)
    @patch('api.real_gnina_scorer.torch')
    @patch.object(RealGNINAScorer, '_extract_molecular_features')
    def test_real_gnina_scoring(self, mock_extract_features, mock_torch):
        """Test real GNINA scoring functionality"""
        # Mock PyTorch components
        mock_model = MagicMock()
        mock_model.return_value = MagicMock()
        mock_model.return_value.item.return_value = -9.2
        
        mock_torch.device.return_value = "cpu"
        mock_torch.FloatTensor.return_value.unsqueeze.return_value.to.return_value = "tensor"
        mock_torch.no_grad.return_value.__enter__ = MagicMock()
        mock_torch.no_grad.return_value.__exit__ = MagicMock()
        
        # Mock feature extraction
        import numpy as np
        mock_extract_features.return_value = np.array([0.5] * 1024, dtype=np.float32)
        
        scorer = RealGNINAScorer()
        scorer.model = mock_model
        
        result = scorer.score_pose(self.test_ligand_sdf, self.test_receptor_pdb, self.test_pose)
        
        self.assertIsInstance(result, dict)
        self.assertIn('gnina_score', result)
        self.assertIn('confidence', result)
        self.assertIn('device', result)
    
    def test_no_mock_fallback_available(self):
        """Test that no mock fallback methods exist"""
        self.assertFalse(hasattr(GNINAScorer, '_rescore_poses_mock'))
        self.assertFalse(hasattr(GNINAScorer, '_get_pose_features_mock'))
        self.assertFalse(hasattr(RealGNINAScorer, '_fallback_scoring'))


class TestRealPocketDetector(TestCase):
    """Test real pocket detection implementation"""
    
    def setUp(self):
        self.test_pdb_content = """
ATOM      1  N   ALA A   1      20.154  16.967  11.897  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  17.829  12.309  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.683  17.119  12.527  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.640  15.908  12.789  1.00 20.00           O  
"""
        self.test_pdb_id = "1ABC"
    
    def test_pocket_detector_requires_real_implementation(self):
        """Test that PocketDetector raises RuntimeError when real implementation not available"""
        with patch('api.real_pocket_detection.RealPocketDetector.is_available', return_value=False):
            with self.assertRaises(RuntimeError) as context:
                PocketDetector.detect_pockets(self.test_pdb_id, self.test_pdb_content)
            
            self.assertIn("Pocket detection libraries are not available", str(context.exception))
    
    @patch('api.real_pocket_detection.POCKET_LIBS_AVAILABLE', True)
    def test_real_pocket_detector_availability(self):
        """Test pocket detector availability check"""
        self.assertTrue(RealPocketDetector.is_available())
    
    @patch('api.real_pocket_detection.POCKET_LIBS_AVAILABLE', False)
    def test_real_pocket_detector_unavailability(self):
        """Test pocket detector unavailability check"""
        self.assertFalse(RealPocketDetector.is_available())
    
    @patch('api.real_pocket_detection.POCKET_LIBS_AVAILABLE', True)
    @patch.object(GeometricPocketDetector, '_parse_pdb_content')
    @patch.object(GeometricPocketDetector, '_extract_coordinates')
    @patch.object(GeometricPocketDetector, '_grid_based_detection')
    @patch.object(GeometricPocketDetector, '_alpha_shape_detection')
    def test_geometric_pocket_detection(self, mock_alpha, mock_grid, mock_coords, mock_parse):
        """Test geometric pocket detection"""
        # Mock structure parsing
        mock_structure = MagicMock()
        mock_parse.return_value = mock_structure
        
        # Mock coordinate extraction
        import numpy as np
        mock_coords.return_value = (
            np.array([[20.0, 17.0, 12.0], [21.0, 18.0, 13.0]]),
            ['CA', 'CB']
        )
        
        # Mock cavity detection methods
        mock_grid.return_value = [
            {
                'center': [20.0, 17.0, 12.0],
                'volume': 800.0,
                'druggability_score': 0.8,
                'method': 'grid_clustering'
            }
        ]
        
        mock_alpha.return_value = [
            {
                'center': [21.0, 18.0, 13.0],
                'volume': 600.0,
                'druggability_score': 0.6,
                'method': 'alpha_shape'
            }
        ]
        
        with GeometricPocketDetector() as detector:
            pockets = detector.detect_cavities(self.test_pdb_content)
            
            self.assertIsInstance(pockets, list)
    
    def test_no_mock_fallback_available(self):
        """Test that no mock fallback methods exist"""
        self.assertFalse(hasattr(PocketDetector, '_detect_pockets_mock'))
        self.assertFalse(hasattr(PocketDetector, '_generate_pocket_residues_mock'))
        self.assertFalse(hasattr(GeometricPocketDetector, '_fallback_detection'))
    
    @patch('api.real_pocket_detection.POCKET_LIBS_AVAILABLE', True)
    def test_pocket_druggability_analysis(self):
        """Test pocket druggability analysis"""
        test_pocket = {
            'center': [20.0, 17.0, 12.0],
            'volume': 800.0,
            'surface_area': 400.0,
            'druggability_score': 0.8,
            'hydrophobicity': 0.6,
            'polarity': 0.3,
            'residues': [
                {'residue_name': 'PHE', 'distance_to_center': 3.5},
                {'residue_name': 'ARG', 'distance_to_center': 4.2},
                {'residue_name': 'ALA', 'distance_to_center': 2.8}
            ]
        }
        
        analysis = PocketDetector.analyze_pocket_druggability(test_pocket)
        
        self.assertIsInstance(analysis, dict)
        self.assertIn('druggability_class', analysis)
        self.assertIn('volume_category', analysis)
        self.assertIn('shape_analysis', analysis)
        self.assertIn('chemical_environment', analysis)
        self.assertIn('recommendations', analysis)
        
        # Check that it classified as high druggability
        self.assertEqual(analysis['druggability_class'], 'high')


class TestEnhancedBindingPocketModel(TestCase):
    """Test enhanced BindingPocket model"""
    
    def setUp(self):
        self.pocket_data = {
            'protein_pdb_id': '1ABC',
            'center_x': 20.0,
            'center_y': 17.0,
            'center_z': 12.0,
            'volume': 800.0,
            'surface_area': 400.0,
            'druggability_score': 0.8,
            'confidence_score': 0.9,
            'pocket_rank': 1,
            'detection_method': 'geometric_analysis',
            'detection_software': 'custom_geometric_detector',
            'residues': [{'residue_name': 'PHE', 'distance': 3.5}]
        }
    
    def test_create_binding_pocket(self):
        """Test creating a BindingPocket with enhanced fields"""
        pocket = BindingPocket.objects.create(**self.pocket_data)
        
        self.assertEqual(pocket.protein_pdb_id, '1ABC')
        self.assertEqual(pocket.center_x, 20.0)
        self.assertEqual(pocket.volume, 800.0)
        self.assertEqual(pocket.druggability_score, 0.8)
        self.assertEqual(pocket.detection_method, 'geometric_analysis')
        self.assertEqual(pocket.pocket_rank, 1)
    
    def test_pocket_methods(self):
        """Test BindingPocket model methods"""
        pocket = BindingPocket.objects.create(**self.pocket_data)
        
        # Test get_center_coordinates
        center = pocket.get_center_coordinates()
        self.assertEqual(center, [20.0, 17.0, 12.0])
        
        # Test distance_to_point
        distance = pocket.distance_to_point(21.0, 18.0, 13.0)
        self.assertAlmostEqual(distance, 1.732, places=2)
        
        # Test is_druggable
        self.assertTrue(pocket.is_druggable())
        
        # Test get_druggability_class
        self.assertEqual(pocket.get_druggability_class(), 'high')
    
    def test_pocket_ordering(self):
        """Test that pockets are ordered by druggability score"""
        pocket1 = BindingPocket.objects.create(
            protein_pdb_id='1ABC',
            center_x=20.0, center_y=17.0, center_z=12.0,
            volume=500.0, druggability_score=0.6, confidence_score=0.8,
            residues=[]
        )
        pocket2 = BindingPocket.objects.create(
            protein_pdb_id='1ABC',
            center_x=21.0, center_y=18.0, center_z=13.0,
            volume=800.0, druggability_score=0.9, confidence_score=0.7,
            residues=[]
        )
        
        pockets = list(BindingPocket.objects.all())
        self.assertEqual(pockets[0], pocket2)  # Higher druggability score first
        self.assertEqual(pockets[1], pocket1)


class TestPocketAPIEndpoints(TestCase):
    """Test API endpoints for pocket functionality"""
    
    def setUp(self):
        self.client = Client()
        self.test_pocket = BindingPocket.objects.create(
            protein_pdb_id='1ABC',
            center_x=20.0,
            center_y=17.0,
            center_z=12.0,
            volume=800.0,
            druggability_score=0.8,
            confidence_score=0.9,
            pocket_rank=1,
            detection_method='geometric_analysis',
            residues=[{'residue_name': 'PHE', 'distance': 3.5}]
        )
    
    def test_get_pocket_suggestions_endpoint(self):
        """Test GET /api/pockets/suggestions endpoint"""
        response = self.client.get('/api/pockets/suggestions?pdb_id=1ABC')
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        
        self.assertEqual(data['pdb_id'], '1ABC')
        self.assertIn('pocket_suggestions', data)
        self.assertIn('summary', data)
        self.assertEqual(len(data['pocket_suggestions']), 1)
        
        pocket = data['pocket_suggestions'][0]
        self.assertEqual(pocket['center_x'], 20.0)
        self.assertEqual(pocket['druggability_score'], 0.8)
    
    def test_get_pocket_suggestions_with_filters(self):
        """Test pocket suggestions with druggability filter"""
        response = self.client.get('/api/pockets/suggestions?pdb_id=1ABC&min_druggability=0.9&max_results=5')
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        
        # Should return no results since our pocket has 0.8 druggability
        self.assertEqual(len(data['pocket_suggestions']), 0)
    
    def test_analyze_docking_site_endpoint(self):
        """Test POST /api/pockets/analyze-site endpoint"""
        request_data = {
            'pdb_id': '1ABC',
            'docking_params': {
                'center_x': 20.5,
                'center_y': 17.2,
                'center_z': 12.1
            }
        }
        
        with patch('api.real_pocket_detection.RealPocketDetector.analyze_site_vs_pockets') as mock_analyze:
            mock_analyze.return_value = {
                'match_found': True,
                'distance': 0.6,
                'match_quality': 0.9,
                'recommendation': 'Good site selection'
            }
            
            response = self.client.post(
                '/api/pockets/analyze-site',
                json.dumps(request_data),
                content_type='application/json'
            )
            
            self.assertEqual(response.status_code, 200)
            data = json.loads(response.content)
            
            self.assertEqual(data['pdb_id'], '1ABC')
            self.assertIn('analysis', data)
            self.assertTrue(data['analysis']['match_found'])
    
    def test_analyze_docking_site_no_pockets(self):
        """Test analyze docking site when no pockets exist"""
        request_data = {
            'pdb_id': '2XYZ',
            'docking_params': {
                'center_x': 20.0,
                'center_y': 17.0,
                'center_z': 12.0
            }
        }
        
        response = self.client.post(
            '/api/pockets/analyze-site',
            json.dumps(request_data),
            content_type='application/json'
        )
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        
        self.assertFalse(data['analysis']['match_found'])
        self.assertIn('No pockets found', data['analysis']['message'])
    
    @patch('api.real_pocket_detection.RealPocketDetector.detect_pockets')
    @patch('api.chem_utils.ChemUtils.fetch_protein_structure')
    def test_detect_binding_pockets_enhanced(self, mock_fetch, mock_detect):
        """Test enhanced pocket detection endpoint"""
        mock_fetch.return_value = "PDB content"
        mock_detect.return_value = [
            {
                'center': [20.0, 17.0, 12.0],
                'volume': 800.0,
                'druggability_score': 0.8,
                'rank': 1,
                'confidence': 0.9,
                'detection_method': 'geometric_analysis',
                'software': 'custom_geometric_detector'
            }
        ]
        
        response = self.client.post(
            '/api/pockets/detect',
            json.dumps({'pdb_id': '1XYZ'}),
            content_type='application/json'
        )
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        
        self.assertEqual(data['pdb_id'], '1XYZ')
        self.assertIn('pockets', data)
        self.assertEqual(len(data['pockets']), 1)
        
        # Check that pocket was saved to database
        saved_pocket = BindingPocket.objects.get(protein_pdb_id='1XYZ')
        self.assertEqual(saved_pocket.center_x, 20.0)
        self.assertEqual(saved_pocket.druggability_score, 0.8)


class TestProductionReadiness(TestCase):
    """Test that Phase 2 is production-ready"""
    
    def test_no_debug_or_mock_code(self):
        """Test that no mock implementation methods exist in production files"""
        # Test GNINA scorer - check that mock methods don't exist
        self.assertFalse(hasattr(GNINAScorer, '_rescore_poses_mock'))
        self.assertFalse(hasattr(GNINAScorer, '_get_pose_features_mock'))
        self.assertFalse(hasattr(RealGNINAScorer, '_fallback_scoring'))
        
        # Test pocket detector - check that mock methods don't exist
        self.assertFalse(hasattr(PocketDetector, '_detect_pockets_mock'))
        self.assertFalse(hasattr(PocketDetector, '_generate_pocket_residues_mock'))
        self.assertFalse(hasattr(GeometricPocketDetector, '_fallback_detection'))
    
    def test_error_handling_without_libraries(self):
        """Test proper error handling when scientific libraries are missing"""
        with patch('api.real_gnina_scorer.GNINA_LIBS_AVAILABLE', False):
            with self.assertRaises(RuntimeError):
                GNINAScorer.rescore_poses([], "", "")
        
        with patch('api.real_pocket_detection.POCKET_LIBS_AVAILABLE', False):
            with self.assertRaises(RuntimeError):
                PocketDetector.detect_pockets("1ABC", "")
    
    def test_docker_dependencies_installation(self):
        """Test that Docker containers include necessary dependencies"""
        import os
        from pathlib import Path
        
        project_root = Path(__file__).parent.parent
        
        # Read Dockerfile.cuda
        cuda_dockerfile_path = project_root / 'Dockerfile.cuda'
        if cuda_dockerfile_path.exists():
            with open(cuda_dockerfile_path, 'r') as f:
                cuda_dockerfile = f.read()
            
            # Check for PyTorch installation
            self.assertIn('pytorch', cuda_dockerfile)
            self.assertIn('pytorch-cuda', cuda_dockerfile)
            
            # Check for GNINA installation
            self.assertIn('gnina', cuda_dockerfile)
        
        # Read Dockerfile.cpu
        cpu_dockerfile_path = project_root / 'Dockerfile.cpu'
        if cpu_dockerfile_path.exists():
            with open(cpu_dockerfile_path, 'r') as f:
                cpu_dockerfile = f.read()
            
            # Check for PyTorch CPU installation
            self.assertIn('pytorch', cpu_dockerfile)
            self.assertIn('cpuonly', cpu_dockerfile)
            
            # Check for GNINA installation
            self.assertIn('gnina', cpu_dockerfile)


# Import inspection for production readiness test
import inspect
