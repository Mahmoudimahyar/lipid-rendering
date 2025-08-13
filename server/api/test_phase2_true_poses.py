"""
Tests for Phase 2: True Docked Poses to Frontend

This module tests SDF extraction from docking poses and spatial differences.
"""

import pytest
import json
import re
from django.test import TestCase, Client
from django.urls import reverse
from unittest.mock import patch, MagicMock

from api.models import DockingJob
from api.docking_utils import DockingEngine


class TestSDFPoseExtraction(TestCase):
    """Test SDF extraction from docking poses"""
    
    def setUp(self):
        self.test_params = {
            'ligand_smiles': 'CCO',  # Ethanol
            'receptor_pdb_id': '1CRN',
            'center_x': 0.0,
            'center_y': 0.0,
            'center_z': 0.0,
            'size_x': 20.0,
            'size_y': 20.0,
            'size_z': 20.0,
            'exhaustiveness': 8,
            'num_modes': 3,
            'seed': 42
        }
    
    def test_mock_docking_includes_valid_sdf(self):
        """Test that mock docking includes valid SDF content for each pose"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result = DockingEngine.run_production_docking(self.test_params)
                
                self.assertTrue(result.get('success'))
                poses = result.get('poses', [])
                self.assertGreater(len(poses), 0)
                
                for i, pose in enumerate(poses):
                    with self.subTest(pose=i):
                        # Check that SDF content exists
                        self.assertIn('sdf', pose)
                        sdf_content = pose['sdf']
                        self.assertIsInstance(sdf_content, str)
                        self.assertGreater(len(sdf_content), 0)
                        
                        # Validate basic SDF structure
                        self.assertIn('$$$$', sdf_content)  # SDF terminator
                        self.assertIn('M  END', sdf_content)  # Molecule end marker
                        self.assertIn('V2000', sdf_content)  # MOL format version
                        
                        # Check pose-specific content
                        self.assertIn(f'Pose {i + 1}', sdf_content)
    
    def test_sdf_content_has_correct_coordinates(self):
        """Test that SDF content includes the correct spatial coordinates"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result = DockingEngine.run_production_docking(self.test_params)
                
                poses = result.get('poses', [])
                
                for pose in poses:
                    sdf_content = pose['sdf']
                    
                    # Extract coordinates from SDF
                    coordinate_lines = []
                    in_atom_block = False
                    
                    for line in sdf_content.split('\n'):
                        # Check if we're in the atom block (after the count line with V2000)
                        if 'V2000' in line:
                            in_atom_block = True
                            continue
                        
                        # End of atom block
                        if line.strip().startswith('M  END'):
                            break
                            
                        # Look for atom coordinate lines (in atom block)
                        if (in_atom_block and line.strip() and
                            len(line.split()) >= 4):
                            parts = line.split()
                            try:
                                # Try to parse first 3 parts as coordinates
                                float(parts[0])
                                float(parts[1]) 
                                float(parts[2])
                                coordinate_lines.append(line)
                            except (ValueError, IndexError):
                                continue
                    
                    # Should have at least some coordinate lines
                    self.assertGreater(len(coordinate_lines), 0)
                    
                    # Check that coordinates are near the docking center
                    for coord_line in coordinate_lines:
                        parts = coord_line.split()
                        if len(parts) >= 3:
                            try:
                                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                                # Coordinates should be within reasonable range of center
                                self.assertGreater(x, -50.0)
                                self.assertLess(x, 50.0)
                                self.assertGreater(y, -50.0)
                                self.assertLess(y, 50.0)
                                self.assertGreater(z, -50.0)
                                self.assertLess(z, 50.0)
                            except ValueError:
                                pass  # Skip non-coordinate lines
    
    def test_poses_have_spatial_differences(self):
        """Test that different poses have different spatial coordinates"""
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result = DockingEngine.run_production_docking(self.test_params)
                
                poses = result.get('poses', [])
                self.assertGreater(len(poses), 1, "Need at least 2 poses to compare")
                
                # Extract first coordinates from each pose
                pose_coords = []
                for pose in poses:
                    sdf_content = pose['sdf']
                    lines = sdf_content.split('\n')
                    in_atom_block = False
                    
                    # Find first coordinate line
                    for line in lines:
                        # Check if we're in the atom block
                        if 'V2000' in line:
                            in_atom_block = True
                            continue
                        
                        # End of atom block
                        if line.strip().startswith('M  END'):
                            break
                            
                        # Look for first atom coordinate line
                        if (in_atom_block and line.strip() and
                            len(line.split()) >= 4):
                            parts = line.split()
                            try:
                                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                                pose_coords.append((x, y, z))
                                break
                            except (ValueError, IndexError):
                                continue
                
                # Check that poses have different coordinates
                self.assertGreater(len(pose_coords), 1)
                for i in range(1, len(pose_coords)):
                    coord1 = pose_coords[0]
                    coord2 = pose_coords[i]
                    # At least one coordinate should be different
                    different = (coord1[0] != coord2[0] or 
                               coord1[1] != coord2[1] or 
                               coord1[2] != coord2[2])
                    self.assertTrue(different, f"Pose 1 and pose {i+1} should have different coordinates")
    
    def test_deterministic_sdf_with_seed(self):
        """Test that same seed produces identical SDF content"""
        seed = 42
        params1 = self.test_params.copy()
        params1['seed'] = seed
        params2 = self.test_params.copy()
        params2['seed'] = seed
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result1 = DockingEngine.run_production_docking(params1)
                result2 = DockingEngine.run_production_docking(params2)
                
                poses1 = result1.get('poses', [])
                poses2 = result2.get('poses', [])
                
                self.assertEqual(len(poses1), len(poses2))
                
                for i, (pose1, pose2) in enumerate(zip(poses1, poses2)):
                    with self.subTest(pose=i):
                        # SDF content should be identical for same seed
                        self.assertEqual(pose1['sdf'], pose2['sdf'])


class TestSDFValidation(TestCase):
    """Test SDF format validation using RDKit"""
    
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
            'num_modes': 3,
            'seed': 42
        }
    
    def test_sdf_parseable_by_rdkit(self):
        """Test that generated SDF content can be parsed by RDKit"""
        try:
            from rdkit import Chem
            rdkit_available = True
        except ImportError:
            rdkit_available = False
        
        if not rdkit_available:
            self.skipTest("RDKit not available for SDF validation")
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result = DockingEngine.run_production_docking(self.test_params)
                
                poses = result.get('poses', [])
                
                for i, pose in enumerate(poses):
                    with self.subTest(pose=i):
                        sdf_content = pose['sdf']
                        
                        # Try to parse SDF with RDKit
                        try:
                            mol = Chem.MolFromMolBlock(sdf_content)
                            # Should successfully create a molecule object
                            self.assertIsNotNone(mol, f"RDKit failed to parse SDF for pose {i+1}")
                            
                            # Check that molecule has atoms
                            if mol:
                                num_atoms = mol.GetNumAtoms()
                                self.assertGreater(num_atoms, 0, f"Pose {i+1} SDF has no atoms")
                                
                        except Exception as e:
                            self.fail(f"RDKit parsing failed for pose {i+1}: {e}")
    
    def test_sdf_has_3d_coordinates(self):
        """Test that SDF content includes 3D coordinates"""
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolTransforms
            rdkit_available = True
        except ImportError:
            rdkit_available = False
        
        if not rdkit_available:
            self.skipTest("RDKit not available for 3D coordinate validation")
        
        with patch.object(DockingEngine, 'is_real_docking_available', return_value=False):
            with patch.dict('django.conf.settings.__dict__', {'DOCKING_ALLOW_MOCK': True}):
                result = DockingEngine.run_production_docking(self.test_params)
                
                poses = result.get('poses', [])
                
                for i, pose in enumerate(poses):
                    with self.subTest(pose=i):
                        sdf_content = pose['sdf']
                        mol = Chem.MolFromMolBlock(sdf_content)
                        
                        if mol and mol.GetNumAtoms() > 0:
                            # Check if molecule has 3D coordinates
                            conf = mol.GetConformer()
                            self.assertIsNotNone(conf, f"Pose {i+1} has no conformer")
                            
                            # Check that coordinates are not all zero
                            positions = conf.GetPositions()
                            non_zero_coords = any(abs(pos[0]) > 0.01 or abs(pos[1]) > 0.01 or abs(pos[2]) > 0.01 
                                                for pos in positions)
                            self.assertTrue(non_zero_coords, f"Pose {i+1} has all zero coordinates")


class TestAPISDFContent(TestCase):
    """Test that API returns poses with SDF content"""
    
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
            'num_modes': 3,
            'seed': 42
        }
    
    def test_docking_api_returns_sdf_in_poses(self):
        """Test that docking API returns SDF content in pose results"""
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
                        data=json.dumps(self.test_data),
                        content_type='application/json'
                    )
                    
                    self.assertEqual(response.status_code, 200)
                    result = json.loads(response.content)
                    job_id = result['job_id']
                    
                    # Get job status to check results
                    status_url = reverse('get-docking-status', args=[job_id])
                    status_response = self.client.get(status_url)
                    
                    self.assertEqual(status_response.status_code, 200)
                    status_data = json.loads(status_response.content)
                    
                    # Check that job completed successfully
                    self.assertEqual(status_data['status'], 'completed')
                    self.assertIn('results', status_data)
                    
                    results = status_data['results']
                    self.assertIn('poses', results)
                    
                    poses = results['poses']
                    self.assertGreater(len(poses), 0)
                    
                    # Check each pose has SDF content
                    for i, pose in enumerate(poses):
                        with self.subTest(pose=i):
                            self.assertIn('sdf', pose)
                            self.assertIsInstance(pose['sdf'], str)
                            self.assertGreater(len(pose['sdf']), 0)
                            self.assertIn('$$$$', pose['sdf'])
    
    def test_job_status_preserves_sdf_content(self):
        """Test that job status API preserves SDF content in database"""
        # Create a completed job with SDF content
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
            num_modes=3,
            seed=42,
            status='completed',
            engine='mock',
            is_mock=True,
            results_json={
                'success': True,
                'poses': [
                    {
                        'mode': 1,
                        'affinity': -8.5,
                        'rmsd_lb': 0.0,
                        'rmsd_ub': 1.2,
                        'center_x': 0.5,
                        'center_y': 0.5,
                        'center_z': 0.5,
                        'sdf': """Test Pose 1
  Mock SDF
  
  3  2  0  0  0  0  0  0  0  0999 V2000
    0.5000    0.5000    0.5000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.5000    0.5000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000    0.5000    0.5000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$"""
                    }
                ]
            }
        )
        
        status_url = reverse('get-docking-status', args=[str(job.job_id)])
        response = self.client.get(status_url)
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        
        # Check that SDF content is preserved
        poses = data['results']['poses']
        self.assertEqual(len(poses), 1)
        
        sdf_content = poses[0]['sdf']
        self.assertIn('Test Pose 1', sdf_content)
        self.assertIn('$$$$', sdf_content)
        self.assertIn('M  END', sdf_content)


class TestRealVinaSDFExtraction(TestCase):
    """Test SDF extraction with real Vina (mocked)"""
    
    def test_real_vina_sdf_extraction_method_exists(self):
        """Test that real Vina implementation has SDF extraction methods"""
        from api.real_docking_engine import RealDockingEngine
        
        engine = RealDockingEngine()
        
        # Check that the new SDF extraction methods exist
        self.assertTrue(hasattr(engine, '_extract_pose_as_sdf'))
        self.assertTrue(hasattr(engine, '_extract_single_pose_from_pdbqt'))
        self.assertTrue(hasattr(engine, '_convert_pdbqt_to_sdf'))
        self.assertTrue(hasattr(engine, '_create_sdf_from_atoms'))
        self.assertTrue(hasattr(engine, '_create_fallback_sdf'))
    
    def test_fallback_sdf_creation(self):
        """Test fallback SDF creation for error cases"""
        from api.real_docking_engine import RealDockingEngine
        
        engine = RealDockingEngine()
        
        # Test fallback SDF
        fallback_sdf = engine._create_fallback_sdf(1)
        
        self.assertIsInstance(fallback_sdf, str)
        self.assertIn('Fallback', fallback_sdf)
        self.assertIn('$$$$', fallback_sdf)
        self.assertIn('M  END', fallback_sdf)
    
    def test_sdf_from_atoms_creation(self):
        """Test SDF creation from atom list"""
        from api.real_docking_engine import RealDockingEngine
        
        engine = RealDockingEngine()
        
        # Test atoms
        atoms = [
            {'element': 'C', 'x': 0.0, 'y': 0.0, 'z': 0.0},
            {'element': 'O', 'x': 1.0, 'y': 0.0, 'z': 0.0},
            {'element': 'H', 'x': -1.0, 'y': 0.0, 'z': 0.0}
        ]
        
        sdf_content = engine._create_sdf_from_atoms(atoms, pose_index=0)
        
        self.assertIsInstance(sdf_content, str)
        self.assertIn('$$$$', sdf_content)
        self.assertIn('M  END', sdf_content)
        self.assertIn('V2000', sdf_content)
        
        # Check atom count in header
        self.assertIn('  3  0', sdf_content)  # 3 atoms, 0 bonds
        
        # Check coordinates are present
        self.assertIn('0.0000', sdf_content)
        self.assertIn('1.0000', sdf_content)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
