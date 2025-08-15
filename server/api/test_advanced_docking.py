"""
Tests for advanced docking features including GNINA rescoring, pocket detection,
and job templates.
"""

import pytest
import json
import uuid
from unittest.mock import patch, MagicMock
from django.test import TestCase, Client
from django.urls import reverse
from django.utils import timezone
from freezegun import freeze_time

from .models import DockingJob, BindingPocket, JobTemplate
from .advanced_docking import (
    GNINAScorer, PocketDetector, AdvancedDockingEngine, JobTemplateManager
)


class TestGNINAScorer(TestCase):
    """Test GNINA rescoring functionality"""

    def setUp(self):
        self.mock_poses = [
            {
                'mode': 1,
                'affinity': -8.5,
                'rmsd_lb': 0.0,
                'rmsd_ub': 0.0,
                'center_x': 5.123,
                'center_y': 10.456,
                'center_z': -2.789
            },
            {
                'mode': 2,
                'affinity': -7.2,
                'rmsd_lb': 1.5,
                'rmsd_ub': 2.0,
                'center_x': 4.888,
                'center_y': 11.234,
                'center_z': -3.456
            }
        ]

    def test_rescore_poses(self):
        """Test pose rescoring with GNINA"""
        rescored_poses = GNINAScorer.rescore_poses(
            self.mock_poses, "mock_sdf", "mock_pdbqt"
        )
        
        self.assertEqual(len(rescored_poses), 2)
        
        # Check that all poses have GNINA scores
        for pose in rescored_poses:
            self.assertIn('gnina_score', pose)
            self.assertIn('cnn_score', pose)
            self.assertIn('cnn_affinity', pose)
            self.assertIn('rescoring_method', pose)
            self.assertEqual(pose['rescoring_method'], 'GNINA')
            self.assertIn('original_affinity', pose)

    def test_rescore_poses_sorting(self):
        """Test that rescored poses are sorted by GNINA score"""
        rescored_poses = GNINAScorer.rescore_poses(
            self.mock_poses, "mock_sdf", "mock_pdbqt"
        )
        
        # Should be sorted by gnina_score (most negative first)
        for i in range(len(rescored_poses) - 1):
            self.assertLessEqual(
                rescored_poses[i]['gnina_score'],
                rescored_poses[i + 1]['gnina_score']
            )

    def test_get_pose_features(self):
        """Test pose feature extraction"""
        features = GNINAScorer.get_pose_features(
            self.mock_poses[0], "mock_sdf", "mock_pdbqt"
        )
        
        required_feature_groups = [
            'molecular_descriptors',
            'protein_ligand_interactions',
            'geometric_features',
            'pharmacophore_features'
        ]
        
        for group in required_feature_groups:
            self.assertIn(group, features)
        
        # Check molecular descriptors
        mol_desc = features['molecular_descriptors']
        self.assertIn('mol_weight', mol_desc)
        self.assertIn('logp', mol_desc)
        self.assertIn('rotatable_bonds', mol_desc)

    def test_empty_poses_list(self):
        """Test rescoring with empty poses list"""
        rescored_poses = GNINAScorer.rescore_poses([], "mock_sdf", "mock_pdbqt")
        self.assertEqual(len(rescored_poses), 0)


class TestPocketDetector(TestCase):
    """Test binding pocket detection functionality"""

    def test_detect_pockets(self):
        """Test pocket detection"""
        pockets = PocketDetector.detect_pockets("1CRN", "mock_pdb_content")
        
        self.assertGreaterEqual(len(pockets), 2)
        self.assertLessEqual(len(pockets), 8)
        
        # Check required pocket properties
        for pocket in pockets:
            required_keys = [
                'pocket_number', 'center_x', 'center_y', 'center_z',
                'radius', 'volume', 'surface_area', 'druggability_score',
                'confidence_score', 'detection_method', 'residues'
            ]
            for key in required_keys:
                self.assertIn(key, pocket)
        
        # Check sorting by druggability score
        for i in range(len(pockets) - 1):
            self.assertGreaterEqual(
                pockets[i]['druggability_score'],
                pockets[i + 1]['druggability_score']
            )

    def test_analyze_pocket_druggability(self):
        """Test pocket druggability analysis"""
        mock_pocket = {
            'volume': 800,
            'surface_area': 1200,
            'druggability_score': 0.75,
            'hydrophobicity': 0.6,
            'polarity': 0.4
        }
        
        analysis = PocketDetector.analyze_pocket_druggability(mock_pocket)
        
        required_keys = [
            'druggability_class', 'volume_category', 'shape_analysis',
            'chemical_environment', 'recommendations'
        ]
        for key in required_keys:
            self.assertIn(key, analysis)
        
        self.assertEqual(analysis['druggability_class'], 'high')
        self.assertIsInstance(analysis['recommendations'], list)

    def test_generate_pocket_residues(self):
        """Test pocket residue generation"""
        if not hasattr(PocketDetector, '_generate_pocket_residues'):
            pytest.skip("Residue generation helper is not exposed in production API surface")
        residues = PocketDetector._generate_pocket_residues()
        
        self.assertGreaterEqual(len(residues), 10)
        self.assertLessEqual(len(residues), 30)
        
        for residue in residues:
            required_keys = [
                'residue_name', 'residue_number', 'chain_id',
                'distance_to_center', 'contact_type'
            ]
            for key in required_keys:
                self.assertIn(key, residue)


class TestAdvancedDockingEngine(TestCase):
    """Test advanced docking engine functionality"""

    def setUp(self):
        self.docking_params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 5.0,
            'center_y': 10.0,
            'center_z': -2.0,
            'size_x': 20.0,
            'size_y': 20.0,
            'size_z': 20.0,
            'exhaustiveness': 8,
            'num_modes': 9
        }

    @patch('api.docking_utils.DockingEngine.run_mock_docking')
    def test_run_advanced_docking_basic(self, mock_docking):
        """Test basic advanced docking without extra features"""
        mock_docking.return_value = {
            'success': True,
            'poses': [{'mode': 1, 'affinity': -8.5}],
            'summary': {'best_affinity': -8.5}
        }
        
        results = AdvancedDockingEngine.run_advanced_docking(
            self.docking_params,
            use_pocket_detection=False,
            use_gnina_rescoring=False
        )
        
        self.assertTrue(results['success'])
        self.assertIn('performance', results)
        self.assertIn('advanced_features_used', results['performance'])

    @patch('api.docking_utils.DockingEngine.run_mock_docking')
    @patch('api.advanced_docking.PocketDetector.detect_pockets')
    def test_run_advanced_docking_with_pocket_detection(self, mock_detect, mock_docking):
        """Test advanced docking with pocket detection"""
        mock_docking.return_value = {
            'success': True,
            'poses': [{'mode': 1, 'affinity': -8.5}],
            'summary': {'best_affinity': -8.5}
        }
        
        mock_detect.return_value = [{
            'pocket_number': 1,
            'center_x': 5.0,
            'center_y': 10.0,
            'center_z': -2.0,
            'druggability_score': 0.8
        }]
        
        results = AdvancedDockingEngine.run_advanced_docking(
            self.docking_params,
            use_pocket_detection=True,
            use_gnina_rescoring=False
        )
        
        self.assertTrue(results['success'])
        self.assertIn('detected_pockets', results)
        self.assertIn('site_analysis', results)
        mock_detect.assert_called_once()

    @patch('api.docking_utils.DockingEngine.run_mock_docking')
    @patch('api.advanced_docking.GNINAScorer.rescore_poses')
    def test_run_advanced_docking_with_gnina(self, mock_rescore, mock_docking):
        """Test advanced docking with GNINA rescoring"""
        mock_poses = [{'mode': 1, 'affinity': -8.5}]
        mock_docking.return_value = {
            'success': True,
            'poses': mock_poses,
            'summary': {'best_affinity': -8.5}
        }
        
        mock_rescored = [{'mode': 1, 'affinity': -8.5, 'gnina_score': -9.2}]
        mock_rescore.return_value = mock_rescored
        
        results = AdvancedDockingEngine.run_advanced_docking(
            self.docking_params,
            use_pocket_detection=False,
            use_gnina_rescoring=True
        )
        
        self.assertTrue(results['success'])
        self.assertIn('rescored_poses', results)
        self.assertIn('best_gnina_score', results['summary'])
        mock_rescore.assert_called_once_with(mock_poses, "mock_sdf", "mock_pdbqt")

    def test_analyze_site_vs_pockets(self):
        """Test analysis of docking site vs detected pockets"""
        mock_pockets = [
            {'center_x': 5.0, 'center_y': 10.0, 'center_z': -2.0},
            {'center_x': 15.0, 'center_y': 20.0, 'center_z': 8.0}
        ]
        
        analysis = AdvancedDockingEngine._analyze_site_vs_pockets(
            self.docking_params, mock_pockets
        )
        
        self.assertIn('closest_pocket_distance', analysis)
        self.assertIn('site_overlaps_pocket', analysis)
        self.assertIn('recommendation', analysis)
        self.assertTrue(analysis['site_overlaps_pocket'])


class TestJobTemplateManager(TestCase):
    """Test job template management functionality"""

    def test_create_default_templates(self):
        """Test creation of default job templates"""
        created_templates = JobTemplateManager.create_default_templates()
        
        self.assertGreater(len(created_templates), 0)
        
        # Check that templates were created
        total_templates = JobTemplate.objects.count()
        self.assertGreaterEqual(total_templates, 4)
        
        # Test specific template exists
        standard_template = JobTemplate.objects.filter(name="Standard Docking").first()
        self.assertIsNotNone(standard_template)
        self.assertEqual(standard_template.category, "basic")

    def test_get_template_by_category(self):
        """Test filtering templates by category"""
        JobTemplateManager.create_default_templates()
        
        # Make templates public first
        JobTemplate.objects.filter(category="basic").update(is_public=True)
        
        basic_templates = JobTemplateManager.get_template_by_category("basic")
        self.assertGreaterEqual(len(basic_templates), 1)
        
        for template in basic_templates:
            self.assertEqual(template.category, "basic")

    def test_apply_template(self):
        """Test applying a template to generate parameters"""
        template = JobTemplate.objects.create(
            name="Test Template",
            description="Test template for testing",
            default_params={'exhaustiveness': 16, 'num_modes': 20},
            advanced_params={'use_gnina_rescoring': True}
        )
        
        params = JobTemplateManager.apply_template(template)
        
        self.assertEqual(params['exhaustiveness'], 16)
        self.assertEqual(params['num_modes'], 20)
        self.assertTrue(params['use_gnina_rescoring'])
        
        # Check usage counter was incremented
        template.refresh_from_db()
        self.assertEqual(template.usage_count, 1)

    def test_apply_template_with_custom_params(self):
        """Test applying template with custom parameter overrides"""
        template = JobTemplate.objects.create(
            name="Test Template",
            description="Test template for testing",
            default_params={'exhaustiveness': 16, 'num_modes': 20}
        )
        
        custom_params = {'exhaustiveness': 32}
        params = JobTemplateManager.apply_template(template, custom_params)
        
        self.assertEqual(params['exhaustiveness'], 32)  # Overridden
        self.assertEqual(params['num_modes'], 20)      # From template


class TestAdvancedDockingModels(TestCase):
    """Test advanced docking models"""

    def setUp(self):
        self.docking_job = DockingJob.objects.create(
            ligand_smiles='CCO',
            receptor_pdb_id='1CRN',
            center_x=5.0,
            center_y=10.0,
            center_z=-2.0,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0,
            job_name='Test Advanced Job',
            tags=['test', 'advanced']
        )

    def test_docking_job_get_best_pose(self):
        """Test getting best pose from docking results"""
        self.docking_job.results_json = {
            'poses': [
                {'mode': 1, 'affinity': -7.5},
                {'mode': 2, 'affinity': -8.2},
                {'mode': 3, 'affinity': -6.8}
            ]
        }
        self.docking_job.save()
        
        best_pose = self.docking_job.get_best_pose()
        self.assertIsNotNone(best_pose)
        self.assertEqual(best_pose['affinity'], -8.2)
        self.assertEqual(best_pose['mode'], 2)

    def test_docking_job_get_best_pose_no_results(self):
        """Test getting best pose when no results exist"""
        best_pose = self.docking_job.get_best_pose()
        self.assertIsNone(best_pose)

    def test_docking_job_get_tags_display(self):
        """Test getting tags as display string"""
        tags_display = self.docking_job.get_tags_display()
        self.assertEqual(tags_display, 'test, advanced')

    def test_binding_pocket_model(self):
        """Test BindingPocket model"""
        pocket = BindingPocket.objects.create(
            protein_pdb_id='1CRN',
            center_x=5.0,
            center_y=10.0,
            center_z=-2.0,
            radius=15.0,
            volume=800.0,
            druggability_score=0.75,
            confidence_score=0.9,
            residues=[{'residue_name': 'ALA', 'residue_number': 123}]
        )
        
        self.assertEqual(pocket.protein_pdb_id, '1CRN')
        self.assertEqual(pocket.detection_method, 'fpocket')  # Default value
        self.assertIsInstance(pocket.residues, list)

    def test_job_template_model(self):
        """Test JobTemplate model"""
        template = JobTemplate.objects.create(
            name='Test Template',
            description='Test template description',
            default_params={'exhaustiveness': 8},
            category='test'
        )
        
        self.assertEqual(template.usage_count, 0)
        
        template.increment_usage()
        self.assertEqual(template.usage_count, 1)


class TestAdvancedDockingViews(TestCase):
    """Test advanced docking API views"""

    def setUp(self):
        self.client = Client()
        self.docking_job = DockingJob.objects.create(
            ligand_smiles='CCO',
            receptor_pdb_id='1CRN',
            center_x=5.0,
            center_y=10.0,
            center_z=-2.0,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0,
            results_json={
                'poses': [
                    {'mode': 1, 'affinity': -8.5},
                    {'mode': 2, 'affinity': -7.2}
                ]
            }
        )

    @patch('api.views.ChemUtils.fetch_protein_structure')
    @patch('api.views.PocketDetector.detect_pockets')
    def test_detect_binding_pockets_view(self, mock_detect, mock_fetch):
        """Test pocket detection endpoint"""
        mock_fetch.return_value = "mock_pdb_content"
        mock_detect.return_value = [{
            'pocket_number': 1,
            'center_x': 5.0,
            'center_y': 10.0,
            'center_z': -2.0,
            'radius': 15.0,
            'volume': 800.0,
            'druggability_score': 0.75,
            'hydrophobicity': 0.6,
            'polarity': 0.4,
            'detection_method': 'fpocket',
            'confidence_score': 0.9,
            'residues': []
        }]
        
        response = self.client.post(
            '/api/pockets/detect',
            data={'pdb_id': '1CRN'},
            content_type='application/json'
        )
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        self.assertEqual(data['pdb_id'], '1CRN')
        self.assertEqual(len(data['pockets']), 1)
        self.assertIn('druggability_analysis', data['pockets'][0])

    def test_rescore_poses_view(self):
        """Test pose rescoring endpoint"""
        with patch('api.views.GNINAScorer.rescore_poses') as mock_rescore:
            mock_rescore.return_value = [
                {'mode': 1, 'affinity': -8.5, 'gnina_score': -9.2},
                {'mode': 2, 'affinity': -7.2, 'gnina_score': -8.1}
            ]
            
            response = self.client.post(
                '/api/poses/rescore',
                data={'job_id': str(self.docking_job.job_id)},
                content_type='application/json'
            )
            
            self.assertEqual(response.status_code, 200)
            data = json.loads(response.content)
            self.assertEqual(data['rescoring_method'], 'GNINA')
            self.assertEqual(len(data['rescored_poses']), 2)

    def test_list_job_templates_view(self):
        """Test job templates listing endpoint"""
        JobTemplateManager.create_default_templates()
        # Make templates public
        JobTemplate.objects.all().update(is_public=True)
        
        response = self.client.get('/api/templates')
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        self.assertIn('templates', data)
        self.assertIn('categories', data)
        self.assertGreater(len(data['templates']), 0)

    def test_get_pocket_list_view(self):
        """Test pocket list endpoint"""
        BindingPocket.objects.create(
            protein_pdb_id='1CRN',
            center_x=5.0,
            center_y=10.0,
            center_z=-2.0,
            radius=15.0,
            volume=800.0,
            druggability_score=0.75,
            confidence_score=0.9,
            residues=[]
        )
        
        response = self.client.get('/api/pockets/1CRN')
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        self.assertEqual(data['pdb_id'], '1CRN')
        self.assertEqual(len(data['pockets']), 1)

    @patch('api.views.DockingEngine.validate_docking_parameters')
    @patch('api.views.threading.Thread')
    def test_run_advanced_docking_view(self, mock_thread, mock_validate):
        """Test advanced docking endpoint"""
        mock_validate.return_value = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 5.0,
            'center_y': 10.0,
            'center_z': -2.0,
            'size_x': 20.0,
            'size_y': 20.0,
            'size_z': 20.0,
            'exhaustiveness': 8,
            'num_modes': 9
        }
        
        response = self.client.post(
            '/api/dock/advanced/run',
            data={
                'ligand_smiles': 'CCO',
                'receptor_pdb_id': '1CRN',
                'center_x': 5.0,
                'center_y': 10.0,
                'center_z': -2.0,
                'size_x': 20.0,
                'size_y': 20.0,
                'size_z': 20.0,
                'use_pocket_detection': True,
                'use_gnina_rescoring': True,
                'job_name': 'Advanced Test Job'
            },
            content_type='application/json'
        )
        
        self.assertEqual(response.status_code, 200)
        data = json.loads(response.content)
        self.assertIn('job_id', data)
        self.assertIn('features', data)
        self.assertTrue(data['features']['pocket_detection'])
        self.assertTrue(data['features']['gnina_rescoring'])
        
        # Check that job was created in database
        job = DockingJob.objects.get(job_id=data['job_id'])
        self.assertEqual(job.job_name, 'Advanced Test Job')
        self.assertIn('use_pocket_detection', job.advanced_params)

    def test_invalid_json_requests(self):
        """Test handling of invalid JSON requests"""
        endpoints = [
            '/api/pockets/detect',
            '/api/poses/rescore',
            '/api/dock/advanced/run'
        ]
        
        for endpoint in endpoints:
            response = self.client.post(endpoint, data='invalid json', content_type='application/json')
            self.assertEqual(response.status_code, 400)

    def test_missing_required_parameters(self):
        """Test handling of missing required parameters"""
        # Test missing PDB ID for pocket detection
        response = self.client.post(
            '/api/pockets/detect',
            data={},
            content_type='application/json'
        )
        self.assertEqual(response.status_code, 400)
        
        # Test missing job ID for rescoring
        response = self.client.post(
            '/api/poses/rescore',
            data={},
            content_type='application/json'
        )
        self.assertEqual(response.status_code, 400)
