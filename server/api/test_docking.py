import pytest
import json
import uuid
from django.urls import reverse
from django.test import TestCase
from django.utils import timezone
from unittest.mock import patch, MagicMock
from freezegun import freeze_time

from .models import DockingJob
from .docking_utils import DockingEngine


class TestDockingModels(TestCase):
    """Test DockingJob model"""
    
    def test_create_docking_job(self):
        """Test creating a docking job"""
        job = DockingJob.objects.create(
            ligand_smiles="CCO",
            receptor_pdb_id="1CRN",
            center_x=5.0,
            center_y=5.0,
            center_z=5.0,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0,
            exhaustiveness=8,
            num_modes=9
        )
        
        assert job.job_id is not None
        assert job.status == 'pending'
        assert job.ligand_smiles == "CCO"
        assert job.receptor_pdb_id == "1CRN"
        assert job.created_at is not None
        assert not job.is_finished()
        
    def test_job_duration(self):
        """Test job duration calculation"""
        with freeze_time("2023-01-01 10:00:00") as frozen_time:
            job = DockingJob.objects.create(
                ligand_smiles="CCO",
                receptor_pdb_id="1CRN",
                center_x=0, center_y=0, center_z=0,
                size_x=20, size_y=20, size_z=20,
                started_at=timezone.now()
            )
            
            # No duration if not completed
            assert job.duration() is None
            
            # Move time forward and complete job
            frozen_time.tick(delta=120)  # 2 minutes
            job.completed_at = timezone.now()
            job.save()
            
            duration = job.duration()
            assert duration == 120.0
            
    def test_job_status_transitions(self):
        """Test job status transitions"""
        job = DockingJob.objects.create(
            ligand_smiles="CCO",
            receptor_pdb_id="1CRN",
            center_x=0, center_y=0, center_z=0,
            size_x=20, size_y=20, size_z=20
        )
        
        # Initial state
        assert not job.is_finished()
        
        # Running state
        job.status = 'running'
        assert not job.is_finished()
        
        # Terminal states
        for status in ['completed', 'failed', 'cancelled']:
            job.status = status
            assert job.is_finished()


class TestDockingEngine(TestCase):
    """Test DockingEngine utility functions"""
    
    def test_validate_docking_parameters_success(self):
        """Test successful parameter validation"""
        params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 5.0,
            'center_y': 5.0,
            'center_z': 5.0,
            'size_x': 20.0,
            'size_y': 20.0,
            'size_z': 20.0,
            'exhaustiveness': 8,
            'num_modes': 9
        }
        
        validated = DockingEngine.validate_docking_parameters(params)
        
        assert validated['ligand_smiles'] == 'CCO'
        assert validated['receptor_pdb_id'] == '1CRN'
        assert validated['center_x'] == 5.0
        assert validated['exhaustiveness'] == 8
        assert validated['num_modes'] == 9
        
    def test_validate_docking_parameters_missing_required(self):
        """Test validation with missing required parameters"""
        params = {
            'ligand_smiles': 'CCO'
            # Missing other required fields
        }
        
        with pytest.raises(ValueError, match="Missing required parameter"):
            DockingEngine.validate_docking_parameters(params)
            
    def test_validate_docking_parameters_invalid_smiles(self):
        """Test validation with empty SMILES"""
        params = {
            'ligand_smiles': '',
            'receptor_pdb_id': '1CRN',
            'center_x': 5.0, 'center_y': 5.0, 'center_z': 5.0,
            'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0
        }
        
        with pytest.raises(ValueError, match="Ligand SMILES cannot be empty"):
            DockingEngine.validate_docking_parameters(params)
            
    def test_validate_docking_parameters_invalid_pdb_id(self):
        """Test validation with invalid PDB ID"""
        params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': 'INVALID_PDB',
            'center_x': 5.0, 'center_y': 5.0, 'center_z': 5.0,
            'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0
        }
        
        with pytest.raises(ValueError, match="Receptor PDB ID must be 4 characters"):
            DockingEngine.validate_docking_parameters(params)
            
    def test_validate_docking_parameters_invalid_coordinates(self):
        """Test validation with invalid coordinates"""
        params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 'invalid',
            'center_y': 5.0, 'center_z': 5.0,
            'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0
        }
        
        with pytest.raises(ValueError, match="Invalid center_x: must be a number"):
            DockingEngine.validate_docking_parameters(params)
            
    def test_validate_docking_parameters_negative_size(self):
        """Test validation with negative box size"""
        params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 5.0, 'center_y': 5.0, 'center_z': 5.0,
            'size_x': -10.0, 'size_y': 20.0, 'size_z': 20.0
        }
        
        with pytest.raises(ValueError, match="Invalid size_x: must be"):
            DockingEngine.validate_docking_parameters(params)
    
    @patch('api.docking_utils.ChemUtils.prepare_ligand_for_docking')
    @patch('api.docking_utils.ChemUtils.prepare_receptor_for_docking')
    def test_prepare_docking_inputs(self, mock_prepare_receptor, mock_prepare_ligand):
        """Test docking input preparation"""
        mock_prepare_ligand.return_value = {'smiles': 'CCO', 'sdf_content': 'MOCK_SDF'}
        mock_prepare_receptor.return_value = {'pdb_id': '1CRN', 'pdb_content': 'MOCK_PDB'}
        
        result = DockingEngine.prepare_docking_inputs('CCO', '1CRN')
        
        assert 'ligand' in result
        assert 'receptor' in result
        assert result['ligand']['smiles'] == 'CCO'
        assert result['receptor']['pdb_id'] == '1CRN'
        
        mock_prepare_ligand.assert_called_once_with('CCO')
        mock_prepare_receptor.assert_called_once_with('1CRN')
    
    def test_run_mock_docking(self):
        """Test mock docking execution"""
        params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 5.0, 'center_y': 5.0, 'center_z': 5.0,
            'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0,
            'exhaustiveness': 8,
            'num_modes': 5
        }
        
        results = DockingEngine.run_mock_docking(params)
        
        assert results['success'] is True
        assert 'poses' in results
        assert 'summary' in results
        assert len(results['poses']) == 5
        
        # Check pose structure
        pose = results['poses'][0]
        assert 'mode' in pose
        assert 'affinity' in pose
        assert 'rmsd_lb' in pose
        assert 'rmsd_ub' in pose
        assert 'coordinates' in pose
        
        # Check summary
        summary = results['summary']
        assert summary['num_poses'] == 5
        assert 'best_affinity' in summary
        assert 'calculation_time' in summary
        
    @patch('api.docking_utils.ChemUtils.get_protein_info')
    def test_estimate_binding_site_auto(self, mock_get_info):
        """Test automatic binding site estimation"""
        mock_get_info.return_value = {
            'title': 'Test Protein',
            'molecular_weight': 10000
        }
        
        result = DockingEngine.estimate_binding_site_auto('1CRN')
        
        assert 'center_x' in result
        assert 'center_y' in result
        assert 'center_z' in result
        assert 'size_x' in result
        assert 'size_y' in result
        assert 'size_z' in result
        assert 'confidence' in result
        assert 'method' in result
        assert 'protein_info' in result
        
        mock_get_info.assert_called_once_with('1CRN')
        
    @patch('api.docking_utils.ChemUtils.get_protein_info')
    def test_estimate_binding_site_with_ligand(self, mock_get_info):
        """Test binding site estimation with ligand"""
        mock_get_info.return_value = {
            'title': 'Test Protein',
            'molecular_weight': 10000
        }
        
        result = DockingEngine.estimate_binding_site_auto('1CRN', 'CCO')
        
        assert result['center_x'] == 5.0  # Expected values for ligand-based estimation
        assert result['size_x'] == 25.0
        assert 'ligand CCO' in result['method']
        
    @patch('api.docking_utils.ChemUtils.get_protein_info')
    def test_estimate_binding_site_error_handling(self, mock_get_info):
        """Test binding site estimation error handling"""
        mock_get_info.side_effect = Exception("API Error")
        
        result = DockingEngine.estimate_binding_site_auto('INVALID')
        
        assert result['confidence'] == 'low'
        assert 'error' in result['method']


@pytest.mark.django_db
class TestDockingAPIViews:
    """Test docking API endpoints"""
    
    def test_estimate_binding_site_endpoint(self, client):
        """Test binding site estimation endpoint"""
        with patch('api.views.DockingEngine.estimate_binding_site_auto') as mock_estimate:
            mock_estimate.return_value = {
                'center_x': 5.0, 'center_y': 5.0, 'center_z': 5.0,
                'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0,
                'confidence': 'medium'
            }
            
            url = reverse('estimate-binding-site')
            payload = {'pdb_id': '1CRN', 'ligand_smiles': 'CCO'}
            response = client.post(url, data=json.dumps(payload), content_type='application/json')
            
            assert response.status_code == 200
            data = response.json()
            assert data['center_x'] == 5.0
            assert data['confidence'] == 'medium'
            mock_estimate.assert_called_once_with('1CRN', 'CCO')
    
    def test_estimate_binding_site_missing_pdb_id(self, client):
        """Test binding site estimation with missing PDB ID"""
        url = reverse('estimate-binding-site')
        response = client.post(url, data=json.dumps({}), content_type='application/json')
        
        assert response.status_code == 400
        assert "PDB ID is required" in response.content.decode()
    
    def test_run_docking_endpoint(self, client):
        """Test docking execution endpoint"""
        with patch('api.views.DockingEngine.validate_docking_parameters') as mock_validate, \
             patch('api.views.threading.Thread') as mock_thread:
            
            mock_validate.return_value = {
                'ligand_smiles': 'CCO', 'receptor_pdb_id': '1CRN',
                'center_x': 5.0, 'center_y': 5.0, 'center_z': 5.0,
                'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0,
                'exhaustiveness': 8, 'num_modes': 9
            }
            
            url = reverse('run-docking')
            payload = {
                'ligand_smiles': 'CCO', 'receptor_pdb_id': '1CRN',
                'center_x': 5.0, 'center_y': 5.0, 'center_z': 5.0,
                'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0
            }
            response = client.post(url, data=json.dumps(payload), content_type='application/json')
            
            assert response.status_code == 200
            data = response.json()
            assert 'job_id' in data
            assert data['status'] == 'pending'
            assert 'Docking job started' in data['message']
            
            # Verify job was created
            job_id = data['job_id']
            job = DockingJob.objects.get(job_id=job_id)
            assert job.ligand_smiles == 'CCO'
            assert job.receptor_pdb_id == '1CRN'
            
            # Verify background thread was started
            mock_thread.assert_called_once()
    
    def test_run_docking_invalid_parameters(self, client):
        """Test docking with invalid parameters"""
        with patch('api.views.DockingEngine.validate_docking_parameters') as mock_validate:
            mock_validate.side_effect = ValueError("Invalid SMILES")
            
            url = reverse('run-docking')
            payload = {'ligand_smiles': 'INVALID'}
            response = client.post(url, data=json.dumps(payload), content_type='application/json')
            
            assert response.status_code == 400
            assert "Invalid SMILES" in response.content.decode()
    
    def test_get_docking_status_endpoint(self, client):
        """Test getting docking job status"""
        # Create a test job
        job = DockingJob.objects.create(
            ligand_smiles="CCO",
            receptor_pdb_id="1CRN",
            center_x=5.0, center_y=5.0, center_z=5.0,
            size_x=20.0, size_y=20.0, size_z=20.0,
            status='completed',
            results_json={'poses': [], 'summary': {}}
        )
        
        url = reverse('get-docking-status', kwargs={'job_id': str(job.job_id)})
        response = client.get(url)
        
        assert response.status_code == 200
        data = response.json()
        assert data['job_id'] == str(job.job_id)
        assert data['status'] == 'completed'
        assert data['ligand_smiles'] == 'CCO'
        assert 'results' in data
    
    def test_get_docking_status_not_found(self, client):
        """Test getting status for non-existent job"""
        fake_job_id = str(uuid.uuid4())
        url = reverse('get-docking-status', kwargs={'job_id': fake_job_id})
        response = client.get(url)
        
        assert response.status_code == 400
        assert "not found" in response.content.decode()
    
    def test_list_docking_jobs_endpoint(self, client):
        """Test listing docking jobs"""
        # Create some test jobs
        for i in range(3):
            DockingJob.objects.create(
                ligand_smiles=f"SMILES_{i}",
                receptor_pdb_id="1CRN",
                center_x=5.0, center_y=5.0, center_z=5.0,
                size_x=20.0, size_y=20.0, size_z=20.0
            )
        
        url = reverse('list-docking-jobs')
        response = client.get(url)
        
        assert response.status_code == 200
        data = response.json()
        assert 'jobs' in data
        assert data['total_count'] == 3
        assert len(data['jobs']) == 3
        
        # Check job structure
        job = data['jobs'][0]
        assert 'job_id' in job
        assert 'status' in job
        assert 'ligand_smiles' in job
        assert 'receptor_pdb_id' in job
        assert 'created_at' in job
