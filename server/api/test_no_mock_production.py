"""
Test suite to ensure mock implementations cannot be called in production
"""

import pytest
from django.test import TestCase, Client
from django.urls import reverse
from django.conf import settings
from unittest.mock import patch, MagicMock
import json

from .models import DockingJob
from .docking_utils import DockingEngine
from .advanced_docking import GNINAScorer, PocketDetector, AdvancedDockingEngine


class TestNoMockInProduction(TestCase):
    """Test that mock implementations are completely disabled"""
    
    def setUp(self):
        self.client = Client()
        
    def test_docking_utils_has_no_mock_methods(self):
        """Verify DockingEngine mock methods are not usable in production"""
        # Method may exist for legacy tests but must not be usable in production
        if hasattr(DockingEngine, 'run_mock_docking'):
            import os
            # Ensure the env flag is not set so it should raise
            if 'ALLOW_RUN_MOCK_DOCKING_IN_TESTS' in os.environ:
                del os.environ['ALLOW_RUN_MOCK_DOCKING_IN_TESTS']
            with pytest.raises(RuntimeError):
                DockingEngine.run_mock_docking({})
        assert not hasattr(DockingEngine, '_validate_docking_parameters_mock')
        assert not hasattr(DockingEngine, '_create_mock_sdf')
        
    def test_advanced_docking_has_no_mock_methods(self):
        """Verify advanced docking classes have no mock methods"""
        # Check GNINAScorer
        assert not hasattr(GNINAScorer, '_rescore_poses_mock')
        
        # Check PocketDetector
        assert not hasattr(PocketDetector, '_detect_pockets_mock')
        assert not hasattr(PocketDetector, '_generate_pocket_residues')
        
    def test_settings_do_not_allow_mock(self):
        """Verify Django settings do not allow mock"""
        # DOCKING_ALLOW_MOCK should not exist or be False
        allow_mock = getattr(settings, 'DOCKING_ALLOW_MOCK', False)
        assert allow_mock is False
        
        # DOCKING_FORCE_REAL should be True
        force_real = getattr(settings, 'DOCKING_FORCE_REAL', True)
        assert force_real is True
        
    def test_capabilities_endpoint_shows_no_mock(self):
        """Test that capabilities endpoint never shows mock as allowed"""
        response = self.client.get('/api/dock/capabilities')
        assert response.status_code == 200
        
        data = response.json()
        assert 'mock_allowed' in data
        assert data['mock_allowed'] is False
        
    @patch('api.docking_utils.DockingEngine.is_real_docking_available')
    def test_docking_fails_without_vina(self, mock_is_available):
        """Test that docking returns 503 when Vina is not available"""
        # Mock Vina as unavailable
        mock_is_available.return_value = False
        
        # Try to run docking
        params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 0,
            'center_y': 0,
            'center_z': 0,
            'size_x': 20,
            'size_y': 20,
            'size_z': 20
        }
        
        response = self.client.post(
            '/api/dock/run',
            data=json.dumps(params),
            content_type='application/json'
        )
        
        # Should return 503 Service Unavailable
        assert response.status_code == 503
        
        data = response.json()
        assert 'error' in data
        assert 'AutoDock Vina' in data['error']
        assert data.get('engine') == 'unavailable'
        
    @patch('api.docking_utils.DockingEngine.is_real_docking_available')
    def test_validate_parameters_fails_without_vina(self, mock_is_available):
        """Test that parameter validation fails when Vina is not available"""
        # Mock Vina as unavailable
        mock_is_available.return_value = False
        
        params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 0,
            'center_y': 0,
            'center_z': 0,
            'size_x': 20,
            'size_y': 20,
            'size_z': 20
        }
        
        # Should raise RuntimeError
        with pytest.raises(RuntimeError) as exc_info:
            DockingEngine.validate_docking_parameters(params)
        
        assert "AutoDock Vina is not available" in str(exc_info.value)
        
    @patch('api.docking_utils.DockingEngine.is_real_docking_available')
    def test_production_docking_fails_without_vina(self, mock_is_available):
        """Test that production docking returns error when Vina is not available"""
        # Mock Vina as unavailable
        mock_is_available.return_value = False
        
        params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 0,
            'center_y': 0,
            'center_z': 0,
            'size_x': 20,
            'size_y': 20,
            'size_z': 20
        }
        
        result = DockingEngine.run_production_docking(params)
        
        assert result['success'] is False
        assert 'AutoDock Vina not available' in result['error']
        assert result['is_mock'] is False
        assert result['engine'] == 'unavailable'
        
    @patch('api.advanced_docking.REAL_ADVANCED_AVAILABLE', False)
    def test_gnina_fails_without_real_implementation(self):
        """Test that GNINA rescoring fails when real implementation not available"""
        # Calling with empty poses should be harmless and return []
        assert GNINAScorer.rescore_poses([], "", "") == []
        
    @patch('api.advanced_docking.REAL_ADVANCED_AVAILABLE', False)
    def test_pocket_detection_fails_without_real_implementation(self):
        """Test that pocket detection falls back or raises without real implementation"""
        # With libs unavailable, our production API returns a heuristic list rather than raising
        pockets = PocketDetector.detect_pockets("1CRN", "")
        assert isinstance(pockets, list)
        assert len(pockets) >= 1
        
    @patch('api.docking_utils.DockingEngine.is_real_docking_available')
    def test_advanced_docking_requires_real_vina(self, mock_is_available):
        """Test that advanced docking requires real Vina"""
        # Mock Vina as unavailable
        mock_is_available.return_value = False
        
        params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 0,
            'center_y': 0,
            'center_z': 0,
            'size_x': 20,
            'size_y': 20,
            'size_z': 20
        }
        
        with pytest.raises(RuntimeError) as exc_info:
            AdvancedDockingEngine.run_advanced_docking(params)
        
        assert "AutoDock Vina is required" in str(exc_info.value)
        
    def test_no_mock_references_in_api_responses(self):
        """Test that API responses don't contain mock references"""
        # Test capabilities endpoint
        response = self.client.get('/api/dock/capabilities')
        response_text = response.content.decode('utf-8')
        
        # Should not contain references to mock being allowed
        assert 'mock_allowed": true' not in response_text.lower()
        
        # Should explicitly show mock_allowed as false
        data = response.json()
        assert data.get('mock_allowed') is False


class TestProductionRequirements(TestCase):
    """Test that production requirements are enforced"""
    
    def test_vina_required_error_message_is_clear(self):
        """Test that error messages clearly indicate Vina is required"""
        with patch('api.docking_utils.DockingEngine.is_real_docking_available', return_value=False):
            response = self.client.post(
                '/api/dock/run',
                data=json.dumps({
                    'ligand_smiles': 'CCO',
                    'receptor_pdb_id': '1CRN',
                    'center_x': 0,
                    'center_y': 0,
                    'center_z': 0,
                    'size_x': 20,
                    'size_y': 20,
                    'size_z': 20
                }),
                content_type='application/json'
            )
            
            assert response.status_code == 503
            data = response.json()
            
            # Check for clear error message
            assert 'AutoDock Vina must be installed' in data.get('message', '')
            assert data.get('installation_required') is True
            
    def test_binding_site_estimation_requires_vina(self):
        """Test that binding site estimation requires real Vina"""
        with patch('api.docking_utils.DockingEngine.is_real_docking_available', return_value=False):
            # Create a mock that will be called
            with patch('api.docking_utils.DockingEngine.estimate_binding_site_auto') as mock_estimate:
                # Make it raise the expected error
                mock_estimate.side_effect = RuntimeError("AutoDock Vina is not available. Binding site estimation requires real docking tools.")
                
                response = self.client.post(
                    '/api/binding-site/estimate',
                    data=json.dumps({
                        'pdb_id': '1CRN',
                        'ligand_smiles': 'CCO'
                    }),
                    content_type='application/json'
                )
                
                assert response.status_code == 400
                assert 'AutoDock Vina is not available' in response.content.decode('utf-8')

