"""
Tests for chemical utilities
"""

import pytest
from unittest.mock import patch, MagicMock
from api.chem_utils import ChemUtils


class TestChemUtils:
    """Test cases for ChemUtils class"""
    
    def test_validate_smiles_valid(self):
        """Test SMILES validation with valid input"""
        with patch('requests.get') as mock_get:
            # Mock PubChem API response
            mock_response = MagicMock()
            mock_response.status_code = 200
            mock_response.json.return_value = {
                'PropertyTable': {
                    'Properties': [{
                        'MolecularFormula': 'C2H6O',
                        'MolecularWeight': 46.07
                    }]
                }
            }
            mock_get.return_value = mock_response
            
            result = ChemUtils.validate_smiles('CCO')
            
            assert result['valid'] is True
            assert result['molecular_formula'] == 'C2H6O'
            assert result['molecular_weight'] == 46.07
            assert result['method'] == 'PubChem API'
    
    def test_validate_smiles_invalid(self):
        """Test SMILES validation with invalid input"""
        with patch('requests.get') as mock_get:
            # Mock PubChem API response for invalid SMILES
            mock_response = MagicMock()
            mock_response.status_code = 404
            mock_get.return_value = mock_response
            
            result = ChemUtils.validate_smiles('INVALID')
            
            assert result['valid'] is False
            assert 'error' in result
            assert result['method'] == 'PubChem API'
    
    def test_smiles_to_sdf_success(self):
        """Test SMILES to SDF conversion"""
        with patch('requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.status_code = 200
            mock_response.text = "SDF_CONTENT_HERE"
            mock_get.return_value = mock_response
            
            result = ChemUtils.smiles_to_sdf('CCO')
            
            assert result == "SDF_CONTENT_HERE"
    
    def test_smiles_to_sdf_failure(self):
        """Test SMILES to SDF conversion failure"""
        with patch('requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.status_code = 404
            mock_get.return_value = mock_response
            
            result = ChemUtils.smiles_to_sdf('INVALID')
            
            assert result is None
    
    def test_generate_3d_coordinates(self):
        """Test 3D coordinate generation"""
        with patch('requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.status_code = 200
            mock_response.text = "ORIGINAL_SDF"
            mock_get.return_value = mock_response
            
            result = ChemUtils.generate_3d_coordinates('CCO')
            
            assert result is not None
            assert "3D coordinates" in result
    
    def test_fetch_protein_structure_success(self):
        """Test successful protein structure fetching"""
        with patch('requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.status_code = 200
            mock_response.text = "PDB_CONTENT_HERE"
            mock_get.return_value = mock_response
            
            result = ChemUtils.fetch_protein_structure('1CRN')
            
            assert result == "PDB_CONTENT_HERE"
    
    def test_fetch_protein_structure_invalid_id(self):
        """Test protein structure fetching with invalid PDB ID"""
        result = ChemUtils.fetch_protein_structure('ABC')  # Too short
        assert result is None
        
        result = ChemUtils.fetch_protein_structure('')  # Empty
        assert result is None
    
    def test_get_protein_info_success(self):
        """Test successful protein info retrieval"""
        with patch('requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.status_code = 200
            mock_response.json.return_value = {
                'struct': {'title': 'Test Protein'},
                'rcsb_entry_info': {'resolution_combined': [2.0]},
                'exptl': [{'method': 'X-RAY DIFFRACTION'}],
                'rcsb_entity_source_organism': [{'organism_names': ['E. coli']}],
                'rcsb_accession_info': {'deposit_date': '2023-01-01'},
                'rcsb_entry_info': {'molecular_weight': 12345}
            }
            mock_get.return_value = mock_response
            
            result = ChemUtils.get_protein_info('1CRN')
            
            assert result is not None
            assert result['pdb_id'] == '1CRN'
            assert result['title'] == 'Test Protein'
            assert result['method'] == 'X-RAY DIFFRACTION'
    
    def test_prepare_ligand_for_docking_success(self):
        """Test successful ligand preparation"""
        with patch.object(ChemUtils, 'validate_smiles') as mock_validate, \
             patch.object(ChemUtils, 'generate_3d_coordinates') as mock_3d:
            
            mock_validate.return_value = {
                'valid': True,
                'molecular_formula': 'C2H6O',
                'molecular_weight': 46.07
            }
            mock_3d.return_value = "SDF_WITH_3D_COORDS"
            
            result = ChemUtils.prepare_ligand_for_docking('CCO')
            
            assert result['success'] is True
            assert result['smiles'] == 'CCO'
            assert result['molecular_formula'] == 'C2H6O'
            assert result['sdf_content'] == "SDF_WITH_3D_COORDS"
    
    def test_prepare_ligand_for_docking_invalid_smiles(self):
        """Test ligand preparation with invalid SMILES"""
        with patch.object(ChemUtils, 'validate_smiles') as mock_validate:
            mock_validate.return_value = {
                'valid': False,
                'error': 'Invalid SMILES'
            }
            
            result = ChemUtils.prepare_ligand_for_docking('INVALID')
            
            assert result['success'] is False
            assert 'error' in result
    
    def test_prepare_receptor_for_docking_success(self):
        """Test successful receptor preparation"""
        with patch.object(ChemUtils, 'get_protein_info') as mock_info, \
             patch.object(ChemUtils, 'fetch_protein_structure') as mock_fetch:
            
            mock_info.return_value = {
                'pdb_id': '1CRN',
                'title': 'Test Protein'
            }
            mock_fetch.return_value = "PDB_CONTENT"
            
            result = ChemUtils.prepare_receptor_for_docking('1CRN')
            
            assert result['success'] is True
            assert result['pdb_id'] == '1CRN'
            assert result['pdb_content'] == "PDB_CONTENT"
    
    def test_prepare_receptor_for_docking_not_found(self):
        """Test receptor preparation with non-existent PDB ID"""
        with patch.object(ChemUtils, 'get_protein_info') as mock_info:
            mock_info.return_value = None
            
            result = ChemUtils.prepare_receptor_for_docking('XXXX')
            
            assert result['success'] is False
            assert 'not found' in result['error']
    
    def test_estimate_binding_site(self):
        """Test binding site estimation"""
        result = ChemUtils.estimate_binding_site("MOCK_PDB_CONTENT")
        
        assert result['success'] is True
        assert 'binding_sites' in result
        assert len(result['binding_sites']) > 0
        assert 'center' in result['binding_sites'][0]
        assert 'volume' in result['binding_sites'][0]
