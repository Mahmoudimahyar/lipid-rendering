import pytest
from django.urls import reverse
from unittest.mock import patch, MagicMock
import json


@pytest.mark.django_db
def test_healthz(client):
    """Test health check endpoint"""
    url = reverse('healthz')
    response = client.get(url)
    assert response.status_code == 200
    assert response.json() == {"ok": True}


@pytest.mark.django_db
@patch('api.chem_utils.ChemUtils.get_protein_info')
@patch('api.chem_utils.ChemUtils.fetch_protein_structure')
@patch('api.chem_utils.ChemUtils.estimate_binding_site')
def test_pdb_proxy_success(mock_binding, mock_fetch, mock_info, client):
    """Test successful PDB proxy request"""
    mock_info.return_value = {
        'title': 'Test Protein',
        'organism': ['E. coli'],
        'method': 'X-RAY DIFFRACTION',
        'resolution': [2.0],
        'molecular_weight': 12345,
        'deposited': '2023-01-01'
    }
    mock_fetch.return_value = "PDB_CONTENT"
    mock_binding.return_value = {
        'binding_sites': [{'site_id': 1, 'center': {'x': 0, 'y': 0, 'z': 0}}]
    }
    
    url = reverse('pdb-proxy', kwargs={'pdb_id': '1CRN'})
    response = client.get(url)
    assert response.status_code == 200
    
    data = response.json()
    assert data['pdb_id'] == '1CRN'
    assert data['title'] == 'Test Protein'
    assert 'binding_sites' in data


@pytest.mark.django_db
def test_pdb_proxy_invalid_id(client):
    """Test PDB proxy with invalid PDB ID"""
    url = reverse('pdb-proxy', kwargs={'pdb_id': 'ABC'})
    response = client.get(url)
    assert response.status_code == 400
    
    data = response.json()
    assert 'error' in data


@pytest.mark.django_db
@patch('api.chem_utils.ChemUtils.prepare_ligand_for_docking')
def test_prepare_ligand_success(mock_prepare, client):
    """Test successful ligand preparation"""
    mock_prepare.return_value = {
        'success': True,
        'smiles': 'CCO',
        'molecular_formula': 'C2H6O',
        'molecular_weight': 46.07,
        'sdf_content': 'SDF_CONTENT',
        'preparation_method': 'ChemUtils'
    }
    
    url = reverse('prepare-ligand')
    payload = {'smiles': 'CCO'}
    response = client.post(url, data=json.dumps(payload), content_type='application/json')
    assert response.status_code == 200
    
    data = response.json()
    assert data['smiles'] == 'CCO'
    assert data['molecular_formula'] == 'C2H6O'
    assert data['preparation_status'] == 'success'


@pytest.mark.django_db
def test_prepare_ligand_missing_smiles(client):
    """Test ligand preparation with missing SMILES"""
    url = reverse('prepare-ligand')
    response = client.post(url, data=json.dumps({}), content_type='application/json')
    assert response.status_code == 400
    
    data = response.json()
    assert 'error' in data


@pytest.mark.django_db
@patch('api.chem_utils.ChemUtils.prepare_ligand_for_docking')
def test_prepare_ligand_failure(mock_prepare, client):
    """Test ligand preparation failure"""
    mock_prepare.return_value = {
        'success': False,
        'error': 'Invalid SMILES'
    }
    
    url = reverse('prepare-ligand')
    payload = {'smiles': 'INVALID'}
    response = client.post(url, data=json.dumps(payload), content_type='application/json')
    assert response.status_code == 400
    
    data = response.json()
    assert 'error' in data


@pytest.mark.django_db
@patch('api.chem_utils.ChemUtils.prepare_receptor_for_docking')
def test_prepare_receptor_success(mock_prepare, client):
    """Test successful receptor preparation"""
    mock_prepare.return_value = {
        'success': True,
        'pdb_id': '1CRN',
        'protein_info': {'title': 'Test Protein'},
        'pdb_content': 'PDB_CONTENT',
        'preparation_method': 'ChemUtils'
    }
    
    url = reverse('prepare-receptor')
    payload = {'pdb_id': '1CRN'}
    response = client.post(url, data=json.dumps(payload), content_type='application/json')
    assert response.status_code == 200
    
    data = response.json()
    assert data['pdb_id'] == '1CRN'
    assert data['preparation_status'] == 'success'


@pytest.mark.django_db
def test_prepare_receptor_invalid_id(client):
    """Test receptor preparation with invalid PDB ID"""
    url = reverse('prepare-receptor')
    payload = {'pdb_id': 'ABC'}
    response = client.post(url, data=json.dumps(payload), content_type='application/json')
    assert response.status_code == 400
    
    data = response.json()
    assert 'error' in data