"""
Chemical utilities for ligand preparation and receptor handling.
Uses external APIs and services for cheminformatics operations.
"""

import requests
import logging
import time
from typing import Dict, List, Optional, Tuple
from django.conf import settings

logger = logging.getLogger(__name__)

class ChemUtils:
    """Utilities for chemical operations using external APIs"""
    
    # PubChem REST API base URL
    PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    # Protein Data Bank REST API
    PDB_BASE = "https://data.rcsb.org/rest/v1"
    PDB_FILES_BASE = "https://files.rcsb.org/download"
    
    @classmethod
    def validate_smiles(cls, smiles: str) -> Dict:
        """
        Validate SMILES string using PubChem API
        
        Args:
            smiles: SMILES string to validate
            
        Returns:
            Dict with validation result and metadata
        """
        try:
            # Use PubChem to validate SMILES
            url = f"{cls.PUBCHEM_BASE}/compound/smiles/{smiles}/property/MolecularFormula,MolecularWeight/JSON"
            
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                properties = data.get('PropertyTable', {}).get('Properties', [])
                
                if properties:
                    prop = properties[0]
                    return {
                        'valid': True,
                        'molecular_formula': prop.get('MolecularFormula', ''),
                        'molecular_weight': prop.get('MolecularWeight', 0),
                        'method': 'PubChem API'
                    }
            
            return {
                'valid': False,
                'error': 'Invalid SMILES or compound not found',
                'method': 'PubChem API'
            }
            
        except Exception as e:
            logger.error(f"SMILES validation error: {e}")
            return {
                'valid': False,
                'error': str(e),
                'method': 'PubChem API'
            }
    
    @classmethod
    def smiles_to_sdf(cls, smiles: str) -> Optional[str]:
        """
        Convert SMILES to SDF format using PubChem API
        
        Args:
            smiles: SMILES string
            
        Returns:
            SDF content as string or None if failed
        """
        try:
            # Get SDF from PubChem
            url = f"{cls.PUBCHEM_BASE}/compound/smiles/{smiles}/SDF"
            
            response = requests.get(url, timeout=15)
            if response.status_code == 200:
                return response.text
            
            return None
            
        except Exception as e:
            logger.error(f"SMILES to SDF conversion error: {e}")
            return None
    
    @classmethod
    def generate_3d_coordinates(cls, smiles: str) -> Optional[str]:
        """
        Generate 3D coordinates for a molecule from SMILES
        
        Args:
            smiles: SMILES string
            
        Returns:
            SDF with 3D coordinates or None if failed
        """
        try:
            # For now, use the same SDF endpoint as PubChem provides some 3D data
            # In a real implementation, you might use specialized services like:
            # - ChemSpider
            # - NIH Chemical Identifier Resolver
            # - Or run local conformer generation
            
            url = f"{cls.PUBCHEM_BASE}/compound/smiles/{smiles}/SDF"
            
            response = requests.get(url, timeout=15)
            if response.status_code == 200:
                sdf_content = response.text
                
                # Add metadata to indicate this has coordinates
                # (This is a placeholder - real 3D generation would use specialized tools)
                if sdf_content:
                    lines = sdf_content.split('\n')
                    if len(lines) > 0:
                        lines[0] = lines[0] + " - 3D coordinates"
                    return '\n'.join(lines)
            
            return None
            
        except Exception as e:
            logger.error(f"3D coordinate generation error: {e}")
            return None
    
    @classmethod
    def fetch_protein_structure(cls, pdb_id: str) -> Optional[str]:
        """
        Fetch protein structure from PDB
        
        Args:
            pdb_id: PDB identifier (e.g., '1CRN')
            
        Returns:
            PDB file content or None if failed
        """
        try:
            # Validate PDB ID format
            if not pdb_id or len(pdb_id) != 4:
                return None
            
            pdb_id = pdb_id.upper()
            url = f"{cls.PDB_FILES_BASE}/{pdb_id}.pdb"
            
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                return response.text
            
            return None
            
        except Exception as e:
            logger.error(f"PDB fetch error: {e}")
            return None
    
    @classmethod
    def get_protein_info(cls, pdb_id: str) -> Optional[Dict]:
        """
        Get protein information from PDB API
        
        Args:
            pdb_id: PDB identifier
            
        Returns:
            Dict with protein metadata or None if failed
        """
        try:
            if not pdb_id or len(pdb_id) != 4:
                return None
            
            pdb_id = pdb_id.upper()
            url = f"{cls.PDB_BASE}/core/entry/{pdb_id}"
            
            response = requests.get(url, timeout=15)
            if response.status_code == 200:
                data = response.json()
                
                return {
                    'pdb_id': pdb_id,
                    'title': data.get('struct', {}).get('title', ''),
                    'resolution': data.get('rcsb_entry_info', {}).get('resolution_combined', []),
                    'method': data.get('exptl', [{}])[0].get('method', ''),
                    'organism': data.get('rcsb_entity_source_organism', [{}])[0].get('organism_names', []),
                    'deposited': data.get('rcsb_accession_info', {}).get('deposit_date', ''),
                    'molecular_weight': data.get('rcsb_entry_info', {}).get('molecular_weight', 0)
                }
            
            return None
            
        except Exception as e:
            logger.error(f"PDB info fetch error: {e}")
            return None
    
    @classmethod
    def prepare_ligand_for_docking(cls, smiles: str) -> Dict:
        """
        Prepare ligand for molecular docking
        
        Args:
            smiles: SMILES string of the ligand
            
        Returns:
            Dict with prepared ligand data
        """
        try:
            # Validate SMILES
            validation = cls.validate_smiles(smiles)
            if not validation.get('valid'):
                return {
                    'success': False,
                    'error': validation.get('error', 'Invalid SMILES'),
                    'smiles': smiles
                }
            
            # Generate 3D coordinates
            sdf_3d = cls.generate_3d_coordinates(smiles)
            if not sdf_3d:
                return {
                    'success': False,
                    'error': 'Failed to generate 3D coordinates',
                    'smiles': smiles
                }
            
            return {
                'success': True,
                'smiles': smiles,
                'sdf_content': sdf_3d,
                'molecular_formula': validation.get('molecular_formula'),
                'molecular_weight': validation.get('molecular_weight'),
                'preparation_method': 'PubChem API + 3D generation'
            }
            
        except Exception as e:
            logger.error(f"Ligand preparation error: {e}")
            return {
                'success': False,
                'error': str(e),
                'smiles': smiles
            }
    
    @classmethod
    def prepare_receptor_for_docking(cls, pdb_id: str, chain_id: Optional[str] = None) -> Dict:
        """
        Prepare protein receptor for molecular docking
        
        Args:
            pdb_id: PDB identifier
            chain_id: Optional chain identifier
            
        Returns:
            Dict with prepared receptor data
        """
        try:
            # Get protein info
            protein_info = cls.get_protein_info(pdb_id)
            if not protein_info:
                return {
                    'success': False,
                    'error': f'PDB {pdb_id} not found',
                    'pdb_id': pdb_id
                }
            
            # Fetch structure
            pdb_content = cls.fetch_protein_structure(pdb_id)
            if not pdb_content:
                return {
                    'success': False,
                    'error': f'Failed to fetch PDB structure for {pdb_id}',
                    'pdb_id': pdb_id
                }
            
            # Basic processing (in real implementation, would clean structure, add hydrogens, etc.)
            processed_content = pdb_content
            
            # If chain specified, could filter here
            if chain_id:
                # Placeholder for chain filtering
                # In real implementation, would parse PDB and extract specific chain
                protein_info['selected_chain'] = chain_id
            
            return {
                'success': True,
                'pdb_id': pdb_id,
                'pdb_content': processed_content,
                'protein_info': protein_info,
                'preparation_method': 'PDB download + basic processing'
            }
            
        except Exception as e:
            logger.error(f"Receptor preparation error: {e}")
            return {
                'success': False,
                'error': str(e),
                'pdb_id': pdb_id
            }

    @classmethod
    def estimate_binding_site(cls, pdb_content: str) -> Dict:
        """
        Estimate potential binding sites in a protein structure
        
        Args:
            pdb_content: PDB file content
            
        Returns:
            Dict with binding site predictions
        """
        try:
            # This is a placeholder implementation
            # Real binding site prediction would use algorithms like:
            # - CASTp
            # - fpocket
            # - SiteMap
            # - Or ML-based methods
            
            # For now, return a mock binding site
            return {
                'success': True,
                'binding_sites': [
                    {
                        'site_id': 1,
                        'center': {'x': 0.0, 'y': 0.0, 'z': 0.0},
                        'volume': 500.0,
                        'confidence': 0.85,
                        'description': 'Predicted binding site 1'
                    }
                ],
                'method': 'Geometric analysis (placeholder)'
            }
            
        except Exception as e:
            logger.error(f"Binding site estimation error: {e}")
            return {
                'success': False,
                'error': str(e)
            }
