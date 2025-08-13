import json
import time
import random
import math
import logging
from typing import Dict, List, Any
from .chem_utils import ChemUtils

# Import real implementations
try:
    from .real_docking_engine import RealDockingUtils
    REAL_DOCKING_AVAILABLE = True
except ImportError:
    REAL_DOCKING_AVAILABLE = False

logger = logging.getLogger(__name__)

class DockingEngine:
    """Molecular docking engine with real and mock implementations"""
    
    @staticmethod
    def is_real_docking_available() -> bool:
        """Check if real docking software is available"""
        if REAL_DOCKING_AVAILABLE:
            from .real_docking_engine import RealDockingEngine
            return RealDockingEngine.is_available()
        return False
    
    @staticmethod
    def validate_docking_parameters(params: Dict[str, Any]) -> Dict[str, Any]:
        """Validate docking parameters using real or mock implementation"""
        if DockingEngine.is_real_docking_available():
            logger.info("Using real docking parameter validation")
            return RealDockingUtils.validate_docking_parameters(params)
        else:
            logger.warning("Real docking not available, using mock validation")
            return DockingEngine._validate_docking_parameters_mock(params)
    
    @staticmethod
    def _validate_docking_parameters_mock(params: Dict[str, Any]) -> Dict[str, Any]:
        """Mock validation for fallback"""
        validated = {}
        
        # Required parameters
        required_fields = ['ligand_smiles', 'receptor_pdb_id', 'center_x', 'center_y', 'center_z', 
                          'size_x', 'size_y', 'size_z']
        
        for field in required_fields:
            if field not in params:
                raise ValueError(f"Missing required parameter: {field}")
        
        # Validate SMILES
        ligand_smiles = params['ligand_smiles'].strip()
        if not ligand_smiles:
            raise ValueError("Ligand SMILES cannot be empty")
        validated['ligand_smiles'] = ligand_smiles
        
        # Validate PDB ID
        receptor_pdb_id = params['receptor_pdb_id'].strip().upper()
        if not receptor_pdb_id or len(receptor_pdb_id) != 4:
            raise ValueError("Receptor PDB ID must be 4 characters")
        validated['receptor_pdb_id'] = receptor_pdb_id
        
        # Validate coordinates and sizes
        for coord in ['center_x', 'center_y', 'center_z']:
            try:
                validated[coord] = float(params[coord])
            except (ValueError, TypeError):
                raise ValueError(f"Invalid {coord}: must be a number")
        
        for size in ['size_x', 'size_y', 'size_z']:
            try:
                size_val = float(params[size])
                if size_val <= 0:
                    raise ValueError(f"Invalid {size}: must be positive")
                validated[size] = size_val
            except (ValueError, TypeError):
                raise ValueError(f"Invalid {size}: must be a positive number")
        
        # Optional parameters with defaults
        validated['exhaustiveness'] = max(1, min(32, int(params.get('exhaustiveness', 8))))
        validated['num_modes'] = max(1, min(20, int(params.get('num_modes', 9))))
        
        # Seed parameter for reproducibility
        seed = params.get('seed')
        if seed is not None:
            try:
                seed_val = int(seed)
            except (ValueError, TypeError):
                raise ValueError("Seed must be an integer")
            
            if seed_val < 0:
                raise ValueError("Seed must be a non-negative integer")
            validated['seed'] = seed_val
        else:
            # Default seed behavior - can be None for auto-random or a fixed value
            validated['seed'] = None  # None means auto-random
        
        return validated
    
    @staticmethod
    def run_production_docking(validated_params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run production docking using ONLY real AutoDock Vina - NO MOCK ALLOWED
        """
        from django.conf import settings
        
        # Check for force real setting
        force_real = getattr(settings, 'DOCKING_FORCE_REAL', False)
        cuda_enabled = getattr(settings, 'DOCKING_CUDA_ENABLED', False)
        
        # Add CUDA preference to parameters
        validated_params['cuda_enabled'] = cuda_enabled
        
        if DockingEngine.is_real_docking_available():
            logger.info(f"Running real AutoDock Vina docking (CUDA: {cuda_enabled})")
            result = RealDockingUtils.run_production_docking(validated_params)
            # Add engine metadata
            result['engine'] = 'vina'
            result['is_mock'] = False
            result['cuda_enabled'] = cuda_enabled
            result['runtime'] = 'server-side'
            return result
        else:
            # Real docking required but not available
            logger.error("AutoDock Vina not available - REAL DOCKING REQUIRED")
            return {
                'success': False,
                'error': 'AutoDock Vina not available - real docking required',
                'engine': 'unavailable',
                'is_mock': False,
                'cuda_enabled': False,
                'runtime': 'server-side',
                'message': 'AutoDock Vina must be installed. Mock docking is disabled.',
                'installation_required': True
            }
    
    @staticmethod
    def prepare_docking_inputs(ligand_smiles: str, receptor_pdb_id: str) -> Dict[str, Any]:
        """Prepare ligand and receptor for docking"""
        
        # Prepare ligand using existing ChemUtils
        ligand_data = ChemUtils.prepare_ligand_for_docking(ligand_smiles)
        
        # Prepare receptor using existing ChemUtils
        receptor_data = ChemUtils.prepare_receptor_for_docking(receptor_pdb_id)
        
        return {
            'ligand': ligand_data,
            'receptor': receptor_data
        }
    
    @staticmethod
    def run_mock_docking(validated_params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Mock implementation of AutoDock Vina docking.
        In a real implementation, this would call actual docking software.
        """
        
        # Set random seed if provided for reproducible mock results
        seed = validated_params.get('seed')
        if seed is not None:
            random.seed(seed)
            logger.info(f"Mock docking using deterministic seed: {seed}")
        else:
            logger.info("Mock docking using random seed (nondeterministic)")
        
        # Simulate docking calculation time (2-5 seconds)
        calc_time = random.uniform(2.0, 5.0)
        time.sleep(calc_time)
        
        # Generate mock docking poses
        num_modes = validated_params['num_modes']
        poses = []
        
        center_x = validated_params['center_x']
        center_y = validated_params['center_y'] 
        center_z = validated_params['center_z']
        
        for i in range(num_modes):
            # Generate mock binding affinity (kcal/mol) - lower is better
            affinity = random.uniform(-12.0, -5.0)
            
            # Generate pose coordinates near the binding site center
            pose_x = center_x + random.uniform(-2.0, 2.0)
            pose_y = center_y + random.uniform(-2.0, 2.0)
            pose_z = center_z + random.uniform(-2.0, 2.0)
            
            # Generate RMSD values
            rmsd_lb = random.uniform(0.0, 2.0)  # RMSD lower bound
            rmsd_ub = rmsd_lb + random.uniform(0.0, 1.0)  # RMSD upper bound
            
            poses.append({
                'mode': i + 1,
                'affinity': round(affinity, 3),
                'rmsd_lb': round(rmsd_lb, 3),
                'rmsd_ub': round(rmsd_ub, 3),
                'center_x': round(pose_x, 3),
                'center_y': round(pose_y, 3),
                'center_z': round(pose_z, 3),
                # In real implementation, this would be actual molecular coordinates
                'coordinates': f"MOCK_POSE_{i+1}_COORDINATES",
                'sdf': DockingEngine._create_mock_sdf(i + 1, pose_x, pose_y, pose_z)
            })
        
        # Sort by affinity (best first)
        poses.sort(key=lambda x: x['affinity'])
        
        # Calculate some summary statistics
        best_affinity = poses[0]['affinity'] if poses else 0
        mean_affinity = sum(p['affinity'] for p in poses) / len(poses) if poses else 0
        
        results = {
            'success': True,
            'poses': poses,
            'summary': {
                'num_poses': len(poses),
                'best_affinity': round(best_affinity, 3),
                'mean_affinity': round(mean_affinity, 3),
                'calculation_time': round(calc_time, 2),
                'software': 'AutoDock Vina (Mock)',
                'version': '1.2.0 (simulated)'
            },
            'parameters': validated_params,
            'timestamp': time.time(),
            'seed_used': seed
        }
        
        return results
    
    @staticmethod
    def _create_mock_sdf(pose_number: int, center_x: float, center_y: float, center_z: float) -> str:
        """Create a mock SDF for testing purposes"""
        # Create a simple mock molecule (e.g., ethanol) positioned at the given coordinates
        return f"""Mock Docked Pose {pose_number}
  Generated by Mock AutoDock Vina

  8  7  0  0  0  0  0  0  0  0999 V2000
{center_x:10.4f}{center_y:10.4f}{center_z:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0
{center_x+1.0:10.4f}{center_y:10.4f}{center_z:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0
{center_x+2.0:10.4f}{center_y:10.4f}{center_z:10.4f} O   0  0  0  0  0  0  0  0  0  0  0  0
{center_x:10.4f}{center_y+1.0:10.4f}{center_z:10.4f} H   0  0  0  0  0  0  0  0  0  0  0  0
{center_x:10.4f}{center_y-1.0:10.4f}{center_z:10.4f} H   0  0  0  0  0  0  0  0  0  0  0  0
{center_x+1.0:10.4f}{center_y+1.0:10.4f}{center_z:10.4f} H   0  0  0  0  0  0  0  0  0  0  0  0
{center_x+1.0:10.4f}{center_y-1.0:10.4f}{center_z:10.4f} H   0  0  0  0  0  0  0  0  0  0  0  0
{center_x+2.0:10.4f}{center_y+1.0:10.4f}{center_z:10.4f} H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  2  6  1  0  0  0  0
  2  7  1  0  0  0  0
  3  8  1  0  0  0  0
M  END
$$$$"""
    
    @staticmethod
    def estimate_binding_site_auto(receptor_pdb_id: str, ligand_smiles: str = None) -> Dict[str, Any]:
        """
        Automatically estimate binding site for receptor.
        This is a simplified version - real implementation would use cavity detection.
        """
        try:
            # Get protein metadata to estimate size
            metadata = ChemUtils.get_protein_info(receptor_pdb_id)
            
            # For mock purposes, create a reasonable binding site
            # In reality, this would use algorithms like fpocket, CASTp, etc.
            
            if ligand_smiles:
                # If ligand provided, create a site that accommodates it
                center_x, center_y, center_z = 5.0, 5.0, 5.0
                size_x, size_y, size_z = 25.0, 25.0, 25.0
                method = f"Estimated based on ligand {ligand_smiles}"
            else:
                # Default site for the protein
                center_x, center_y, center_z = 0.0, 0.0, 0.0
                size_x, size_y, size_z = 20.0, 20.0, 20.0
                method = "Default geometric center"
            
            return {
                'center_x': center_x,
                'center_y': center_y,
                'center_z': center_z,
                'size_x': size_x,
                'size_y': size_y,
                'size_z': size_z,
                'confidence': 'medium',
                'method': method,
                'protein_info': {
                    'pdb_id': receptor_pdb_id,
                    'title': metadata.get('title', 'Unknown'),
                    'molecular_weight': metadata.get('molecular_weight', 0)
                }
            }
            
        except Exception as e:
            # Fallback to default values
            return {
                'center_x': 0.0,
                'center_y': 0.0,
                'center_z': 0.0,
                'size_x': 20.0,
                'size_y': 20.0,
                'size_z': 20.0,
                'confidence': 'low',
                'method': f"Default (error: {str(e)})",
                'protein_info': {
                    'pdb_id': receptor_pdb_id,
                    'title': 'Unknown',
                    'molecular_weight': 0
                }
            }
