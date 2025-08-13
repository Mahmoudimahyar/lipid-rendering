# Reproducibility in Molecular Docking

This document explains how Lipid Rendering ensures reproducible molecular docking results through seeded calculations, deterministic algorithms, and comprehensive parameter logging.

## Overview

Reproducibility is critical for scientific computing and regression testing. Lipid Rendering implements multiple layers of determinism to ensure identical results across runs, environments, and time.

## Seeding Strategy

### Random Seed Propagation

The system uses a single random seed that propagates through the entire docking pipeline:

```
User Input Seed → Ligand Preparation → Receptor Preparation → Docking Calculation → Result Generation
```

### Seed Sources

1. **Explicit User Seed**: Provided via API parameter `seed`
2. **Default System Seed**: `42` used when no seed specified
3. **CI/Testing Seed**: Fixed seed for smoke tests and benchmarks

```json
{
  "ligand_smiles": "CCO",
  "receptor_pdb_id": "1CRN",
  "seed": 42,
  "exhaustiveness": 8,
  "num_modes": 9
}
```

## Deterministic Pipeline Components

### 1. Ligand Preparation

**SMILES → 3D Coordinates → PDBQT**

```python
def prepare_ligand(self, smiles: str, random_seed: int = 42):
    """
    Deterministic ligand preparation pipeline
    
    Steps:
    1. Parse SMILES to canonical form
    2. Add hydrogens at pH 7.4 (deterministic protonation)
    3. Generate 3D coordinates with ETKDG (seeded)
    4. Optimize geometry with MMFF94
    5. Convert to PDBQT format
    """
    # Set RDKit random seed for ETKDG
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    
    # Deterministic 3D coordinate generation
    AllChem.EmbedMolecule(mol, params)
```

**Deterministic Elements:**
- **Canonical SMILES**: Standardized molecular representation
- **Protonation State**: Fixed pH 7.4 for consistent ionization
- **3D Embedding**: ETKDG algorithm with explicit random seed
- **Force Field**: MMFF94 optimization (deterministic for same starting geometry)
- **PDBQT Conversion**: Reproducible format conversion

### 2. Receptor Preparation

**PDB → Cleaned Structure → PDBQT**

```python
def prepare_receptor(self, pdb_id: str):
    """
    Deterministic receptor preparation
    
    Policy:
    - First model only (if multiple models)
    - Highest occupancy alternative locations
    - Remove HETATM except essential cofactors
    - Add hydrogens at pH 7.4
    - Standardized PDBQT conversion
    """
```

**Deterministic Elements:**
- **Model Selection**: Always use first model from PDB
- **Occupancy Selection**: Highest occupancy for alternative conformations  
- **Cofactor Policy**: Consistent retention of essential cofactors (FAD, NAD, ATP, etc.)
- **Protonation**: Standardized pH 7.4 hydrogen addition
- **Processing Order**: Fixed sequence of cleaning operations

### 3. AutoDock Vina Execution

**Seeded Docking Calculation**

```python
def run_vina_docking(self, ligand_pdbqt: str, receptor_pdbqt: str, 
                    center: tuple, size: tuple, seed: int = 42):
    """
    Execute AutoDock Vina with explicit seed
    """
    vina_cmd = [
        'vina',
        '--ligand', ligand_pdbqt,
        '--receptor', receptor_pdbqt,
        '--center_x', str(center[0]),
        '--center_y', str(center[1]), 
        '--center_z', str(center[2]),
        '--size_x', str(size[0]),
        '--size_y', str(size[1]),
        '--size_z', str(size[2]),
        '--seed', str(seed),  # Explicit seed for Vina
        '--exhaustiveness', str(exhaustiveness),
        '--num_modes', str(num_modes),
        '--out', output_pdbqt
    ]
```

**Deterministic Elements:**
- **Random Seed**: Explicit `--seed` parameter to Vina
- **Search Parameters**: Fixed exhaustiveness and mode count
- **Binding Site**: Exact center coordinates and search box size
- **Algorithm**: Vina's internal random number generator seeded consistently

### 4. Result Processing

**PDBQT → SDF Conversion**

```python
def _extract_pose_as_sdf(self, pdbqt_content: str, pose_index: int):
    """
    Deterministic pose extraction and SDF conversion
    
    Process:
    1. Parse PDBQT pose by pose index
    2. Extract coordinates in fixed order
    3. Convert to SDF with consistent atom ordering
    4. Add metadata in standardized format
    """
```

**Deterministic Elements:**
- **Pose Ordering**: Vina outputs poses sorted by affinity (deterministic)
- **Coordinate Extraction**: Fixed parsing order from PDBQT
- **SDF Generation**: Consistent atom ordering and metadata
- **Metadata Addition**: Standardized property embedding

## Parameter Logging

### Comprehensive Tracking

All parameters affecting reproducibility are logged and stored:

```json
{
  "job_metadata": {
    "timestamp": "2025-01-15T10:30:45Z",
    "git_commit": "abc123def456",
    "version": "1.0.0"
  },
  "input_parameters": {
    "ligand_smiles": "CCO",
    "receptor_pdb_id": "1CRN",
    "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
    "size_x": 20.0, "size_y": 20.0, "size_z": 20.0,
    "exhaustiveness": 8,
    "num_modes": 9,
    "seed": 42
  },
  "preparation_metadata": {
    "ligand_preparation": {
      "smiles_canonical": "CCO",
      "protonation_ph": 7.4,
      "embedding_seed": 42,
      "mmff_optimization": true,
      "optimization_converged": true,
      "num_conformers": 1
    },
    "receptor_preparation": {
      "pdb_source": "RCSB",
      "model_selected": 1,
      "chains_retained": ["A"],
      "cofactors_retained": [],
      "hydrogen_addition": true
    }
  },
  "engine_parameters": {
    "engine": "autodock_vina",
    "version": "1.2.5",
    "algorithm_seed": 42,
    "search_algorithm": "monte_carlo_simulated_annealing"
  }
}
```

### Structured Logging

Runtime logging captures all reproducibility-critical events:

```python
logger.info("Ligand preparation started", extra={
    'job_id': job_id,
    'smiles': smiles,
    'preparation_seed': seed,
    'step': 'ligand_preparation'
})

logger.info("3D embedding completed", extra={
    'job_id': job_id,
    'embedding_seed': seed,
    'conformers_generated': num_conformers,
    'optimization_converged': converged
})

logger.info("Vina docking initiated", extra={
    'job_id': job_id,
    'vina_seed': seed,
    'exhaustiveness': exhaustiveness,
    'search_space_volume': volume
})
```

## Reproducibility Validation

### Testing Framework

The system includes multiple validation mechanisms:

#### 1. Smoke Tests
```python
def test_deterministic_docking():
    """Verify identical results with same seed"""
    params = {
        'ligand_smiles': 'CCO',
        'receptor_pdb_id': '1CRN',
        'seed': 42,
        'exhaustiveness': 4,
        'num_modes': 3
    }
    
    result1 = run_docking(params)
    result2 = run_docking(params)
    
    # Should be identical
    assert result1['poses'][0]['affinity'] == result2['poses'][0]['affinity']
    assert result1['poses'][0]['sdf'] == result2['poses'][0]['sdf']
```

#### 2. Benchmark Validation
```python
def test_benchmark_reproducibility():
    """Test benchmark results are consistent"""
    benchmark_id = "1STP"
    
    # Run same benchmark multiple times
    results = [run_benchmark(benchmark_id) for _ in range(3)]
    
    # RMSD and affinities should be identical
    rmsds = [r['rmsd_vs_crystal'] for r in results]
    assert len(set(rmsds)) == 1  # All identical
```

#### 3. Cross-Environment Validation
```python
def test_cross_platform_reproducibility():
    """Verify results across different environments"""
    # Same results on Linux, macOS, Windows (with same library versions)
    assert linux_result == macos_result == windows_result
```

### Reproducibility Hash

Each job generates a reproducibility hash for comparison:

```python
def calculate_reproducibility_hash(results: dict) -> str:
    """Calculate hash of reproducibility-critical components"""
    hash_components = {
        'affinities': [pose['affinity'] for pose in results['poses']],
        'num_poses': len(results['poses']),
        'parameters': {
            'seed': results['docking_parameters']['seed'],
            'exhaustiveness': results['docking_parameters']['exhaustiveness'],
            'num_modes': results['docking_parameters']['num_modes']
        },
        'engine': results['engine_metadata']['engine']
    }
    
    hash_string = json.dumps(hash_components, sort_keys=True)
    return hashlib.sha256(hash_string.encode()).hexdigest()[:16]
```

## Common Reproducibility Issues

### 1. Library Version Differences

**Problem**: Different versions of scientific libraries produce different results
```
RDKit 2023.09.1 vs 2024.03.1 → Different 3D embeddings
Vina 1.2.3 vs 1.2.5 → Algorithm improvements
```

**Solution**: Pin exact library versions in requirements
```python
# requirements.txt
rdkit==2024.03.1
autodock-vina==1.2.5
meeko==0.5.0
```

### 2. Platform Differences

**Problem**: Floating-point arithmetic varies across platforms
**Solution**: Use consistent compilation flags and library builds

### 3. Threading Non-determinism

**Problem**: Multi-threaded Vina calculations may vary
**Solution**: Control threading and use consistent random seeds

### 4. Time-based Variations

**Problem**: Timeout-based early termination
**Solution**: Use deterministic termination criteria instead of wall-clock time

## Best Practices

### For Developers

1. **Always Use Seeds**: Provide explicit seeds for all random operations
2. **Log Parameters**: Record all input parameters and intermediate states
3. **Validate Reproducibility**: Include reproducibility tests in CI/CD
4. **Pin Dependencies**: Lock scientific library versions
5. **Document Changes**: Note any changes that might affect reproducibility

### For Users

1. **Specify Seeds**: Use explicit seed parameters for important calculations
2. **Record Parameters**: Save complete parameter sets for important results
3. **Environment Consistency**: Use same library versions for comparison studies
4. **Validate Results**: Compare reproducibility hashes when rerunning calculations

### For CI/CD

1. **Fixed Seeds**: Use consistent seeds across test runs
2. **Environment Control**: Use containerized environments for testing
3. **Regression Detection**: Compare results against known baselines
4. **Artifact Storage**: Store reproducibility artifacts for comparison

## Environment Setup

### Conda Environment

```yaml
# environment.yml for reproducible setup
name: lipid-rendering-reproducible
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11.7
  - rdkit=2024.03.1
  - autodock-vina=1.2.5
  - openbabel=3.1.1
  - pip:
    - meeko==0.5.0
    - django==5.2.5
```

### Docker Configuration

```dockerfile
# Use specific base image for reproducibility
FROM mambaorg/micromamba:1.5.8

# Install exact versions of scientific libraries
RUN micromamba install -y -n base -c conda-forge \
    python=3.11.7 \
    rdkit=2024.03.1 \
    autodock-vina=1.2.5 \
    openbabel=3.1.1

# Set deterministic Python hash seed
ENV PYTHONHASHSEED=0
```

## Validation Commands

### Check Reproducibility

```bash
# Run same job twice and compare hashes
python manage.py shell -c "
from api.test_smoke_seeded import SmokeTestRunner
runner = SmokeTestRunner()
result1 = runner.run_full_smoke_test()
result2 = runner.run_full_smoke_test()
print(f'Hash 1: {result1[\"reproducibility_hash\"]}')
print(f'Hash 2: {result2[\"reproducibility_hash\"]}')
print(f'Identical: {result1[\"reproducibility_hash\"] == result2[\"reproducibility_hash\"]}')
"
```

### Environment Validation

```bash
# Check library versions
python -c "
import rdkit; print(f'RDKit: {rdkit.__version__}')
import subprocess; print(f'Vina: {subprocess.check_output([\"vina\", \"--version\"]).decode()}')
"
```

## Related Documentation

- [DOCKING_MODES.md](DOCKING_MODES.md) - Real vs mock docking behavior
- [CI.md](CI.md) - Continuous integration and testing
- [VISUALIZATION.md](VISUALIZATION.md) - Result visualization and validation
