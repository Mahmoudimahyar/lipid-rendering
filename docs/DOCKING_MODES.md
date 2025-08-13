# Docking Modes: Real vs Mock

This document explains the docking engine modes available in Lipid Rendering, their policies, default behavior, and how to identify which mode is being used.

## Overview

Lipid Rendering supports two primary docking modes:

1. **Real Docking** - Uses actual AutoDock Vina for production-quality molecular docking
2. **Mock Docking** - Provides simulated results for development, testing, and environments without scientific libraries

## Docking Engine Policies

### Default Behavior

- **Production**: Real docking is enforced by default (`DOCKING_ALLOW_MOCK = False`)
- **Development**: Mock docking is allowed when scientific libraries are unavailable
- **CI/Testing**: Mock docking is enabled for reliable testing (`DOCKING_ALLOW_MOCK = True`)

### Mode Selection Logic

The system automatically selects the appropriate docking mode based on:

1. **Scientific Library Availability**: Checks for AutoDock Vina installation
2. **Configuration Setting**: Respects `DOCKING_ALLOW_MOCK` in Django settings
3. **Fallback Policy**: Falls back to mock only if explicitly allowed

```python
# Mode selection pseudocode
if vina_available:
    use_real_docking()
elif DOCKING_ALLOW_MOCK:
    use_mock_docking()
else:
    return_error("Vina unavailable and mock disabled")
```

## Real Docking Mode

### Requirements

- **AutoDock Vina**: Must be installed and accessible via PATH
- **RDKit**: Required for ligand preparation and validation
- **Meeko**: Preferred for PDBQT preparation (OpenBabel fallback)
- **Scientific Environment**: Full conda/pip scientific stack

### Features

- **Production Quality**: Actual molecular docking calculations
- **Deterministic Results**: Seeded calculations for reproducibility
- **Complete Pipeline**: SMILES → 3D → PDBQT → Dock → SDF
- **RMSD Validation**: True structural validation against crystal structures
- **Performance Metrics**: Real calculation times and resource usage

### Capabilities

```json
{
  "vina_available": true,
  "engine_default": "vina",
  "vina_version": "1.2.5",
  "scientific_libs": {
    "rdkit": "2024.03.6",
    "meeko": "0.5.0",
    "openbabel": "3.1.1"
  }
}
```

### Configuration

```python
# settings.py
DOCKING_ALLOW_MOCK = False  # Enforce real docking
```

## Mock Docking Mode

### Purpose

Mock mode is designed for:

- **Development**: Local development without complex scientific setup
- **Testing**: Reliable CI/CD pipeline execution
- **Demonstration**: Quick system validation and UI testing
- **Fallback**: Graceful degradation when real libraries unavailable

### Features

- **Deterministic Results**: Seeded random generation for consistent testing
- **Realistic Output**: Scientifically plausible affinity values and pose counts
- **Fast Execution**: Sub-second docking simulation
- **Complete API**: Full compatibility with real docking API
- **SDF Generation**: Valid SDF structures for visualization testing

### Mock Result Generation

```python
def generate_mock_results(seed: int, num_modes: int):
    random.seed(seed)
    poses = []
    for mode in range(1, num_modes + 1):
        poses.append({
            'mode': mode,
            'affinity': random.uniform(-12.0, -4.0),  # Realistic range
            'rmsd_lb': random.uniform(0.0, 5.0),
            'rmsd_ub': random.uniform(0.0, 8.0),
            'sdf': generate_mock_sdf(ligand_smiles, mode)
        })
    return sorted(poses, key=lambda p: p['affinity'])  # Best first
```

### Configuration

```python
# settings.py
DOCKING_ALLOW_MOCK = True  # Allow fallback to mock
```

## Engine Identification

### API Indicators

All docking responses include engine metadata to identify the mode used:

```json
{
  "job_id": "dock_abc123",
  "status": "completed",
  "results": { ... },
  "engine_metadata": {
    "engine": "vina",           // "vina" or "mock"
    "is_mock": false,           // true for mock mode
    "version": "1.2.5",         // vina version or "mock"
    "method": "autodock_vina",  // "autodock_vina" or "simulation"
    "calculation_time": 45.2    // actual or simulated time
  }
}
```

### Frontend Indicators

The UI displays clear indicators of the docking mode:

- **Real Mode**: Green "REAL" badge with Vina version
- **Mock Mode**: Red "MOCK" badge with warning indicator
- **Capabilities**: API capabilities endpoint shows current configuration

### Logging

Structured logging captures engine selection:

```python
logger.info("Docking job started", extra={
    'job_id': job_id,
    'engine': 'vina' if vina_available else 'mock',
    'is_mock': not vina_available,
    'docking_params': params
})
```

## Configuration Management

### Environment Variables

```bash
# Force real docking only
export DOCKING_ALLOW_MOCK=false

# Allow mock fallback (development/testing)
export DOCKING_ALLOW_MOCK=true
```

### Django Settings

```python
# Production settings
DOCKING_ALLOW_MOCK = False

# Development/CI settings  
DOCKING_ALLOW_MOCK = True
```

### Runtime Detection

```python
from api.real_docking_engine import DockingEngine

# Check current capabilities
engine = DockingEngine()
if engine.is_real_docking_available():
    print("Real docking available")
else:
    print("Using mock docking")
```

## Best Practices

### Development

1. **Local Setup**: Install scientific libraries for full functionality
2. **Mock Testing**: Use mock mode for rapid iteration and UI development
3. **Configuration**: Set `DOCKING_ALLOW_MOCK = True` for development
4. **Validation**: Test both real and mock modes before deployment

### Production

1. **Real Only**: Set `DOCKING_ALLOW_MOCK = False` for production
2. **Monitoring**: Monitor engine metadata in API responses
3. **Validation**: Verify Vina installation during deployment
4. **Fallback**: Have clear error handling when real docking fails

### Testing

1. **CI Configuration**: Use mock mode for reliable CI/CD
2. **Integration Tests**: Test both modes when possible
3. **Smoke Tests**: Include engine capability validation
4. **Regression Tests**: Validate deterministic behavior with seeds

## Troubleshooting

### Common Issues

**Error: "AutoDock Vina not available and mock docking is disabled"**
- **Cause**: Vina not installed and `DOCKING_ALLOW_MOCK = False`
- **Solution**: Install Vina or enable mock mode for development

**Mock mode in production**
- **Detection**: Check `engine_metadata.is_mock` in API responses
- **Fix**: Install scientific libraries and verify Vina availability

**Inconsistent results**
- **Cause**: Switching between real and mock modes
- **Solution**: Use consistent mode and verify engine metadata

### Verification Commands

```bash
# Check Vina installation
vina --version

# Check Python scientific libraries
python -c "import rdkit; print(rdkit.__version__)"
python -c "import meeko; print(meeko.__version__)"

# Test API capabilities
curl http://localhost:8000/api/dock/capabilities
```

## Migration Guide

### From Mock to Real

1. **Install Dependencies**: Set up scientific environment
2. **Update Configuration**: Set `DOCKING_ALLOW_MOCK = False`
3. **Verify Installation**: Check capabilities endpoint
4. **Test Pipeline**: Run smoke test with real docking
5. **Monitor Results**: Verify `is_mock: false` in responses

### From Real to Mock (Development)

1. **Update Configuration**: Set `DOCKING_ALLOW_MOCK = True`
2. **Remove Dependencies**: Optional for lightweight development
3. **Test Fallback**: Verify mock mode activation
4. **Update Tests**: Adapt tests for mock behavior expectations

## Related Documentation

- [REPRODUCIBILITY.md](REPRODUCIBILITY.md) - Deterministic behavior and seeding
- [CI.md](CI.md) - Testing and deployment pipelines
- [CONTRIBUTING.md](CONTRIBUTING.md) - Development guidelines and testing requirements
