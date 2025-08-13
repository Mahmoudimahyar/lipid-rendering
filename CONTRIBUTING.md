# Contributing to Lipid Rendering

Welcome to Lipid Rendering! This guide explains how to contribute to the project, including development workflow, testing requirements, and quality standards.

## Table of Contents

- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Testing Requirements](#testing-requirements)
- [Code Quality Standards](#code-quality-standards)
- [Scientific Validation](#scientific-validation)
- [Documentation Standards](#documentation-standards)
- [Pull Request Process](#pull-request-process)

## Getting Started

### Prerequisites

**Required Tools:**
- Python 3.11+
- Node.js 18+
- Git
- Docker (optional, for container testing)

**Scientific Libraries (for full functionality):**
- AutoDock Vina
- RDKit
- OpenBabel
- Meeko (optional)

### Development Environment Setup

**1. Clone the Repository**
```bash
git clone https://github.com/your-org/lipid-rendering.git
cd lipid-rendering
```

**2. Backend Setup**
```bash
cd server

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install development dependencies
pip install pytest pytest-django coverage factory-boy freezegun black flake8

# Run migrations
python manage.py migrate

# Verify setup
python manage.py test
```

**3. Frontend Setup**
```bash
cd lipid_viewer

# Install dependencies
npm ci

# Verify setup
npm test

# Start development server
npm start
```

**4. Scientific Libraries Setup (Optional)**

For full docking functionality, install scientific libraries:

```bash
# Using conda/mamba (recommended)
conda create -n lipid-rendering python=3.11
conda activate lipid-rendering
conda install -c conda-forge rdkit autodock-vina openbabel
pip install meeko

# Verify installation
vina --version
python -c "import rdkit; print(rdkit.__version__)"
```

## Development Workflow

### Branch Strategy

**Main Branches:**
- `main` - Production-ready code
- `develop` - Integration branch for features
- `feature/*` - Feature development branches
- `bugfix/*` - Bug fix branches
- `hotfix/*` - Critical production fixes

**Workflow:**
1. Create feature branch from `develop`
2. Implement feature with tests
3. Run local validation
4. Submit pull request to `develop`
5. After review and CI, merge to `develop`
6. Release branches merge `develop` to `main`

### Local Development Commands

**Backend Development:**
```bash
# Run development server
python manage.py runserver

# Run specific tests
pytest api/test_docking.py -v

# Run with coverage
pytest --cov=api --cov-report=term-missing

# Code formatting
black api/
flake8 api/

# Type checking (if enabled)
mypy api/
```

**Frontend Development:**
```bash
# Development server with hot reload
npm start

# Run tests in watch mode
npm test

# Run all tests once
npm test -- --watchAll=false

# Linting and formatting
npm run lint
npm run format

# Build for production
npm run build
```

## Testing Requirements

### Mandatory Testing Rules

**🚨 All features MUST include tests before merge**

1. **New Features**: Comprehensive test coverage for all functionality
2. **Bug Fixes**: Regression tests to prevent reoccurrence
3. **API Changes**: Integration tests for all modified endpoints
4. **UI Changes**: Component tests and accessibility validation
5. **Scientific Features**: Validation against known benchmarks

### Test Categories

**Backend Tests (pytest + Django):**

```python
# 1. Unit Tests - Individual function/method testing
def test_validate_docking_parameters():
    params = {'exhaustiveness': 8, 'num_modes': 9}
    result = validate_docking_parameters(params)
    assert result['valid'] == True

# 2. Integration Tests - API endpoint testing  
def test_docking_api_endpoint():
    response = client.post('/api/dock/run', data=valid_params)
    assert response.status_code == 200
    assert 'job_id' in response.json()

# 3. Model Tests - Database functionality
def test_docking_job_model():
    job = DockingJob.objects.create(ligand_smiles='CCO')
    assert job.status == 'pending'
    assert job.created_at is not None

# 4. Mock Tests - Behavior without scientific libraries
@patch('api.real_docking_engine.vina', None)
def test_mock_docking_fallback():
    result = run_docking(params)
    assert result['engine_metadata']['is_mock'] == True
```

**Frontend Tests (Jest + React Testing Library):**

```javascript
// 1. Component Tests - UI component behavior
test('renders docking form correctly', () => {
  render(<DockingForm onSubmit={mockSubmit} />);
  expect(screen.getByLabelText(/ligand smiles/i)).toBeInTheDocument();
});

// 2. Integration Tests - Component interaction
test('submits docking job successfully', async () => {
  render(<DockingPage />);
  fireEvent.change(screen.getByLabelText(/smiles/i), {target: {value: 'CCO'}});
  fireEvent.click(screen.getByText(/submit/i));
  await waitFor(() => expect(mockApiCall).toHaveBeenCalled());
});

// 3. Accessibility Tests - ARIA compliance
test('form is accessible', () => {
  const { container } = render(<DockingForm />);
  expect(container).toHaveNoViolations();
});

// 4. Snapshot Tests - Visual regression prevention
test('matches snapshot', () => {
  const tree = renderer.create(<MoleculeViewer />).toJSON();
  expect(tree).toMatchSnapshot();
});
```

### Testing Standards

**Coverage Requirements:**
- **Backend**: Minimum 80% line coverage for new code
- **Frontend**: Minimum 75% line coverage for new components
- **Critical Paths**: 95% coverage for docking pipeline and API endpoints

**Performance Requirements:**
- **Unit Tests**: < 100ms per test
- **Integration Tests**: < 5 seconds per test
- **Full Test Suite**: < 10 minutes total execution time

### Smoke Test Compliance

**🚨 All PRs must pass smoke tests**

The smoke test validates the complete docking pipeline:

```bash
# Run local smoke test before PR submission
cd server
python test_smoke_local.py
```

**Smoke Test Validation:**
- ✅ API health and capabilities
- ✅ Complete docking pipeline (SMILES → SDF)
- ✅ Deterministic results with fixed seed
- ✅ Engine metadata validation
- ✅ Result structure compliance

## Code Quality Standards

### Python Code Standards

**Style Guide:**
- Follow PEP 8 with Black formatting
- Line length: 88 characters (Black default)
- Use type hints for all public functions
- Comprehensive docstrings for all modules/classes/functions

**Example:**
```python
from typing import Dict, List, Optional

def validate_docking_parameters(params: Dict[str, Any]) -> Dict[str, bool]:
    """
    Validate docking parameters for scientific accuracy and API compliance.
    
    Args:
        params: Dictionary containing docking parameters
        
    Returns:
        Validation result with 'valid' boolean and 'errors' list
        
    Raises:
        ValueError: If critical parameters are missing or invalid
    """
    errors = []
    
    # Validate required parameters
    required_fields = ['ligand_smiles', 'receptor_pdb_id']
    for field in required_fields:
        if field not in params:
            errors.append(f"Missing required field: {field}")
    
    return {
        'valid': len(errors) == 0,
        'errors': errors
    }
```

**Quality Tools:**
```bash
# Formatting
black api/ --line-length 88

# Linting  
flake8 api/ --max-line-length=88 --ignore=E203,W503

# Import sorting
isort api/ --profile black

# Type checking (if configured)
mypy api/ --strict
```

### JavaScript Code Standards

**Style Guide:**
- Use ES6+ features and modern React patterns
- Prefer function components with hooks
- Use TypeScript for complex components
- ESLint + Prettier configuration

**Example:**
```javascript
import React, { useState, useCallback, useEffect } from 'react';
import PropTypes from 'prop-types';

/**
 * Enhanced molecular viewer component with 3D visualization
 * 
 * @param {Object} props - Component props
 * @param {Array} props.poses - Array of docked poses with SDF content
 * @param {Function} props.onPoseSelect - Callback when pose is selected
 * @param {Object} props.viewerConfig - 3Dmol.js configuration options
 */
const Enhanced3DViewer = ({ poses, onPoseSelect, viewerConfig = {} }) => {
  const [currentPose, setCurrentPose] = useState(0);
  const [viewer, setViewer] = useState(null);

  const handlePoseChange = useCallback((poseIndex) => {
    setCurrentPose(poseIndex);
    onPoseSelect?.(poses[poseIndex]);
  }, [poses, onPoseSelect]);

  // Component implementation...
  
  return (
    <div className="enhanced-3d-viewer">
      {/* Viewer implementation */}
    </div>
  );
};

Enhanced3DViewer.propTypes = {
  poses: PropTypes.arrayOf(PropTypes.shape({
    mode: PropTypes.number.isRequired,
    affinity: PropTypes.number.isRequired,
    sdf: PropTypes.string.isRequired
  })).isRequired,
  onPoseSelect: PropTypes.func,
  viewerConfig: PropTypes.object
};
```

### Scientific Code Standards

**Reproducibility Requirements:**
- All random operations must accept explicit seeds
- Document all algorithms and parameter choices
- Include scientific references for methods used
- Validate against known benchmarks

**Example:**
```python
def prepare_ligand_deterministic(smiles: str, random_seed: int = 42) -> str:
    """
    Prepare ligand for docking with deterministic 3D coordinate generation.
    
    Implements the ETKDG algorithm (Riniker & Landrum, 2015) with explicit
    random seeding for reproducible conformer generation.
    
    Args:
        smiles: SMILES string of ligand molecule
        random_seed: Random seed for coordinate generation
        
    Returns:
        PDBQT format string ready for AutoDock Vina
        
    References:
        Riniker, S. & Landrum, G.A. J Chem Inf Model 55, 2562-2574 (2015)
    """
    # Implementation with explicit seeding...
```

## Scientific Validation

### Benchmark Compliance

**Required Validation:**
- All docking algorithms must pass benchmark validation
- RMSD calculations must align with crystal structures
- Performance metrics must meet established thresholds

**Benchmark Execution:**
```bash
# Run benchmark subset (for PR validation)
cd benchmarks/scripts
python benchmark_runner.py --mode pr

# Run full benchmark suite (for comprehensive validation)
python benchmark_runner.py --mode full
```

### Scientific Accuracy Standards

**RMSD Validation:**
- Top pose RMSD ≤ 2.0 Å for successful docking
- Mean RMSD ≤ 3.0 Å across all poses
- Structural alignment using Maximum Common Substructure

**Affinity Validation:**
- Binding affinities within realistic ranges (-15 to -2 kcal/mol)
- Correlation with experimental data where available
- Consistency across multiple runs with same seed

### Documentation of Scientific Methods

**Required Documentation:**
- Algorithm descriptions with scientific references
- Parameter choices with justification
- Validation methodology and results
- Known limitations and assumptions

## Documentation Standards

### Code Documentation

**Python Docstrings:**
```python
def run_production_docking(ligand_smiles: str, receptor_pdb: str, **kwargs) -> Dict:
    """
    Execute production-quality molecular docking using AutoDock Vina.
    
    This function implements the complete docking pipeline including ligand
    preparation, receptor processing, and Vina execution with deterministic
    seeding for reproducible results.
    
    Args:
        ligand_smiles (str): SMILES representation of ligand molecule
        receptor_pdb (str): PDB ID or file path for receptor structure
        **kwargs: Additional docking parameters (center, size, exhaustiveness, etc.)
        
    Returns:
        Dict containing:
            - poses: List of docked poses with affinities and SDF content
            - metadata: Engine information and timing data
            - parameters: Complete parameter set used for docking
            
    Raises:
        ValidationError: If input parameters are invalid
        DockingError: If docking calculation fails
        
    Example:
        >>> result = run_production_docking(
        ...     ligand_smiles='CCO',
        ...     receptor_pdb='1CRN',
        ...     center_x=0.0, center_y=0.0, center_z=0.0,
        ...     size_x=20.0, size_y=20.0, size_z=20.0,
        ...     seed=42
        ... )
        >>> print(f"Found {len(result['poses'])} poses")
    """
```

**JavaScript JSDoc:**
```javascript
/**
 * Renders interactive 3D molecular visualization using 3Dmol.js
 * 
 * @component
 * @param {Object} props - Component properties
 * @param {Pose[]} props.poses - Array of molecular poses to visualize
 * @param {Function} [props.onPoseSelect] - Callback when user selects a pose
 * @param {ViewerConfig} [props.config] - 3Dmol.js viewer configuration
 * @param {string} [props.className] - Additional CSS classes
 * 
 * @example
 * <Enhanced3DViewer
 *   poses={dockingResults.poses}
 *   onPoseSelect={(pose) => console.log('Selected:', pose)}
 *   config={{ backgroundColor: 'white' }}
 * />
 * 
 * @returns {JSX.Element} The rendered 3D viewer component
 */
```

### API Documentation

**All API endpoints must include:**
- OpenAPI/Swagger specifications
- Request/response examples
- Error code documentation
- Authentication requirements
- Rate limiting information

### Feature Documentation

**New features require:**
- User guide documentation in `docs/`
- Developer documentation for APIs
- Configuration examples
- Migration guides (if applicable)

## Pull Request Process

### PR Requirements Checklist

**Before Submitting:**
- [ ] All tests pass locally (`pytest` and `npm test`)
- [ ] Smoke test passes (`python test_smoke_local.py`)
- [ ] Code follows style guidelines (Black, ESLint)
- [ ] New code has appropriate test coverage
- [ ] Documentation updated for new features
- [ ] No sensitive information in commit history

**PR Description Must Include:**
- [ ] Clear description of changes and motivation
- [ ] List of new/modified functionality
- [ ] Testing strategy and coverage
- [ ] Screenshots for UI changes
- [ ] Breaking changes and migration notes
- [ ] Related issue numbers

### Review Process

**Code Review Requirements:**
1. **At least one reviewer** must approve
2. **All CI checks must pass** (tests, benchmarks, quality gates)
3. **No unresolved review comments**
4. **Up-to-date with target branch**

**Review Focus Areas:**
- Code quality and adherence to standards
- Test coverage and quality
- Security considerations
- Performance implications
- Scientific accuracy
- Documentation completeness

### CI/CD Requirements

**Automated Checks:**
- ✅ Unit tests (backend + frontend)
- ✅ Integration tests
- ✅ Smoke tests with deterministic validation
- ✅ Code quality (linting, formatting)
- ✅ Security scanning
- ✅ Benchmark subset validation (for relevant changes)

**Manual Verification:**
- Scientific accuracy review
- User experience validation
- Performance impact assessment

## Common Development Tasks

### Adding New Docking Features

**1. Implement Backend Logic**
```python
# 1. Add to docking_utils.py
def new_docking_feature(params):
    # Implementation with proper error handling
    pass

# 2. Add API endpoint in views.py
@require_POST 
def new_feature_endpoint(request):
    # API endpoint implementation
    pass

# 3. Add URL routing in urls.py
path('dock/new-feature', new_feature_endpoint, name='new-feature')
```

**2. Add Frontend Integration**
```javascript
// 1. Add API function in dockingApi.js
export const callNewFeature = async (params) => {
  const response = await fetch('/api/dock/new-feature', {
    method: 'POST',
    body: JSON.stringify(params)
  });
  return response.json();
};

// 2. Integrate in UI components
const handleNewFeature = async () => {
  const result = await callNewFeature(params);
  setResults(result);
};
```

**3. Add Tests**
```python
# Backend tests
def test_new_docking_feature():
    result = new_docking_feature(valid_params)
    assert result['success'] == True

def test_new_feature_api():
    response = client.post('/api/dock/new-feature', data=params)
    assert response.status_code == 200
```

```javascript
// Frontend tests
test('new feature integration', async () => {
  render(<ComponentWithNewFeature />);
  await userEvent.click(screen.getByText('New Feature'));
  expect(await screen.findByText('Success')).toBeInTheDocument();
});
```

### Adding New Tests

**Test File Organization:**
```
server/api/
├── test_phase1_seeded_vina.py      # Phase-specific tests
├── test_phase2_true_poses.py
├── test_docking.py                 # Core docking functionality
├── test_models.py                  # Database model tests
└── test_views.py                   # API endpoint tests

lipid_viewer/src/
├── components/
│   ├── MoleculeViewer.test.jsx
│   └── DockingForm.test.jsx
└── utils/
    └── dockingApi.test.js
```

**Test Naming Convention:**
```python
# Descriptive test names
def test_seeded_docking_produces_identical_results():
def test_mock_mode_fallback_when_vina_unavailable():
def test_api_returns_proper_error_for_invalid_smiles():
```

## Getting Help

### Resources

- **Documentation**: See `docs/` directory for comprehensive guides
- **API Reference**: Available at `/api/docs/` when server is running
- **Examples**: Check `examples/` directory for usage patterns
- **Tests**: Existing tests provide implementation examples

### Community

- **Issues**: Report bugs and request features via GitHub Issues
- **Discussions**: Use GitHub Discussions for questions and ideas
- **Code Review**: Tag relevant team members for specialized reviews

### Development Setup Issues

**Common Problems:**

1. **Scientific Library Installation**: See conda environment setup
2. **Test Failures**: Ensure all dependencies installed correctly
3. **Performance Issues**: Check system resources and library versions
4. **Docker Issues**: Verify Docker setup and permissions

**Getting Support:**
- Check existing GitHub Issues for similar problems
- Include full error messages and environment details
- Provide steps to reproduce the issue

## License and Legal

By contributing to Lipid Rendering, you agree that your contributions will be licensed under the same license as the project. Ensure you have the right to submit your contributions and that they don't violate any third-party rights.

---

Thank you for contributing to Lipid Rendering! Your efforts help advance open science and molecular modeling accessibility.
