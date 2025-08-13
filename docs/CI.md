# Continuous Integration and Testing

This document explains the CI/CD pipeline for Lipid Rendering, including GitHub Actions workflows, smoke testing, nightly benchmarks, and quality assurance processes.

## Overview

Lipid Rendering uses a comprehensive CI/CD pipeline that ensures code quality, functional correctness, and scientific accuracy through automated testing, benchmarks, and deployment validation.

## CI/CD Architecture

```
Pull Request → Unit Tests → Smoke Tests → Benchmark Subset → Merge
     ↓
   Main Branch → Full Tests → Nightly Benchmarks → Deployment → Monitoring
```

## GitHub Actions Workflows

### 1. Main CI Pipeline (`.github/workflows/ci.yml`)

**Triggers:**
- Push to `main` or `develop` branches
- Pull requests to `main` or `develop`

**Jobs:**
1. **Backend Testing** - Python/Django with scientific libraries
2. **Frontend Testing** - React/JavaScript with coverage
3. **Smoke Testing** - End-to-end seeded docking validation
4. **Docker Build** - Container validation (optional)

**Pipeline Flow:**

```yaml
backend-test:
  - Setup micromamba environment
  - Install scientific dependencies (RDKit, Vina, OpenBabel)
  - Run Django system checks
  - Execute unit tests with coverage
  - Upload test results and coverage

frontend-test:
  - Setup Node.js environment
  - Run ESLint and formatting checks
  - Execute Jest unit tests with coverage
  - Upload test results

smoke-test:
  - Start Django development server
  - Run seeded smoke test with real/mock docking
  - Validate deterministic results
  - Upload artifacts (JSON, SDF, logs)

docker-build:
  - Build production Docker container
  - Test container health endpoints
  - Validate deployment configuration
```

### 2. Nightly Benchmark Suite (`.github/workflows/nightly-benchmarks.yml`)

**Triggers:**
- Scheduled: Daily at 2 AM UTC
- Manual dispatch with configurable parameters

**Features:**
- **Complete Benchmark Suite**: All 5 protein-ligand complexes
- **Regression Detection**: Automated threshold monitoring
- **Performance Tracking**: Historical comparison and trending
- **Issue Creation**: Automatic GitHub issues on regression

**Benchmark Execution:**

```yaml
nightly-benchmarks:
  - Setup complete scientific environment
  - Run full benchmark suite (1STP, 1HVR, 4DFR, 3PTB, 1ERE)
  - Calculate RMSD vs crystal structures
  - Generate performance reports
  - Compare against historical baselines
  - Create GitHub issues on regression detection
```

## Testing Strategy

### 1. Unit Testing

**Backend Tests (pytest + Django)**

```python
# Test categories
- Model tests: Database schema and relationships
- View tests: API endpoints and business logic  
- Integration tests: End-to-end pipeline validation
- Mock tests: Behavior without scientific libraries
- Phase tests: Feature-specific validation suites

# Test execution
pytest -v --cov=api --cov-report=xml --maxfail=5
```

**Frontend Tests (Jest + React Testing Library)**

```javascript
// Test categories
- Component tests: UI component behavior
- Integration tests: Component interaction
- Utility tests: Helper function validation
- Snapshot tests: Visual regression detection
- Accessibility tests: ARIA compliance

// Test execution  
npm test -- --coverage --watchAll=false
```

### 2. Smoke Testing

**Purpose**: Fast, deterministic validation of core docking pipeline

**Test Configuration:**
```python
SMOKE_TEST_PARAMS = {
    'ligand_smiles': 'CCO',        # Simple ethanol
    'receptor_pdb_id': '1CRN',     # Small crambin protein
    'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0,
    'size_x': 15.0, 'size_y': 15.0, 'size_z': 15.0,
    'exhaustiveness': 4,           # Fast for CI
    'num_modes': 3,               # Minimal for speed
    'seed': 42                    # Fixed for reproducibility
}
```

**Validation Criteria:**
- ✅ API health and capabilities check
- ✅ Docking job submission and completion
- ✅ Result structure validation (poses, affinities, SDF)
- ✅ Engine metadata verification (real vs mock)
- ✅ Reproducibility hash consistency

**Execution Time:** Target < 5 minutes for complete smoke test

### 3. Benchmark Testing

**PR Benchmarks (Lightweight Subset):**
- **1STP** (Streptavidin-Biotin): Strong binding validation
- **3PTB** (Trypsin-Benzamidine): Simple system verification
- **Time Limit**: 10 minutes maximum
- **Success Criteria**: RMSD ≤ 2.0Å, completion without errors

**Nightly Benchmarks (Complete Suite):**
- **All 5 Complexes**: Comprehensive scientific validation
- **RMSD Calculation**: vs crystal structures using RDKit alignment
- **Success Thresholds**: 60% overall success, 50% RMSD success
- **Time Limit**: 2 hours maximum

## Scientific Library Management

### Conda Environment Setup

**Micromamba Integration:**
```yaml
- name: Setup micromamba
  uses: mamba-org/setup-micromamba@v1
  with:
    micromamba-version: '1.5.8'
    environment-name: lipid-rendering-ci
    cache-environment: true
    cache-downloads: true
```

**Scientific Dependencies:**
```yaml
- name: Install scientific dependencies
  run: |
    micromamba install -c conda-forge -y \
      rdkit \
      openbabel \
      autodock-vina \
      numpy \
      scipy
```

**Fallback Strategy:**
- **Real Docking**: When scientific libraries available
- **Mock Docking**: Graceful fallback for environments without libraries
- **CI Configuration**: `DOCKING_ALLOW_MOCK = True` for testing

## Configuration Management

### Environment-Specific Settings

**CI Settings (`core/settings_ci.py`):**
```python
# Enable mock docking for testing
DOCKING_ALLOW_MOCK = True

# Use in-memory database for speed
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': ':memory:',
    }
}

# Optimized logging for CI
LOGGING = {
    'handlers': {
        'console': {'class': 'logging.StreamHandler'},
        'file': {'class': 'logging.FileHandler', 'filename': 'ci_test.log'}
    }
}
```

**Production Settings:**
```python
# Enforce real docking only
DOCKING_ALLOW_MOCK = False

# Production database configuration
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        # ... production config
    }
}
```

## Artifact Management

### Test Artifacts

**Backend Artifacts:**
- Unit test results (XML format)
- Coverage reports (XML for Codecov integration)
- Performance metrics and timing data

**Frontend Artifacts:**
- Jest test results and coverage
- ESLint reports
- Build artifacts and bundle analysis

**Smoke Test Artifacts:**
- JSON results with full metadata
- SDF files for each docked pose
- Structured logs with execution details
- Reproducibility hashes for comparison

**Benchmark Artifacts:**
- Complete benchmark results (JSON)
- RMSD calculations and statistics
- PDB structures and derived files
- Performance metrics and timing data

### Retention Policies

```yaml
retention-days: 30   # Standard test artifacts
retention-days: 90   # Benchmark results and CI summaries
retention-days: 365  # Production deployment artifacts
```

## Quality Gates

### Required Checks for PR Merge

1. **✅ All Unit Tests Pass**: Backend and frontend test suites
2. **✅ Smoke Test Passes**: Core docking pipeline validation
3. **✅ PR Benchmarks Pass**: Subset benchmark validation
4. **✅ Code Coverage**: Maintain minimum coverage thresholds
5. **✅ Code Quality**: ESLint, formatting, and style checks

### Nightly Quality Validation

1. **✅ Full Benchmark Suite**: All 5 protein-ligand complexes
2. **✅ Regression Detection**: Performance within historical ranges
3. **✅ Scientific Accuracy**: RMSD validation against crystal structures
4. **✅ Performance Monitoring**: Execution time and resource usage

## Monitoring and Alerting

### Regression Detection

**Automated Threshold Monitoring:**
```python
REGRESSION_THRESHOLDS = {
    'min_overall_success_rate': 0.6,  # 60% minimum
    'min_rmsd_success_rate': 0.5,     # 50% minimum  
    'max_mean_rmsd': 3.0,             # 3.0Å maximum
    'max_execution_time': 7200,       # 2 hours maximum
}
```

**GitHub Issue Creation:**
```yaml
- name: Create issue on regression
  if: failure()
  uses: actions/github-script@v7
  with:
    script: |
      github.rest.issues.create({
        owner: context.repo.owner,
        repo: context.repo.repo,
        title: `Nightly Benchmark Regression - ${new Date().toISOString().split('T')[0]}`,
        body: `Benchmark regression detected...`,
        labels: ['regression', 'benchmarks', 'automated']
      });
```

### Performance Tracking

**Metrics Collection:**
- Test execution times
- Benchmark RMSD statistics
- Docking calculation performance
- Resource utilization (memory, CPU)

**Historical Comparison:**
- Compare against previous runs
- Detect performance trends
- Alert on significant changes

## Deployment Validation

### Docker Container Testing

**Container Health Checks:**
```yaml
- name: Test Docker container
  run: |
    docker run -d --name lipid-test -p 8000:8000 lipid-rendering:ci
    sleep 15
    curl -f http://localhost:8000/api/healthz
    docker stop lipid-test && docker rm lipid-test
```

**Multi-Stage Build Validation:**
```dockerfile
# Production build stage
FROM python:3.11-slim as production
COPY requirements.txt .
RUN pip install -r requirements.txt

# Test stage
FROM production as test  
COPY . .
RUN python manage.py test

# Final production stage
FROM production as final
COPY . .
EXPOSE 8000
CMD ["gunicorn", "core.wsgi:application"]
```

## Local Development Integration

### Pre-commit Hooks

**Setup:**
```bash
# Install pre-commit
pip install pre-commit

# Setup hooks
pre-commit install

# Manual run
pre-commit run --all-files
```

**Hook Configuration (`.pre-commit-config.yaml`):**
```yaml
repos:
  - repo: https://github.com/psf/black
    rev: 22.3.0
    hooks:
      - id: black
        language_version: python3.11

  - repo: https://github.com/pycqa/flake8
    rev: 4.0.1
    hooks:
      - id: flake8

  - repo: https://github.com/pre-commit/mirrors-eslint
    rev: v8.44.0
    hooks:
      - id: eslint
        files: \.(js|jsx)$
```

### Local Testing Commands

**Backend Testing:**
```bash
# Quick test run
pytest -q --maxfail=1 --disable-warnings

# Full test suite with coverage
pytest -v --cov=api --cov-report=term-missing

# Specific test phases
pytest api/test_phase1_seeded_vina.py -v
```

**Frontend Testing:**
```bash
# Unit tests
npm test

# Coverage report
npm test -- --coverage

# Linting
npm run lint
```

**Smoke Test:**
```bash
# Local smoke test validation
python server/test_smoke_local.py

# Manual smoke test run
cd server && python -c "
from api.test_smoke_seeded import SmokeTestRunner
runner = SmokeTestRunner()
result = runner.run_full_smoke_test()
print(f'Status: {result[\"status\"]}')
"
```

## CI Performance Optimization

### Parallel Execution

**Job Parallelization:**
- Backend and frontend tests run simultaneously
- Independent artifact upload and processing
- Concurrent Docker builds when applicable

**Test Parallelization:**
```yaml
# Backend testing with parallel execution
strategy:
  matrix:
    test-group: [unit, integration, phase-tests]
    
steps:
  - name: Run test group
    run: pytest tests/${{ matrix.test-group }} -v
```

### Caching Strategy

**Dependency Caching:**
```yaml
# Conda environment caching
cache-environment: true
cache-downloads: true

# Node.js dependency caching
cache: 'npm'
cache-dependency-path: lipid_viewer/package-lock.json
```

**Docker Layer Caching:**
```yaml
- name: Set up Docker Buildx
  uses: docker/setup-buildx-action@v3
  with:
    buildkitd-flags: --allow-insecure-entitlement security.insecure
```

## Troubleshooting

### Common CI Issues

**Scientific Library Installation Failures:**
```yaml
# Retry strategy for conda installs
- name: Install dependencies with retry
  run: |
    for i in {1..3}; do
      micromamba install -c conda-forge rdkit && break
      echo "Attempt $i failed, retrying..."
      sleep 10
    done
```

**Test Timeouts:**
```yaml
# Appropriate timeouts for different job types
timeout-minutes: 20   # Frontend tests
timeout-minutes: 30   # Backend tests  
timeout-minutes: 120  # Nightly benchmarks
```

**Flaky Test Handling:**
```python
# Retry decorator for flaky tests
@pytest.mark.flaky(reruns=3, reruns_delay=2)
def test_network_dependent_feature():
    # Test that might fail due to network issues
    pass
```

### Debug Strategies

**Enhanced Logging:**
```yaml
- name: Enable debug logging
  run: |
    export DJANGO_LOG_LEVEL=DEBUG
    export PYTEST_VERBOSE=1
    pytest -v -s --tb=long
```

**Artifact Collection for Debugging:**
```yaml
- name: Upload debug artifacts
  if: failure()
  uses: actions/upload-artifact@v4
  with:
    name: debug-logs
    path: |
      server/*.log
      server/db.sqlite3
      /tmp/pytest-*
```

## Related Documentation

- [DOCKING_MODES.md](DOCKING_MODES.md) - Understanding testing with real vs mock docking
- [REPRODUCIBILITY.md](REPRODUCIBILITY.md) - Deterministic testing and validation
- [CONTRIBUTING.md](CONTRIBUTING.md) - Development workflow and testing requirements
