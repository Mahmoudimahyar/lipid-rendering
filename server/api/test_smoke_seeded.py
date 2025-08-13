"""
Seeded Smoke Test for CI Pipeline

This module implements a deterministic smoke test that runs a complete
docking pipeline with fixed parameters to ensure reproducibility and
detect regressions across CI runs.
"""

import os
import json
import time
import logging
import requests
from datetime import datetime, timezone
from typing import Dict, Any, List
import hashlib


class SmokeTestRunner:
    """
    Runs deterministic smoke tests for the docking pipeline
    
    Features:
    - Fixed seed for reproducible results
    - Fast parameters for CI efficiency  
    - Complete pipeline validation
    - Artifact generation for comparison
    """
    
    def __init__(self, base_url: str = "http://localhost:8000"):
        self.base_url = base_url
        self.api_url = f"{base_url}/api"
        self.test_timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        
        # Fixed test parameters for deterministic results
        self.test_params = {
            'ligand_smiles': 'CCO',  # Simple ethanol for fast testing
            'receptor_pdb_id': '1CRN',  # Small crambin protein
            'center_x': 0.0,
            'center_y': 0.0, 
            'center_z': 0.0,
            'size_x': 15.0,  # Smaller search space for speed
            'size_y': 15.0,
            'size_z': 15.0,
            'exhaustiveness': 4,  # Reduced for speed
            'num_modes': 3,  # Minimal modes for speed
            'seed': 42  # Fixed seed for reproducibility
        }
        
        # Setup logging
        self.setup_logging()
    
    def setup_logging(self):
        """Setup structured logging for smoke test"""
        log_filename = f"smoke_test_{self.test_timestamp}.log"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_filename),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        self.log_file = log_filename
    
    def test_api_health(self) -> bool:
        """Test that the API is responding"""
        try:
            response = requests.get(f"{self.api_url}/healthz", timeout=10)
            if response.status_code == 200:
                self.logger.info("[PASS] API health check passed")
                return True
            else:
                self.logger.error(f"[FAIL] API health check failed: {response.status_code}")
                return False
        except Exception as e:
            self.logger.error(f"[FAIL] API health check exception: {e}")
            return False
    
    def test_capabilities(self) -> Dict[str, Any]:
        """Test docking capabilities endpoint"""
        try:
            response = requests.get(f"{self.api_url}/dock/capabilities", timeout=10)
            if response.status_code == 200:
                data = response.json()
                self.logger.info("[PASS] Capabilities check passed")
                self.logger.info(f"   Vina available: {data.get('vina_available', False)}")
                self.logger.info(f"   Default engine: {data.get('engine_default', 'unknown')}")
                return data
            else:
                self.logger.error(f"[FAIL] Capabilities check failed: {response.status_code}")
                return {}
        except Exception as e:
            self.logger.error(f"[FAIL] Capabilities check exception: {e}")
            return {}
    
    def run_seeded_docking(self) -> Dict[str, Any]:
        """Run the main seeded docking test"""
        self.logger.info("[START] Starting seeded docking test")
        self.logger.info(f"   Parameters: {json.dumps(self.test_params, indent=2)}")
        
        start_time = time.time()
        
        try:
            # Submit docking job
            response = requests.post(
                f"{self.api_url}/dock/run",
                json=self.test_params,
                timeout=30
            )
            
            # If docking fails due to no Vina, try enabling mock mode
            if response.status_code == 503:
                response_data = response.json()
                if "AutoDock Vina not available" in response_data.get('error', ''):
                    self.logger.warning("[WARN] Vina not available, attempting with mock mode")
                    
                    # Try to enable mock mode temporarily by calling capabilities to check settings
                    caps_response = requests.get(f"{self.api_url}/dock/capabilities", timeout=10)
                    if caps_response.status_code == 200:
                        caps_data = caps_response.json()
                        if not caps_data.get('mock_allowed', False):
                            self.logger.warning("[WARN] Mock mode not allowed in settings, smoke test will use mock anyway for CI")
                    
                    # Note: In a real implementation, we might need to modify settings or use a test-specific endpoint
                    # For now, we'll document this as expected behavior for local testing
                    raise Exception(f"Real Vina not available and mock mode disabled. This is expected for local testing without scientific libraries. In CI, mock mode should be enabled. Error: {response.text}")
            
            if response.status_code != 200:
                raise Exception(f"Docking submission failed: {response.status_code} - {response.text}")
            
            job_data = response.json()
            job_id = job_data['job_id']
            self.logger.info(f"[PASS] Docking job submitted: {job_id}")
            
            # Poll for completion
            max_wait_time = 300  # 5 minutes max for smoke test
            poll_interval = 5  # Poll every 5 seconds
            waited_time = 0
            
            while waited_time < max_wait_time:
                time.sleep(poll_interval)
                waited_time += poll_interval
                
                status_response = requests.get(
                    f"{self.api_url}/dock/status/{job_id}",
                    timeout=10
                )
                
                if status_response.status_code != 200:
                    raise Exception(f"Status check failed: {status_response.status_code}")
                
                status_data = status_response.json()
                current_status = status_data['status']
                
                self.logger.info(f"   Status after {waited_time}s: {current_status}")
                
                if current_status == 'completed':
                    elapsed_time = time.time() - start_time
                    self.logger.info(f"[PASS] Docking completed in {elapsed_time:.1f}s")
                    return status_data
                elif current_status == 'failed':
                    error_msg = status_data.get('error_message', 'Unknown error')
                    raise Exception(f"Docking job failed: {error_msg}")
            
            raise Exception(f"Docking job timeout after {max_wait_time}s")
            
        except Exception as e:
            self.logger.error(f"[FAIL] Seeded docking test failed: {e}")
            raise
    
    def validate_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Validate docking results for expected structure and content"""
        validation_report = {
            'passed': True,
            'checks': [],
            'warnings': [],
            'errors': []
        }
        
        def add_check(name: str, passed: bool, message: str = ""):
            validation_report['checks'].append({
                'name': name,
                'passed': passed,
                'message': message
            })
            if not passed:
                validation_report['passed'] = False
                validation_report['errors'].append(f"{name}: {message}")
        
        def add_warning(message: str):
            validation_report['warnings'].append(message)
        
        # Check basic structure
        add_check(
            "Has results",
            'results' in results and results['results'] is not None,
            "Results section missing"
        )
        
        if 'results' in results and results['results']:
            result_data = results['results']
            
            # Check poses
            poses = result_data.get('poses', [])
            add_check(
                "Has poses",
                len(poses) > 0,
                f"Expected poses, got {len(poses)}"
            )
            
            add_check(
                "Correct number of poses",
                len(poses) == self.test_params['num_modes'],
                f"Expected {self.test_params['num_modes']} poses, got {len(poses)}"
            )
            
            # Check pose structure
            if poses:
                first_pose = poses[0]
                required_fields = ['mode', 'affinity', 'sdf']
                for field in required_fields:
                    add_check(
                        f"Pose has {field}",
                        field in first_pose,
                        f"Missing {field} in pose"
                    )
                
                # Check SDF content
                if 'sdf' in first_pose:
                    sdf_content = first_pose['sdf']
                    add_check(
                        "SDF not empty",
                        len(sdf_content.strip()) > 0,
                        "SDF content is empty"
                    )
                    
                    add_check(
                        "SDF has molecule end marker",
                        '$$$$' in sdf_content,
                        "SDF missing molecule end marker"
                    )
        
        # Check engine metadata
        engine_meta = results.get('engine_metadata', {})
        add_check(
            "Has engine metadata",
            bool(engine_meta),
            "Engine metadata missing"
        )
        
        if engine_meta:
            add_check(
                "Engine specified",
                'engine' in engine_meta,
                "Engine type not specified"
            )
            
            add_check(
                "Mock status specified",
                'is_mock' in engine_meta,
                "Mock status not specified"
            )
        
        # Check parameters preservation
        params = results.get('docking_parameters', {})
        add_check(
            "Parameters preserved",
            params.get('seed') == self.test_params['seed'],
            f"Seed mismatch: expected {self.test_params['seed']}, got {params.get('seed')}"
        )
        
        # Log validation results
        if validation_report['passed']:
            self.logger.info("[PASS] Results validation passed")
        else:
            self.logger.error("[FAIL] Results validation failed")
            for error in validation_report['errors']:
                self.logger.error(f"   {error}")
        
        for warning in validation_report['warnings']:
            self.logger.warning(f"   {warning}")
        
        return validation_report
    
    def calculate_result_hash(self, results: Dict[str, Any]) -> str:
        """Calculate reproducibility hash from key result components"""
        # Extract deterministic components for hashing
        hash_components = {}
        
        if 'results' in results and results['results']:
            result_data = results['results']
            poses = result_data.get('poses', [])
            
            # Hash pose affinities (should be deterministic with same seed)
            if poses:
                hash_components['affinities'] = [pose.get('affinity') for pose in poses]
                hash_components['num_poses'] = len(poses)
        
        # Hash key parameters
        params = results.get('docking_parameters', {})
        hash_components['parameters'] = {
            'seed': params.get('seed'),
            'exhaustiveness': params.get('exhaustiveness'),
            'num_modes': params.get('num_modes')
        }
        
        # Hash engine info
        engine = results.get('engine_metadata', {})
        hash_components['engine'] = {
            'engine': engine.get('engine'),
            'is_mock': engine.get('is_mock')
        }
        
        # Create hash
        hash_string = json.dumps(hash_components, sort_keys=True)
        result_hash = hashlib.sha256(hash_string.encode()).hexdigest()[:16]
        
        self.logger.info(f"[HASH] Result hash: {result_hash}")
        return result_hash
    
    def save_artifacts(self, results: Dict[str, Any], validation: Dict[str, Any]) -> str:
        """Save test artifacts for analysis and comparison"""
        artifact_data = {
            'smoke_test_metadata': {
                'timestamp': self.test_timestamp,
                'test_parameters': self.test_params,
                'git_commit': os.environ.get('GITHUB_SHA', 'unknown'),
                'git_ref': os.environ.get('GITHUB_REF_NAME', 'unknown'),
                'build_number': os.environ.get('GITHUB_RUN_NUMBER', 'unknown'),
            },
            'docking_results': results,
            'validation_report': validation,
            'reproducibility_hash': self.calculate_result_hash(results)
        }
        
        # Save main results file
        results_filename = f"smoke_test_results_{self.test_timestamp}.json"
        with open(results_filename, 'w') as f:
            json.dump(artifact_data, f, indent=2, default=str)
        
        self.logger.info(f"[SAVE] Results saved to {results_filename}")
        
        # Save individual SDF files for poses
        if 'results' in results and results['results']:
            poses = results['results'].get('poses', [])
            for i, pose in enumerate(poses):
                if 'sdf' in pose:
                    sdf_filename = f"smoke_test_pose_{i+1}_{self.test_timestamp}.sdf"
                    with open(sdf_filename, 'w') as f:
                        f.write(pose['sdf'])
                    self.logger.info(f"[SAVE] Pose {i+1} SDF saved to {sdf_filename}")
        
        return results_filename
    
    def run_full_smoke_test(self) -> Dict[str, Any]:
        """Run the complete smoke test pipeline"""
        self.logger.info("[START] Starting CI Smoke Test Pipeline")
        self.logger.info(f"   Timestamp: {self.test_timestamp}")
        
        try:
            # 1. Test API health
            if not self.test_api_health():
                raise Exception("API health check failed")
            
            # 2. Test capabilities
            capabilities = self.test_capabilities()
            
            # 3. Run seeded docking
            results = self.run_seeded_docking()
            
            # 4. Validate results
            validation = self.validate_results(results)
            
            # 5. Save artifacts
            output_file = self.save_artifacts(results, validation)
            
            # 6. Final status
            if validation['passed']:
                self.logger.info("[SUCCESS] Smoke test PASSED - All checks successful")
                return {
                    'status': 'PASSED',
                    'output_file': output_file,
                    'validation': validation,
                    'capabilities': capabilities
                }
            else:
                self.logger.error("[FAIL] Smoke test FAILED - Validation errors")
                return {
                    'status': 'FAILED',
                    'output_file': output_file,
                    'validation': validation,
                    'capabilities': capabilities
                }
                
        except Exception as e:
            self.logger.error(f"[CRASH] Smoke test CRASHED: {e}")
            return {
                'status': 'CRASHED',
                'error': str(e),
                'timestamp': self.test_timestamp
            }


def main():
    """Main entry point for standalone execution"""
    import sys
    
    runner = SmokeTestRunner()
    result = runner.run_full_smoke_test()
    
    if result['status'] != 'PASSED':
        sys.exit(1)
    
    print("Smoke test completed successfully!")


if __name__ == '__main__':
    main()
