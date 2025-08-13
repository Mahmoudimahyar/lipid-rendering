#!/usr/bin/env python3
"""
Benchmark Suite Runner for Lipid Rendering

Runs validation benchmarks against known protein-ligand complexes to detect
regressions and validate docking accuracy. Calculates RMSD vs crystal structures
and tracks binding affinity predictions.
"""

import os
import sys
import json
import yaml
import time
import logging
import requests
import tempfile
import subprocess
from pathlib import Path
from datetime import datetime, timezone
from typing import Dict, List, Any, Optional, Tuple
import argparse

# Add server path for Django imports
benchmark_dir = Path(__file__).parent.parent
server_dir = benchmark_dir.parent / "server"
sys.path.insert(0, str(server_dir))

# Setup Django environment
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'core.settings')
import django
django.setup()


class BenchmarkRunner:
    """
    Runs docking benchmarks against known protein-ligand complexes
    
    Features:
    - Automated PDB download and preparation
    - Ligand extraction from crystal structures
    - Automated docking with standardized parameters
    - RMSD calculation vs native poses
    - Performance tracking and regression detection
    """
    
    def __init__(self, config_path: str = None, server_url: str = "http://localhost:8000"):
        self.server_url = server_url
        self.api_url = f"{server_url}/api"
        
        # Load configuration
        if config_path is None:
            config_path = benchmark_dir / "benchmark_config.yaml"
        
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        # Setup directories
        self.assets_dir = benchmark_dir / "assets"
        self.results_dir = benchmark_dir / "results"
        self.scripts_dir = benchmark_dir / "scripts"
        
        # Ensure directories exist
        self.assets_dir.mkdir(exist_ok=True)
        self.results_dir.mkdir(exist_ok=True)
        
        # Setup logging
        self.setup_logging()
        
        # Benchmark settings
        self.global_settings = self.config['global_settings']
        self.benchmarks = self.config['benchmarks']
        self.thresholds = self.global_settings['thresholds']
        
        # Results storage
        self.results = {}
        
    def setup_logging(self):
        """Setup logging for benchmark runs"""
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        log_file = self.results_dir / f"benchmark_run_{timestamp}.log"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        self.log_file = log_file
    
    def download_pdb(self, pdb_id: str) -> str:
        """Download PDB file from RCSB"""
        pdb_file = self.assets_dir / f"{pdb_id.lower()}.pdb"
        
        if pdb_file.exists():
            self.logger.info(f"PDB {pdb_id} already exists, skipping download")
            return str(pdb_file)
        
        self.logger.info(f"Downloading PDB {pdb_id}")
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            with open(pdb_file, 'w') as f:
                f.write(response.text)
            
            self.logger.info(f"Downloaded PDB {pdb_id} to {pdb_file}")
            return str(pdb_file)
            
        except Exception as e:
            self.logger.error(f"Failed to download PDB {pdb_id}: {e}")
            raise
    
    def extract_ligand_from_pdb(self, pdb_file: str, ligand_resname: str) -> Optional[str]:
        """Extract ligand from PDB file and return SMILES if possible"""
        try:
            # Try to use RDKit if available
            try:
                from rdkit import Chem
                from rdkit.Chem import AllChem
                
                # Read PDB and extract ligand
                mol = Chem.MolFromPDBFile(pdb_file)
                if mol is not None:
                    # Try to get SMILES
                    smiles = Chem.MolToSmiles(mol)
                    self.logger.info(f"Extracted ligand SMILES: {smiles}")
                    return smiles
                    
            except ImportError:
                self.logger.warning("RDKit not available for ligand extraction")
            
            # Fallback: manual PDB parsing
            with open(pdb_file, 'r') as f:
                ligand_lines = []
                for line in f:
                    if line.startswith('HETATM') and ligand_resname in line:
                        ligand_lines.append(line)
            
            if ligand_lines:
                self.logger.info(f"Found {len(ligand_lines)} ligand atoms for {ligand_resname}")
                return None  # Would need complex parsing for SMILES
            else:
                self.logger.warning(f"No ligand {ligand_resname} found in PDB")
                return None
                
        except Exception as e:
            self.logger.error(f"Error extracting ligand: {e}")
            return None
    
    def calculate_rmsd(self, docked_sdf: str, reference_pdb: str, ligand_resname: str) -> Optional[float]:
        """Calculate RMSD between docked pose and crystal structure"""
        try:
            from .rmsd_calculator import RMSDCalculator
            
            calculator = RMSDCalculator()
            rmsd = calculator.calculate_rmsd(docked_sdf, reference_pdb, ligand_resname)
            
            if rmsd is not None:
                self.logger.info(f"Calculated RMSD: {rmsd:.3f} Å")
            else:
                self.logger.warning("RMSD calculation failed")
            
            return rmsd
            
        except Exception as e:
            self.logger.error(f"Error calculating RMSD: {e}")
            return None
    
    def run_docking_benchmark(self, benchmark_id: str, benchmark_config: Dict) -> Dict[str, Any]:
        """Run docking benchmark for a single complex"""
        self.logger.info(f"Running benchmark: {benchmark_id} - {benchmark_config['name']}")
        
        start_time = time.time()
        result = {
            'benchmark_id': benchmark_id,
            'name': benchmark_config['name'],
            'timestamp': datetime.now(timezone.utc).isoformat(),
            'status': 'running',
            'errors': []
        }
        
        try:
            # 1. Download PDB
            pdb_file = self.download_pdb(benchmark_config['pdb_id'])
            
            # 2. Prepare docking parameters
            docking_params = {
                'ligand_smiles': benchmark_config['ligand']['smiles'],
                'receptor_pdb_id': benchmark_config['pdb_id'],
                **benchmark_config['binding_site'],
                **self.global_settings['docking_params']
            }
            
            self.logger.info(f"Docking parameters: {docking_params}")
            
            # 3. Submit docking job
            response = requests.post(
                f"{self.api_url}/dock/run",
                json=docking_params,
                timeout=60
            )
            
            if response.status_code != 200:
                raise Exception(f"Docking submission failed: {response.status_code} - {response.text}")
            
            job_data = response.json()
            job_id = job_data['job_id']
            self.logger.info(f"Submitted docking job: {job_id}")
            
            # 4. Poll for completion
            max_wait = 30 * 60  # 30 minutes
            poll_interval = 10  # 10 seconds
            waited = 0
            
            while waited < max_wait:
                time.sleep(poll_interval)
                waited += poll_interval
                
                status_response = requests.get(f"{self.api_url}/dock/status/{job_id}", timeout=10)
                if status_response.status_code != 200:
                    raise Exception(f"Status check failed: {status_response.status_code}")
                
                status_data = status_response.json()
                current_status = status_data['status']
                
                self.logger.info(f"Status after {waited}s: {current_status}")
                
                if current_status == 'completed':
                    break
                elif current_status == 'failed':
                    error_msg = status_data.get('error_message', 'Unknown error')
                    raise Exception(f"Docking failed: {error_msg}")
            
            if current_status != 'completed':
                raise Exception(f"Docking timeout after {max_wait}s")
            
            # 5. Analyze results
            docking_results = status_data['results']
            poses = docking_results.get('poses', [])
            
            if not poses:
                raise Exception("No poses generated")
            
            # Get best pose (first one, sorted by affinity)
            best_pose = poses[0]
            
            # 6. Calculate RMSD vs crystal
            rmsd = self.calculate_rmsd(
                best_pose.get('sdf', ''),
                pdb_file,
                benchmark_config['ligand']['resname']
            )
            
            # 7. Compile results
            elapsed_time = time.time() - start_time
            
            result.update({
                'status': 'completed',
                'job_id': job_id,
                'elapsed_time': elapsed_time,
                'num_poses': len(poses),
                'best_affinity': best_pose.get('affinity'),
                'rmsd_vs_crystal': rmsd,
                'all_affinities': [pose.get('affinity') for pose in poses],
                'engine_metadata': status_data.get('engine_metadata', {}),
                'docking_parameters': docking_params
            })
            
            # 8. Evaluate success
            success_criteria = self.evaluate_benchmark_success(result, benchmark_config)
            result['success_criteria'] = success_criteria
            result['overall_success'] = all(success_criteria.values())
            
            self.logger.info(f"Benchmark {benchmark_id} completed: Success={result['overall_success']}")
            
        except Exception as e:
            result['status'] = 'failed'
            result['error'] = str(e)
            result['elapsed_time'] = time.time() - start_time
            self.logger.error(f"Benchmark {benchmark_id} failed: {e}")
        
        return result
    
    def evaluate_benchmark_success(self, result: Dict, benchmark_config: Dict) -> Dict[str, bool]:
        """Evaluate if benchmark meets success criteria"""
        criteria = {}
        
        # RMSD criteria
        if result.get('rmsd_vs_crystal') is not None:
            rmsd = result['rmsd_vs_crystal']
            criteria['rmsd_success'] = rmsd <= self.thresholds['rmsd_success']
            criteria['rmsd_excellent'] = rmsd <= self.thresholds['rmsd_excellent']
        else:
            criteria['rmsd_success'] = False
            criteria['rmsd_excellent'] = False
        
        # Pose count criteria
        criteria['sufficient_poses'] = result.get('num_poses', 0) >= self.thresholds['min_poses']
        
        # Affinity reasonableness (compared to expected)
        expected_affinity = benchmark_config.get('expected_results', {}).get('typical_affinity')
        if expected_affinity and result.get('best_affinity'):
            affinity_diff = abs(result['best_affinity'] - expected_affinity)
            criteria['affinity_reasonable'] = affinity_diff <= self.thresholds['affinity_variance'] * 2
        else:
            criteria['affinity_reasonable'] = True  # No baseline to compare
        
        # Completion criteria
        criteria['completed_successfully'] = result['status'] == 'completed'
        
        return criteria
    
    def run_benchmark_suite(self, benchmark_ids: List[str] = None, mode: str = "full") -> Dict[str, Any]:
        """Run complete benchmark suite"""
        if benchmark_ids is None:
            if mode == "pr":
                benchmark_ids = self.config['global_settings']['ci_settings']['pr_benchmarks']
            else:
                benchmark_ids = list(self.benchmarks.keys())
        
        self.logger.info(f"Running benchmark suite: {mode} mode")
        self.logger.info(f"Benchmarks: {benchmark_ids}")
        
        suite_start = time.time()
        suite_results = {
            'suite_info': {
                'mode': mode,
                'benchmark_ids': benchmark_ids,
                'timestamp': datetime.now(timezone.utc).isoformat(),
                'config_version': self.config.get('version', '1.0')
            },
            'individual_results': {},
            'summary': {}
        }
        
        # Run each benchmark
        for benchmark_id in benchmark_ids:
            if benchmark_id not in self.benchmarks:
                self.logger.error(f"Unknown benchmark: {benchmark_id}")
                continue
            
            benchmark_config = self.benchmarks[benchmark_id]
            result = self.run_docking_benchmark(benchmark_id, benchmark_config)
            suite_results['individual_results'][benchmark_id] = result
        
        # Calculate summary statistics
        suite_results['summary'] = self.calculate_suite_summary(suite_results['individual_results'])
        suite_results['suite_info']['elapsed_time'] = time.time() - suite_start
        
        # Save results
        self.save_results(suite_results, mode)
        
        return suite_results
    
    def calculate_suite_summary(self, results: Dict[str, Dict]) -> Dict[str, Any]:
        """Calculate summary statistics for benchmark suite"""
        completed_results = [r for r in results.values() if r['status'] == 'completed']
        failed_results = [r for r in results.values() if r['status'] == 'failed']
        
        if not completed_results:
            return {
                'total_benchmarks': len(results),
                'completed': 0,
                'failed': len(failed_results),
                'success_rate': 0.0
            }
        
        # Success rates
        overall_successes = [r for r in completed_results if r.get('overall_success', False)]
        rmsd_successes = [r for r in completed_results 
                         if r.get('success_criteria', {}).get('rmsd_success', False)]
        
        # RMSD statistics
        rmsds = [r['rmsd_vs_crystal'] for r in completed_results 
                if r.get('rmsd_vs_crystal') is not None]
        
        # Affinity statistics
        affinities = [r['best_affinity'] for r in completed_results 
                     if r.get('best_affinity') is not None]
        
        summary = {
            'total_benchmarks': len(results),
            'completed': len(completed_results),
            'failed': len(failed_results),
            'overall_success_rate': len(overall_successes) / len(completed_results),
            'rmsd_success_rate': len(rmsd_successes) / len(completed_results) if completed_results else 0,
            'rmsd_statistics': {
                'mean': sum(rmsds) / len(rmsds) if rmsds else None,
                'median': sorted(rmsds)[len(rmsds)//2] if rmsds else None,
                'min': min(rmsds) if rmsds else None,
                'max': max(rmsds) if rmsds else None,
            } if rmsds else None,
            'affinity_statistics': {
                'mean': sum(affinities) / len(affinities) if affinities else None,
                'min': min(affinities) if affinities else None,
                'max': max(affinities) if affinities else None,
            } if affinities else None
        }
        
        return summary
    
    def save_results(self, results: Dict, mode: str):
        """Save benchmark results to file"""
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        results_file = self.results_dir / f"benchmark_results_{mode}_{timestamp}.json"
        
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        self.logger.info(f"Results saved to {results_file}")
        
        # Also save as latest for comparison
        latest_file = self.results_dir / f"benchmark_results_{mode}_latest.json"
        with open(latest_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
    
    def check_api_health(self) -> bool:
        """Check if the docking API is available"""
        try:
            response = requests.get(f"{self.api_url}/healthz", timeout=10)
            return response.status_code == 200
        except:
            return False


def main():
    """Main entry point for benchmark runner"""
    parser = argparse.ArgumentParser(description="Run docking benchmarks")
    parser.add_argument('--mode', choices=['full', 'pr', 'custom'], default='full',
                       help='Benchmark mode: full suite, PR subset, or custom')
    parser.add_argument('--benchmarks', nargs='+', 
                       help='Specific benchmarks to run (for custom mode)')
    parser.add_argument('--server-url', default='http://localhost:8000',
                       help='URL of the docking server')
    parser.add_argument('--config', help='Path to benchmark configuration file')
    
    args = parser.parse_args()
    
    # Create runner
    runner = BenchmarkRunner(config_path=args.config, server_url=args.server_url)
    
    # Check API health
    if not runner.check_api_health():
        runner.logger.error("Docking API is not available")
        sys.exit(1)
    
    # Determine benchmarks to run
    if args.mode == 'custom' and args.benchmarks:
        benchmark_ids = args.benchmarks
    else:
        benchmark_ids = None
    
    # Run benchmarks
    try:
        results = runner.run_benchmark_suite(benchmark_ids, args.mode)
        
        # Print summary
        summary = results['summary']
        print(f"\nBenchmark Suite Summary ({args.mode} mode):")
        print(f"Total benchmarks: {summary['total_benchmarks']}")
        print(f"Completed: {summary['completed']}")
        print(f"Failed: {summary['failed']}")
        print(f"Overall success rate: {summary['overall_success_rate']:.2%}")
        print(f"RMSD success rate: {summary['rmsd_success_rate']:.2%}")
        
        if summary['rmsd_statistics']:
            rmsd_stats = summary['rmsd_statistics']
            print(f"RMSD stats: mean={rmsd_stats['mean']:.2f}Å, "
                  f"median={rmsd_stats['median']:.2f}Å, "
                  f"range={rmsd_stats['min']:.2f}-{rmsd_stats['max']:.2f}Å")
        
        # Exit code based on success
        if summary['overall_success_rate'] < 0.5:  # Less than 50% success
            sys.exit(1)
        
    except Exception as e:
        runner.logger.error(f"Benchmark suite failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
