"""
Tests for Phase 6: Benchmark Suite & Regression Tests

This module tests the benchmark infrastructure including RMSD calculation,
benchmark runner functionality, and CI integration.
"""

import pytest
import json
import tempfile
import os
import sys
import yaml
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
from django.test import TestCase

# Add benchmarks path for imports
benchmark_dir = Path(__file__).parent.parent.parent / "benchmarks"


class TestBenchmarkConfiguration(TestCase):
    """Test benchmark configuration loading and validation"""
    
    def setUp(self):
        self.config_path = benchmark_dir / "benchmark_config.yaml"
    
    def test_benchmark_config_exists(self):
        """Test that benchmark configuration file exists"""
        self.assertTrue(self.config_path.exists(), "benchmark_config.yaml not found")
    
    def test_benchmark_config_valid_yaml(self):
        """Test that benchmark configuration is valid YAML"""
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        self.assertIsInstance(config, dict)
        self.assertIn('global_settings', config)
        self.assertIn('benchmarks', config)
    
    def test_benchmark_config_has_required_benchmarks(self):
        """Test that all required benchmarks are defined"""
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        benchmarks = config['benchmarks']
        required_benchmarks = ["1STP", "1HVR", "4DFR", "3PTB", "1ERE"]
        
        for benchmark_id in required_benchmarks:
            self.assertIn(benchmark_id, benchmarks, f"Missing benchmark: {benchmark_id}")
    
    def test_benchmark_entries_have_required_fields(self):
        """Test that benchmark entries have all required fields"""
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        required_fields = [
            'name', 'description', 'pdb_id', 'ligand', 'protein', 
            'binding_site', 'expected_results'
        ]
        
        for benchmark_id, benchmark in config['benchmarks'].items():
            for field in required_fields:
                self.assertIn(field, benchmark, 
                            f"Benchmark {benchmark_id} missing field: {field}")
            
            # Test ligand fields
            self.assertIn('smiles', benchmark['ligand'])
            self.assertIn('resname', benchmark['ligand'])
            
            # Test binding site fields
            binding_site = benchmark['binding_site']
            for coord in ['center_x', 'center_y', 'center_z', 'size_x', 'size_y', 'size_z']:
                self.assertIn(coord, binding_site)
                self.assertIsInstance(binding_site[coord], (int, float))
    
    def test_global_settings_structure(self):
        """Test that global settings have correct structure"""
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        global_settings = config['global_settings']
        
        # Test docking parameters
        self.assertIn('docking_params', global_settings)
        docking_params = global_settings['docking_params']
        self.assertIn('exhaustiveness', docking_params)
        self.assertIn('num_modes', docking_params)
        self.assertIn('seed', docking_params)
        
        # Test thresholds
        self.assertIn('thresholds', global_settings)
        thresholds = global_settings['thresholds']
        self.assertIn('rmsd_success', thresholds)
        self.assertIn('affinity_variance', thresholds)
        
        # Test CI settings
        self.assertIn('ci_settings', global_settings)
        ci_settings = global_settings['ci_settings']
        self.assertIn('pr_benchmarks', ci_settings)
        self.assertIsInstance(ci_settings['pr_benchmarks'], list)


class TestRMSDCalculator(TestCase):
    """Test RMSD calculation functionality"""
    
    def test_rmsd_calculator_concept(self):
        """Test RMSD calculator concept and structure"""
        # Test that RMSD calculator file exists
        rmsd_calc_file = benchmark_dir / "scripts" / "rmsd_calculator.py"
        self.assertTrue(rmsd_calc_file.exists(), "rmsd_calculator.py not found")
        
        # Test that the file contains expected classes/functions
        with open(rmsd_calc_file, 'r') as f:
            content = f.read()
        
        self.assertIn('class RMSDCalculator', content)
        self.assertIn('def calculate_rmsd', content)
        self.assertIn('def calculate_detailed_analysis', content)
    
    def test_rmsd_fallback_function(self):
        """Test that RMSD calculation has fallback functionality"""
        rmsd_calc_file = benchmark_dir / "scripts" / "rmsd_calculator.py"
        
        with open(rmsd_calc_file, 'r') as f:
            content = f.read()
        
        # Should have fallback for when RDKit is not available
        self.assertIn('_fallback_rmsd_calculation', content)
        self.assertIn('rdkit_available', content)
    
    def test_rmsd_alignment_methods(self):
        """Test that RMSD calculator supports multiple alignment methods"""
        rmsd_calc_file = benchmark_dir / "scripts" / "rmsd_calculator.py"
        
        with open(rmsd_calc_file, 'r') as f:
            content = f.read()
        
        # Should support multiple alignment strategies
        self.assertIn('_mcs_based_alignment', content)
        self.assertIn('_symmetry_aware_alignment', content)
        self.assertIn('_calculate_aligned_rmsd', content)


class TestBenchmarkRunner(TestCase):
    """Test benchmark runner functionality"""
    
    def test_benchmark_runner_structure(self):
        """Test benchmark runner file structure and key components"""
        runner_file = benchmark_dir / "scripts" / "benchmark_runner.py"
        self.assertTrue(runner_file.exists(), "benchmark_runner.py not found")
        
        with open(runner_file, 'r') as f:
            content = f.read()
        
        # Test key classes and methods exist
        self.assertIn('class BenchmarkRunner', content)
        self.assertIn('def run_benchmark_suite', content)
        self.assertIn('def run_docking_benchmark', content)
        self.assertIn('def calculate_rmsd', content)
        self.assertIn('def evaluate_benchmark_success', content)
    
    def test_benchmark_runner_imports_rmsd(self):
        """Test that benchmark runner imports RMSD calculator"""
        runner_file = benchmark_dir / "scripts" / "benchmark_runner.py"
        
        with open(runner_file, 'r') as f:
            content = f.read()
        
        # Should import and use RMSD calculator
        self.assertIn('rmsd_calculator', content)
        self.assertIn('RMSDCalculator', content)
    
    def test_benchmark_runner_config_handling(self):
        """Test that benchmark runner handles configuration properly"""
        runner_file = benchmark_dir / "scripts" / "benchmark_runner.py"
        
        with open(runner_file, 'r') as f:
            content = f.read()
        
        # Should handle YAML configuration
        self.assertIn('yaml', content)
        self.assertIn('config', content)
        self.assertIn('thresholds', content)
    
    def test_benchmark_evaluation_logic(self):
        """Test benchmark evaluation criteria structure"""
        runner_file = benchmark_dir / "scripts" / "benchmark_runner.py"
        
        with open(runner_file, 'r') as f:
            content = f.read()
        
        # Should have evaluation criteria
        self.assertIn('rmsd_success', content)
        self.assertIn('rmsd_excellent', content)
        self.assertIn('affinity_reasonable', content)
        self.assertIn('sufficient_poses', content)


class TestBenchmarkCIIntegration(TestCase):
    """Test CI integration for benchmarks"""
    
    def test_benchmark_pr_subset_defined(self):
        """Test that PR benchmark subset is properly defined"""
        config_path = benchmark_dir / "benchmark_config.yaml"
        
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        pr_benchmarks = config['global_settings']['ci_settings']['pr_benchmarks']
        
        self.assertIsInstance(pr_benchmarks, list)
        self.assertGreater(len(pr_benchmarks), 0)
        
        # Verify PR benchmarks exist in main benchmark list
        all_benchmarks = config['benchmarks'].keys()
        for benchmark_id in pr_benchmarks:
            self.assertIn(benchmark_id, all_benchmarks, 
                        f"PR benchmark {benchmark_id} not found in main benchmark list")
    
    def test_benchmark_timeouts_reasonable(self):
        """Test that benchmark timeouts are reasonable for CI"""
        config_path = benchmark_dir / "benchmark_config.yaml"
        
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        ci_settings = config['global_settings']['ci_settings']
        
        # Timeouts should be reasonable for CI
        self.assertIn('timeout_minutes', ci_settings)
        timeout = ci_settings['timeout_minutes']
        self.assertGreaterEqual(timeout, 5)   # At least 5 minutes
        self.assertLessEqual(timeout, 60)     # No more than 1 hour
    
    def test_performance_targets_defined(self):
        """Test that performance targets are defined"""
        config_path = benchmark_dir / "benchmark_config.yaml"
        
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        if 'performance_targets' in config:
            targets = config['performance_targets']
            
            # Check time targets
            if 'pr_benchmark_time' in targets:
                self.assertLessEqual(targets['pr_benchmark_time'], 600)  # 10 minutes max
            
            if 'nightly_benchmark_time' in targets:
                self.assertLessEqual(targets['nightly_benchmark_time'], 3600)  # 1 hour max


class TestBenchmarkDirectoryStructure(TestCase):
    """Test benchmark directory structure and files"""
    
    def test_benchmark_directories_exist(self):
        """Test that benchmark directories exist"""
        directories = ['assets', 'scripts', 'results']
        
        for dir_name in directories:
            dir_path = benchmark_dir / dir_name
            self.assertTrue(dir_path.exists(), f"Directory {dir_name} does not exist")
            self.assertTrue(dir_path.is_dir(), f"{dir_name} is not a directory")
    
    def test_benchmark_scripts_exist(self):
        """Test that required benchmark scripts exist"""
        scripts = ['benchmark_runner.py', 'rmsd_calculator.py']
        
        for script_name in scripts:
            script_path = benchmark_dir / "scripts" / script_name
            self.assertTrue(script_path.exists(), f"Script {script_name} does not exist")
            self.assertTrue(script_path.is_file(), f"{script_name} is not a file")
    
    def test_benchmark_config_accessible(self):
        """Test that benchmark configuration is accessible"""
        config_path = benchmark_dir / "benchmark_config.yaml"
        
        self.assertTrue(config_path.exists(), "benchmark_config.yaml does not exist")
        self.assertTrue(config_path.is_file(), "benchmark_config.yaml is not a file")
        
        # Test that it's readable
        with open(config_path, 'r') as f:
            content = f.read()
            self.assertGreater(len(content), 0, "benchmark_config.yaml is empty")


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
