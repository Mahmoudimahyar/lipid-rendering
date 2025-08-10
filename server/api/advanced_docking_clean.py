"""
Advanced docking utilities including GNINA rescoring and pocket detection
"""

import json
import time
import random
import math
import subprocess
import tempfile
import os
from typing import Dict, List, Any, Optional, Tuple
from django.conf import settings
from .models import BindingPocket, JobTemplate


class GNINAScorer:
    """GNINA neural network-based scoring and rescoring utilities"""
    
    @staticmethod
    def rescore_poses(poses: List[Dict], ligand_sdf: str, receptor_pdbqt: str) -> List[Dict]:
        """
        Rescore docking poses using GNINA neural network scoring
        This is a mock implementation - real version would call GNINA
        """
        rescored_poses = []
        
        for pose in poses:
            # Mock GNINA rescoring - in reality this would:
            # 1. Write pose to SDF file
            # 2. Call GNINA with --score_only flag
            # 3. Parse CNN scores and other metrics
            
            # Simulate neural network scoring with slight variations
            original_affinity = pose.get('affinity', 0)
            
            # GNINA tends to be more accurate, so add some realistic variation
            gnina_score = original_affinity + random.uniform(-0.5, 0.5)
            cnn_score = random.uniform(0.1, 0.9)  # CNN probability
            cnn_affinity = random.uniform(-12.0, -4.0)  # CNN-predicted affinity
            
            # Additional GNINA-specific scores
            gnina_metrics = {
                'gnina_score': round(gnina_score, 3),
                'cnn_score': round(cnn_score, 4),
                'cnn_affinity': round(cnn_affinity, 3),
                'cnn_vs_refined_aff_diff': round(abs(cnn_affinity - original_affinity), 3),
                'cnn_vs_refined_aff_rank_diff': random.randint(-2, 2),
                'gauss1': random.uniform(-0.5, 0.1),
                'gauss2': random.uniform(-0.3, 0.2),
                'repulsion': random.uniform(0.0, 0.8),
                'hydrophobic': random.uniform(-0.8, 0.0),
                'hydrogen': random.uniform(-0.6, 0.0)
            }
            
            # Create rescored pose
            rescored_pose = pose.copy()
            rescored_pose.update(gnina_metrics)
            rescored_pose['rescoring_method'] = 'GNINA'
            rescored_pose['original_affinity'] = original_affinity
            
            rescored_poses.append(rescored_pose)
        
        # Sort by GNINA score (more negative = better)
        rescored_poses.sort(key=lambda x: x.get('gnina_score', 0))
        
        return rescored_poses
    
    @staticmethod
    def get_pose_features(pose: Dict, ligand_sdf: str, receptor_pdbqt: str) -> Dict:
        """
        Extract detailed molecular features for a pose using GNINA
        """
        # Mock feature extraction - real implementation would use GNINA's
        # feature extraction capabilities
        
        features = {
            'molecular_descriptors': {
                'mol_weight': random.uniform(200, 600),
                'logp': random.uniform(-2, 6),
                'rotatable_bonds': random.randint(0, 15),
                'hbd': random.randint(0, 8),  # Hydrogen bond donors
                'hba': random.randint(0, 12),  # Hydrogen bond acceptors
                'tpsa': random.uniform(20, 200)  # Topological polar surface area
            },
            'protein_ligand_interactions': {
                'hydrogen_bonds': random.randint(0, 6),
                'hydrophobic_contacts': random.randint(2, 15),
                'pi_pi_stacking': random.randint(0, 3),
                'salt_bridges': random.randint(0, 2),
                'van_der_waals': random.randint(5, 25)
            },
            'geometric_features': {
                'buried_surface_area': random.uniform(300, 800),
                'ligand_efficiency': random.uniform(0.2, 0.6),
                'contact_surface': random.uniform(400, 1000),
                'shape_complementarity': random.uniform(0.4, 0.8)
            },
            'pharmacophore_features': {
                'aromatic_rings': random.randint(0, 4),
                'basic_groups': random.randint(0, 3),
                'acidic_groups': random.randint(0, 2),
                'hydrophobic_centers': random.randint(1, 6)
            }
        }
        
        return features


class PocketDetector:
    """Binding pocket detection and analysis utilities"""
    
    @staticmethod
    def detect_pockets(pdb_id: str, pdb_content: str) -> List[Dict]:
        """
        Detect binding pockets in a protein structure
        This is a mock implementation - real version would use fpocket, CASTp, etc.
        """
        # Mock pocket detection - in reality this would:
        # 1. Write PDB to temporary file
        # 2. Run fpocket or similar tool
        # 3. Parse pocket detection results
        
        num_pockets = random.randint(2, 8)
        pockets = []
        
        for i in range(num_pockets):
            # Generate realistic pocket properties
            center_x = random.uniform(-20, 20)
            center_y = random.uniform(-20, 20)
            center_z = random.uniform(-20, 20)
            
            pocket = {
                'pocket_number': i + 1,
                'center_x': round(center_x, 3),
                'center_y': round(center_y, 3),
                'center_z': round(center_z, 3),
                'radius': round(random.uniform(8, 25), 2),
                'volume': round(random.uniform(200, 2000), 1),
                'surface_area': round(random.uniform(400, 3000), 1),
                'druggability_score': round(random.uniform(0.3, 1.0), 3),
                'hydrophobicity': round(random.uniform(0.2, 0.8), 3),
                'polarity': round(random.uniform(0.1, 0.7), 3),
                'confidence_score': round(random.uniform(0.6, 0.95), 3),
                'detection_method': 'fpocket',
                'residues': PocketDetector._generate_pocket_residues(),
                'properties': {
                    'is_druggable': random.choice([True, False]),
                    'cavity_type': random.choice(['deep', 'shallow', 'tunnel', 'cleft']),
                    'accessibility': random.choice(['buried', 'surface', 'partially_buried']),
                    'conservation_score': round(random.uniform(0.0, 1.0), 3)
                }
            }
            
            pockets.append(pocket)
        
        # Sort by druggability score (highest first)
        pockets.sort(key=lambda x: x['druggability_score'], reverse=True)
        
        return pockets
    
    @staticmethod
    def _generate_pocket_residues() -> List[Dict]:
        """Generate mock residue list for pocket"""
        amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 
                      'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                      'THR', 'TRP', 'TYR', 'VAL']
        
        num_residues = random.randint(10, 30)
        residues = []
        
        for i in range(num_residues):
            residue = {
                'residue_name': random.choice(amino_acids),
                'residue_number': random.randint(50, 300),
                'chain_id': random.choice(['A', 'B']),
                'distance_to_center': round(random.uniform(3.0, 15.0), 2),
                'contact_type': random.choice(['direct', 'water_mediated', 'van_der_waals'])
            }
            residues.append(residue)
        
        return residues
    
    @staticmethod
    def analyze_pocket_druggability(pocket: Dict) -> Dict:
        """
        Analyze pocket druggability using various descriptors
        """
        volume = pocket.get('volume', 0)
        surface_area = pocket.get('surface_area', 0)
        
        # Calculate druggability metrics
        volume_sa_ratio = volume / surface_area if surface_area > 0 else 0
        
        analysis = {
            'druggability_class': 'high' if pocket.get('druggability_score', 0) > 0.7 else 
                                 'medium' if pocket.get('druggability_score', 0) > 0.4 else 'low',
            'volume_category': 'large' if volume > 1000 else 'medium' if volume > 500 else 'small',
            'shape_analysis': {
                'volume_surface_ratio': round(volume_sa_ratio, 3),
                'sphericity': round(random.uniform(0.3, 0.9), 3),
                'elongation': round(random.uniform(0.1, 0.8), 3)
            },
            'chemical_environment': {
                'hydrophobic_percentage': round(pocket.get('hydrophobicity', 0) * 100, 1),
                'polar_percentage': round(pocket.get('polarity', 0) * 100, 1),
                'charged_residues': random.randint(2, 8),
                'aromatic_residues': random.randint(1, 5)
            },
            'recommendations': PocketDetector._get_druggability_recommendations(pocket)
        }
        
        return analysis
    
    @staticmethod
    def _get_druggability_recommendations(pocket: Dict) -> List[str]:
        """Generate druggability recommendations based on pocket properties"""
        recommendations = []
        
        volume = pocket.get('volume', 0)
        druggability = pocket.get('druggability_score', 0)
        
        if druggability > 0.8:
            recommendations.append("Excellent druggability - high priority target")
        elif druggability > 0.6:
            recommendations.append("Good druggability - suitable for drug development")
        else:
            recommendations.append("Low druggability - consider fragment-based approaches")
        
        if volume < 300:
            recommendations.append("Small pocket - consider fragment screening")
        elif volume > 1500:
            recommendations.append("Large pocket - may require larger molecules or protein-protein interaction inhibitors")
        
        if pocket.get('hydrophobicity', 0) > 0.7:
            recommendations.append("Highly hydrophobic pocket - lipophilic compounds preferred")
        
        return recommendations


class AdvancedDockingEngine:
    """Enhanced docking engine with advanced features"""
    
    @staticmethod
    def run_advanced_docking(params: Dict, use_pocket_detection: bool = False, 
                           use_gnina_rescoring: bool = False) -> Dict:
        """
        Run docking with advanced features
        """
        # Run basic docking first
        from .docking_utils import DockingEngine
        basic_results = DockingEngine.run_mock_docking(params)
        
        if not basic_results.get('success'):
            return basic_results
        
        advanced_results = basic_results.copy()
        
        # Add pocket detection if requested
        if use_pocket_detection:
            pdb_id = params.get('receptor_pdb_id', '')
            if pdb_id:
                pockets = PocketDetector.detect_pockets(pdb_id, "mock_pdb_content")
                advanced_results['detected_pockets'] = pockets
                
                # Analyze how well docking site matches detected pockets
                site_analysis = AdvancedDockingEngine._analyze_site_vs_pockets(params, pockets)
                advanced_results['site_analysis'] = site_analysis
        
        # Add GNINA rescoring if requested  
        if use_gnina_rescoring:
            poses = advanced_results.get('poses', [])
            if poses:
                rescored_poses = GNINAScorer.rescore_poses(poses, "mock_sdf", "mock_pdbqt")
                advanced_results['rescored_poses'] = rescored_poses
                
                # Add pose features for best poses
                for pose in rescored_poses[:3]:  # Top 3 poses
                    pose['molecular_features'] = GNINAScorer.get_pose_features(
                        pose, "mock_sdf", "mock_pdbqt"
                    )
                
                # Update summary with rescoring info
                if rescored_poses:
                    best_rescored = min(rescored_poses, key=lambda x: x.get('gnina_score', 0))
                    advanced_results['summary']['best_gnina_score'] = best_rescored.get('gnina_score')
                    advanced_results['summary']['rescoring_method'] = 'GNINA'
        
        # Add computational performance metrics
        advanced_results['performance'] = {
            'cpu_time': round(random.uniform(30, 300), 2),
            'memory_usage': random.randint(100, 500),
            'advanced_features_used': {
                'pocket_detection': use_pocket_detection,
                'gnina_rescoring': use_gnina_rescoring
            }
        }
        
        return advanced_results
    
    @staticmethod
    def _analyze_site_vs_pockets(docking_params: Dict, detected_pockets: List[Dict]) -> Dict:
        """
        Analyze how well the docking site matches detected pockets
        """
        site_center = (
            docking_params.get('center_x', 0),
            docking_params.get('center_y', 0), 
            docking_params.get('center_z', 0)
        )
        
        closest_pocket = None
        min_distance = float('inf')
        
        for pocket in detected_pockets:
            pocket_center = (pocket['center_x'], pocket['center_y'], pocket['center_z'])
            distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(site_center, pocket_center)))
            
            if distance < min_distance:
                min_distance = distance
                closest_pocket = pocket
        
        analysis = {
            'closest_pocket_distance': round(min_distance, 2),
            'site_overlaps_pocket': min_distance < 10.0,
            'closest_pocket': closest_pocket,
            'recommendation': (
                "Docking site well aligned with detected pocket" if min_distance < 5.0 else
                "Docking site partially overlaps pocket - consider adjustment" if min_distance < 10.0 else
                "Docking site may not correspond to natural binding pocket"
            )
        }
        
        return analysis


class JobTemplateManager:
    """Manager for job templates and presets"""
    
    @staticmethod
    def create_default_templates():
        """Create default job templates"""
        templates = [
            {
                'name': 'Standard Docking',
                'description': 'Standard molecular docking with default parameters',
                'category': 'basic',
                'is_public': True,
                'default_params': {
                    'exhaustiveness': 8,
                    'num_modes': 9,
                    'energy_range': 3.0
                },
                'advanced_params': {
                    'use_gnina_rescoring': False,
                    'use_pocket_detection': False
                }
            },
            {
                'name': 'High Precision Docking',
                'description': 'Thorough docking with increased exhaustiveness',
                'category': 'precision',
                'is_public': True,
                'default_params': {
                    'exhaustiveness': 16,
                    'num_modes': 20,
                    'energy_range': 5.0
                },
                'advanced_params': {
                    'use_gnina_rescoring': True,
                    'use_pocket_detection': True
                }
            },
            {
                'name': 'Fragment Screening',
                'description': 'Optimized for small molecule fragments',
                'category': 'fragments',
                'is_public': True,
                'default_params': {
                    'exhaustiveness': 12,
                    'num_modes': 15,
                    'energy_range': 4.0
                },
                'advanced_params': {
                    'use_gnina_rescoring': True,
                    'use_pocket_detection': True,
                    'fragment_mode': True
                }
            },
            {
                'name': 'Virtual Screening',
                'description': 'Fast docking for large compound libraries',
                'category': 'screening',
                'is_public': True,
                'default_params': {
                    'exhaustiveness': 4,
                    'num_modes': 5,
                    'energy_range': 2.0
                },
                'advanced_params': {
                    'use_gnina_rescoring': False,
                    'use_pocket_detection': False,
                    'fast_mode': True
                }
            }
        ]
        
        created_templates = []
        for template_data in templates:
            template, created = JobTemplate.objects.get_or_create(
                name=template_data['name'],
                defaults=template_data
            )
            if created:
                created_templates.append(template)
        
        return created_templates
    
    @staticmethod
    def get_template_by_category(category: str) -> List[JobTemplate]:
        """Get templates by category"""
        return JobTemplate.objects.filter(category=category, is_public=True)
    
    @staticmethod
    def apply_template(template: JobTemplate, custom_params: Dict = None) -> Dict:
        """Apply template to generate docking parameters"""
        params = template.default_params.copy()
        
        if template.advanced_params:
            params.update(template.advanced_params)
        
        if custom_params:
            params.update(custom_params)
        
        # Increment usage counter
        template.increment_usage()
        
        return params
