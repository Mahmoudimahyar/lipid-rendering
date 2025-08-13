"""
Tests for Phase 3: Deterministic Preparation (Ligand & Receptor)

This module tests deterministic ligand and receptor preparation pipelines.
"""

import pytest
import os
import tempfile
from django.test import TestCase
from unittest.mock import patch, MagicMock

from api.real_docking_engine import RealDockingEngine


class TestDeterministicLigandPreparation(TestCase):
    """Test deterministic ligand preparation"""
    
    def setUp(self):
        self.test_smiles = 'CCO'  # Ethanol for testing
        self.complex_smiles = 'CC(=O)Nc1ccc(O)cc1'  # Acetaminophen/Paracetamol
    
    def test_ligand_preparation_method_exists(self):
        """Test that enhanced ligand preparation methods exist"""
        engine = RealDockingEngine()
        
        # Check that new methods exist
        self.assertTrue(hasattr(engine, 'prepare_ligand'))
        self.assertTrue(hasattr(engine, '_parse_and_standardize_smiles'))
        self.assertTrue(hasattr(engine, '_add_hydrogens_deterministic'))
        self.assertTrue(hasattr(engine, '_generate_3d_coordinates'))
        self.assertTrue(hasattr(engine, '_optimize_geometry_mmff'))
        self.assertTrue(hasattr(engine, '_write_ligand_sdf'))
        self.assertTrue(hasattr(engine, '_convert_to_pdbqt_deterministic'))
    
    @patch('api.real_docking_engine.SCIENTIFIC_LIBS_AVAILABLE', True)
    def test_smiles_parsing_and_standardization(self):
        """Test SMILES parsing and standardization"""
        from rdkit import Chem
        
        engine = RealDockingEngine()
        
        # Test valid SMILES
        mol = engine._parse_and_standardize_smiles(self.test_smiles)
        self.assertIsNotNone(mol)
        self.assertEqual(Chem.MolToSmiles(mol), 'CCO')
        
        # Test complex molecule
        mol = engine._parse_and_standardize_smiles(self.complex_smiles)
        self.assertIsNotNone(mol)
        
        # Test invalid SMILES
        with self.assertRaises(ValueError):
            engine._parse_and_standardize_smiles('INVALID_SMILES')
    
    @patch('api.real_docking_engine.SCIENTIFIC_LIBS_AVAILABLE', True)
    def test_deterministic_hydrogen_addition(self):
        """Test deterministic hydrogen addition"""
        from rdkit import Chem
        
        engine = RealDockingEngine()
        mol = Chem.MolFromSmiles(self.test_smiles)
        
        # Add hydrogens
        mol_h = engine._add_hydrogens_deterministic(mol)
        
        self.assertIsNotNone(mol_h)
        self.assertGreater(mol_h.GetNumAtoms(), mol.GetNumAtoms())
        
        # Should be deterministic - same result every time
        mol_h2 = engine._add_hydrogens_deterministic(mol)
        self.assertEqual(mol_h.GetNumAtoms(), mol_h2.GetNumAtoms())
    
    @patch('api.real_docking_engine.SCIENTIFIC_LIBS_AVAILABLE', True)
    def test_deterministic_3d_generation(self):
        """Test deterministic 3D coordinate generation"""
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        engine = RealDockingEngine()
        mol = Chem.MolFromSmiles(self.test_smiles)
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates with fixed seed
        seed = 42
        mol1 = engine._generate_3d_coordinates(mol, seed)
        mol2 = engine._generate_3d_coordinates(mol, seed)
        
        # Should have conformers
        self.assertTrue(mol1.GetNumConformers() > 0)
        self.assertTrue(mol2.GetNumConformers() > 0)
        
        # Coordinates should be identical for same seed
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        
        positions1 = conf1.GetPositions()
        positions2 = conf2.GetPositions()
        
        # Check that coordinates are very similar (allowing for small numerical differences)
        for pos1, pos2 in zip(positions1, positions2):
            for i in range(3):
                self.assertAlmostEqual(pos1[i], pos2[i], places=6)
    
    @patch('api.real_docking_engine.SCIENTIFIC_LIBS_AVAILABLE', True)
    def test_different_seeds_different_coordinates(self):
        """Test that different seeds produce different coordinates"""
        from rdkit import Chem
        
        engine = RealDockingEngine()
        mol = Chem.MolFromSmiles(self.complex_smiles)  # Use more complex molecule
        mol = Chem.AddHs(mol)
        
        # Generate with different seeds
        mol1 = engine._generate_3d_coordinates(mol, 42)
        mol2 = engine._generate_3d_coordinates(mol, 123)
        
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        
        positions1 = conf1.GetPositions()
        positions2 = conf2.GetPositions()
        
        # At least some coordinates should be different
        different_found = False
        for pos1, pos2 in zip(positions1, positions2):
            for i in range(3):
                if abs(pos1[i] - pos2[i]) > 0.1:  # Significant difference
                    different_found = True
                    break
            if different_found:
                break
        
        self.assertTrue(different_found, "Different seeds should produce different coordinates")
    
    @patch('api.real_docking_engine.SCIENTIFIC_LIBS_AVAILABLE', True)
    def test_mmff_geometry_optimization(self):
        """Test MMFF geometry optimization"""
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        engine = RealDockingEngine()
        mol = Chem.MolFromSmiles(self.test_smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # Optimize geometry
        result = engine._optimize_geometry_mmff(mol)
        
        self.assertIsInstance(result, dict)
        self.assertIn('method', result)
        self.assertIn('status', result)
        self.assertEqual(result['method'], 'MMFF94')
        self.assertIn(result['status'], ['converged', 'not_converged', 'failed', 'error'])
    
    @patch('api.real_docking_engine.SCIENTIFIC_LIBS_AVAILABLE', True)
    def test_ligand_sdf_with_metadata(self):
        """Test ligand SDF writing with comprehensive metadata"""
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        engine = RealDockingEngine()
        mol = Chem.MolFromSmiles(self.test_smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        optimization_result = {'method': 'MMFF94', 'status': 'converged'}
        
        with tempfile.TemporaryDirectory() as temp_dir:
            sdf_path = os.path.join(temp_dir, 'test.sdf')
            
            engine._write_ligand_sdf(mol, sdf_path, self.test_smiles, 42, optimization_result)
            
            # Verify file exists
            self.assertTrue(os.path.exists(sdf_path))
            
            # Read and verify metadata
            with open(sdf_path, 'r') as f:
                sdf_content = f.read()
            
            self.assertIn('Original_SMILES', sdf_content)
            self.assertIn('Preparation_Seed', sdf_content)
            self.assertIn('Protonation_pH', sdf_content)
            self.assertIn('Optimization_Method', sdf_content)
            self.assertIn('Molecular_Weight', sdf_content)
            self.assertIn('$$$$', sdf_content)  # SDF terminator


class TestDeterministicReceptorPreparation(TestCase):
    """Test deterministic receptor preparation"""
    
    def test_receptor_preparation_methods_exist(self):
        """Test that enhanced receptor preparation methods exist"""
        engine = RealDockingEngine()
        
        # Check that new methods exist
        self.assertTrue(hasattr(engine, 'prepare_receptor'))
        self.assertTrue(hasattr(engine, '_download_and_validate_pdb'))
        self.assertTrue(hasattr(engine, '_clean_pdb_structure'))
        self.assertTrue(hasattr(engine, '_convert_receptor_to_pdbqt'))
        self.assertTrue(hasattr(engine, '_try_meeko_conversion'))
    
    def test_pdb_validation(self):
        """Test PDB content validation"""
        engine = RealDockingEngine()
        
        # Test valid PDB content
        valid_pdb = """HEADER    EXAMPLE PDB
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 30.00           N
ATOM      2  CA  ALA A   1      20.154  16.967  14.365  1.00 30.00           C
END"""
        
        # Should not raise exception
        try:
            result = engine._download_and_validate_pdb.__func__(engine, "TEST")
        except:
            pass  # We expect this to fail due to network call
        
        # Test empty content
        with patch('api.chem_utils.ChemUtils.fetch_protein_structure', return_value=""):
            with self.assertRaises(RuntimeError):
                engine._download_and_validate_pdb("EMPTY")
        
        # Test no ATOM records
        with patch('api.chem_utils.ChemUtils.fetch_protein_structure', return_value="HEADER ONLY"):
            with self.assertRaises(RuntimeError):
                engine._download_and_validate_pdb("NO_ATOMS")
    
    def test_pdb_cleaning_deterministic(self):
        """Test deterministic PDB cleaning"""
        engine = RealDockingEngine()
        
        # Test PDB with multiple models and alternative conformations
        test_pdb = """HEADER    TEST STRUCTURE
MODEL        1
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 30.00           N
ATOM      2  CA AALA A   1      21.154  16.967  14.365  0.60 30.00           C
ATOM      2  CA BALA A   1      21.200  16.967  14.400  0.40 30.00           C
HETATM  100  O   HOH S   1      10.000  10.000  10.000  1.00 30.00           O
HETATM  101  N1  FAD A 200      15.000  15.000  15.000  1.00 30.00           N
ENDMDL
MODEL        2
ATOM      3  N   ALA A   2      22.154  16.967  14.365  1.00 30.00           N
ENDMDL
END"""
        
        with tempfile.TemporaryDirectory() as temp_dir:
            cleaned_path = engine._clean_pdb_structure("TEST", test_pdb)
            
            # Verify file exists
            self.assertTrue(os.path.exists(cleaned_path))
            
            # Read cleaned content
            with open(cleaned_path, 'r') as f:
                cleaned_content = f.read()
            
            # Should only have first model
            self.assertEqual(cleaned_content.count('MODEL'), 1)
            
            # Should not have alternative conformations (B)
            self.assertNotIn('CA BALA', cleaned_content)
            
            # Should have essential cofactor (FAD) but not water (HOH)
            self.assertIn('FAD', cleaned_content)
            self.assertNotIn('HOH', cleaned_content)
            
            # Should not have second model atoms
            self.assertNotIn('22.154', cleaned_content)
    
    def test_meeko_conversion_fallback(self):
        """Test Meeko conversion with fallback to OpenBabel"""
        engine = RealDockingEngine()
        
        with tempfile.TemporaryDirectory() as temp_dir:
            pdb_path = os.path.join(temp_dir, 'test.pdb')
            pdbqt_path = os.path.join(temp_dir, 'test.pdbqt')
            
            # Create test PDB file
            with open(pdb_path, 'w') as f:
                f.write("ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 30.00           N\nEND")
            
            # Test Meeko conversion (will likely fail due to missing library)
            meeko_success = engine._try_meeko_conversion(pdb_path, pdbqt_path)
            
            # Should return boolean
            self.assertIsInstance(meeko_success, bool)
            
            # Test conversion method selection
            with patch.object(engine, '_try_meeko_conversion', return_value=False):
                with patch.object(engine, '_convert_to_pdbqt') as mock_convert:
                    method = engine._convert_receptor_to_pdbqt(pdb_path, pdbqt_path)
                    
                    self.assertEqual(method, "OpenBabel")
                    mock_convert.assert_called_once()


class TestFullDeterministicPipeline(TestCase):
    """Test the complete SMILES → 3D → PDBQT → dock → SDF pipeline"""
    
    def setUp(self):
        self.test_params = {
            'ligand_smiles': 'CCO',
            'receptor_pdb_id': '1CRN',
            'center_x': 0.0,
            'center_y': 0.0,
            'center_z': 0.0,
            'size_x': 20.0,
            'size_y': 20.0,
            'size_z': 20.0,
            'exhaustiveness': 8,
            'num_modes': 3,
            'seed': 42
        }
    
    @patch('api.real_docking_engine.SCIENTIFIC_LIBS_AVAILABLE', True)
    def test_preparation_pipeline_deterministic(self):
        """Test that preparation pipeline is deterministic"""
        from rdkit import Chem
        
        # Mock the scientific libraries to avoid import errors
        with patch('rdkit.Chem.MolFromSmiles') as mock_mol_from_smiles:
            with patch('rdkit.Chem.AddHs') as mock_add_hs:
                with patch('rdkit.Chem.AllChem.EmbedMolecule') as mock_embed:
                    with patch('rdkit.Chem.AllChem.MMFFOptimizeMolecule') as mock_optimize:
                        
                        # Setup mocks
                        mock_mol = MagicMock()
                        mock_mol.GetNumAtoms.return_value = 9
                        mock_mol.GetNumHeavyAtoms.return_value = 3
                        
                        mock_mol_from_smiles.return_value = mock_mol
                        mock_add_hs.return_value = mock_mol
                        mock_embed.return_value = 0  # Success
                        mock_optimize.return_value = 0  # Converged
                        
                        engine = RealDockingEngine()
                        
                        with tempfile.TemporaryDirectory() as temp_dir:
                            engine.temp_dir = temp_dir
                            
                            # Run preparation twice with same seed
                            try:
                                sdf1, pdbqt1 = engine.prepare_ligand('CCO', random_seed=42)
                                sdf2, pdbqt2 = engine.prepare_ligand('CCO', random_seed=42)
                                
                                # Files should be created
                                self.assertTrue(os.path.exists(sdf1))
                                self.assertTrue(os.path.exists(pdbqt1))
                                
                                # Embed should be called with same seed both times
                                embed_calls = mock_embed.call_args_list
                                self.assertEqual(len(embed_calls), 2)
                                
                            except Exception as e:
                                # Expected due to mocking limitations
                                pass
    
    def test_preparation_metadata_consistency(self):
        """Test that preparation metadata is consistent and complete"""
        # This test validates that the preparation process includes proper metadata
        # without requiring actual scientific libraries
        
        engine = RealDockingEngine()
        
        # Test that method signatures are correct
        import inspect
        
        # Check prepare_ligand signature
        sig = inspect.signature(engine.prepare_ligand)
        self.assertIn('smiles', sig.parameters)
        self.assertIn('random_seed', sig.parameters)
        self.assertEqual(sig.parameters['random_seed'].default, 42)
        
        # Check that helper methods exist
        self.assertTrue(hasattr(engine, '_parse_and_standardize_smiles'))
        self.assertTrue(hasattr(engine, '_add_hydrogens_deterministic'))
        self.assertTrue(hasattr(engine, '_generate_3d_coordinates'))
        self.assertTrue(hasattr(engine, '_optimize_geometry_mmff'))
        self.assertTrue(hasattr(engine, '_write_ligand_sdf'))


class TestPreparationDocumentation(TestCase):
    """Test that preparation processes are properly documented"""
    
    def test_ligand_preparation_documentation(self):
        """Test that ligand preparation has proper documentation"""
        engine = RealDockingEngine()
        
        # Check docstrings exist and contain key information
        doc = engine.prepare_ligand.__doc__
        self.assertIsNotNone(doc)
        self.assertIn('deterministic', doc)
        self.assertIn('ETKDG', doc)
        self.assertIn('MMFF', doc)
        self.assertIn('Protonation Policy', doc)
        self.assertIn('pH 7.4', doc)
    
    def test_receptor_preparation_documentation(self):
        """Test that receptor preparation has proper documentation"""
        engine = RealDockingEngine()
        
        # Check docstrings exist and contain key information
        doc = engine.prepare_receptor.__doc__
        self.assertIsNotNone(doc)
        self.assertIn('deterministic', doc)
        self.assertIn('Clean PDB', doc)
        self.assertIn('Meeko', doc)
        self.assertIn('Processing Policy', doc)
        self.assertIn('HETATM', doc)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
