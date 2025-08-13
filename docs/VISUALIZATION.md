# Molecular Visualization in Lipid Rendering

This document explains the molecular visualization pipeline, including pose extraction from AutoDock Vina, PDBQT processing with Meeko, and 3D rendering with 3Dmol.js.

## Overview

Lipid Rendering provides interactive 3D molecular visualization that transforms docking results into rich, scientifically accurate visual representations. The pipeline handles the complete flow from Vina PDBQT output to browser-based 3D visualization.

## Visualization Pipeline

```
AutoDock Vina → PDBQT Poses → SDF Conversion → 3Dmol.js → Interactive 3D
```

### 1. Pose Extraction from Vina

**PDBQT Output Processing**

AutoDock Vina outputs multiple docked poses in PDBQT format, ranked by binding affinity:

```pdbqt
MODEL 1
REMARK VINA RESULT:    -8.5      0.000      0.000
ATOM      1  C1  LIG A   1      22.500   4.000  24.000  1.00  0.00     0.123 C 
ATOM      2  C2  LIG A   1      23.500   4.000  24.000  1.00  0.00     0.076 C 
...
ENDMDL
MODEL 2
REMARK VINA RESULT:    -7.8      1.245      2.876
...
```

**Extraction Process**

```python
def _extract_pose_as_sdf(self, pdbqt_content: str, pose_index: int):
    """
    Extract individual pose from Vina PDBQT output
    
    Steps:
    1. Parse PDBQT by MODEL blocks
    2. Extract affinity from VINA RESULT comment
    3. Convert coordinates to SDF format
    4. Add metadata (affinity, RMSD, pose index)
    """
    models = re.split(r'MODEL\s+\d+', pdbqt_content)[1:]  # Skip header
    pose_block = models[pose_index]
    
    # Extract Vina result line
    vina_result = re.search(r'REMARK VINA RESULT:\s+([+-]?\d+\.?\d*)', pose_block)
    affinity = float(vina_result.group(1))
    
    # Convert coordinates to SDF format
    sdf_content = self._convert_pdbqt_to_sdf(pose_block)
    
    return {
        'mode': pose_index + 1,
        'affinity': affinity,
        'sdf': sdf_content
    }
```

### 2. PDBQT to SDF Conversion

**Coordinate Transformation**

The system converts PDBQT coordinate format to SDF for better visualization support:

```python
def _convert_pdbqt_to_sdf(self, pdbqt_block: str) -> str:
    """
    Convert PDBQT coordinates to SDF format
    
    PDBQT format:
    ATOM      1  C1  LIG A   1      22.500   4.000  24.000  1.00  0.00     0.123 C
    
    SDF format:
    22.5000    4.0000   24.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    """
    atoms = []
    bonds = []
    
    for line in pdbqt_block.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Parse PDBQT atom line
            x = float(line[30:38])
            y = float(line[38:46]) 
            z = float(line[46:54])
            element = line[77:79].strip()
            
            atoms.append({
                'x': x, 'y': y, 'z': z,
                'element': element,
                'charge': 0
            })
    
    # Generate SDF format
    sdf_content = self._create_sdf_from_atoms(atoms, bonds)
    return sdf_content
```

**Metadata Embedding**

SDF files include comprehensive metadata for each pose:

```sdf
  RDKit          3D

  15 14  0  0  0  0  0  0  0  0999 V2000
   22.5000    4.0000   24.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   23.5000    4.0000   24.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   ...
M  END
> <AFFINITY>
-8.5

> <POSE_INDEX>
1

> <RMSD_LB>
0.000

> <RMSD_UB>
0.000

> <DOCKING_PARAMETERS>
{"seed": 42, "exhaustiveness": 8}

$$$$
```

### 3. Meeko Integration

**Preferred PDBQT Processing**

Meeko is the preferred tool for PDBQT preparation and processing:

```python
def _try_meeko_conversion(self, pdb_path: str, pdbqt_path: str) -> bool:
    """
    Use Meeko for high-quality PDBQT conversion
    """
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        
        # Prepare molecule with Meeko
        preparator = MoleculePreparation()
        preparator.prepare(pdb_path)
        
        # Write PDBQT
        writer = PDBQTWriterLegacy()
        writer.write_pdbqt_string(preparator.mol)
        
        with open(pdbqt_path, 'w') as f:
            f.write(writer.pdbqt_string)
        
        return True
        
    except ImportError:
        logger.warning("Meeko not available, using OpenBabel fallback")
        return False
```

**OpenBabel Fallback**

When Meeko is unavailable, OpenBabel provides compatibility:

```python
def _convert_to_pdbqt(self, input_path: str, output_path: str, is_ligand: bool = True):
    """
    Convert molecular structure to PDBQT using OpenBabel
    """
    cmd = ['obabel', input_path, '-O', output_path]
    
    if is_ligand:
        cmd.extend(['-h'])  # Add hydrogens
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        logger.info(f"OpenBabel conversion successful: {output_path}")
    except subprocess.CalledProcessError as e:
        logger.error(f"OpenBabel conversion failed: {e}")
        raise
```

## 3Dmol.js Rendering

### Frontend Visualization Components

**Enhanced3DViewer Component**

The main visualization component handles multiple poses and interaction modes:

```javascript
import $3Dmol from '3dmol';

const Enhanced3DViewer = ({ poses, receptorData, onPoseSelect }) => {
  const [viewer, setViewer] = useState(null);
  const [currentPose, setCurrentPose] = useState(0);
  
  useEffect(() => {
    // Initialize 3Dmol.js viewer
    const viewerElement = document.getElementById('3dmol-viewer');
    const molviewer = $3Dmol.createViewer(viewerElement, {
      backgroundColor: 'white',
      antialias: true,
      camerax: 0,
      cameray: 0,
      cameraz: 100
    });
    
    setViewer(molviewer);
  }, []);
```

**SDF Loading and Display**

```javascript
const loadPose = useCallback((poseIndex) => {
  if (!viewer || !poses[poseIndex]) return;
  
  // Clear previous pose
  viewer.removeAllModels();
  
  // Load receptor (if available)
  if (receptorData) {
    viewer.addModel(receptorData, 'pdb');
    viewer.setStyle({}, {
      cartoon: { color: 'lightgray', opacity: 0.8 },
      stick: { hidden: true }
    });
  }
  
  // Load ligand pose from SDF
  const sdf = poses[poseIndex].sdf;
  viewer.addModel(sdf, 'sdf');
  
  // Style ligand with pose-specific coloring
  const poseColor = getPoseColor(poseIndex);
  viewer.setStyle({model: -1}, {
    stick: { 
      colorscheme: poseColor,
      radius: 0.15 
    },
    sphere: { 
      colorscheme: poseColor,
      radius: 0.4 
    }
  });
  
  // Render and center
  viewer.render();
  viewer.zoomTo();
}, [viewer, poses, receptorData]);
```

### Visualization Features

**Multi-Pose Display**

Users can cycle through different docked poses:

```javascript
const PoseSelector = ({ poses, currentPose, onPoseChange }) => (
  <div className="pose-selector">
    <h3>Docked Poses</h3>
    {poses.map((pose, index) => (
      <button
        key={index}
        className={`pose-button ${index === currentPose ? 'active' : ''}`}
        onClick={() => onPoseChange(index)}
      >
        <div className="pose-info">
          <span className="pose-number">Pose {pose.mode}</span>
          <span className="affinity">{pose.affinity.toFixed(1)} kcal/mol</span>
        </div>
      </button>
    ))}
  </div>
);
```

**Interactive Controls**

```javascript
const ViewControls = ({ viewer, onReset, onSaveImage }) => (
  <div className="view-controls">
    <button onClick={() => viewer.rotate(90, 'y')}>
      Rotate Y
    </button>
    <button onClick={() => viewer.zoomTo()}>
      Center
    </button>
    <button onClick={onReset}>
      Reset View
    </button>
    <button onClick={() => exportViewer(viewer)}>
      Save Image
    </button>
  </div>
);
```

**Styling and Coloring**

```javascript
const applyVisualizationStyle = (viewer, styleMode) => {
  switch (styleMode) {
    case 'stick':
      viewer.setStyle({}, {
        stick: { colorscheme: 'default', radius: 0.15 }
      });
      break;
      
    case 'sphere':
      viewer.setStyle({}, {
        sphere: { colorscheme: 'default', radius: 0.4 }
      });
      break;
      
    case 'cartoon':
      viewer.setStyle({}, {
        cartoon: { color: 'spectrum' }
      });
      break;
      
    case 'surface':
      viewer.addSurface('VDW', {
        opacity: 0.7,
        colorscheme: 'default'
      });
      break;
  }
  
  viewer.render();
};
```

## Advanced Visualization Features

### Binding Site Highlighting

**Pocket Visualization**

```javascript
const highlightBindingSite = (viewer, bindingSiteResidues) => {
  // Highlight binding site residues
  viewer.setStyle(
    { resi: bindingSiteResidues },
    {
      stick: { colorscheme: 'orange', radius: 0.2 },
      cartoon: { color: 'orange', opacity: 0.9 }
    }
  );
  
  // Add binding site surface
  viewer.addSurface('VDW', {
    opacity: 0.3,
    color: 'orange'
  }, { resi: bindingSiteResidues });
  
  viewer.render();
};
```

### Interaction Analysis

**Hydrogen Bond Display**

```javascript
const showHydrogenBonds = (viewer, ligandSel, receptorSel) => {
  // Find potential hydrogen bonds
  const hbonds = viewer.findHydrogenBonds(ligandSel, receptorSel, {
    maxDistance: 3.5,
    angleThreshold: 30
  });
  
  // Display hydrogen bonds as dashed lines
  hbonds.forEach(bond => {
    viewer.addLine({
      start: bond.donor,
      end: bond.acceptor,
      color: 'yellow',
      linewidth: 2,
      dashed: true
    });
  });
  
  viewer.render();
};
```

### Export Capabilities

**Image Export**

```javascript
const exportViewer = (viewer, options = {}) => {
  const defaultOptions = {
    format: 'png',
    width: 1920,
    height: 1080,
    antialias: true
  };
  
  const finalOptions = { ...defaultOptions, ...options };
  
  viewer.pngURI((uri) => {
    const link = document.createElement('a');
    link.href = uri;
    link.download = `molecular_view_${Date.now()}.png`;
    link.click();
  }, finalOptions);
};
```

**SDF Export**

```javascript
const exportCurrentPose = (poses, currentPose) => {
  const sdf = poses[currentPose].sdf;
  const blob = new Blob([sdf], { type: 'chemical/x-mdl-sdfile' });
  const url = URL.createObjectURL(blob);
  
  const link = document.createElement('a');
  link.href = url;
  link.download = `pose_${currentPose + 1}.sdf`;
  link.click();
  
  URL.revokeObjectURL(url);
};
```

## Performance Optimization

### Efficient Rendering

**Model Management**

```javascript
const optimizeViewer = (viewer) => {
  // Enable performance optimizations
  viewer.setBackgroundColor('white', 0);
  viewer.setProjection('perspective');
  
  // Optimize for large molecules
  viewer.setStyle({}, {
    stick: { 
      radius: 0.1,  // Smaller radius for performance
      quality: 20   // Lower quality for speed
    }
  });
  
  // Use level-of-detail for complex structures
  viewer.enableFog(true);
  viewer.setFogNear(50);
  viewer.setFogFar(100);
};
```

### Memory Management

**Cleanup Strategies**

```javascript
const cleanupViewer = (viewer) => {
  // Clear all models and surfaces
  viewer.removeAllModels();
  viewer.removeAllSurfaces();
  viewer.removeAllLabels();
  viewer.removeAllShapes();
  
  // Force garbage collection
  viewer.render();
};
```

## Testing and Validation

### Visual Regression Testing

**Snapshot Testing**

```javascript
describe('Molecular Visualization', () => {
  test('renders poses correctly', async () => {
    const { container } = render(
      <Enhanced3DViewer poses={mockPoses} />
    );
    
    // Wait for 3Dmol.js to render
    await waitFor(() => {
      expect(container.querySelector('canvas')).toBeInTheDocument();
    });
    
    // Snapshot test for visual regression
    expect(container).toMatchSnapshot();
  });
});
```

### SDF Validation

**Structure Validation**

```python
def test_sdf_validity():
    """Test that generated SDF files are valid"""
    from rdkit import Chem
    
    poses = extract_poses_from_vina_output(pdbqt_content)
    
    for pose in poses:
        sdf = pose['sdf']
        
        # Parse with RDKit to validate
        mol = Chem.MolFromMolBlock(sdf)
        assert mol is not None, "Invalid SDF structure"
        
        # Check basic properties
        assert mol.GetNumAtoms() > 0
        assert mol.GetNumBonds() >= 0
```

## Configuration and Customization

### 3Dmol.js Configuration

```javascript
const viewerConfig = {
  backgroundColor: 'white',
  antialias: true,
  camerax: 0,
  cameray: 0, 
  cameraz: 100,
  enableShiftDragZoom: true,
  enableFog: false,
  fog: {
    near: 50,
    far: 100,
    color: 'white'
  }
};
```

### Styling Presets

```javascript
const stylePresets = {
  publication: {
    stick: { radius: 0.2, quality: 40 },
    sphere: { radius: 0.4 },
    surface: { opacity: 0.8 },
    backgroundColor: 'white'
  },
  
  interactive: {
    stick: { radius: 0.15, quality: 20 },
    sphere: { radius: 0.3 },
    surface: { opacity: 0.6 },
    backgroundColor: '#f0f0f0'
  },
  
  presentation: {
    stick: { radius: 0.25, quality: 50 },
    sphere: { radius: 0.5 },
    surface: { opacity: 0.9 },
    backgroundColor: 'black'
  }
};
```

## Related Documentation

- [DOCKING_MODES.md](DOCKING_MODES.md) - Understanding docking engine behavior
- [REPRODUCIBILITY.md](REPRODUCIBILITY.md) - Ensuring consistent results
- [CONTRIBUTING.md](CONTRIBUTING.md) - Development guidelines for visualization features
