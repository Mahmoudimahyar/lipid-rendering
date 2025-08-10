import React, { useEffect, useRef, useState } from 'react'

/**
 * Specialized 3D viewer for protein-ligand docking visualization
 * Renders both protein structure and multiple ligand poses
 */
function ProteinLigandViewer({
  ligandSmiles,
  receptorPdbId,
  dockingPoses = [],
  selectedPose = 0,
  visiblePoses = new Set([0]),
  colorScheme = 'affinity',
  className = ''
}) {
  const containerRef = useRef(null)
  const viewerRef = useRef(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState(null)
  const [proteinLoaded, setProteinLoaded] = useState(false)
  const [ligandLoaded, setLigandLoaded] = useState(false)

  // Initialize 3Dmol viewer
  useEffect(() => {
    if (!containerRef.current || !window.$3Dmol) {
      setError('3Dmol.js library not available')
      setIsLoading(false)
      return
    }

    try {
      // Create unique container ID
      const containerId = `protein-ligand-viewer-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`
      containerRef.current.id = containerId

      // Initialize viewer
      const viewer = window.$3Dmol.createViewer(containerRef.current, {
        defaultcolors: window.$3Dmol.rasmolElementColors
      })
      
      viewerRef.current = viewer
      
      // Load protein and ligand
      loadProteinAndLigand()
      
    } catch (err) {
      console.error('ProteinLigandViewer initialization error:', err)
      setError(err.message)
      setIsLoading(false)
    }

    // Cleanup
    return () => {
      if (viewerRef.current) {
        try {
          viewerRef.current.clear()
        } catch (e) {
          console.warn('Error during 3Dmol cleanup:', e)
        }
      }
    }
  }, [receptorPdbId])

  // Update poses when visibility changes
  useEffect(() => {
    if (viewerRef.current && proteinLoaded && dockingPoses.length > 0) {
      updateLigandPoses()
    }
  }, [dockingPoses, visiblePoses, colorScheme, proteinLoaded])

  // Focus on selected pose
  useEffect(() => {
    if (viewerRef.current && dockingPoses[selectedPose]) {
      focusOnPose(selectedPose)
    }
  }, [selectedPose, dockingPoses])

  const loadProteinAndLigand = async () => {
    if (!viewerRef.current) return

    setIsLoading(true)
    setError(null)

    try {
      // Load protein structure
      await loadProteinStructure()
      
      // Load ligand poses if available
      if (dockingPoses.length > 0) {
        updateLigandPoses()
      } else if (ligandSmiles) {
        // Load just the ligand molecule if no poses
        await loadLigandMolecule()
      }

      // Render and zoom to fit
      viewerRef.current.zoomTo()
      viewerRef.current.render()
      
      setIsLoading(false)
      
    } catch (err) {
      console.error('Error loading protein and ligand:', err)
      setError(err.message)
      setIsLoading(false)
    }
  }

  const loadProteinStructure = async () => {
    if (!receptorPdbId || !viewerRef.current) return

    try {
      // Fetch protein structure from PDB
      const pdbUrl = `https://files.rcsb.org/download/${receptorPdbId.toUpperCase()}.pdb`
      const response = await fetch(pdbUrl)
      
      if (!response.ok) {
        throw new Error(`Failed to fetch protein ${receptorPdbId}: ${response.statusText}`)
      }
      
      const pdbData = await response.text()
      
      // Add protein model to viewer
      viewerRef.current.addModel(pdbData, 'pdb')
      
      // Style protein as cartoon (apply to all models)
      viewerRef.current.setStyle({}, {
        cartoon: { 
          color: 'lightgray',
          opacity: 0.8
        }
      })
      
      // Add surface representation (optional)
      // viewerRef.current.addSurface(window.$3Dmol.SurfaceType.VDW, {
      //   opacity: 0.3,
      //   color: 'white'
      // }, { model: proteinModel })
      
      setProteinLoaded(true)
      console.log(`Protein ${receptorPdbId} loaded successfully`)
      
    } catch (err) {
      console.error('Error loading protein:', err)
      setError(`Failed to load protein ${receptorPdbId}: ${err.message}`)
    }
  }

  const loadLigandMolecule = async () => {
    if (!ligandSmiles || !viewerRef.current) return

    try {
      // Generate 3D coordinates for ligand
      const response = await fetch('http://localhost:8000/api/ligand/prepare', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          smiles: ligandSmiles,
          num_conformers: 1
        })
      })

      if (!response.ok) {
        throw new Error(`Failed to prepare ligand: ${response.statusText}`)
      }

      const ligandData = await response.json()
      
      // Add ligand model
      viewerRef.current.addModel(ligandData.sdf, 'sdf')
      
      // Style ligand (apply to the last model)
      viewerRef.current.setStyle({}, {
        stick: {
          colorscheme: 'default',
          radius: 0.2
        }
      })
      
      setLigandLoaded(true)
      console.log('Ligand loaded successfully')
      
    } catch (err) {
      console.error('Error loading ligand:', err)
      // Don't set error for ligand issues, just log
    }
  }

  const updateLigandPoses = () => {
    if (!viewerRef.current || !dockingPoses.length) return

    try {
      // Clear all models except protein (model 0)
      // Note: 3Dmol.js doesn't have removeModel, so we clear all and re-add protein
      const proteinAdded = proteinLoaded
      
      // Clear and re-add protein if it was loaded
      if (proteinAdded) {
        // We'll re-add the protein in loadProteinStructure if needed
      }

      // Add visible poses
      dockingPoses.forEach((pose, index) => {
        if (!visiblePoses.has(index)) return

        try {
          // Create ligand SDF for this pose
          const ligandSDF = generatePoseSDF(pose, ligandSmiles)
          
          if (ligandSDF) {
            viewerRef.current.addModel(ligandSDF, 'sdf')
            
            // Get color for this pose
            const color = getPoseColor(pose, index, dockingPoses, colorScheme)
            
            // Style this pose (apply styling to last model added)
            viewerRef.current.setStyle({}, {
              stick: {
                color: color,
                radius: index === selectedPose ? 0.3 : 0.2
              }
            })
          }
        } catch (err) {
          console.warn(`Error adding pose ${index}:`, err)
        }
      })

      // Re-render
      viewerRef.current.render()
      
    } catch (err) {
      console.error('Error updating ligand poses:', err)
    }
  }

  const generatePoseSDF = (pose, smiles) => {
    // Mock SDF generation - in reality, this would use the pose coordinates
    // For now, we'll create a basic SDF structure
    const sdf = `
  Generated from docking pose
  
 13 12  0  0  0  0  0  0  0  0999 V2000
   ${(pose.center_x || 0).toFixed(4)}   ${(pose.center_y || 0).toFixed(4)}   ${(pose.center_z || 0).toFixed(4)} C   0  0  0  0  0  0  0  0  0  0  0  0
   ${(pose.center_x + 1 || 1).toFixed(4)}   ${(pose.center_y || 0).toFixed(4)}   ${(pose.center_z || 0).toFixed(4)} C   0  0  0  0  0  0  0  0  0  0  0  0
   ${(pose.center_x || 0).toFixed(4)}   ${(pose.center_y + 1 || 1).toFixed(4)}   ${(pose.center_z || 0).toFixed(4)} O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
$$$$`
    return sdf
  }

  const getPoseColor = (pose, index, poses, scheme) => {
    switch (scheme) {
      case 'affinity':
        const affinity = Math.abs(pose.affinity || -7.0)
        const maxAffinity = Math.max(...poses.map(p => Math.abs(p.affinity || -7.0)))
        const minAffinity = Math.min(...poses.map(p => Math.abs(p.affinity || -7.0)))
        const ratio = (affinity - minAffinity) / (maxAffinity - minAffinity || 1)
        
        return `rgb(${Math.floor(255 * ratio)}, ${Math.floor(255 * (1 - ratio))}, 0)`
      
      case 'rmsd':
        const rmsd = pose.rmsd_lb || 1.0
        return `rgb(0, ${Math.floor(255 * (1 - Math.min(rmsd / 3, 1)))}, 255)`
      
      case 'rainbow':
      default:
        const colors = ['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet']
        return colors[index % colors.length]
    }
  }

  const focusOnPose = (poseIndex) => {
    if (!viewerRef.current || !dockingPoses[poseIndex]) return

    const pose = dockingPoses[poseIndex]
    
    try {
      // Focus camera on pose using the correct 3Dmol.js API
      viewerRef.current.zoomTo()
      viewerRef.current.center()
      viewerRef.current.render()
      
      console.log(`Focused on pose ${poseIndex}`)
    } catch (err) {
      console.warn('Error focusing on pose:', err)
    }
  }

  if (error) {
    return (
      <div className={`flex items-center justify-center h-96 bg-red-50 border border-red-200 rounded-lg ${className}`}>
        <div className="text-center text-red-700">
          <div className="font-semibold mb-2">Visualization Error</div>
          <div className="text-sm">{error}</div>
        </div>
      </div>
    )
  }

  return (
    <div className={`relative ${className}`}>
      {/* 3D Viewer Container */}
      <div
        ref={containerRef}
        className="w-full h-96 bg-black rounded-lg border border-gray-300"
        style={{ minHeight: '400px' }}
      />
      
      {/* Loading Overlay */}
      {isLoading && (
        <div className="absolute inset-0 flex items-center justify-center bg-black bg-opacity-75 rounded-lg">
          <div className="text-center text-white">
            <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-white mx-auto mb-4"></div>
            <p>Loading protein and ligand...</p>
          </div>
        </div>
      )}
      
      {/* Status Overlay */}
      <div className="absolute bottom-4 left-4 bg-black bg-opacity-75 text-white rounded px-3 py-2 text-sm">
        <div>Protein: {proteinLoaded ? '✅' : '⏳'} {receptorPdbId || 'N/A'}</div>
        <div>Ligand: {ligandLoaded || dockingPoses.length > 0 ? '✅' : '⏳'} {dockingPoses.length || 1} pose(s)</div>
        <div>Visible: {visiblePoses.size} of {dockingPoses.length || 1}</div>
      </div>
    </div>
  )
}

export default ProteinLigandViewer



/**
 * Specialized 3D viewer for protein-ligand docking visualization
 * Renders both protein structure and multiple ligand poses
 */
function ProteinLigandViewer({
  ligandSmiles,
  receptorPdbId,
  dockingPoses = [],
  selectedPose = 0,
  visiblePoses = new Set([0]),
  colorScheme = 'affinity',
  className = ''
}) {
  const containerRef = useRef(null)
  const viewerRef = useRef(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState(null)
  const [proteinLoaded, setProteinLoaded] = useState(false)
  const [ligandLoaded, setLigandLoaded] = useState(false)

  // Initialize 3Dmol viewer
  useEffect(() => {
    if (!containerRef.current || !window.$3Dmol) {
      setError('3Dmol.js library not available')
      setIsLoading(false)
      return
    }

    try {
      // Create unique container ID
      const containerId = `protein-ligand-viewer-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`
      containerRef.current.id = containerId

      // Initialize viewer
      const viewer = window.$3Dmol.createViewer(containerRef.current, {
        defaultcolors: window.$3Dmol.rasmolElementColors
      })
      
      viewerRef.current = viewer
      
      // Load protein and ligand
      loadProteinAndLigand()
      
    } catch (err) {
      console.error('ProteinLigandViewer initialization error:', err)
      setError(err.message)
      setIsLoading(false)
    }

    // Cleanup
    return () => {
      if (viewerRef.current) {
        try {
          viewerRef.current.clear()
        } catch (e) {
          console.warn('Error during 3Dmol cleanup:', e)
        }
      }
    }
  }, [receptorPdbId])

  // Update poses when visibility changes
  useEffect(() => {
    if (viewerRef.current && proteinLoaded && dockingPoses.length > 0) {
      updateLigandPoses()
    }
  }, [dockingPoses, visiblePoses, colorScheme, proteinLoaded])

  // Focus on selected pose
  useEffect(() => {
    if (viewerRef.current && dockingPoses[selectedPose]) {
      focusOnPose(selectedPose)
    }
  }, [selectedPose, dockingPoses])

  const loadProteinAndLigand = async () => {
    if (!viewerRef.current) return

    setIsLoading(true)
    setError(null)

    try {
      // Load protein structure
      await loadProteinStructure()
      
      // Load ligand poses if available
      if (dockingPoses.length > 0) {
        updateLigandPoses()
      } else if (ligandSmiles) {
        // Load just the ligand molecule if no poses
        await loadLigandMolecule()
      }

      // Render and zoom to fit
      viewerRef.current.zoomTo()
      viewerRef.current.render()
      
      setIsLoading(false)
      
    } catch (err) {
      console.error('Error loading protein and ligand:', err)
      setError(err.message)
      setIsLoading(false)
    }
  }

  const loadProteinStructure = async () => {
    if (!receptorPdbId || !viewerRef.current) return

    try {
      // Fetch protein structure from PDB
      const pdbUrl = `https://files.rcsb.org/download/${receptorPdbId.toUpperCase()}.pdb`
      const response = await fetch(pdbUrl)
      
      if (!response.ok) {
        throw new Error(`Failed to fetch protein ${receptorPdbId}: ${response.statusText}`)
      }
      
      const pdbData = await response.text()
      
      // Add protein model to viewer
      viewerRef.current.addModel(pdbData, 'pdb')
      
      // Style protein as cartoon (apply to all models)
      viewerRef.current.setStyle({}, {
        cartoon: { 
          color: 'lightgray',
          opacity: 0.8
        }
      })
      
      // Add surface representation (optional)
      // viewerRef.current.addSurface(window.$3Dmol.SurfaceType.VDW, {
      //   opacity: 0.3,
      //   color: 'white'
      // }, { model: proteinModel })
      
      setProteinLoaded(true)
      console.log(`Protein ${receptorPdbId} loaded successfully`)
      
    } catch (err) {
      console.error('Error loading protein:', err)
      setError(`Failed to load protein ${receptorPdbId}: ${err.message}`)
    }
  }

  const loadLigandMolecule = async () => {
    if (!ligandSmiles || !viewerRef.current) return

    try {
      // Generate 3D coordinates for ligand
      const response = await fetch('http://localhost:8000/api/ligand/prepare', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          smiles: ligandSmiles,
          num_conformers: 1
        })
      })

      if (!response.ok) {
        throw new Error(`Failed to prepare ligand: ${response.statusText}`)
      }

      const ligandData = await response.json()
      
      // Add ligand model
      viewerRef.current.addModel(ligandData.sdf, 'sdf')
      
      // Style ligand (apply to the last model)
      viewerRef.current.setStyle({}, {
        stick: {
          colorscheme: 'default',
          radius: 0.2
        }
      })
      
      setLigandLoaded(true)
      console.log('Ligand loaded successfully')
      
    } catch (err) {
      console.error('Error loading ligand:', err)
      // Don't set error for ligand issues, just log
    }
  }

  const updateLigandPoses = () => {
    if (!viewerRef.current || !dockingPoses.length) return

    try {
      // Clear all models except protein (model 0)
      // Note: 3Dmol.js doesn't have removeModel, so we clear all and re-add protein
      const proteinAdded = proteinLoaded
      
      // Clear and re-add protein if it was loaded
      if (proteinAdded) {
        // We'll re-add the protein in loadProteinStructure if needed
      }

      // Add visible poses
      dockingPoses.forEach((pose, index) => {
        if (!visiblePoses.has(index)) return

        try {
          // Create ligand SDF for this pose
          const ligandSDF = generatePoseSDF(pose, ligandSmiles)
          
          if (ligandSDF) {
            viewerRef.current.addModel(ligandSDF, 'sdf')
            
            // Get color for this pose
            const color = getPoseColor(pose, index, dockingPoses, colorScheme)
            
            // Style this pose (apply styling to last model added)
            viewerRef.current.setStyle({}, {
              stick: {
                color: color,
                radius: index === selectedPose ? 0.3 : 0.2
              }
            })
          }
        } catch (err) {
          console.warn(`Error adding pose ${index}:`, err)
        }
      })

      // Re-render
      viewerRef.current.render()
      
    } catch (err) {
      console.error('Error updating ligand poses:', err)
    }
  }

  const generatePoseSDF = (pose, smiles) => {
    // Mock SDF generation - in reality, this would use the pose coordinates
    // For now, we'll create a basic SDF structure
    const sdf = `
  Generated from docking pose
  
 13 12  0  0  0  0  0  0  0  0999 V2000
   ${(pose.center_x || 0).toFixed(4)}   ${(pose.center_y || 0).toFixed(4)}   ${(pose.center_z || 0).toFixed(4)} C   0  0  0  0  0  0  0  0  0  0  0  0
   ${(pose.center_x + 1 || 1).toFixed(4)}   ${(pose.center_y || 0).toFixed(4)}   ${(pose.center_z || 0).toFixed(4)} C   0  0  0  0  0  0  0  0  0  0  0  0
   ${(pose.center_x || 0).toFixed(4)}   ${(pose.center_y + 1 || 1).toFixed(4)}   ${(pose.center_z || 0).toFixed(4)} O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
$$$$`
    return sdf
  }

  const getPoseColor = (pose, index, poses, scheme) => {
    switch (scheme) {
      case 'affinity':
        const affinity = Math.abs(pose.affinity || -7.0)
        const maxAffinity = Math.max(...poses.map(p => Math.abs(p.affinity || -7.0)))
        const minAffinity = Math.min(...poses.map(p => Math.abs(p.affinity || -7.0)))
        const ratio = (affinity - minAffinity) / (maxAffinity - minAffinity || 1)
        
        return `rgb(${Math.floor(255 * ratio)}, ${Math.floor(255 * (1 - ratio))}, 0)`
      
      case 'rmsd':
        const rmsd = pose.rmsd_lb || 1.0
        return `rgb(0, ${Math.floor(255 * (1 - Math.min(rmsd / 3, 1)))}, 255)`
      
      case 'rainbow':
      default:
        const colors = ['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet']
        return colors[index % colors.length]
    }
  }

  const focusOnPose = (poseIndex) => {
    if (!viewerRef.current || !dockingPoses[poseIndex]) return

    const pose = dockingPoses[poseIndex]
    
    try {
      // Focus camera on pose using the correct 3Dmol.js API
      viewerRef.current.zoomTo()
      viewerRef.current.center()
      viewerRef.current.render()
      
      console.log(`Focused on pose ${poseIndex}`)
    } catch (err) {
      console.warn('Error focusing on pose:', err)
    }
  }

  if (error) {
    return (
      <div className={`flex items-center justify-center h-96 bg-red-50 border border-red-200 rounded-lg ${className}`}>
        <div className="text-center text-red-700">
          <div className="font-semibold mb-2">Visualization Error</div>
          <div className="text-sm">{error}</div>
        </div>
      </div>
    )
  }

  return (
    <div className={`relative ${className}`}>
      {/* 3D Viewer Container */}
      <div
        ref={containerRef}
        className="w-full h-96 bg-black rounded-lg border border-gray-300"
        style={{ minHeight: '400px' }}
      />
      
      {/* Loading Overlay */}
      {isLoading && (
        <div className="absolute inset-0 flex items-center justify-center bg-black bg-opacity-75 rounded-lg">
          <div className="text-center text-white">
            <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-white mx-auto mb-4"></div>
            <p>Loading protein and ligand...</p>
          </div>
        </div>
      )}
      
      {/* Status Overlay */}
      <div className="absolute bottom-4 left-4 bg-black bg-opacity-75 text-white rounded px-3 py-2 text-sm">
        <div>Protein: {proteinLoaded ? '✅' : '⏳'} {receptorPdbId || 'N/A'}</div>
        <div>Ligand: {ligandLoaded || dockingPoses.length > 0 ? '✅' : '⏳'} {dockingPoses.length || 1} pose(s)</div>
        <div>Visible: {visiblePoses.size} of {dockingPoses.length || 1}</div>
      </div>
    </div>
  )
}

export default ProteinLigandViewer


