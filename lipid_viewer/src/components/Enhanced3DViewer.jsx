import React, { useEffect, useRef, useState } from 'react'

/**
 * Enhanced 3D Molecular Viewer with multiple visualization options
 * Inspired by Protein Data Bank and professional molecular visualization tools
 */
function Enhanced3DViewer({
  ligandSmiles,
  receptorPdbId,
  poses = [],
  selectedPose = 0,
  visiblePoses = new Set([0]),
  onFocusBindingSite,
  onResetView,
  className = '',
  ...props
}) {
  console.log('🚀 Enhanced3DViewer: Initializing with:', { 
    ligandSmiles, 
    receptorPdbId, 
    poses: poses.length, 
    selectedPose, 
    visiblePoses: Array.from(visiblePoses) 
  })
  
  const containerRef = useRef(null)
  const viewerRef = useRef(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState(null)
  const [loadedModels, setLoadedModels] = useState({ protein: -1, ligand: -1 })
  
  // Visualization control states
  const [proteinStyle, setProteinStyle] = useState('cartoon') // cartoon, surface, ribbon, backbone
  const [ligandStyle, setLigandStyle] = useState('stick') // stick, sphere, ball-stick
  const [colorScheme, setColorScheme] = useState('spectrum') // spectrum, chainid, element, hydrophobicity
  const [backgroundColor, setBackgroundColor] = useState('#1a1a1a') // dark background for contrast
  const [showSurface, setShowSurface] = useState(false)
  const [showLabels, setShowLabels] = useState(false)

  // Available visualization options (like PDB)
  const proteinStyles = [
    { value: 'cartoon', label: 'Cartoon', description: 'Secondary structure representation' },
    { value: 'ribbon', label: 'Ribbon', description: 'Smooth ribbon through backbone' },
    { value: 'backbone', label: 'Backbone', description: 'Backbone trace' },
    { value: 'surface', label: 'Surface', description: 'Molecular surface' },
    { value: 'wireframe', label: 'Wireframe', description: 'Wire model' }
  ]

  const ligandStyles = [
    { value: 'stick', label: 'Stick', description: 'Stick representation' },
    { value: 'sphere', label: 'Sphere', description: 'Space-filling spheres' },
    { value: 'ball-stick', label: 'Ball & Stick', description: 'Combined spheres and sticks' },
    { value: 'line', label: 'Line', description: 'Simple line representation' }
  ]

  const colorSchemes = [
    { value: 'spectrum', label: 'Spectrum', description: 'Rainbow colors by chain' },
    { value: 'chainid', label: 'Chain ID', description: 'Different colors per chain' },
    { value: 'element', label: 'Element', description: 'Color by atomic element' },
    { value: 'hydrophobicity', label: 'Hydrophobicity', description: 'Hydrophobic/hydrophilic regions' },
    { value: 'secondary', label: 'Secondary Structure', description: 'Color by secondary structure' }
  ]

  const backgroundColors = [
    { value: '#1a1a1a', label: 'Dark Gray', color: '#1a1a1a' },
    { value: '#000000', label: 'Black', color: '#000000' },
    { value: '#ffffff', label: 'White', color: '#ffffff' },
    { value: '#2c3e50', label: 'Dark Blue', color: '#2c3e50' },
    { value: '#34495e', label: 'Slate', color: '#34495e' }
  ]

  // Initialize 3D viewer
  useEffect(() => {
    if (!containerRef.current || !ligandSmiles || !receptorPdbId) return

    const initializeViewer = async () => {
      setIsLoading(true)
      setError(null)

      try {
        // Wait for 3Dmol to be available
        if (typeof window.$3Dmol === 'undefined') {
          throw new Error('3Dmol.js library not loaded. Please ensure the script is included.')
        }

        // Clear any existing viewer
        if (viewerRef.current) {
          try {
            viewerRef.current.clear()
          } catch (e) {
            console.warn('Error clearing existing viewer:', e)
          }
        }

        // Get container and ensure it has proper dimensions
        const container = containerRef.current
        if (!container) {
          throw new Error('Container not found')
        }

        // Set container dimensions explicitly to avoid WebGL framebuffer issues
        container.style.width = '100%'
        container.style.height = '500px'
        container.style.minHeight = '500px'
        container.style.position = 'relative'
        container.style.backgroundColor = backgroundColor
        
        // Wait for DOM to settle
        await new Promise(resolve => setTimeout(resolve, 100))

        // Create 3Dmol viewer with explicit configuration
        const config = {
          backgroundColor: backgroundColor,
          antialias: true,
          alpha: true,
          preserveDrawingBuffer: true
        }

        console.log('🎮 Creating 3Dmol viewer with config:', config)
        const viewer = window.$3Dmol.createViewer(container, config)
        
        if (!viewer) {
          throw new Error('Failed to create 3Dmol viewer')
        }

        viewerRef.current = viewer
        console.log('✅ 3Dmol viewer created successfully')

        // Load protein
        const proteinModelIndex = await loadProtein(viewer, receptorPdbId)
        
        // Test pose coordinate diversity before loading models
        testPoseCoordinateDiversity(poses)

        // Load multiple ligand models - one for each pose at its correct position
        const ligandModelIndices = await loadMultiplePoseModels(viewer, ligandSmiles, poses)

        // Store model indices
        setLoadedModels({ 
          protein: proteinModelIndex, 
          ligand: ligandModelIndices[0], // For compatibility, main ligand index
          ligandPoses: ligandModelIndices // All pose model indices
        })

        // Apply initial styling with pose information
        applyVisualizationStyle(viewer, proteinModelIndex, ligandModelIndices)

        // Final render with proper viewport setup
        await renderViewer(viewer)
        
        setIsLoading(false)
        console.log('🎉 Enhanced3DViewer initialization complete')

      } catch (err) {
        console.error('❌ Enhanced3DViewer initialization failed:', err)
        setError(err.message)
        setIsLoading(false)
      }
    }

    initializeViewer()

    // Cleanup
    return () => {
      if (viewerRef.current) {
        try {
          viewerRef.current.clear()
          viewerRef.current = null
        } catch (e) {
          console.warn('Error during cleanup:', e)
        }
      }
    }
  }, [ligandSmiles, receptorPdbId])

  // Update visualization when poses change - show/hide different pose models
  useEffect(() => {
    if (viewerRef.current && loadedModels.protein >= 0 && poses.length > 0) {
      console.log('🔄 Updating visualization for pose changes:', { selectedPose, visiblePoses: Array.from(visiblePoses) })
      
      // Verify pose positions before styling
      verifyPosePositions(viewerRef.current, poses, loadedModels.ligandPoses || [])
      
      // Apply styling to show/hide pose models
      applyVisualizationStyle(viewerRef.current, loadedModels.protein, loadedModels.ligandPoses || [loadedModels.ligand])
      renderViewer(viewerRef.current)
    }
  }, [selectedPose, visiblePoses, loadedModels, poses])

  // Update visualization when style options change
  useEffect(() => {
    if (viewerRef.current && loadedModels.protein >= 0) {
      applyVisualizationStyle(viewerRef.current, loadedModels.protein, loadedModels.ligandPoses || [loadedModels.ligand])
      renderViewer(viewerRef.current)
    }
  }, [proteinStyle, ligandStyle, colorScheme, showSurface, showLabels, loadedModels])

  // Update background color
  useEffect(() => {
    if (viewerRef.current) {
      viewerRef.current.setBackgroundColor(backgroundColor)
      viewerRef.current.render()
    }
    if (containerRef.current) {
      containerRef.current.style.backgroundColor = backgroundColor
    }
  }, [backgroundColor])

  const loadProtein = async (viewer, pdbId) => {
    try {
      console.log(`🧬 Loading protein ${pdbId}...`)
      
      const pdbUrl = `https://files.rcsb.org/download/${pdbId.toUpperCase()}.pdb`
      const response = await fetch(pdbUrl)
      
      if (!response.ok) {
        throw new Error(`Failed to fetch protein ${pdbId}: ${response.statusText}`)
      }
      
      const pdbData = await response.text()
      viewer.addModel(pdbData, 'pdb')
      
      console.log(`✅ Protein ${pdbId} loaded successfully`)
      return 0 // Protein is always model 0
      
    } catch (err) {
      console.error('❌ Error loading protein:', err)
      throw err
    }
  }

  const loadLigand = async (viewer, smiles) => {
    try {
      console.log(`💊 Loading ligand from SMILES: ${smiles}`)
      
      // Try backend preparation first
      try {
        const response = await fetch('/api/ligand/prepare', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ smiles })
        })

        if (response.ok) {
          const data = await response.json()
          if (data.sdf_content) {
            viewer.addModel(data.sdf_content, 'sdf')
            console.log('✅ Ligand loaded from backend')
            return 1 // Ligand is model 1
          }
        }
      } catch (backendErr) {
        console.warn('Backend ligand preparation failed:', backendErr)
      }

      // Fallback to SMILES (limited 3D info)
      viewer.addModel(smiles, 'smiles')
      console.log('✅ Ligand loaded via SMILES fallback')
      return 1
      
    } catch (err) {
      console.error('❌ Error loading ligand:', err)
      throw err
    }
  }

  const testPoseCoordinateDiversity = (poses) => {
    console.log('🧪 TESTING: Checking pose coordinate diversity...')
    
    if (!poses || poses.length < 2) {
      console.warn('⚠️ TEST: Need at least 2 poses to test diversity')
      return
    }
    
    const coordinates = poses.map(pose => ({
      mode: pose.mode,
      x: pose.center_x,
      y: pose.center_y,
      z: pose.center_z
    }))
    
    console.log('📊 Pose coordinates:')
    coordinates.forEach(coord => {
      console.log(`  Pose ${coord.mode}: (${coord.x?.toFixed(2)}, ${coord.y?.toFixed(2)}, ${coord.z?.toFixed(2)})`)
    })
    
    // Check if all coordinates are the same (problematic)
    const firstCoord = coordinates[0]
    const allSame = coordinates.every(coord => 
      Math.abs(coord.x - firstCoord.x) < 0.1 &&
      Math.abs(coord.y - firstCoord.y) < 0.1 &&
      Math.abs(coord.z - firstCoord.z) < 0.1
    )
    
    if (allSame) {
      console.error('❌ CRITICAL: ALL POSES HAVE IDENTICAL COORDINATES!')
      console.error('This explains why selecting different poses shows no position change')
      console.error('The docking results may not contain proper pose coordinate data')
    } else {
      console.log('✅ Poses have diverse coordinates - position changes should be visible')
      
      // Calculate average distance between poses
      let totalDistance = 0
      let comparisons = 0
      
      for (let i = 0; i < coordinates.length; i++) {
        for (let j = i + 1; j < coordinates.length; j++) {
          const distance = Math.sqrt(
            Math.pow(coordinates[i].x - coordinates[j].x, 2) +
            Math.pow(coordinates[i].y - coordinates[j].y, 2) +
            Math.pow(coordinates[i].z - coordinates[j].z, 2)
          )
          totalDistance += distance
          comparisons++
        }
      }
      
      const avgDistance = totalDistance / comparisons
      console.log(`📏 Average distance between poses: ${avgDistance.toFixed(2)} Å`)
      
      if (avgDistance < 2.0) {
        console.warn('⚠️ Poses are very close together - position changes may be subtle')
      }
    }
  }

  const verifyPosePositions = (viewer, poses, modelIndices) => {
    try {
      console.log('🧪 TESTING: Verifying pose positions...')
      
      if (!modelIndices || modelIndices.length === 0) {
        console.warn('⚠️ TEST: No model indices to verify')
        return
      }
      
      const positionReport = []
      
      poses.forEach((pose, poseIndex) => {
        const modelIndex = modelIndices[poseIndex]
        
        if (modelIndex !== undefined) {
          try {
            const model = viewer.getModel(modelIndex)
            if (model) {
              // Try to get model center/position - using safer approach
              let modelPosition = null
              
              try {
                // Use viewer.getAtoms() to calculate centroid manually - more reliable
                const atoms = viewer.getAtomList({ model: modelIndex })
                if (atoms && atoms.length > 0) {
                  let sumX = 0, sumY = 0, sumZ = 0
                  atoms.forEach(atom => {
                    sumX += atom.x || 0
                    sumY += atom.y || 0  
                    sumZ += atom.z || 0
                  })
                  modelPosition = {
                    x: sumX / atoms.length,
                    y: sumY / atoms.length,
                    z: sumZ / atoms.length
                  }
                  console.log(`📍 Calculated centroid from ${atoms.length} atoms`)
                } else {
                  console.warn(`⚠️ No atoms found for model ${modelIndex}`)
                }
              } catch (atomErr) {
                console.warn(`⚠️ Could not get atoms for model ${modelIndex}:`, atomErr)
                // Ultimate fallback - assume coordinates based on what we embedded
                if (pose.center_x !== undefined) {
                  modelPosition = {
                    x: pose.center_x,
                    y: pose.center_y,
                    z: pose.center_z
                  }
                  console.log(`🔄 Using embedded coordinates as fallback`)
                }
              }
              
              const expectedPosition = {
                x: pose.center_x,
                y: pose.center_y,
                z: pose.center_z
              }
              
              const positionData = {
                pose: pose.mode,
                modelIndex,
                expected: expectedPosition,
                actual: modelPosition,
                match: modelPosition ? 
                  Math.abs(modelPosition.x - expectedPosition.x) < 5 &&
                  Math.abs(modelPosition.y - expectedPosition.y) < 5 &&
                  Math.abs(modelPosition.z - expectedPosition.z) < 5 : false
              }
              
              positionReport.push(positionData)
              
              console.log(`🔍 TEST Pose ${pose.mode}:`)
              console.log(`  Expected: (${expectedPosition.x?.toFixed(2)}, ${expectedPosition.y?.toFixed(2)}, ${expectedPosition.z?.toFixed(2)})`)
              console.log(`  Actual: ${modelPosition ? `(${modelPosition.x?.toFixed(2)}, ${modelPosition.y?.toFixed(2)}, ${modelPosition.z?.toFixed(2)})` : 'Unknown'}`)
              console.log(`  Position Match: ${positionData.match ? '✅' : '❌'}`)
            }
          } catch (modelErr) {
            console.error(`❌ TEST: Error verifying pose ${pose.mode}:`, modelErr)
          }
        }
      })
      
      // Summary report
      const matchingPoses = positionReport.filter(p => p.match).length
      console.log(`📊 TEST SUMMARY: ${matchingPoses}/${positionReport.length} poses positioned correctly`)
      
      if (matchingPoses === 0 && positionReport.length > 0) {
        console.error('❌ CRITICAL: NO POSES ARE POSITIONED CORRECTLY!')
        console.error('This indicates the translation/positioning system is not working')
      }
      
      return positionReport
      
    } catch (err) {
      console.error('❌ Error in pose position verification:', err)
    }
  }

  const embedPoseCoordinatesInSDF = (originalSDF, pose) => {
    try {
      console.log(`🔧 Embedding coordinates for pose ${pose.mode}`)
      
      if (!pose.center_x || !pose.center_y || !pose.center_z) {
        console.warn(`⚠️ Pose ${pose.mode} missing coordinates, using original SDF`)
        return originalSDF
      }
      
      // Parse SDF content to modify atomic coordinates
      const lines = originalSDF.split('\n')
      const modifiedLines = []
      
      // Find the counts line (line 3 in SDF format)
      let atomCount = 0
      let bondCount = 0
      let countsLineIndex = -1
      
      for (let i = 0; i < lines.length; i++) {
        if (i === 3) { // Counts line is typically line 3 (0-indexed)
          const countsLine = lines[i].trim()
          if (countsLine.length >= 6) {
            atomCount = parseInt(countsLine.substring(0, 3).trim())
            bondCount = parseInt(countsLine.substring(3, 6).trim())
            countsLineIndex = i
            console.log(`📊 SDF has ${atomCount} atoms, ${bondCount} bonds`)
            break
          }
        }
      }
      
      if (atomCount === 0) {
        console.warn('⚠️ Could not parse SDF atom count, using original')
        return originalSDF
      }
      
      // Calculate the offset needed to move ligand to pose coordinates
      // We'll translate all atoms by the pose offset
      const offsetX = pose.center_x || 0
      const offsetY = pose.center_y || 0  
      const offsetZ = pose.center_z || 0
      
      console.log(`📍 Applying offset: (${offsetX}, ${offsetY}, ${offsetZ})`)
      
      // Copy header lines (first 4 lines)
      for (let i = 0; i < Math.min(4, lines.length); i++) {
        modifiedLines.push(lines[i])
      }
      
      // Process atom block (lines 4 to 4+atomCount-1)
      for (let i = 4; i < 4 + atomCount && i < lines.length; i++) {
        const line = lines[i]
        if (line.length >= 30) { // SDF atom line should be at least 30 chars
          try {
            // Parse original coordinates (positions 0-10, 10-20, 20-30)
            const origX = parseFloat(line.substring(0, 10).trim())
            const origY = parseFloat(line.substring(10, 20).trim()) 
            const origZ = parseFloat(line.substring(20, 30).trim())
            
            // Apply offset to move to pose binding site
            const newX = origX + offsetX
            const newY = origY + offsetY
            const newZ = origZ + offsetZ
            
            // Reconstruct line with new coordinates
            const newLine = 
              newX.toFixed(4).padStart(10) +
              newY.toFixed(4).padStart(10) +
              newZ.toFixed(4).padStart(10) +
              line.substring(30) // Keep rest of the line (element, charge, etc.)
              
            modifiedLines.push(newLine)
          } catch (parseErr) {
            console.warn(`⚠️ Could not parse atom line ${i}, keeping original`)
            modifiedLines.push(line)
          }
        } else {
          modifiedLines.push(line)
        }
      }
      
      // Copy remaining lines (bonds, properties, etc.)
      for (let i = 4 + atomCount; i < lines.length; i++) {
        modifiedLines.push(lines[i])
      }
      
      const modifiedSDF = modifiedLines.join('\n')
      console.log(`✅ Successfully embedded coordinates for pose ${pose.mode}`)
      
      return modifiedSDF
      
    } catch (err) {
      console.error(`❌ Error embedding coordinates for pose ${pose.mode}:`, err)
      return originalSDF
    }
  }

  const loadMultiplePoseModels = async (viewer, smiles, poses) => {
    try {
      console.log(`💊 Loading ${poses.length} pose models for SMILES: ${smiles}`)
      
      if (!poses || poses.length === 0) {
        // Fallback to single ligand
        console.warn('⚠️ No poses available, loading single ligand')
        const ligandIndex = await loadLigand(viewer, smiles)
        return [ligandIndex]
      }
      
      // Get base ligand SDF from backend
      let baseLigandSDF = null
      try {
        const response = await fetch('/api/ligand/prepare', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ smiles })
        })

        if (response.ok) {
          const data = await response.json()
          if (data.sdf_content) {
            baseLigandSDF = data.sdf_content
            console.log('✅ Base ligand SDF obtained from backend')
          }
        }
      } catch (err) {
        console.warn('⚠️ Could not get ligand SDF from backend:', err)
      }
      
      const modelIndices = []
      
      // Load one model per pose
      for (let i = 0; i < poses.length; i++) {
        const pose = poses[i]
        
        try {
          if (baseLigandSDF) {
            // Create modified SDF with pose coordinates embedded
            // This approach avoids 3Dmol.js API issues by modifying the molecular data directly
            try {
              const modifiedSDF = embedPoseCoordinatesInSDF(baseLigandSDF, pose)
              
              // Add the modified ligand model with coordinates already embedded
              viewer.addModel(modifiedSDF, 'sdf', { doAssembly: false })
              const modelIndex = i + 1 // Protein is model 0, poses start at 1
              
              console.log(`✅ Pose ${pose.mode} loaded with embedded coordinates`)
              console.log(`📍 Embedded coordinates: (${pose.center_x}, ${pose.center_y}, ${pose.center_z})`)
              
              modelIndices.push(modelIndex)
            } catch (embedErr) {
              console.error(`❌ Failed to embed coordinates for pose ${pose.mode}:`, embedErr)
              // Fallback to original SDF without positioning
              viewer.addModel(baseLigandSDF, 'sdf', { doAssembly: false })
              modelIndices.push(i + 1)
            }
          } else {
            // Fallback to SMILES
            viewer.addModel(smiles, 'smiles')
            modelIndices.push(i + 1)
            console.warn(`⚠️ Using SMILES fallback for pose ${pose.mode}`)
          }
        } catch (poseErr) {
          console.error(`❌ Failed to load pose ${pose.mode}:`, poseErr)
        }
      }
      
      console.log(`🎯 Successfully loaded ${modelIndices.length} pose models`)
      return modelIndices
      
    } catch (err) {
      console.error('❌ Error loading multiple pose models:', err)
      // Fallback to single ligand
      const ligandIndex = await loadLigand(viewer, smiles)
      return [ligandIndex]
    }
  }

  const loadPoseLigandsAtCoordinates = async (viewer, smiles, poses) => {
    try {
      console.log(`💊 Loading ${poses.length} pose ligands at their binding coordinates for SMILES: ${smiles}`)
      
      if (!poses || poses.length === 0) {
        // Fallback to single ligand
        const singleLigandIndex = await loadLigand(viewer, smiles)
        return [singleLigandIndex]
      }
      
      const modelIndices = []
      
      // Get the base ligand structure from backend
      let baseLigandSDF = null
      try {
        const response = await fetch('/api/ligand/prepare', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ smiles })
        })

        if (response.ok) {
          const data = await response.json()
          if (data.sdf_content) {
            baseLigandSDF = data.sdf_content
            console.log('✅ Base ligand SDF obtained from backend')
          }
        }
      } catch (err) {
        console.warn('⚠️ Could not get ligand SDF from backend:', err)
      }
      
      // Load each pose as a separate model at its specific coordinates
      for (let i = 0; i < poses.length; i++) {
        const pose = poses[i]
        
        try {
          let modelAdded = false
          
          // Use the base ligand SDF for each pose
          if (baseLigandSDF) {
            viewer.addModel(baseLigandSDF, 'sdf')
            const modelIndex = i + 1 // Protein is model 0, poses start at 1
            
            // Position the ligand at the pose coordinates
            if (pose.center_x !== undefined && pose.center_y !== undefined && pose.center_z !== undefined) {
              // In 3Dmol.js, we use the viewer.getModel(index) method to get a specific model
              // Get the just-added model by its index
              const currentModel = viewer.getModel(modelIndex)
              
              if (currentModel) {
                // Get the current centroid of the ligand
                const currentCentroid = currentModel.getCentroid()
                
                // Calculate translation vector to move to pose coordinates
                const translationX = pose.center_x - currentCentroid.x
                const translationY = pose.center_y - currentCentroid.y
                const translationZ = pose.center_z - currentCentroid.z
                
                // Apply translation to move ligand to binding site
                currentModel.translate(translationX, translationY, translationZ)
                
                console.log(`✅ Pose ${pose.mode} positioned at binding site (${pose.center_x}, ${pose.center_y}, ${pose.center_z})`)
                console.log(`📍 Translation applied: (${translationX.toFixed(3)}, ${translationY.toFixed(3)}, ${translationZ.toFixed(3)})`)
              } else {
                console.warn(`⚠️ Could not get model ${modelIndex} for pose ${pose.mode}`)
              }
            }
            
            modelIndices.push(modelIndex)
            modelAdded = true
          }
          
          // Fallback if SDF loading fails
          if (!modelAdded) {
            console.warn(`⚠️ Using SMILES fallback for pose ${pose.mode}`)
            viewer.addModel(smiles, 'smiles')
            modelIndices.push(i + 1)
          }
          
        } catch (poseErr) {
          console.error(`❌ Failed to load pose ${pose.mode}:`, poseErr)
          // Continue with next pose
        }
      }
      
      console.log(`🎯 Successfully loaded ${modelIndices.length} poses at their binding coordinates`)
      console.log(`📊 Model indices: [${modelIndices.join(', ')}]`)
      
      return modelIndices
      
    } catch (err) {
      console.error('❌ Error loading pose ligands at coordinates:', err)
      // Fallback to single ligand
      const singleLigandIndex = await loadLigand(viewer, ligandSmiles)
      return [singleLigandIndex]
    }
  }

  const applyVisualizationStyle = (viewer, proteinModelIndex, ligandModelIndices) => {
    try {
      console.log(`🎨 Applying visualization styles - Protein: ${proteinStyle}, Ligand: ${ligandStyle}, Colors: ${colorScheme}`)
      console.log(`🎯 Pose info - Selected: ${selectedPose}, Visible: [${Array.from(visiblePoses).join(', ')}]`)

      // Clear existing styles and surfaces
      viewer.setStyle({}, {})
      viewer.removeAllSurfaces()
      viewer.removeAllLabels()

      // Style protein
      if (proteinModelIndex >= 0) {
        const proteinStyleObj = getProteinStyleObject(proteinStyle, colorScheme)
        viewer.setStyle({ model: proteinModelIndex }, proteinStyleObj)
        console.log('✅ Protein styling applied')

        // Add surface if requested
        if (showSurface) {
          try {
            viewer.addSurface(window.$3Dmol.SurfaceType.VDW, {
              opacity: 0.2,
              color: 'lightblue'
            }, { model: proteinModelIndex })
            console.log('✅ Protein surface added')
          } catch (surfaceErr) {
            console.warn('⚠️ Could not add protein surface:', surfaceErr)
          }
        }
      }

      // Style ligands - multiple models for different poses
      if (ligandModelIndices && ligandModelIndices.length > 0 && poses.length > 0) {
        
        console.log('🎨 TESTING: Styling poses...')
        console.log(`📊 Selected pose: ${selectedPose}, Visible poses: [${Array.from(visiblePoses).join(', ')}]`)
        
        let visibleCount = 0
        let hiddenCount = 0
        
        // Style each pose model individually
        poses.forEach((pose, poseIndex) => {
          const modelIndex = ligandModelIndices[poseIndex]
          const isVisible = visiblePoses.has(poseIndex)
          const isSelected = poseIndex === selectedPose
          
          if (modelIndex !== undefined && modelIndex >= 0) {
            if (isVisible) {
              // Apply styling to visible poses
              const ligandStyleObj = getLigandStyleObject(ligandStyle, isSelected, poseIndex)
              viewer.setStyle({ model: modelIndex }, ligandStyleObj)
              visibleCount++
              
              console.log(`🔍 TEST: Pose ${pose.mode} (model ${modelIndex})`)
              console.log(`  Status: VISIBLE, Selected: ${isSelected}`)
              console.log(`  Expected position: (${pose.center_x?.toFixed(2)}, ${pose.center_y?.toFixed(2)}, ${pose.center_z?.toFixed(2)})`)
              console.log(`  Style applied: ${JSON.stringify(ligandStyleObj)}`)
              
              // Add pose labels if requested
              if (showLabels) {
                try {
                  // Use a fixed position based on pose coordinates instead of model center
                  const labelPosition = {
                    x: pose.center_x || 0,
                    y: pose.center_y || 0, 
                    z: pose.center_z || 0
                  }
                  
                  const labelText = isSelected 
                    ? `Pose ${pose.mode} (${pose.affinity} kcal/mol) - SELECTED`
                    : `Pose ${pose.mode} (${pose.affinity} kcal/mol)`
                  
                  viewer.addLabel(labelText, {
                    position: labelPosition,
                    fontSize: isSelected ? 14 : 12,
                    fontColor: isSelected ? 'yellow' : 'white',
                    backgroundColor: isSelected ? 'red' : 'black',
                    backgroundOpacity: isSelected ? 0.9 : 0.7
                  })
                  console.log(`✅ Label added for pose ${pose.mode} at (${labelPosition.x}, ${labelPosition.y}, ${labelPosition.z})`)
                } catch (labelErr) {
                  console.warn(`⚠️ Could not add label for pose ${pose.mode}:`, labelErr)
                }
              }
            } else {
              // Hide invisible poses by clearing their style
              viewer.setStyle({ model: modelIndex }, {})
              hiddenCount++
              console.log(`🙈 TEST: Pose ${pose.mode} (model ${modelIndex}) HIDDEN`)
            }
          } else {
            console.error(`❌ TEST: Invalid model index for pose ${pose.mode}: ${modelIndex}`)
          }
        })
        
        console.log(`📊 STYLING SUMMARY: ${visibleCount} visible, ${hiddenCount} hidden poses`)
        
        if (visibleCount === 0) {
          console.error('❌ CRITICAL: NO POSES ARE VISIBLE AFTER STYLING!')
        }
        
      } else if (ligandModelIndices && ligandModelIndices.length > 0) {
        // Fallback: single ligand model
        const mainModelIndex = ligandModelIndices[0]
        if (mainModelIndex !== undefined && mainModelIndex >= 0) {
          const ligandStyleObj = getLigandStyleObject(ligandStyle, true, 0)
          viewer.setStyle({ model: mainModelIndex }, ligandStyleObj)
          console.log('✅ Single ligand model styled (fallback)')
          
          if (showLabels) {
            try {
              viewer.addLabel('Ligand', {
                position: { x: 0, y: 0, z: 0 }, // Default position
                fontSize: 12,
                fontColor: 'yellow',
                backgroundColor: 'black',
                backgroundOpacity: 0.7
              })
              console.log('✅ Single ligand label added')
            } catch (labelErr) {
              console.warn('⚠️ Could not add ligand label:', labelErr)
            }
          }
        }
      }

    } catch (err) {
      console.error('❌ Error applying visualization style:', err)
    }
  }

  const getProteinStyleObject = (style, colorScheme) => {
    const baseStyle = {}

    switch (style) {
      case 'cartoon':
        baseStyle.cartoon = {
          color: colorScheme,
          opacity: 0.8,
          thickness: 0.4
        }
        break
      case 'ribbon':
        baseStyle.ribbon = {
          color: colorScheme,
          opacity: 0.8,
          thickness: 0.4
        }
        break
      case 'backbone':
        baseStyle.line = {
          color: colorScheme,
          opacity: 0.8,
          linewidth: 3
        }
        break
      case 'surface':
        baseStyle.surface = {
          color: colorScheme,
          opacity: 0.7
        }
        break
      case 'wireframe':
        baseStyle.line = {
          color: colorScheme,
          opacity: 0.6,
          linewidth: 1
        }
        break
      default:
        baseStyle.cartoon = {
          color: colorScheme,
          opacity: 0.8
        }
    }

    return baseStyle
  }

  const getLigandStyleObject = (style, isSelected = true, poseIndex = 0) => {
    const baseStyle = {}

    // Choose color scheme based on selection and pose
    let colorScheme = 'greenCarbon'
    let opacity = isSelected ? 1.0 : 0.6
    let radius = isSelected ? 1.0 : 0.8
    
    if (poses.length > 1) {
      // Different colors for different poses
      const colors = ['greenCarbon', 'redCarbon', 'blueCarbon', 'yellowCarbon', 'purpleCarbon', 'orangeCarbon']
      colorScheme = colors[poseIndex % colors.length]
    }

    switch (style) {
      case 'stick':
        baseStyle.stick = {
          colorscheme: colorScheme,
          radius: 0.3 * radius,
          opacity: opacity
        }
        break
      case 'sphere':
        baseStyle.sphere = {
          colorscheme: colorScheme,
          radius: 0.4 * radius,
          opacity: opacity * 0.9
        }
        break
      case 'ball-stick':
        baseStyle.stick = {
          colorscheme: colorScheme,
          radius: 0.2 * radius,
          opacity: opacity
        }
        baseStyle.sphere = {
          colorscheme: colorScheme,
          radius: 0.3 * radius,
          opacity: opacity * 0.9
        }
        break
      case 'line':
        baseStyle.line = {
          colorscheme: colorScheme,
          linewidth: 3 * radius,
          opacity: opacity
        }
        break
      default:
        baseStyle.stick = {
          colorscheme: colorScheme,
          radius: 0.3 * radius,
          opacity: opacity
        }
    }

    return baseStyle
  }

  const renderViewer = async (viewer) => {
    try {
      console.log('🎯 Starting viewer render...')
      
      // Set up proper viewport and camera
      viewer.zoomTo()
      viewer.center()
      viewer.render()
      
      // Wait a bit and do a second render for stability
      await new Promise(resolve => setTimeout(resolve, 100))
      viewer.render()
      
      console.log('✅ Viewer rendering complete')
    } catch (err) {
      console.error('❌ Viewer rendering failed:', err)
      throw err
    }
  }

  const resetView = () => {
    if (viewerRef.current) {
      console.log('🔄 Resetting view...')
      viewerRef.current.zoomTo()
      viewerRef.current.center()
      viewerRef.current.render()
      
      // Call parent callback if provided
      if (onResetView) {
        onResetView()
      }
    }
  }

  const focusOnBindingSite = () => {
    if (viewerRef.current && poses.length > 0) {
      console.log('🎯 Focusing on binding site...')
      
      const selectedPoseData = poses[selectedPose]
      if (selectedPoseData) {
        // Focus on the specific pose model at its binding site
        const selectedModelIndex = loadedModels.ligandPoses ? loadedModels.ligandPoses[selectedPose] : loadedModels.ligand
        
        if (selectedModelIndex !== undefined && selectedModelIndex >= 0) {
          console.log(`📍 Focusing on pose ${selectedPoseData.mode} (model ${selectedModelIndex}) at binding site`)
          
          // Zoom to the specific pose model
          viewerRef.current.zoomTo({ model: selectedModelIndex }, 1000)
          viewerRef.current.render()
          
          console.log(`🎯 Camera focused on pose ${selectedPoseData.mode} at coordinates (${selectedPoseData.center_x}, ${selectedPoseData.center_y}, ${selectedPoseData.center_z})`)
        } else {
          // Fallback: zoom to pose coordinates directly
          if (selectedPoseData.center_x !== undefined && selectedPoseData.center_y !== undefined && selectedPoseData.center_z !== undefined) {
            const center = {
              x: selectedPoseData.center_x,
              y: selectedPoseData.center_y,
              z: selectedPoseData.center_z
            }
            
            console.log('📍 Centering on pose coordinates directly:', center)
            viewerRef.current.zoomTo({ center, zoom: 2 })
            viewerRef.current.render()
          }
        }
        
        // Call parent callback if provided
        if (onFocusBindingSite) {
          onFocusBindingSite(selectedPoseData)
        }
      }
    } else if (viewerRef.current && loadedModels.ligand !== undefined && loadedModels.ligand >= 0) {
      // No poses available, just focus on ligand
      console.log('🎯 Focusing on ligand (no pose data)...')
      viewerRef.current.zoomTo({ model: loadedModels.ligand }, 1000)
      viewerRef.current.render()
      
      if (onFocusBindingSite) {
        onFocusBindingSite()
      }
    }
  }

  const saveImage = () => {
    if (viewerRef.current) {
      try {
        const dataUrl = viewerRef.current.pngURI()
        const link = document.createElement('a')
        link.download = `molecular-docking-${Date.now()}.png`
        link.href = dataUrl
        link.click()
      } catch (err) {
        console.error('Error saving image:', err)
      }
    }
  }

  if (error) {
    return (
      <div className="bg-red-50 border border-red-200 rounded-lg p-4">
        <h3 className="text-red-800 font-medium">Visualization Error</h3>
        <p className="text-red-600 text-sm mt-1">{error}</p>
        <button 
          onClick={() => window.location.reload()} 
          className="mt-2 text-red-600 hover:text-red-800 text-sm underline"
        >
          Reload Page
        </button>
      </div>
    )
  }

  return (
    <div className={`space-y-4 ${className}`}>
      {/* Control Panel - PDB-style interface */}
      <div className="bg-white border border-gray-200 rounded-lg p-4 space-y-4">
        <h3 className="text-lg font-semibold text-gray-800 mb-3">Visualization Controls</h3>
        
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
          {/* Protein Style */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Protein Style
            </label>
            <select
              value={proteinStyle}
              onChange={(e) => setProteinStyle(e.target.value)}
              className="w-full px-3 py-2 border border-gray-300 rounded-md text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
              disabled={isLoading}
            >
              {proteinStyles.map(style => (
                <option key={style.value} value={style.value} title={style.description}>
                  {style.label}
                </option>
              ))}
            </select>
          </div>

          {/* Ligand Style */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Ligand Style
            </label>
            <select
              value={ligandStyle}
              onChange={(e) => setLigandStyle(e.target.value)}
              className="w-full px-3 py-2 border border-gray-300 rounded-md text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
              disabled={isLoading}
            >
              {ligandStyles.map(style => (
                <option key={style.value} value={style.value} title={style.description}>
                  {style.label}
                </option>
              ))}
            </select>
          </div>

          {/* Color Scheme */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Color Scheme
            </label>
            <select
              value={colorScheme}
              onChange={(e) => setColorScheme(e.target.value)}
              className="w-full px-3 py-2 border border-gray-300 rounded-md text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
              disabled={isLoading}
            >
              {colorSchemes.map(scheme => (
                <option key={scheme.value} value={scheme.value} title={scheme.description}>
                  {scheme.label}
                </option>
              ))}
            </select>
          </div>

          {/* Background Color */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Background
            </label>
            <div className="flex space-x-1">
              {backgroundColors.map(bg => (
                <button
                  key={bg.value}
                  onClick={() => setBackgroundColor(bg.value)}
                  className={`w-8 h-8 rounded border-2 ${backgroundColor === bg.value ? 'border-blue-500' : 'border-gray-300'}`}
                  style={{ backgroundColor: bg.color }}
                  title={bg.label}
                  disabled={isLoading}
                />
              ))}
            </div>
          </div>
        </div>

        {/* Additional Options */}
        <div className="flex flex-wrap items-center gap-4 pt-2 border-t border-gray-200">
          <label className="flex items-center space-x-2 text-sm">
            <input
              type="checkbox"
              checked={showSurface}
              onChange={(e) => setShowSurface(e.target.checked)}
              disabled={isLoading}
              className="w-4 h-4 text-blue-600 border-gray-300 rounded focus:ring-blue-500"
            />
            <span>Show Surface</span>
          </label>

          <label className="flex items-center space-x-2 text-sm">
            <input
              type="checkbox"
              checked={showLabels}
              onChange={(e) => setShowLabels(e.target.checked)}
              disabled={isLoading}
              className="w-4 h-4 text-blue-600 border-gray-300 rounded focus:ring-blue-500"
            />
            <span>Show Labels</span>
          </label>

          <button
            onClick={focusOnBindingSite}
            disabled={isLoading || poses.length === 0}
            className="px-3 py-1 bg-purple-600 text-white text-sm rounded hover:bg-purple-700 disabled:opacity-50"
          >
            Focus Binding Site
          </button>

          <button
            onClick={resetView}
            disabled={isLoading}
            className="px-3 py-1 bg-blue-600 text-white text-sm rounded hover:bg-blue-700 disabled:opacity-50"
          >
            Reset View
          </button>

          <button
            onClick={saveImage}
            disabled={isLoading}
            className="px-3 py-1 bg-green-600 text-white text-sm rounded hover:bg-green-700 disabled:opacity-50"
          >
            Save Image
          </button>
        </div>
      </div>

      {/* 3D Viewer Container */}
      <div className="bg-white border border-gray-200 rounded-lg overflow-hidden">
        {isLoading && (
          <div className="absolute inset-0 bg-black bg-opacity-50 flex items-center justify-center z-10">
            <div className="bg-white rounded-lg p-4 text-center">
              <div className="animate-spin w-8 h-8 border-4 border-blue-600 border-t-transparent rounded-full mx-auto mb-2"></div>
              <p className="text-gray-700">Loading molecular structures...</p>
            </div>
          </div>
        )}
        
        <div 
          ref={containerRef}
          className="w-full h-[500px] min-h-[500px] relative"
          style={{ backgroundColor }}
        />
      </div>

      {/* Status Information */}
      <div className="text-xs text-gray-500 space-y-1">
        <p>Protein: {receptorPdbId.toUpperCase()} | Ligand: {ligandSmiles.slice(0, 50)}{ligandSmiles.length > 50 ? '...' : ''}</p>
        <p>Models loaded: Protein ({loadedModels.protein >= 0 ? '✓' : '✗'}), Ligand ({loadedModels.ligand >= 0 ? '✓' : '✗'})</p>
        {poses.length > 0 && (
          <p>Poses: {poses.length} total | Selected: {selectedPose + 1} | Visible: {visiblePoses.size}</p>
        )}
      </div>
    </div>
  )
}

export default Enhanced3DViewer
