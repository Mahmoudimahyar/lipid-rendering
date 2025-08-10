import React, { useEffect, useRef, useState } from 'react'

/**
 * Extended molecule viewer specifically for docking visualization
 * Renders both protein structure and ligand poses
 */
function DockingMoleculeViewer({
  ligandSmiles,
  receptorPdbId,
  className = '',
  ...props
}) {
  // Debug logging for component lifecycle
  console.log('🚀 DockingMoleculeViewer: Component rendered with props:', { ligandSmiles, receptorPdbId })
  
  const containerRef = useRef(null)
  const viewerRef = useRef(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState(null)
  const [svgContent, setSvgContent] = useState('')

  // Initialize and load protein + ligand
  useEffect(() => {
    if (!containerRef.current) return

    const loadDockingVisualization = async () => {
      setIsLoading(true)
      setError(null)

      try {
        // Check if 3Dmol is available
        if (typeof window.$3Dmol === 'undefined') {
          setError('3Dmol.js library not available')
          setIsLoading(false)
          return
        }

        // Create a container for 3Dmol with specific dimensions to avoid WebGL issues
        const containerId = `docking-viewer-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`
        const container3D = `<div id="${containerId}" style="width: 100%; height: 400px; min-height: 400px; border-radius: 8px; background: #000; position: relative; overflow: hidden;"></div>`
        
        console.log('📦 Creating container with ID:', containerId)
        console.log('📏 Container HTML:', container3D)
        
        setSvgContent(container3D)

        // Wait longer for DOM to update and ensure container has proper dimensions
        setTimeout(async () => {
          const element = document.getElementById(containerId)
          if (!element) {
            setError('Could not create 3D container')
            setIsLoading(false)
            return
          }

          // Ensure the container has proper dimensions before creating viewer
          const rect = element.getBoundingClientRect()
          console.log('📐 Container dimensions:', {
            width: rect.width,
            height: rect.height,
            top: rect.top,
            left: rect.left,
            bottom: rect.bottom,
            right: rect.right
          })
          
          if (rect.width === 0 || rect.height === 0) {
            console.warn('❌ Container has zero dimensions, retrying...', rect)
            setTimeout(() => loadDockingVisualization(), 200)
            return
          } else {
            console.log('✅ Container has valid dimensions')
          }

          try {
            console.log('🧬 DockingMoleculeViewer: Starting protein+ligand visualization')
            
            // Create 3Dmol viewer with WebGL error handling
            const viewer = window.$3Dmol.createViewer(element, {
              defaultcolors: window.$3Dmol.rasmolElementColors,
              antialias: false,  // Disable antialiasing to avoid WebGL issues
              preserveDrawingBuffer: true
            })
            
            if (!viewer) {
              throw new Error('Failed to create 3Dmol viewer')
            }
            
            viewerRef.current = viewer
            console.log('✅ 3Dmol viewer created successfully')

            let proteinModelIndex = -1
            let ligandModelIndex = -1

            // Load protein if provided
            if (receptorPdbId) {
              proteinModelIndex = await loadProtein(viewer, receptorPdbId)
            }

            // Load ligand if provided
            if (ligandSmiles) {
              ligandModelIndex = await loadLigand(viewer, ligandSmiles)
            }

            // Apply distinct styling to protein and ligand
            applyMolecularStyling(viewer, proteinModelIndex, ligandModelIndex)

            // Center and render with extensive debugging
            try {
              console.log('🎯 Starting render process...')
              
              // Check viewer state before rendering
              console.log('🔍 Viewer state:', {
                hasViewer: !!viewer,
                modelCount: viewer.getModels ? viewer.getModels().length : 'unknown',
                element: element,
                elementDimensions: element.getBoundingClientRect()
              })
              
              // Check if viewer has models
              try {
                const models = viewer.getModels ? viewer.getModels() : []
                console.log('📊 Models in viewer:', models.length)
                models.forEach((model, index) => {
                  try {
                    const atoms = model.selectedAtoms ? model.selectedAtoms({}) : []
                    console.log(`   Model ${index}: ${atoms.length} atoms`)
                  } catch (e) {
                    console.log(`   Model ${index}: Could not count atoms`)
                  }
                })
              } catch (modelErr) {
                console.warn('⚠️ Could not check models:', modelErr)
              }
              
              // Try each render step individually with error handling
              console.log('🎯 Step 1: zoomTo()')
              viewer.zoomTo()
              console.log('✅ zoomTo completed')
              
              console.log('🎯 Step 2: center()')
              viewer.center()
              console.log('✅ center completed')
              
              console.log('🎯 Step 3: render()')
              viewer.render()
              console.log('✅ Initial render completed')
              
              // Check if canvas exists and has content
              setTimeout(() => {
                const canvas = element.querySelector('canvas')
                if (canvas) {
                  console.log('🖼️ Canvas found:', {
                    width: canvas.width,
                    height: canvas.height,
                    clientWidth: canvas.clientWidth,
                    clientHeight: canvas.clientHeight,
                    style: canvas.style.cssText
                  })
                  
                  // Try to get WebGL context info
                  try {
                    const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl')
                    if (gl) {
                      console.log('🎮 WebGL context found:', {
                        vendor: gl.getParameter(gl.VENDOR),
                        renderer: gl.getParameter(gl.RENDERER),
                        version: gl.getParameter(gl.VERSION),
                        maxTextureSize: gl.getParameter(gl.MAX_TEXTURE_SIZE),
                        viewportDims: gl.getParameter(gl.VIEWPORT)
                      })
                    } else {
                      console.error('❌ No WebGL context available!')
                    }
                  } catch (glErr) {
                    console.error('❌ WebGL context check failed:', glErr)
                  }
                } else {
                  console.error('❌ No canvas found in element!')
                }
                
                // Try secondary render
                if (viewerRef.current) {
                  try {
                    console.log('🎯 Secondary render...')
                    viewerRef.current.render()
                    console.log('✅ Secondary render completed')
                  } catch (renderErr) {
                    console.warn('⚠️ Secondary render failed:', renderErr)
                  }
                }
              }, 500)
              
            } catch (renderErr) {
              console.error('❌ Render failed:', renderErr)
              setError(`Rendering failed: ${renderErr.message}`)
            }

            setIsLoading(false)

          } catch (err) {
            console.error('Error creating docking visualization:', err)
            setError(err.message)
            setIsLoading(false)
          }
        }, 300)

      } catch (err) {
        console.error('Error initializing docking viewer:', err)
        setError(err.message)
        setIsLoading(false)
      }
    }

    loadDockingVisualization()

    // Cleanup
    return () => {
      if (viewerRef.current) {
        try {
          viewerRef.current.clear()
        } catch (e) {
          console.warn('Error during cleanup:', e)
        }
      }
    }
  }, [ligandSmiles, receptorPdbId])

  const loadProtein = async (viewer, pdbId) => {
    try {
      console.log(`Loading protein ${pdbId}...`)
      
      // Fetch protein structure from PDB
      const pdbUrl = `https://files.rcsb.org/download/${pdbId.toUpperCase()}.pdb`
      const response = await fetch(pdbUrl)
      
      if (!response.ok) {
        throw new Error(`Failed to fetch protein ${pdbId}: ${response.statusText}`)
      }
      
      const pdbData = await response.text()
      
      // Add protein model and return its index
      viewer.addModel(pdbData, 'pdb')
      const modelIndex = 0 // Protein is always the first model (index 0)
      
      console.log(`Protein ${pdbId} loaded successfully as model ${modelIndex}`)
      return modelIndex
      
    } catch (err) {
      console.error('Error loading protein:', err)
      return -1 // Return -1 to indicate failure
    }
  }

  const loadLigand = async (viewer, smiles) => {
    try {
      console.log(`Loading ligand from SMILES: ${smiles}`)
      
      // Try to use the ligand preparation endpoint
      try {
        const response = await fetch('/api/ligand/prepare', {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            smiles: smiles,
            num_conformers: 1
          })
        })

        if (response.ok) {
          const ligandData = await response.json()
          
          // Add ligand model from SDF and return its index
          viewer.addModel(ligandData.sdf, 'sdf')
          const modelIndex = 1 // Ligand is the second model (index 1)
          
          console.log(`Ligand loaded from backend successfully as model ${modelIndex}`)
          return modelIndex
        }
      } catch (backendErr) {
        console.warn('Backend ligand preparation failed, using fallback:', backendErr)
      }

      // Fallback: try adding SMILES directly (limited support)
      try {
        viewer.addModel(smiles, 'smiles')
        const modelIndex = 1 // Ligand is the second model (index 1)
        
        console.log(`Ligand loaded via SMILES fallback as model ${modelIndex}`)
        return modelIndex
      } catch (smilesErr) {
        console.error('Failed to load ligand via SMILES:', smilesErr)
        return -1
      }
      
    } catch (err) {
      console.error('Error loading ligand:', err)
      return -1 // Return -1 to indicate failure
    }
  }

  // Apply professional molecular styling based on industry standards
  const applyMolecularStyling = (viewer, proteinModelIndex, ligandModelIndex) => {
    try {
      console.log(`Applying styling - Protein: ${proteinModelIndex}, Ligand: ${ligandModelIndex}`)

      // Style protein (if loaded) - Professional approach: cartoon + surface
      if (proteinModelIndex >= 0) {
        console.log(`🧬 Styling protein model ${proteinModelIndex}`)
        
        // Primary protein representation: cartoon (secondary structure)
        viewer.setStyle({ model: proteinModelIndex }, {
          cartoon: {
            color: 'spectrum',  // Color by secondary structure
            opacity: 0.8
          }
        })

        console.log('✅ Applied protein cartoon styling')

        // Add transparent surface to show binding cavity (optional - can be heavy)
        try {
          viewer.addSurface(window.$3Dmol.SurfaceType.VDW, {
            opacity: 0.15,
            color: 'lightblue'
          }, { model: proteinModelIndex })
          console.log('✅ Applied protein surface')
        } catch (surfaceErr) {
          console.warn('⚠️ Could not add protein surface:', surfaceErr)
        }
      }

      // Style ligand (if loaded) - Professional approach: bright colored sticks + spheres
      if (ligandModelIndex >= 0) {
        console.log(`💊 Styling ligand model ${ligandModelIndex}`)
        
        viewer.setStyle({ model: ligandModelIndex }, {
          stick: {
            colorscheme: 'greenCarbon',  // Green carbons, standard heteroatom colors
            radius: 0.4,  // Make thicker for visibility
            opacity: 1.0
          },
          sphere: {
            colorscheme: 'greenCarbon',
            radius: 0.3,  // Make larger for visibility
            opacity: 0.9
          }
        })

        console.log('✅ Applied ligand stick+sphere styling')

        // Add a glow effect around the ligand for emphasis
        try {
          viewer.addSurface(window.$3Dmol.SurfaceType.VDW, {
            opacity: 0.4,
            color: 'yellow'
          }, { model: ligandModelIndex })
          console.log('✅ Applied ligand glow surface')
        } catch (surfaceErr) {
          console.warn('⚠️ Could not add ligand surface:', surfaceErr)
        }
      }

      // Alternative styling if models are not tracked properly
      if (proteinModelIndex < 0 && ligandModelIndex < 0) {
        console.warn('Model indices not available, applying fallback styling')
        
        // Try to style by model index (protein loaded first = 0, ligand = 1)
        try {
          // Assume protein is model 0 and ligand is model 1
          viewer.setStyle({ model: 0 }, {
            cartoon: {
              color: 'spectrum',
              opacity: 0.8
            }
          })
          
          viewer.setStyle({ model: 1 }, {
            stick: {
              colorscheme: 'greenCarbon',
              radius: 0.3
            },
            sphere: {
              colorscheme: 'greenCarbon',
              radius: 0.25,
              opacity: 0.8
            }
          })
          
          console.log('Applied fallback styling: model 0 as protein, model 1 as ligand')
        } catch (fallbackErr) {
          console.warn('Fallback styling failed:', fallbackErr)
        }
      }

    } catch (err) {
      console.error('Error applying molecular styling:', err)
      
      // Fallback: basic styling for all models
      viewer.setStyle({}, {
        stick: { radius: 0.2 },
        sphere: { radius: 0.3 }
      })
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
      {/* Content container */}
      <div 
        ref={containerRef}
        className="w-full h-96"
        dangerouslySetInnerHTML={{ __html: svgContent }}
      />
      
      {/* Loading overlay */}
      {isLoading && (
        <div className="absolute inset-0 flex items-center justify-center bg-black bg-opacity-75 rounded-lg">
          <div className="text-center text-white">
            <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-white mx-auto mb-4"></div>
            <p>Loading docking visualization...</p>
            <p className="text-sm opacity-75 mt-1">
              {receptorPdbId && `Protein: ${receptorPdbId}`}
              {receptorPdbId && ligandSmiles && ' + '}
              {ligandSmiles && `Ligand: ${ligandSmiles.substring(0, 20)}...`}
            </p>
          </div>
        </div>
      )}
      
      {/* Status overlay */}
      <div className="absolute bottom-4 right-4 bg-black bg-opacity-75 text-white rounded px-3 py-2 text-sm">
        <div>🧬 {receptorPdbId ? `Protein: ${receptorPdbId}` : 'No protein'}</div>
        <div>💊 {ligandSmiles ? 'Ligand loaded' : 'No ligand'}</div>
      </div>
    </div>
  )
}

export default DockingMoleculeViewer


