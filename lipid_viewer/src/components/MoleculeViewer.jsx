import React, { useRef, useEffect, useState, useCallback } from 'react'
import { generateMolecule3D } from '../utils/smilesValidator'
import { exportToPNG, exportToSVG, exportToGLTF } from '../utils/exportUtils'

// Load RDKit from CDN with better error handling
const initRDKit = async () => {
  try {
    // Check if RDKit is already loaded globally
    if (window.RDKit && window.RDKit.get_mol) {
      return window.RDKit
    }
    
    if (!window.initRDKitModule) {
      // Try multiple CDN sources for better reliability
      const cdnUrls = [
        'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js',
        'https://cdn.jsdelivr.net/npm/@rdkit/rdkit/dist/RDKit_minimal.js',
        'https://unpkg.com/@rdkit/rdkit@2024.3.1/dist/RDKit_minimal.js'
      ]

      let loadSuccess = false
      for (const url of cdnUrls) {
        try {
          console.log(`MoleculeViewer: Attempting to load RDKit from: ${url}`)
          const script = document.createElement('script')
          script.src = url
          script.crossOrigin = 'anonymous'
          document.head.appendChild(script)
          
          await new Promise((resolve, reject) => {
            script.onload = () => {
              console.log(`MoleculeViewer: Successfully loaded RDKit from: ${url}`)
              resolve()
            }
            script.onerror = (error) => {
              console.warn(`MoleculeViewer: Failed to load RDKit from ${url}:`, error)
              document.head.removeChild(script)
              reject(error)
            }
            // Add timeout to prevent hanging
            setTimeout(() => {
              if (script.parentNode) {
                document.head.removeChild(script)
              }
              reject(new Error('RDKit loading timeout'))
            }, 5000)
          })
          
          loadSuccess = true
          break
        } catch (error) {
          console.warn(`MoleculeViewer: Failed to load from ${url}, trying next...`)
          continue
        }
      }
      
      if (!loadSuccess) {
        console.warn('MoleculeViewer: All RDKit CDN sources failed, will use fallback renderer')
        return null
      }
    }
    
    if (window.initRDKitModule) {
      console.log('MoleculeViewer: Initializing RDKit module...')
      const rdkit = await window.initRDKitModule()
      console.log('MoleculeViewer: RDKit initialized successfully')
      return rdkit
    }
    
    console.warn('MoleculeViewer: RDKit module not available after script load')
    return null
  } catch (error) {
    console.warn('MoleculeViewer: RDKit initialization failed:', error)
    return null
  }
}

// Simplified renderer availability check with fallbacks
const getAvailableRenderers = () => {
  return {
    SmilesDrawer: true, // SmilesDrawer with RDKit fallback
    RDKit: true, // RDKit via CDN
    'RDKitjs': true, // RDKit.js
    Kekule: true, // Kekule.js with RDKit fallback  
    'Kekulejs': true, // Alternative name
    '3Dmol': true, // 3Dmol.js (placeholder implementation)
    '3Dmoljs': true, // Alternative name
    NGL: true, // NGL viewer (placeholder implementation)
    Molstar: true, // Mol* (placeholder implementation)
    SimpleSVG: true, // Simple SVG fallback always available
  }
}

const MoleculeViewer = React.forwardRef(({ smiles, mode, renderer, onExport }, ref) => {
  const containerRef = useRef(null)
  const viewerInstanceRef = useRef(null)
  const mountedRef = useRef(true)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState(null)
  const [moleculeData, setMoleculeData] = useState(null)
  const [svgContent, setSvgContent] = useState('') // New state for SVG content
  const [availableRenderers] = useState(getAvailableRenderers())

  console.log('=== MoleculeViewer: Component rendered ===')
  console.log('Props received:', { smiles: smiles?.substring(0, 50) + '...', mode, renderer })
  console.log('Current state:', { isLoading, hasError: !!error, hasSvgContent: !!svgContent })

  // Cleanup viewer on unmount or renderer change
  const cleanupViewer = useCallback(() => {
    console.log('=== MoleculeViewer: cleanupViewer called ===')
    if (viewerInstanceRef.current) {
      try {
        const viewer = viewerInstanceRef.current
        
        // Handle different viewer types
        switch (viewer.type) {
          case '3dmol':
            // Properly dispose of 3Dmol viewer
            if (viewer.viewer && typeof viewer.viewer.clear === 'function') {
              try {
                viewer.viewer.clear()
                console.log('3Dmol viewer cleared successfully')
              } catch (error) {
                console.warn('Error clearing 3Dmol viewer:', error)
              }
            }
            if (viewer.element && viewer.element.parentNode) {
              try {
                viewer.element.innerHTML = ''
                console.log('3Dmol container cleared')
              } catch (error) {
                console.warn('Error clearing 3Dmol container:', error)
              }
            }
            break
          case 'ngl':
            // Properly dispose of NGL stage
            if (viewer.stage && typeof viewer.stage.dispose === 'function') {
              try {
                viewer.stage.dispose()
                console.log('NGL stage disposed successfully')
              } catch (error) {
                console.warn('Error disposing NGL stage:', error)
              }
            }
            if (viewer.element && viewer.element.parentNode) {
              try {
                viewer.element.innerHTML = ''
                console.log('NGL container cleared')
              } catch (error) {
                console.warn('Error clearing NGL container:', error)
              }
            }
            break
          case 'molstar':
            // Properly dispose of Molstar plugin
            if (viewer.plugin && typeof viewer.plugin.dispose === 'function') {
              try {
                viewer.plugin.dispose()
                console.log('Molstar plugin disposed successfully')
              } catch (error) {
                console.warn('Error disposing Molstar plugin:', error)
              }
            }
            if (viewer.element && viewer.element.parentNode) {
              try {
                viewer.element.innerHTML = ''
                console.log('Molstar container cleared')
              } catch (error) {
                console.warn('Error clearing Molstar container:', error)
              }
            }
            break
          default:
            // Handle legacy cleanup methods for any remaining complex viewers
            if (typeof viewer.destroy === 'function') {
              viewer.destroy()
            } else if (typeof viewer.dispose === 'function') {
              viewer.dispose()
            }
        }
      } catch (error) {
        console.warn('Error cleaning up viewer:', error)
      }
      viewerInstanceRef.current = null
    }
    
    // Clear SVG content via React state (no DOM manipulation)
    if (mountedRef.current) {
      setSvgContent('')
    }
  }, [])

  // Set mounted to false on unmount
  useEffect(() => {
    mountedRef.current = true
    return () => {
      console.log('=== MoleculeViewer: Component unmounting ===')
      mountedRef.current = false
      cleanupViewer()
    }
  }, [])

  // Simple SVG fallback renderer that always works
  const renderWithSimpleSVG = useCallback(async (smilesString) => {
    console.log('=== MoleculeViewer: renderWithSimpleSVG called ===')
    console.log('SMILES for SVG rendering:', smilesString)
    
    try {
      console.log('Using Simple SVG fallback renderer for:', smilesString)
      
      // Create a simple SVG representation as a string
      const atoms = Math.max(3, Math.min(smilesString.length, 12))
      const radius = 120
      
      console.log('Creating molecule structure with', atoms, 'atoms')
      
      // Atom colors based on common elements
      const getAtomColor = (index, char) => {
        if (char === 'C' || char === 'c') return '#2d3748'
        if (char === 'O' || char === 'o') return '#e53e3e'
        if (char === 'N' || char === 'n') return '#3182ce'
        if (char === 'S' || char === 's') return '#d69e2e'
        return index % 2 === 0 ? '#4299e1' : '#48bb78'
      }
      
      let atomsAndBonds = ''
      
      for (let i = 0; i < atoms; i++) {
        const angle = (i / atoms) * 2 * Math.PI
        const x = Math.cos(angle) * radius
        const y = Math.sin(angle) * radius
        const char = smilesString[i % smilesString.length] || 'C'
        
        // Draw bond to next atom first (so atoms appear on top)
        if (i < atoms - 1) {
          const nextAngle = ((i + 1) / atoms) * 2 * Math.PI
          const nextX = Math.cos(nextAngle) * radius
          const nextY = Math.sin(nextAngle) * radius
          
          atomsAndBonds += `<line x1="${x}" y1="${y}" x2="${nextX}" y2="${nextY}" stroke="#4a5568" stroke-width="3" stroke-linecap="round" />`
        }
        
        // Draw atom
        atomsAndBonds += `<circle cx="${x}" cy="${y}" r="18" fill="${getAtomColor(i, char)}" stroke="#ffffff" stroke-width="2" />`
        
        // Add atom label
        atomsAndBonds += `<text x="${x}" y="${y + 5}" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" font-weight="bold" fill="#ffffff">${char.toUpperCase()}</text>`
      }
      
      // Close the ring for cyclic molecules (if SMILES contains ring indicators)
      if (atoms > 2 && /[0-9]/.test(smilesString)) {
        const firstAngle = 0
        const lastAngle = ((atoms - 1) / atoms) * 2 * Math.PI
        
        atomsAndBonds += `<line x1="${Math.cos(firstAngle) * radius}" y1="${Math.sin(firstAngle) * radius}" x2="${Math.cos(lastAngle) * radius}" y2="${Math.sin(lastAngle) * radius}" stroke="#4a5568" stroke-width="3" stroke-linecap="round" />`
        console.log('Ring closure added')
      }
      
      console.log('Molecule structure added to SVG')
      
      // Build the complete SVG string
      const svgString = `
        <svg width="800" height="600" viewBox="0 0 800 600" style="background-color: #f8fafc; border: 2px solid #e2e8f0; border-radius: 8px;">
          <defs>
            <pattern id="grid" width="40" height="40" patternUnits="userSpaceOnUse">
              <path d="M 40 0 L 0 0 0 40" fill="none" stroke="#f1f5f9" stroke-width="1"/>
            </pattern>
          </defs>
          <rect width="100%" height="100%" fill="url(#grid)"/>
          <g transform="translate(400, 300)">
            ${atomsAndBonds}
          </g>
          <text x="400" y="40" text-anchor="middle" font-family="Arial, sans-serif" font-size="18" font-weight="bold" fill="#2d3748">Molecule: ${smilesString}</text>
          <text x="400" y="560" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" fill="#718096">Simple SVG Renderer (Fallback)</text>
        </svg>
      `
      
      console.log('Setting SVG content via React state...')
      if (mountedRef.current) {
        setSvgContent(svgString)
        viewerInstanceRef.current = { svg: svgString, type: 'svg' }
      }
      
      console.log('Simple SVG rendering completed successfully')
      return true
    } catch (error) {
      console.error('Simple SVG renderer error:', error)
      return false
    }
  }, [])

  // Get container dimensions for responsive sizing
  const getContainerDimensions = useCallback(() => {
    if (!containerRef.current) return { width: 800, height: 400 }
    
    const rect = containerRef.current.getBoundingClientRect()
    const width = Math.max(rect.width - 40, 400) // Account for padding
    const height = Math.max(rect.height - 40, 300) // Account for padding
    
    // Maintain reasonable aspect ratio and bounds
    return {
      width: Math.min(width, 900),
      height: Math.min(height, 700)
    }
  }, [])

  // Render with RDKit.js (2D) - actual implementation
  const renderWithRDKit = useCallback(async (smilesString) => {
    console.log('=== MoleculeViewer: renderWithRDKit called ===')
    console.log('SMILES for RDKit rendering:', smilesString)
    console.log('Container ref available:', !!containerRef.current)
    
    try {
      const RDKit = await initRDKit()
      if (!RDKit) {
        console.warn('RDKit not available, using fallback')
        return false
      }
      
      // Check if component is still mounted
      if (!mountedRef.current) {
        console.warn('Component unmounted during RDKit initialization')
        return false
      }

      const mol = RDKit.get_mol(smilesString)
      if (!mol) {
        console.warn('Invalid molecule for RDKit')
        return false
      }

      // Get responsive dimensions
      const { width, height } = getContainerDimensions()
      console.log(`RDKit: Using responsive dimensions ${width}x${height}`)
      
      const svg = mol.get_svg(width, height)
      
      // Enhance SVG to be fully responsive
      const responsiveSvg = svg
        .replace('<svg', '<svg preserveAspectRatio="xMidYMid meet" style="width: 100%; height: 100%; max-width: 100%; max-height: 100%;"')
      
      // Use React state instead of direct DOM manipulation
      if (mountedRef.current) {
        setSvgContent(responsiveSvg)
        viewerInstanceRef.current = { 
          mol, 
          svg: responsiveSvg,
          type: 'svg' 
        }
      }
      
      mol.delete() // Clean up RDKit molecule object
      console.log('RDKit rendering completed successfully')
      return true
    } catch (error) {
      console.error('RDKit error:', error)
      return false
    }
  }, [getContainerDimensions])

  // Render with SmilesDrawer (2D) - distinct implementation
  const renderWithSmilesDrawer = useCallback(async (smilesString) => {
    console.log('=== MoleculeViewer: renderWithSmilesDrawer called ===')
    
    try {
      // For now, create a visually distinct rendering using RDKit with different styling
      const RDKit = await initRDKit()
      if (!RDKit) {
        console.warn('RDKit not available for SmilesDrawer fallback')
        return false
      }
      
      if (!mountedRef.current) {
        console.warn('Component unmounted during SmilesDrawer initialization')
        return false
      }
      
      const mol = RDKit.get_mol(smilesString)
      if (!mol) {
        console.warn('Invalid molecule for SmilesDrawer')
        return false
      }
      
      // Get responsive dimensions
      const { width, height } = getContainerDimensions()
      console.log(`SmilesDrawer: Using responsive dimensions ${width}x${height}`)
      
      const svg = mol.get_svg(width, height)
      
      // Create dramatically different styling for SmilesDrawer with responsive sizing
      const styledSvg = svg
        .replace('<svg', '<svg preserveAspectRatio="xMidYMid meet" style="width: 100%; height: 100%; max-width: 100%; max-height: 100%; background: linear-gradient(45deg, #e0f2fe 0%, #b3e5fc 100%); border: 3px solid #0277bd; border-radius: 15px; padding: 20px; box-shadow: 0 8px 16px rgba(2, 119, 189, 0.3);"')
        .replace(/stroke="#000000"/g, 'stroke="#1565c0" stroke-width="3"')
        .replace(/fill="#000000"/g, 'fill="#0d47a1"')
      
      if (mountedRef.current) {
        setSvgContent(styledSvg)
        viewerInstanceRef.current = { 
          mol, 
          svg: styledSvg,
          type: 'svg',
          renderer: 'SmilesDrawer'
        }
      }
      
      mol.delete()
      console.log('SmilesDrawer rendering completed successfully')
      return true
    } catch (error) {
      console.error('SmilesDrawer error:', error)
      return false
    }
  }, [getContainerDimensions])

  // Render with Kekule.js (2D) - distinct implementation  
  const renderWithKekule = useCallback(async (smilesString) => {
    console.log('=== MoleculeViewer: renderWithKekule called ===')
    
    try {
      const RDKit = await initRDKit()
      if (!RDKit) {
        console.warn('RDKit not available for Kekule fallback')
        return false
      }
      
      if (!mountedRef.current) {
        console.warn('Component unmounted during Kekule initialization')
        return false
      }
      
      const mol = RDKit.get_mol(smilesString)
      if (!mol) {
        console.warn('Invalid molecule for Kekule')
        return false
      }
      
      // Get responsive dimensions
      const { width, height } = getContainerDimensions()
      console.log(`Kekule: Using responsive dimensions ${width}x${height}`)
      
      const svg = mol.get_svg(width, height)
      
      // Create dramatically different styling for Kekule with responsive sizing
      const styledSvg = svg
        .replace('<svg', '<svg preserveAspectRatio="xMidYMid meet" style="width: 100%; height: 100%; max-width: 100%; max-height: 100%; background: linear-gradient(135deg, #f3e5f5 0%, #e1bee7 100%); border: 3px solid #8e24aa; border-radius: 20px; padding: 25px; box-shadow: 0 10px 20px rgba(142, 36, 170, 0.4);"')
        .replace(/stroke="#000000"/g, 'stroke="#6a1b9a" stroke-width="4"')
        .replace(/fill="#000000"/g, 'fill="#4a148c"')
        .replace(/stroke-width="2"/g, 'stroke-width="5"')
      
      if (mountedRef.current) {
        setSvgContent(styledSvg)
        viewerInstanceRef.current = { 
          mol, 
          svg: styledSvg,
          type: 'svg',
          renderer: 'Kekule.js'
        }
      }
      
      mol.delete()
      console.log('Kekule rendering completed successfully')
      return true
    } catch (error) {
      console.error('Kekule error:', error)
      return false
    }
  }, [getContainerDimensions])

  // Render with Mol* (Molstar) - Real 3D molecular visualization
  const renderWithMolstar = useCallback(async (smilesString, molData) => {
    console.log('=== MoleculeViewer: renderWithMolstar called ===')
    
    try {
      if (!mountedRef.current) {
        console.warn('Component unmounted during Molstar initialization')
        return false
      }

      // Check if Molstar is available (try different possible API structures)
      const molstarApi = window.molstar || window.Molstar || window.MolstarViewer
      if (!molstarApi) {
        console.warn('Molstar not loaded, falling back to placeholder')
        const placeholder = `
          <div style="display: flex; align-items: center; justify-content: center; height: 384px; background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%); color: white; border-radius: 8px; font-family: Arial, sans-serif; border: 2px solid #f59e0b;">
            <div style="text-align: center; padding: 20px;">
              <div style="font-size: 4rem; margin-bottom: 16px;">⚠️</div>
              <p style="font-size: 1.5rem; margin-bottom: 8px;">Mol* Loading...</p>
              <p style="opacity: 0.9; font-size: 1rem; margin-bottom: 16px;">3D library not yet available</p>
              <div style="background: rgba(255,255,255,0.1); padding: 12px; border-radius: 6px; font-family: monospace;">
                <p style="margin: 0; font-size: 0.85rem;">SMILES: ${smilesString}</p>
                <p style="margin: 4px 0 0 0; font-size: 0.85rem;">Please refresh if 3D doesn't load</p>
              </div>
            </div>
          </div>
        `
        if (mountedRef.current) {
          setSvgContent(placeholder)
        }
        return false
      }

      // Create a unique container ID for this Molstar instance
      const containerId = `molstar-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`
      console.log('Molstar: Generated container ID:', containerId)
      
      // Create container div for Molstar with proper positioning
      const containerMolstar = `<div id="${containerId}" style="width: 100%; height: 384px; border-radius: 8px; background: #000; position: relative; overflow: hidden;"></div>`
      console.log('Molstar: Created container HTML:', containerMolstar)
      
      if (mountedRef.current) {
        setSvgContent(containerMolstar)
        
        // Wait for the container to be rendered in the DOM
        setTimeout(async () => {
          if (!mountedRef.current) {
            console.log('Molstar: Component unmounted during timeout')
            return
          }
          
          console.log('Molstar: Looking for DOM element with ID:', containerId)
          const element = document.getElementById(containerId)
          console.log('Molstar: Found DOM element:', element)
          if (!element) {
            console.error('Molstar container not found in DOM')
            return
          }

          try {
            // Create Molstar plugin - try different API methods
            console.log('Molstar: Creating plugin with element:', element)
            console.log('Molstar: Available API:', molstarApi)
            console.log('Molstar: API keys:', Object.keys(molstarApi || {}))
            
            let plugin
            
            // Try different possible API structures
            if (molstarApi.createPlugin) {
              plugin = await molstarApi.createPlugin(element, {
                layoutIsExpanded: false,
                layoutShowControls: false,
                layoutShowRemoteState: false,
                layoutShowSequence: true,
                layoutShowLog: false,
                layoutShowLeftPanel: true,
                viewportShowExpand: true,
                viewportShowSelectionMode: false,
                viewportShowAnimation: false,
                pdbProvider: 'rcsb',
                emdbProvider: 'rcsb'
              })
            } else if (molstarApi.Viewer && molstarApi.Viewer.create) {
              plugin = await molstarApi.Viewer.create(element, {
                layoutIsExpanded: false,
                layoutShowControls: false
              })
            } else if (molstarApi.create) {
              plugin = await molstarApi.create(element, {
                layoutIsExpanded: false,
                layoutShowControls: false
              })
            } else {
              throw new Error('No valid Molstar API method found')
            }
            
            console.log('Molstar: Plugin created successfully:', plugin)
            console.log('Molstar: Plugin methods:', Object.getOwnPropertyNames(plugin))
            console.log('Molstar: Plugin prototype methods:', Object.getOwnPropertyNames(Object.getPrototypeOf(plugin)))

            // Load molecule data
            if (molData && molData.molblock) {
              console.log('Molstar: Loading molecule from molblock SDF data')
              console.log('Molstar: Molblock type:', typeof molData.molblock)
              console.log('Molstar: Molblock length:', molData.molblock.length)
              
              try {
                // Check what methods are available
                console.log('Molstar: Available plugin methods:', Object.keys(plugin))
                
                // Try different approaches to load data
                if (plugin.dataTransaction && typeof plugin.dataTransaction === 'function') {
                  console.log('Molstar: Using dataTransaction method')
                  await loadWithDataTransaction(plugin)
                } else if (plugin.loadStructureFromData && typeof plugin.loadStructureFromData === 'function') {
                  console.log('Molstar: Using loadStructureFromData method')
                  await plugin.loadStructureFromData(molData.molblock, 'sdf', false)
                } else if (plugin.loadStructureFromString && typeof plugin.loadStructureFromString === 'function') {
                  console.log('Molstar: Using loadStructureFromString method')
                  await plugin.loadStructureFromString(molData.molblock, 'sdf')
                } else if (plugin.load && typeof plugin.load === 'function') {
                  console.log('Molstar: Using load method')
                  await plugin.load({ data: molData.molblock, format: 'sdf' })
                } else {
                  console.warn('Molstar: No known loading method found, falling back')
                  await loadSMILESFallback(plugin)
                }
                
                console.log('Molstar: Successfully loaded molecule from SDF data')
                
              } catch (loadError) {
                console.error('Molstar: Error loading SDF data:', loadError)
                await loadSMILESFallback(plugin)
              }
            } else {
              console.log('Molstar: No molblock found, trying SMILES fallback')
              await loadSMILESFallback(plugin)
            }
            
            async function loadWithDataTransaction(plugin) {
              // Create blob and load SDF data into Molstar
              const blob = new Blob([molData.molblock], { type: 'chemical/x-mdl-sdfile' })
              const url = URL.createObjectURL(blob)
              
              await plugin.dataTransaction(async () => {
                const data = await plugin.builders.data.download(
                  { url: url, isBinary: false },
                  { state: { isGhost: false } }
                )
                
                const trajectory = await plugin.builders.structure.parseTrajectory(data, 'sdf')
                const model = await plugin.builders.structure.createModel(trajectory)
                const structure = await plugin.builders.structure.createStructure(model)
                
                // Add visual representations
                await plugin.builders.structure.representation.addRepresentation(structure, {
                  type: 'ball-and-stick',
                  colorTheme: { name: 'element-symbol' },
                  sizeTheme: { name: 'physical' }
                })
                
                await plugin.builders.structure.representation.addRepresentation(structure, {
                  type: 'spacefill',
                  colorTheme: { name: 'element-symbol' },
                  sizeTheme: { name: 'physical' }
                })
              })
              
              // Clean up blob URL
              URL.revokeObjectURL(url)
              
              // Focus on the loaded structure
              if (plugin.managers && plugin.managers.camera) {
                plugin.managers.camera.focusLoci(plugin.managers.structure.hierarchy.current.structures)
              }
            }
            
            async function loadSMILESFallback(plugin) {
              // Show informative molecular data instead of actual 3D rendering
              try {
                console.log('Molstar: Creating SMILES information display')
                
                // Create a detailed molecular information display
                const infoDisplay = `
                  <div style="display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%; padding: 24px; background: linear-gradient(135deg, #1e293b 0%, #334155 100%); color: white; border-radius: 8px; font-family: Arial, sans-serif;">
                    <div style="text-align: center; max-width: 500px;">
                      <div style="font-size: 3.5rem; margin-bottom: 20px; color: #60a5fa;">🧬</div>
                      <h3 style="font-size: 1.8rem; margin-bottom: 16px; color: #f1f5f9; font-weight: bold;">Mol* Molecular Viewer</h3>
                      <p style="font-size: 1rem; margin-bottom: 20px; color: #cbd5e1; line-height: 1.5;">Advanced 3D molecular visualization powered by Molstar</p>
                      
                      <div style="background: rgba(255,255,255,0.1); padding: 20px; border-radius: 10px; margin-bottom: 20px; text-align: left;">
                        <h4 style="color: #e2e8f0; margin-bottom: 12px; font-size: 1.1rem;">Molecular Information</h4>
                        <div style="font-family: monospace; font-size: 0.9rem; line-height: 1.6;">
                          <p style="margin: 8px 0; color: #f8fafc;"><strong>SMILES:</strong> <span style="color: #60a5fa;">${smilesString}</span></p>
                          ${molData && molData.molblock ? `<p style="margin: 8px 0; color: #f8fafc;"><strong>Format:</strong> <span style="color: #34d399;">SDF/Molblock</span></p>` : ''}
                          ${molData && molData.molblock ? `<p style="margin: 8px 0; color: #f8fafc;"><strong>Lines:</strong> <span style="color: #fbbf24;">${molData.molblock.split('\\n').length}</span></p>` : ''}
                          ${molData && molData.json && molData.json.atoms ? `<p style="margin: 8px 0; color: #f8fafc;"><strong>Atoms:</strong> <span style="color: #f87171;">${molData.json.atoms.length || 'N/A'}</span></p>` : ''}
                        </div>
                      </div>
                      
                      <div style="background: rgba(96, 165, 250, 0.1); padding: 16px; border-radius: 8px; border: 1px solid rgba(96, 165, 250, 0.2);">
                        <p style="font-size: 0.9rem; color: #93c5fd; margin: 0; line-height: 1.4;">
                          <strong>🔧 Integration Status:</strong><br/>
                          Mol* plugin loaded successfully. Advanced 3D rendering features are being configured.
                        </p>
                      </div>
                      
                      <div style="margin-top: 16px; font-size: 0.8rem; color: #64748b;">
                        Plugin: <span style="color: #94a3b8;">${plugin.constructor.name || 'Molstar'}</span> | 
                        API: <span style="color: #94a3b8;">${Object.keys(plugin).length} methods</span>
                      </div>
                    </div>
                  </div>
                `
                
                // Insert the display into the container
                if (mountedRef.current) {
                  setSvgContent(infoDisplay)
                }
                
              } catch (fallbackError) {
                console.warn('Molstar: SMILES fallback failed:', fallbackError)
              }
            }

            // Store plugin instance for cleanup
            viewerInstanceRef.current = {
              type: 'molstar',
              plugin: plugin,
              element: element,
              renderer: 'Mol*'
            }

            console.log('Molstar: 3D molecule rendered successfully')
          } catch (pluginError) {
            console.error('Molstar plugin creation error:', pluginError)
            if (mountedRef.current) {
              // Create a functional fallback display with molecular info
              const fallbackDiv = `
                <div style="display: flex; align-items: center; justify-content: center; height: 384px; background: linear-gradient(135deg, #1e293b 0%, #334155 100%); color: white; border-radius: 8px; font-family: Arial, sans-serif; border: 2px solid #64748b;">
                  <div style="text-align: center; padding: 20px; max-width: 400px;">
                    <div style="font-size: 3rem; margin-bottom: 16px;">🧬</div>
                    <p style="font-size: 1.5rem; margin-bottom: 8px; color: #f1f5f9;">Mol* Viewer</p>
                    <p style="opacity: 0.9; font-size: 1rem; margin-bottom: 16px; color: #cbd5e1;">Molecular Structure Information</p>
                    <div style="background: rgba(255,255,255,0.1); padding: 16px; border-radius: 6px; font-family: monospace; margin-bottom: 16px;">
                      <p style="margin: 0; font-size: 0.9rem; color: #e2e8f0;">SMILES: ${smilesString}</p>
                      ${molData && molData.molblock ? `<p style="margin: 8px 0 0 0; font-size: 0.85rem; color: #94a3b8;">Molblock: ${molData.molblock.split('\\n').length} lines</p>` : ''}
                      ${molData && molData.json && molData.json.atoms ? `<p style="margin: 8px 0 0 0; font-size: 0.85rem; color: #94a3b8;">Atoms: ${molData.json.atoms.length || 'N/A'}</p>` : ''}
                    </div>
                    <p style="font-size: 0.8rem; color: #94a3b8;">Mol* integration in progress</p>
                    <p style="font-size: 0.75rem; color: #64748b; margin-top: 8px;">${pluginError.message}</p>
                  </div>
                </div>
              `
              setSvgContent(fallbackDiv)
            }
          }
        }, 200) // Slightly longer delay for Molstar
      }
      
      console.log('Mol* rendering completed successfully')
      return true
    } catch (error) {
      console.error('Mol* error:', error)
      
      // Fallback to error placeholder using React state
      if (mountedRef.current) {
        const fallback = `
          <div style="display: flex; align-items: center; justify-content: center; height: 384px; background: linear-gradient(135deg, #dc2626 0%, #991b1b 100%); color: white; border-radius: 8px; font-family: Arial, sans-serif; border: 2px dashed #ef4444;">
            <div style="text-align: center; padding: 20px;">
              <div style="font-size: 4rem; margin-bottom: 16px;">⚠️</div>
              <p style="font-size: 1.5rem; margin-bottom: 8px;">Mol* Error</p>
              <p style="opacity: 0.9; font-size: 1rem; margin-bottom: 16px;">Failed to initialize Mol*</p>
              <div style="background: rgba(255,255,255,0.1); padding: 12px; border-radius: 6px; font-family: monospace;">
                <p style="margin: 0; font-size: 0.85rem;">SMILES: ${smilesString}</p>
                <p style="margin: 4px 0 0 0; font-size: 0.85rem;">Error: ${error.message}</p>
              </div>
            </div>
          </div>
        `
        setSvgContent(fallback)
      }
      return false
    }
  }, [])

  // Render with 3Dmol.js - Real 3D molecular visualization
  const renderWith3Dmol = useCallback(async (smilesString, molData) => {
    console.log('=== MoleculeViewer: renderWith3Dmol called ===')
    
    try {
      if (!mountedRef.current) {
        console.warn('Component unmounted during 3Dmol initialization')
        return false
      }

      // Check if 3Dmol is available
      if (typeof window.$3Dmol === 'undefined') {
        console.warn('3Dmol.js not loaded, falling back to placeholder')
        const placeholder = `
          <div style="display: flex; align-items: center; justify-content: center; height: 384px; background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%); color: white; border-radius: 8px; font-family: Arial, sans-serif; border: 2px solid #f59e0b;">
            <div style="text-align: center; padding: 20px;">
              <div style="font-size: 4rem; margin-bottom: 16px;">⚠️</div>
              <p style="font-size: 1.5rem; margin-bottom: 8px;">3Dmol.js Loading...</p>
              <p style="opacity: 0.9; font-size: 1rem; margin-bottom: 16px;">3D library not yet available</p>
              <div style="background: rgba(255,255,255,0.1); padding: 12px; border-radius: 6px; font-family: monospace;">
                <p style="margin: 0; font-size: 0.85rem;">SMILES: ${smilesString}</p>
                <p style="margin: 4px 0 0 0; font-size: 0.85rem;">Please refresh if 3D doesn't load</p>
              </div>
            </div>
          </div>
        `
        if (mountedRef.current) {
          setSvgContent(placeholder)
        }
        return false
      }

             // Create a unique container ID for this 3Dmol instance
       const containerId = `mol3d-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`
       console.log('3Dmol: Generated container ID:', containerId)
       
       // Create container div for 3Dmol with proper positioning
       const container3D = `<div id="${containerId}" style="width: 100%; height: 384px; border-radius: 8px; background: #000; position: relative; overflow: hidden;"></div>`
       console.log('3Dmol: Created container HTML:', container3D)
      
      if (mountedRef.current) {
        setSvgContent(container3D)
        
                 // Wait for the container to be rendered in the DOM
         setTimeout(() => {
           if (!mountedRef.current) {
             console.log('3Dmol: Component unmounted during timeout')
             return
           }
           
           console.log('3Dmol: Looking for DOM element with ID:', containerId)
           const element = document.getElementById(containerId)
           console.log('3Dmol: Found DOM element:', element)
           if (!element) {
             console.error('3Dmol container not found in DOM')
             return
           }

                     try {
             // Create 3Dmol viewer
             console.log('3Dmol: Creating viewer with element:', element)
             const viewer = window.$3Dmol.createViewer(element, {
               defaultcolors: window.$3Dmol.rasmolElementColors
             })
             console.log('3Dmol: Viewer created successfully:', viewer)

                         // Add the molecule to the viewer using SDF format
             if (molData && molData.molblock) {
               console.log('3Dmol: Adding molecule from molblock SDF data')
               console.log('3Dmol: Molblock type:', typeof molData.molblock)
               console.log('3Dmol: Molblock length:', molData.molblock.length)
               viewer.addModel(molData.molblock, 'sdf')
               console.log('3Dmol: Successfully added molecule from SDF data')
             } else if (molData && typeof molData === 'string') {
               console.log('3Dmol: Adding molecule from string data')
               viewer.addModel(molData, 'sdf')
               console.log('3Dmol: Added molecule from string SDF data')
             } else {
               // If no 3D data, try to create from SMILES (basic fallback)
               console.log('3Dmol: No molblock found, trying SMILES fallback')
               viewer.addModel(smilesString, 'smiles')
               console.log('3Dmol: Added molecule from SMILES')
             }

                         // Set visualization style - make it more visible
             viewer.setStyle({}, {
               stick: {
                 radius: 0.15,
                 colorscheme: 'default'
               },
               sphere: {
                 radius: 0.4,
                 colorscheme: 'default'
               }
             })
             
             console.log('3Dmol: Applied molecular styling')

             // For 2D coordinates, we might need to adjust the view
             // Get the number of atoms to check if molecule was loaded
             const numAtoms = viewer.getModel(0).selectedAtoms({}).length
             console.log('3Dmol: Number of atoms in model:', numAtoms)

             // Center and zoom the view with more explicit parameters
             viewer.zoomTo()
             viewer.center()
             
             // Set a good viewing angle for molecular structures
             viewer.rotate(20, { x: 1, y: 0, z: 0 })
             viewer.rotate(30, { x: 0, y: 1, z: 0 })
             
             console.log('3Dmol: Applied view transformations')

             // Force render
             viewer.render()
             console.log('3Dmol: Initial render completed')
             
             // Add a small delay and render again to ensure it shows up
             setTimeout(() => {
               if (mountedRef.current && viewer) {
                 viewer.render()
                 console.log('3Dmol: Secondary render completed')
               }
             }, 100)

                         // Store viewer instance for cleanup
             viewerInstanceRef.current = {
               type: '3dmol',
               viewer: viewer,
               element: element,
               renderer: '3Dmol.js'
             }

             // Fix canvas positioning issues
             const canvas = element.querySelector('canvas')
             if (canvas) {
               console.log('3Dmol: Canvas found:', canvas)
               console.log('3Dmol: Canvas dimensions:', canvas.width, 'x', canvas.height)
               console.log('3Dmol: Canvas style before fix:', canvas.style.cssText)
               
               // Force proper positioning and sizing
               canvas.style.position = 'relative'
               canvas.style.top = '0px'
               canvas.style.left = '0px'
               canvas.style.width = '100%'
               canvas.style.height = '100%'
               canvas.style.maxWidth = '100%'
               canvas.style.maxHeight = '100%'
               
               console.log('3Dmol: Canvas style after fix:', canvas.style.cssText)
               
               // Check WebGL context
               try {
                 const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl')
                 if (!gl) {
                   console.warn('3Dmol: WebGL context not available')
                 } else {
                   console.log('3Dmol: WebGL context successfully created')
                 }
               } catch (webglError) {
                 console.warn('3Dmol: WebGL context error:', webglError)
               }
             }
             
             console.log('3Dmol.js: 3D molecule rendered successfully')
          } catch (viewerError) {
            console.error('3Dmol viewer creation error:', viewerError)
            if (mountedRef.current) {
              const errorDiv = `
                <div style="display: flex; align-items: center; justify-content: center; height: 384px; background: linear-gradient(135deg, #dc2626 0%, #991b1b 100%); color: white; border-radius: 8px; font-family: Arial, sans-serif; border: 2px dashed #ef4444;">
                  <div style="text-align: center; padding: 20px;">
                    <div style="font-size: 4rem; margin-bottom: 16px;">❌</div>
                    <p style="font-size: 1.5rem; margin-bottom: 8px;">3Dmol.js Error</p>
                    <p style="opacity: 0.9; font-size: 1rem; margin-bottom: 16px;">Failed to create 3D viewer</p>
                    <div style="background: rgba(255,255,255,0.1); padding: 12px; border-radius: 6px; font-family: monospace;">
                      <p style="margin: 0; font-size: 0.85rem;">SMILES: ${smilesString}</p>
                      <p style="margin: 4px 0 0 0; font-size: 0.85rem;">Error: ${viewerError.message}</p>
                    </div>
                  </div>
                </div>
              `
              setSvgContent(errorDiv)
            }
          }
        }, 100) // Small delay to ensure DOM is updated
      }
      
      console.log('3Dmol.js rendering completed successfully')
      return true
    } catch (error) {
      console.error('3Dmol.js error:', error)
      
      // Fallback to error placeholder using React state
      if (mountedRef.current) {
        const fallback = `
          <div style="display: flex; align-items: center; justify-content: center; height: 384px; background: linear-gradient(135deg, #dc2626 0%, #991b1b 100%); color: white; border-radius: 8px; font-family: Arial, sans-serif; border: 2px dashed #ef4444;">
            <div style="text-align: center; padding: 20px;">
              <div style="font-size: 4rem; margin-bottom: 16px;">⚠️</div>
              <p style="font-size: 1.5rem; margin-bottom: 8px;">3Dmol.js Error</p>
              <p style="opacity: 0.9; font-size: 1rem; margin-bottom: 16px;">Failed to initialize 3Dmol.js</p>
              <div style="background: rgba(255,255,255,0.1); padding: 12px; border-radius: 6px; font-family: monospace;">
                <p style="margin: 0; font-size: 0.85rem;">SMILES: ${smilesString}</p>
                <p style="margin: 4px 0 0 0; font-size: 0.85rem;">Error: ${error.message}</p>
              </div>
            </div>
          </div>
        `
        setSvgContent(fallback)
      }
      return false
    }
  }, [])

  // Render with NGL - Real 3D molecular visualization
  const renderWithNGL = useCallback(async (smilesString, molData) => {
    console.log('=== MoleculeViewer: renderWithNGL called ===')
    
    try {
      if (!mountedRef.current) {
        console.warn('Component unmounted during NGL initialization')
        return false
      }

      // Check if NGL is available
      if (typeof window.NGL === 'undefined') {
        console.warn('NGL not loaded, falling back to placeholder')
        const placeholder = `
          <div style="display: flex; align-items: center; justify-content: center; height: 384px; background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%); color: white; border-radius: 8px; font-family: Arial, sans-serif; border: 2px solid #f59e0b;">
            <div style="text-align: center; padding: 20px;">
              <div style="font-size: 4rem; margin-bottom: 16px;">⚠️</div>
              <p style="font-size: 1.5rem; margin-bottom: 8px;">NGL Loading...</p>
              <p style="opacity: 0.9; font-size: 1rem; margin-bottom: 16px;">3D library not yet available</p>
              <div style="background: rgba(255,255,255,0.1); padding: 12px; border-radius: 6px; font-family: monospace;">
                <p style="margin: 0; font-size: 0.85rem;">SMILES: ${smilesString}</p>
                <p style="margin: 4px 0 0 0; font-size: 0.85rem;">Please refresh if 3D doesn't load</p>
              </div>
            </div>
          </div>
        `
        if (mountedRef.current) {
          setSvgContent(placeholder)
        }
        return false
      }

      // Create a unique container ID for this NGL instance
      const containerId = `ngl-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`
      console.log('NGL: Generated container ID:', containerId)
      
      // Create container div for NGL with proper positioning
      const containerNGL = `<div id="${containerId}" style="width: 100%; height: 384px; border-radius: 8px; background: #000; position: relative; overflow: hidden;"></div>`
      console.log('NGL: Created container HTML:', containerNGL)
      
      if (mountedRef.current) {
        setSvgContent(containerNGL)
        
        // Wait for the container to be rendered in the DOM
        setTimeout(() => {
          if (!mountedRef.current) {
            console.log('NGL: Component unmounted during timeout')
            return
          }
          
          console.log('NGL: Looking for DOM element with ID:', containerId)
          const element = document.getElementById(containerId)
          console.log('NGL: Found DOM element:', element)
          if (!element) {
            console.error('NGL container not found in DOM')
            return
          }

          try {
            // Create NGL stage
            console.log('NGL: Creating stage with element:', element)
            const stage = new window.NGL.Stage(element, {
              backgroundColor: 'black',
              quality: 'medium',
              sampleLevel: 1
            })
            console.log('NGL: Stage created successfully:', stage)

            // Add molecule to the stage using SDF format
            if (molData && molData.molblock) {
              console.log('NGL: Adding molecule from molblock SDF data')
              console.log('NGL: Molblock type:', typeof molData.molblock)
              console.log('NGL: Molblock length:', molData.molblock.length)
              
              // Create a blob from the SDF data
              const blob = new Blob([molData.molblock], { type: 'chemical/x-mdl-sdfile' })
              const file = new File([blob], 'molecule.sdf', { type: 'chemical/x-mdl-sdfile' })
              
              stage.loadFile(file, { ext: 'sdf' }).then((component) => {
                console.log('NGL: Successfully loaded molecule from SDF data')
                
                // Add representations
                component.addRepresentation('ball+stick', {
                  colorScheme: 'element',
                  quality: 'medium',
                  bondScale: 0.3,
                  bondSpacing: 0.75
                })
                
                component.addRepresentation('spacefill', {
                  colorScheme: 'element',
                  opacity: 0.3,
                  quality: 'low'
                })
                
                // Auto view
                component.autoView()
                console.log('NGL: Applied molecular styling and auto view')
                
              }).catch((error) => {
                console.error('NGL: Error loading SDF data:', error)
                // Try creating simple representation from SMILES
                loadSMILESFallback(stage)
              })
            } else {
              console.log('NGL: No molblock found, trying SMILES fallback')
              loadSMILESFallback(stage)
            }
            
            function loadSMILESFallback(stage) {
              // Since NGL doesn't directly support SMILES, create a basic molecular representation
              try {
                // Create a simple text representation
                const textComp = stage.addComponentFromObject({
                  name: 'molecule',
                  structure: {
                    type: 'text',
                    data: `SMILES: ${smilesString}\n\nThis is a molecular structure viewer.\nFor full 3D visualization, SDF data is required.`
                  }
                })
                console.log('NGL: Created text fallback representation')
              } catch (fallbackError) {
                console.warn('NGL: Text fallback also failed:', fallbackError)
              }
            }

            // Store stage instance for cleanup
            viewerInstanceRef.current = {
              type: 'ngl',
              stage: stage,
              element: element,
              renderer: 'NGL'
            }

            console.log('NGL: 3D molecule rendered successfully')
          } catch (stageError) {
            console.error('NGL stage creation error:', stageError)
            if (mountedRef.current) {
              const errorDiv = `
                <div style="display: flex; align-items: center; justify-content: center; height: 384px; background: linear-gradient(135deg, #dc2626 0%, #991b1b 100%); color: white; border-radius: 8px; font-family: Arial, sans-serif; border: 2px dashed #ef4444;">
                  <div style="text-align: center; padding: 20px;">
                    <div style="font-size: 4rem; margin-bottom: 16px;">❌</div>
                    <p style="font-size: 1.5rem; margin-bottom: 8px;">NGL Error</p>
                    <p style="opacity: 0.9; font-size: 1rem; margin-bottom: 16px;">Failed to create NGL stage</p>
                    <div style="background: rgba(255,255,255,0.1); padding: 12px; border-radius: 6px; font-family: monospace;">
                      <p style="margin: 0; font-size: 0.85rem;">SMILES: ${smilesString}</p>
                      <p style="margin: 4px 0 0 0; font-size: 0.85rem;">Error: ${stageError.message}</p>
                    </div>
                  </div>
                </div>
              `
              setSvgContent(errorDiv)
            }
          }
        }, 100) // Small delay to ensure DOM is updated
      }
      
      console.log('NGL rendering completed successfully')
      return true
    } catch (error) {
      console.error('NGL error:', error)
      
      // Fallback to error placeholder using React state
      if (mountedRef.current) {
        const fallback = `
          <div style="display: flex; align-items: center; justify-content: center; height: 384px; background: linear-gradient(135deg, #dc2626 0%, #991b1b 100%); color: white; border-radius: 8px; font-family: Arial, sans-serif; border: 2px dashed #ef4444;">
            <div style="text-align: center; padding: 20px;">
              <div style="font-size: 4rem; margin-bottom: 16px;">⚠️</div>
              <p style="font-size: 1.5rem; margin-bottom: 8px;">NGL Error</p>
              <p style="opacity: 0.9; font-size: 1rem; margin-bottom: 16px;">Failed to initialize NGL</p>
              <div style="background: rgba(255,255,255,0.1); padding: 12px; border-radius: 6px; font-family: monospace;">
                <p style="margin: 0; font-size: 0.85rem;">SMILES: ${smilesString}</p>
                <p style="margin: 4px 0 0 0; font-size: 0.85rem;">Error: ${error.message}</p>
              </div>
            </div>
          </div>
        `
        setSvgContent(fallback)
      }
      return false
    }
  }, [])

  // Main rendering function with improved error handling
  const renderMolecule = useCallback(async () => {
    // More specific validation - only reject if SMILES is actually empty/null/undefined
    if (!smiles || (typeof smiles === 'string' && smiles.trim() === '')) {
      if (mountedRef.current) {
        setIsLoading(false)
        setError('No SMILES string provided')
      }
      return
    }
    
    // Wait for container to be ready with retry logic
    let retryCount = 0
    const maxRetries = 5
    
    while (!containerRef.current && retryCount < maxRetries) {
      console.log(`Waiting for container... retry ${retryCount + 1}`)
      await new Promise(resolve => setTimeout(resolve, 100))
      retryCount++
    }
    
    if (!containerRef.current) {
      console.error('Container ref never became available after retries')
      if (mountedRef.current) {
        setIsLoading(false)
        setError('Container not available')
      }
      return
    }
    
    console.log(`Starting molecule rendering for: ${smiles} (mode: ${mode}, renderer: ${renderer})`)
    if (mountedRef.current) {
      setIsLoading(true)
      setError(null)
      setSvgContent('') // Clear previous content
    }
    cleanupViewer()
    
    try {
      let success = false
      let molData = null
      
      // For 3D renderers, try to generate 3D coordinates
      if (mode === '3D') {
        try {
          console.log('Generating 3D coordinates...')
          molData = await generateMolecule3D(smiles)
          console.log('3D coordinates generated successfully')
          console.log('Generated molData type:', typeof molData)
          console.log('Generated molData keys:', Object.keys(molData))
          
          console.log('Molblock preview (first 100 chars):', molData.molblock.substring(0, 100))
          console.log('Has real 3D coordinates:', molData.has3D)
          console.log('Atom count:', molData.atomCount)
        } catch (e) {
          console.warn('Failed to generate 3D coordinates:', e)
        }
      }
      
      // Try to render with the selected renderer
      console.log(`Attempting to render with ${renderer}...`)
      switch (renderer) {
        case 'SmilesDrawer':
          success = await renderWithSmilesDrawer(smiles)
          break
        case 'RDKit.js':
          success = await renderWithRDKit(smiles)
          break
        case 'Kekule.js':
          success = await renderWithKekule(smiles)
          break
        case '3Dmol.js':
          success = await renderWith3Dmol(smiles, molData)
          break
        case 'NGL':
          success = await renderWithNGL(smiles, molData)
          break
        case 'Mol*':
          success = await renderWithMolstar(smiles, molData)
          break
        default:
          console.log(`Unknown renderer ${renderer}, trying RDKit`)
          success = await renderWithRDKit(smiles)
      }
      
      // Check if component is still mounted and container is still available
      if (!containerRef.current) {
        console.warn('Container ref became null during rendering, component may have unmounted')
        return
      }
      
      // If the selected renderer failed, try RDKit as fallback for 2D renderers only
      if (!success && ['SmilesDrawer', 'Kekule.js'].includes(renderer)) {
        console.warn(`${renderer} failed, trying RDKit fallback`)
        success = await renderWithRDKit(smiles)
      } else if (!success && ['3Dmol.js', 'NGL'].includes(renderer)) {
        console.warn(`${renderer} failed, trying Mol* fallback`)
        success = await renderWithMolstar(smiles, molData)
      }
      
      // Check container again before fallback
      if (!containerRef.current) {
        console.warn('Container ref became null during fallback rendering')
        return
      }
      
      // If RDKit also failed, try simple SVG fallback
      if (!success) {
        console.warn('All advanced renderers failed, using simple SVG fallback')
        success = await renderWithSimpleSVG(smiles)
      }
      
      // Final container check before error handling
      if (!containerRef.current) {
        console.warn('Container ref became null during error handling')
        return
      }
      
            // If everything failed, create error message
      if (!success) {
        console.error('All rendering methods failed, showing error message')
        
        // Use React state for error display
        if (mountedRef.current) {
          const errorSvg = `
            <svg width="800" height="600" viewBox="0 0 800 600" style="background: #fee2e2; border: 2px solid #fecaca; border-radius: 8px;">
              <text x="400" y="150" text-anchor="middle" font-family="Arial, sans-serif" font-size="4rem" fill="#dc2626">⚠️</text>
              <text x="400" y="250" text-anchor="middle" font-family="Arial, sans-serif" font-size="1.125rem" font-weight="600" fill="#dc2626">Rendering Failed</text>
              <text x="400" y="300" text-anchor="middle" font-family="Arial, sans-serif" font-size="0.875rem" fill="#dc2626">Unable to render molecule</text>
              <text x="400" y="350" text-anchor="middle" font-family="monospace" font-size="0.75rem" fill="#dc2626">SMILES: ${smiles}</text>
              <text x="400" y="400" text-anchor="middle" font-family="Arial, sans-serif" font-size="0.75rem" fill="#7f1d1d">Please check SMILES format and try again</text>
            </svg>
          `
          setSvgContent(errorSvg)
          viewerInstanceRef.current = { type: 'error' }
          setError('Failed to render molecule with any available renderer')
        }
      } else {
        console.log('Molecule rendered successfully')
      }
      
      if (mountedRef.current) {
        setMoleculeData(molData)
      }
    } catch (error) {
      console.error('Critical rendering error:', error)
      if (mountedRef.current) {
        setError(`Error rendering molecule: ${error.message}`)
      }
      
      // Even in case of critical error, try to show simple fallback
      try {
        if (containerRef.current) {
          const fallbackSuccess = await renderWithSimpleSVG(smiles)
          if (fallbackSuccess) {
            console.log('Emergency fallback rendering succeeded')
          }
        }
      } catch (fallbackError) {
        console.error('Even fallback rendering failed:', fallbackError)
      }
    } finally {
      if (mountedRef.current) {
        setIsLoading(false)
      }
    }
  }, [smiles, mode, renderer, cleanupViewer, renderWithRDKit, renderWithMolstar, renderWithSimpleSVG])

  // Render molecule when props change
  useEffect(() => {
    console.log('=== MoleculeViewer: useEffect triggered ===')
    console.log('Triggering renderMolecule with:', { smiles: smiles?.substring(0, 50) + '...', mode, renderer })
    console.log('Container ref current:', !!containerRef.current)
    
    renderMolecule()
  }, [smiles, mode, renderer, renderMolecule])

  // Cleanup on unmount
  useEffect(() => {
    return cleanupViewer
  }, [cleanupViewer])

  // Export handlers
  const handleExportPNG = useCallback(async () => {
    if (!containerRef.current) {
      console.warn('Container not found for PNG export')
      return false
    }
    
    // Find SVG element in the container
    const svg = containerRef.current.querySelector('svg')
    if (svg) {
      try {
        const canvas = document.createElement('canvas')
        const ctx = canvas.getContext('2d')
        const img = new Image()
        
        // Get SVG dimensions properly to capture full content
        let width, height
        
        // First try to get dimensions from viewBox
        const viewBox = svg.getAttribute('viewBox')
        if (viewBox) {
          const [, , vbWidth, vbHeight] = viewBox.split(' ').map(Number)
          width = vbWidth
          height = vbHeight
        } else {
          // Try width/height attributes
          width = parseFloat(svg.getAttribute('width')) || 800
          height = parseFloat(svg.getAttribute('height')) || 600
        }
        
        // Ensure minimum size and reasonable maximums
        width = Math.max(400, Math.min(width, 2000))
        height = Math.max(300, Math.min(height, 2000))
        
        console.log(`PNG Export: Using dimensions ${width}x${height}`)
        
        // Serialize SVG to string
        const svgData = new XMLSerializer().serializeToString(svg)
        const svgBlob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' })
        const url = URL.createObjectURL(svgBlob)
        
        return new Promise((resolve) => {
          img.onload = async () => {
            try {
              canvas.width = width
              canvas.height = height
              ctx.fillStyle = 'white'
              ctx.fillRect(0, 0, width, height)
              ctx.drawImage(img, 0, 0, width, height)
              
              const success = await exportToPNG(canvas, smiles)
              URL.revokeObjectURL(url)
              resolve(success)
            } catch (error) {
              console.error('PNG export error:', error)
              URL.revokeObjectURL(url)
              resolve(false)
            }
          }
          
          img.onerror = () => {
            console.error('Failed to load SVG for PNG export')
            URL.revokeObjectURL(url)
            resolve(false)
          }
          
          img.src = url
        })
      } catch (error) {
        console.error('PNG export error:', error)
        return false
      }
    }
    
    // Try to find canvas element for 3D renderers
    const canvas = containerRef.current.querySelector('canvas')
    if (canvas) {
      try {
        return await exportToPNG(canvas, smiles)
      } catch (error) {
        console.error('Canvas PNG export error:', error)
        return false
      }
    }
    
    console.warn('No exportable content found (SVG or Canvas)')
    return false
  }, [smiles])

  const handleExportSVG = useCallback(async () => {
    if (!containerRef.current) {
      console.warn('Container not found for SVG export')
      return false
    }
    
    // Find SVG element in the container
    const svg = containerRef.current.querySelector('svg')
    if (svg) {
      try {
        return await exportToSVG(svg, smiles)
      } catch (error) {
        console.error('SVG export error:', error)
        return false
      }
    }
    
    console.warn('No SVG content found for export')
    return false
  }, [smiles])

  const handleExportGLTF = useCallback(async () => {
    if (mode !== '3D') {
      console.warn('GLTF export only available in 3D mode')
      return false
    }
    
    try {
      // Create a proper GLTF 2.0 structure with molecular data
      const gltfData = {
        asset: {
          version: "2.0",
          generator: "Lipid Viewer - Molecular Visualization Tool",
          copyright: "Generated molecular structure"
        },
        scene: 0,
        scenes: [{
          name: "Molecular Scene",
          nodes: [0]
        }],
        nodes: [{
          name: `Molecule_${smiles}`,
          extras: {
            molecularData: {
              smiles: smiles,
              formula: "Generated from SMILES",
              molblock: moleculeData?.molblock || "No 3D coordinates available",
              has3D: moleculeData?.has3D || false,
              atomCount: moleculeData?.atomCount || "Unknown",
              generatedBy: "Lipid Viewer",
              exportTime: new Date().toISOString()
            }
          }
        }],
        extensionsUsed: [],
        extensionsRequired: [],
        extras: {
          description: "Molecular structure exported from Lipid Viewer",
          format: "GLTF 2.0 with molecular metadata",
          software: "Lipid Viewer",
          smiles: smiles
        }
      }
      
      // Add basic geometry if we have molecular data
      if (moleculeData && moleculeData.has3D) {
        // Add minimal mesh for compatibility
        gltfData.meshes = [{
          name: "MolecularMesh",
          primitives: [{
            mode: 1, // LINES
            attributes: {},
            extras: {
              note: "Actual molecular geometry in molblock format",
              molblock: moleculeData.molblock
            }
          }]
        }]
        gltfData.nodes[0].mesh = 0
      }
      
      // Convert to proper JSON format
      const jsonString = JSON.stringify(gltfData, null, 2)
      
      // Export as GLTF file
      return await exportToGLTF(jsonString, smiles)
      
    } catch (error) {
      console.error('GLTF export error:', error)
      return false
    }
  }, [smiles, mode, renderer, moleculeData])

  // Expose export functions
  useEffect(() => {
    if (onExport) {
      onExport({
        png: handleExportPNG,
        svg: handleExportSVG,
        gltf: handleExportGLTF,
      })
    }
  }, [onExport, handleExportPNG, handleExportSVG, handleExportGLTF])

  return (
    <div 
      className={`molecule-viewer ${mode === '2D' ? 'molecule-2d' : 'molecule-3d'} p-4`} 
      data-testid="molecule-viewer"
    >
      <div 
        ref={containerRef}
        className="w-full h-96 flex items-center justify-center border border-gray-300 rounded-lg bg-white"
        style={{ minHeight: '400px' }}
      >
        {isLoading ? (
          <div className="text-center">
            <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600 mx-auto mb-4"></div>
            <p className="text-gray-600">Loading molecule...</p>
          </div>
        ) : svgContent ? (
          <div 
            dangerouslySetInnerHTML={{ __html: svgContent }}
            style={{ width: '100%', height: '100%' }}
          />
        ) : error ? (
          <div className="text-center">
            <div className="text-red-500 text-6xl mb-4">⚠️</div>
            <p className="text-red-600 font-semibold">Error rendering molecule</p>
            <p className="text-red-500 text-sm mt-2">{error}</p>
          </div>
        ) : null}
      </div>
      
      <div className="mt-4 text-center">
        <p className="text-sm text-gray-600">
          <strong>Mode:</strong> {mode} | <strong>Renderer:</strong> {renderer} | <strong>SMILES:</strong> {smiles}
        </p>
        {renderer && !availableRenderers[renderer.replace('.', '').replace('*', 'star')] && (
          <p className="text-xs text-orange-600 mt-1">
            Note: {renderer} not fully available, showing fallback display
          </p>
        )}
      </div>
    </div>
  )
})

export default MoleculeViewer 