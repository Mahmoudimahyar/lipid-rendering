import React, { useState, useEffect, useCallback, useRef } from 'react'
import { Toaster, toast } from 'react-hot-toast'
import { Link } from 'react-router-dom'
import SMILESInput from './components/SMILESInput'
import RendererSelector from './components/RendererSelector'
import MoleculeViewer from './components/MoleculeViewer'
import ViewControls from './components/ViewControls'

function App() {
  // State management
  const [currentSMILES, setCurrentSMILES] = useState('')
  const [mode, setMode] = useState('2D')
  const [renderer, setRenderer] = useState('SmilesDrawer')
  const [isValidSMILES, setIsValidSMILES] = useState(false)
  const [isMoleculeLoaded, setIsMoleculeLoaded] = useState(false)
  const [exportFunctions, setExportFunctions] = useState(null)

  console.log('=== App.jsx: Component rendered ===')
  console.log('Current state:', { currentSMILES, mode, renderer, isValidSMILES, isMoleculeLoaded })

  // Log whenever currentSMILES changes
  useEffect(() => {
    console.log('=== App.jsx: currentSMILES changed ===')
    console.log('New currentSMILES value:', currentSMILES)
    console.log('Will render MoleculeViewer:', !!currentSMILES)
  }, [currentSMILES])
  
  // Refs for managing viewer interactions
  const viewerRef = useRef(null)
  const zoomLevel = useRef(1)

  // Handle mode changes and update default renderer
  const handleModeChange = useCallback((newMode) => {
    console.log('=== App.jsx: handleModeChange called ===')
    console.log('Mode changed from', mode, 'to', newMode)
    setMode(newMode)
    
    // Set default renderer for the mode
    if (newMode === '2D') {
      setRenderer('SmilesDrawer')
    } else {
      setRenderer('3Dmol.js')
    }
  }, [mode])

  // Handle renderer changes
  const handleRendererChange = useCallback((newRenderer) => {
    console.log('=== App.jsx: handleRendererChange called ===')
    console.log('Renderer changed from', renderer, 'to', newRenderer)
    setRenderer(newRenderer)
  }, [renderer])

  // Handle SMILES validation
  const handleSMILESValidation = useCallback((smiles, isValid) => {
    setIsValidSMILES(isValid)
    if (!isValid && smiles.trim()) {
      toast.error('Invalid SMILES string', {
        duration: 3000,
        position: 'top-center',
      })
    }
  }, [])

  // Handle SMILES submission
  const handleSMILESSubmit = useCallback((smiles) => {
    console.log('=== App.jsx: handleSMILESSubmit called ===')
    console.log('Received SMILES:', smiles)
    console.log('Current state before update:', { currentSMILES, mode, renderer })
    
    if (!smiles.trim()) {
      toast.error('Please enter a SMILES string', {
        duration: 4000,
        position: 'top-center',
      })
      return
    }
    
    // The MoleculeViewer will handle rendering validation internally
    setCurrentSMILES(smiles)
    setIsMoleculeLoaded(true)
    console.log('SMILES state updated to:', smiles)
    console.log('isMoleculeLoaded set to true')
    
    // Show success message
    toast.success(`Loading molecule: ${smiles}`)
    
    // Show warning if validation indicated invalid SMILES
    if (!isValidSMILES) {
      toast.error('Warning: SMILES may be invalid, but attempting to render', {
        duration: 5000,
        position: 'top-center',
      })
    }
  }, [isValidSMILES, currentSMILES, mode, renderer])

  // Handle export functions from MoleculeViewer
  const handleExportFunctions = useCallback((functions) => {
    setExportFunctions(functions)
  }, [])

  // View control handlers
  const handleZoomIn = useCallback(() => {
    zoomLevel.current = Math.min(zoomLevel.current * 1.2, 5)
    toast.success(`Zoom: ${Math.round(zoomLevel.current * 100)}%`, {
      duration: 1000,
      position: 'bottom-right',
    })
  }, [])

  const handleZoomOut = useCallback(() => {
    zoomLevel.current = Math.max(zoomLevel.current / 1.2, 0.2)
    toast.success(`Zoom: ${Math.round(zoomLevel.current * 100)}%`, {
      duration: 1000,
      position: 'bottom-right',
    })
  }, [])

  const handleResetView = useCallback(() => {
    zoomLevel.current = 1
    toast.success('View reset', {
      duration: 1000,
      position: 'bottom-right',
    })
  }, [])

  // Export handlers
  const handleExportPNG = useCallback(async () => {
    if (exportFunctions?.png) {
      toast.promise(
        exportFunctions.png(),
        {
          loading: 'Exporting PNG...',
          success: 'PNG exported successfully!',
          error: 'Failed to export PNG',
        }
      )
    }
  }, [exportFunctions])

  const handleExportSVG = useCallback(async () => {
    if (exportFunctions?.svg) {
      toast.promise(
        exportFunctions.svg(),
        {
          loading: 'Exporting SVG...',
          success: 'SVG exported successfully!',
          error: 'Failed to export SVG',
        }
      )
    }
  }, [exportFunctions])

  const handleExportGLTF = useCallback(async () => {
    if (exportFunctions?.gltf) {
      toast.promise(
        exportFunctions.gltf(),
        {
          loading: 'Exporting GLTF...',
          success: 'GLTF exported successfully!',
          error: 'Failed to export GLTF',
        }
      )
    }
  }, [exportFunctions])

  // Keyboard shortcuts
  useEffect(() => {
    const handleKeyPress = (event) => {
      // Don't trigger shortcuts when typing in input fields
      if (event.target.tagName === 'INPUT' || event.target.tagName === 'TEXTAREA') {
        return
      }

      switch (event.key) {
        case '+':
        case '=':
          event.preventDefault()
          handleZoomIn()
          break
        case '-':
        case '_':
          event.preventDefault()
          handleZoomOut()
          break
        case ' ':
          event.preventDefault()
          handleResetView()
          break
        case '2':
          if (event.ctrlKey || event.metaKey) {
            event.preventDefault()
            handleModeChange('2D')
          }
          break
        case '3':
          if (event.ctrlKey || event.metaKey) {
            event.preventDefault()
            handleModeChange('3D')
          }
          break
        default:
          break
      }
    }

    window.addEventListener('keydown', handleKeyPress)
    return () => window.removeEventListener('keydown', handleKeyPress)
  }, [handleZoomIn, handleZoomOut, handleResetView, handleModeChange])

  // Update molecule loaded state when SMILES changes
  useEffect(() => {
    if (!currentSMILES) {
      setIsMoleculeLoaded(false)
    }
  }, [currentSMILES])

  return (
    <div className="min-h-screen bg-gray-50">
      <Toaster />
      
      {/* Header */}
      <header className="bg-white shadow-sm border-b">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-6">
          <div className="flex items-center justify-between">
            <div>
              <h1 className="text-3xl font-bold text-gray-900">Lipid Viewer</h1>
              <p className="mt-1 text-sm text-gray-600">
                Interactive 2D/3D molecular visualization tool
              </p>
            </div>
            <div className="flex items-center space-x-4">
              <div className="text-sm text-gray-500">
                Mode: <span className="font-semibold">{mode}</span>
              </div>
              <div className="text-sm text-gray-500">
                Renderer: <span className="font-semibold">{renderer}</span>
              </div>
              <Link to="/dock" className="ml-4 inline-flex items-center px-3 py-1.5 rounded bg-blue-600 text-white hover:bg-blue-700 text-sm">Docking</Link>
            </div>
          </div>
        </div>
      </header>

      {/* Main Content */}
      <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="space-y-8">
          {/* SMILES Input Section */}
          <SMILESInput
            onSubmit={handleSMILESSubmit}
            onValidation={handleSMILESValidation}
            isValid={isValidSMILES}
          />

          {/* Renderer Selection */}
          <RendererSelector
            mode={mode}
            renderer={renderer}
            onModeChange={handleModeChange}
            onRendererChange={handleRendererChange}
          />

          {/* Molecule Viewer */}
          {currentSMILES && (
            <div className="bg-white rounded-lg shadow-lg p-6">
              <h2 className="text-2xl font-bold text-gray-800 mb-4">Molecule Visualization</h2>
              <MoleculeViewer
                smiles={currentSMILES}
                mode={mode}
                renderer={renderer}
                onExport={handleExportFunctions}
              />
            </div>
          )}

          {/* View Controls */}
          {isMoleculeLoaded && (
            <ViewControls
              mode={mode}
              onZoomIn={handleZoomIn}
              onZoomOut={handleZoomOut}
              onResetView={handleResetView}
              onExportPNG={handleExportPNG}
              onExportSVG={handleExportSVG}
              onExportGLTF={handleExportGLTF}
              isLoaded={isMoleculeLoaded}
            />
          )}

          {/* Instructions */}
          {!currentSMILES && (
            <div className="bg-blue-50 border border-blue-200 rounded-lg p-6">
              <h3 className="text-lg font-semibold text-blue-800 mb-3">
                Welcome to Lipid Viewer
              </h3>
              <div className="text-blue-700 space-y-2">
                <p>• Enter a SMILES string above to visualize molecular structures</p>
                <p>• Choose between 2D and 3D visualization modes</p>
                <p>• Select different rendering engines for optimal display</p>
                <p>• Export your visualizations as PNG, SVG, or GLTF files</p>
                <p>• Use keyboard shortcuts: +/- for zoom, Space to reset view</p>
              </div>
            </div>
          )}

          {/* Framework Preparation Notice */}
          <div className="bg-yellow-50 border border-yellow-200 rounded-lg p-4">
            <h4 className="text-sm font-semibold text-yellow-800 mb-2">
              🚧 Framework Preparation
            </h4>
            <p className="text-yellow-700 text-sm">
              Electrostatic surface visualization feature is prepared for future implementation.
              The current framework supports multiple rendering engines and can be extended
              to include surface rendering capabilities.
            </p>
          </div>
        </div>
      </main>

      {/* Footer */}
      <footer className="bg-white border-t mt-12">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-6">
          <div className="text-center text-sm text-gray-600">
            <p>
              Powered by{' '}
              <span className="font-semibold">SmilesDrawer</span>,{' '}
              <span className="font-semibold">RDKit.js</span>,{' '}
              <span className="font-semibold">Kekule.js</span>,{' '}
              <span className="font-semibold">3Dmol.js</span>,{' '}
              <span className="font-semibold">Mol*</span>, and{' '}
              <span className="font-semibold">NGL</span>
            </p>
            <p className="mt-1">
              Built with React 18, Vite, and Tailwind CSS
            </p>
          </div>
        </div>
      </footer>
    </div>
  )
}

export default App 