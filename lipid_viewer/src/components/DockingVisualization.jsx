import React, { useState, useEffect, useRef } from 'react'
import DockingMoleculeViewer from './DockingMoleculeViewer'
import Enhanced3DViewer from './Enhanced3DViewer'

/**
 * Enhanced docking visualization component with pose overlay and scoring
 * This component wraps the existing MoleculeViewer with docking-specific features
 */
function DockingVisualization({ 
  ligandSmiles,
  receptorPdbId,
  dockingResults,
  mode = '3D',
  renderer = '3Dmol.js',
  className = ''
}) {
  const [selectedPose, setSelectedPose] = useState(0)
  const [visiblePoses, setVisiblePoses] = useState(new Set([0]))
  const [colorScheme, setColorScheme] = useState('affinity') // 'affinity', 'rmsd', 'rainbow'
  const [showLegend, setShowLegend] = useState(true)
  const [cameraFocused, setCameraFocused] = useState(false)
  const [useEnhancedViewer, setUseEnhancedViewer] = useState(true)
  const viewerRef = useRef(null)

  // Color schemes for different visualization modes
  const colorSchemes = {
    affinity: {
      name: 'Binding Affinity',
      getColor: (pose, index, poses) => {
        const affinity = Math.abs(pose.affinity)
        const maxAffinity = Math.max(...poses.map(p => Math.abs(p.affinity)))
        const minAffinity = Math.min(...poses.map(p => Math.abs(p.affinity)))
        const ratio = (affinity - minAffinity) / (maxAffinity - minAffinity)
        
        // Green (best) to Red (worst) gradient
        const r = Math.floor(255 * ratio)
        const g = Math.floor(255 * (1 - ratio))
        return `rgb(${r}, ${g}, 0)`
      },
      description: 'Green = Better affinity, Red = Worse affinity'
    },
    rmsd: {
      name: 'RMSD',
      getColor: (pose, index, poses) => {
        const rmsd = pose.rmsd_lb
        const maxRmsd = Math.max(...poses.map(p => p.rmsd_lb))
        const minRmsd = Math.min(...poses.map(p => p.rmsd_lb))
        const ratio = (rmsd - minRmsd) / (maxRmsd - minRmsd)
        
        // Blue (low RMSD) to Yellow (high RMSD) gradient
        const r = Math.floor(255 * ratio)
        const g = Math.floor(255 * ratio)
        const b = Math.floor(255 * (1 - ratio))
        return `rgb(${r}, ${g}, ${b})`
      },
      description: 'Blue = Low RMSD, Yellow = High RMSD'
    },
    rainbow: {
      name: 'Pose Order',
      getColor: (pose, index, poses) => {
        const hue = (index / poses.length) * 300 // 0-300 degrees (avoid red overlap)
        return `hsl(${hue}, 70%, 50%)`
      },
      description: 'Different colors for each pose'
    }
  }

  // Get sorted poses for display
  const poses = dockingResults?.poses || []
  const sortedPoses = [...poses].sort((a, b) => a.affinity - b.affinity) // Best first

  // Handle pose selection
  const handlePoseSelect = (poseIndex) => {
    setSelectedPose(poseIndex)
    setVisiblePoses(new Set([poseIndex]))
  }

  // Handle pose visibility toggle
  const togglePoseVisibility = (poseIndex) => {
    const newVisible = new Set(visiblePoses)
    if (newVisible.has(poseIndex)) {
      newVisible.delete(poseIndex)
    } else {
      newVisible.add(poseIndex)
    }
    setVisiblePoses(newVisible)
  }

  // Show/hide all poses
  const showAllPoses = () => {
    setVisiblePoses(new Set(poses.map((_, index) => index)))
  }

  const hideAllPoses = () => {
    setVisiblePoses(new Set())
  }

  // Focus camera on binding site
  const focusOnBindingSite = (poseData) => {
    console.log('🎯 DockingVisualization: Focus on binding site requested', poseData)
    setCameraFocused(true)
    
    // Additional logic can be added here if needed
    if (poses.length > 0) {
      const selectedPoseData = poseData || poses[selectedPose]
      if (selectedPoseData) {
        console.log('📍 Focusing on pose:', selectedPoseData)
      }
    }
  }

  // Reset camera view
  const resetCamera = () => {
    console.log('🔄 DockingVisualization: Reset camera view requested')
    setCameraFocused(false)
    
    // Additional logic can be added here if needed
  }

  return (
    <div className={`docking-visualization ${className}`}>
      {/* Main Viewer Area */}
      <div className="relative">
        {/* Viewer Toggle Controls - moved outside viewer */}
        <div className="flex justify-between items-center mb-2">
          <h3 className="text-lg font-semibold text-gray-800">3D Molecular Visualization</h3>
          <div className="flex items-center space-x-2">
            <span className="text-gray-600 text-sm">Viewer Mode:</span>
            <button
              onClick={() => setUseEnhancedViewer(!useEnhancedViewer)}
              className={`px-3 py-1 text-sm rounded transition-colors ${
                useEnhancedViewer 
                  ? 'bg-blue-600 text-white' 
                  : 'bg-gray-200 text-gray-700 hover:bg-gray-300'
              }`}
            >
              {useEnhancedViewer ? 'Enhanced' : 'Basic'}
            </button>
          </div>
        </div>

        {/* 3D Viewer Container */}
        <div 
          ref={viewerRef}
          className="w-full bg-black rounded-lg border-2 border-gray-300 relative overflow-hidden"
        >

          {/* Conditional Viewer Rendering */}
          {useEnhancedViewer ? (
            <Enhanced3DViewer
              ligandSmiles={ligandSmiles}
              receptorPdbId={receptorPdbId}
              poses={poses}
              selectedPose={selectedPose}
              visiblePoses={visiblePoses}
              onFocusBindingSite={focusOnBindingSite}
              onResetView={resetCamera}
              className="w-full"
              key={`enhanced-viewer-${ligandSmiles}-${receptorPdbId}`}
            />
          ) : (
            <div className="h-96">
              <DockingMoleculeViewer
                ligandSmiles={ligandSmiles}
                receptorPdbId={receptorPdbId}
                className="w-full h-full"
                key={`basic-viewer-${ligandSmiles}-${receptorPdbId}`}
              />
            </div>
          )}
          
          {/* Info Overlay - only show for basic viewer */}
          {!useEnhancedViewer && (
            <div className="absolute bottom-4 left-4 bg-black bg-opacity-75 text-white rounded px-3 py-2 text-sm">
              <div>Ligand: {ligandSmiles || 'N/A'}</div>
              <div>Receptor: {receptorPdbId || 'N/A'}</div>
              <div>Visible poses: {visiblePoses.size} of {poses.length}</div>
            </div>
          )}

          {/* Overlay Controls - only show for basic viewer */}
          {!useEnhancedViewer && (
            <div className="absolute top-4 right-4 space-y-2">
              <button
                onClick={focusOnBindingSite}
                className={`px-3 py-1 text-xs rounded ${cameraFocused ? 'bg-blue-600 text-white' : 'bg-white text-gray-700'} shadow`}
              >
                Focus Site
              </button>
              <button
                onClick={resetCamera}
                className="block px-3 py-1 text-xs bg-white text-gray-700 rounded shadow"
              >
                Reset View
              </button>
            </div>
          )}

          {/* Color Scheme Legend - only show for basic viewer */}
          {!useEnhancedViewer && showLegend && (
            <div className="absolute top-4 left-4 bg-white rounded-lg shadow-lg p-3 max-w-xs">
              <div className="text-sm font-semibold mb-2">
                Color Scheme: {colorSchemes[colorScheme].name}
              </div>
              <div className="text-xs text-gray-600 mb-2">
                {colorSchemes[colorScheme].description}
              </div>
              <div className="flex items-center space-x-2 text-xs">
                <button
                  onClick={() => setShowLegend(false)}
                  className="text-gray-400 hover:text-gray-600"
                >
                  ✕
                </button>
              </div>
            </div>
          )}
        </div>

        {/* Control Panel */}
        <div className="mt-4 grid grid-cols-1 lg:grid-cols-2 gap-4">
          {/* Pose Selection Panel */}
          <div className="bg-white rounded-lg shadow p-4">
            <div className="flex items-center justify-between mb-3">
              <h3 className="text-lg font-semibold">Poses</h3>
              <div className="space-x-2">
                <button
                  onClick={showAllPoses}
                  className="px-2 py-1 text-xs bg-blue-100 text-blue-700 rounded hover:bg-blue-200"
                >
                  Show All
                </button>
                <button
                  onClick={hideAllPoses}
                  className="px-2 py-1 text-xs bg-gray-100 text-gray-700 rounded hover:bg-gray-200"
                >
                  Hide All
                </button>
              </div>
            </div>
            
            <div className="max-h-48 overflow-y-auto space-y-2">
              {sortedPoses.map((pose, index) => {
                const originalIndex = poses.findIndex(p => p.mode === pose.mode)
                const isSelected = selectedPose === originalIndex
                const isVisible = visiblePoses.has(originalIndex)
                const poseColor = colorSchemes[colorScheme].getColor(pose, originalIndex, poses)
                
                return (
                  <div
                    key={pose.mode}
                    className={`flex items-center p-2 rounded border ${
                      isSelected ? 'border-blue-500 bg-blue-50' : 'border-gray-200'
                    }`}
                  >
                    {/* Color indicator */}
                    <div
                      className="w-4 h-4 rounded-full border mr-3"
                      style={{ backgroundColor: poseColor }}
                    />
                    
                    {/* Pose info */}
                    <div 
                      className="flex-1 cursor-pointer"
                      onClick={() => handlePoseSelect(originalIndex)}
                    >
                      <div className="flex justify-between items-center">
                        <span className="font-medium text-sm">Mode {pose.mode}</span>
                        <span className="text-sm font-bold text-green-600">
                          {pose.affinity} kcal/mol
                        </span>
                      </div>
                      <div className="text-xs text-gray-500">
                        RMSD: {pose.rmsd_lb} - {pose.rmsd_ub} Å
                      </div>
                    </div>
                    
                    {/* Visibility toggle */}
                    <button
                      onClick={() => togglePoseVisibility(originalIndex)}
                      className={`ml-2 px-2 py-1 text-xs rounded ${
                        isVisible 
                          ? 'bg-green-100 text-green-700 hover:bg-green-200' 
                          : 'bg-gray-100 text-gray-500 hover:bg-gray-200'
                      }`}
                    >
                      {isVisible ? '👁️' : '👁️‍🗨️'}
                    </button>
                  </div>
                )
              })}
            </div>
          </div>

          {/* Visualization Controls Panel */}
          <div className="bg-white rounded-lg shadow p-4">
            <h3 className="text-lg font-semibold mb-3">Visualization</h3>
            
            {/* Color Scheme Selection */}
            <div className="mb-4">
              <label className="block text-sm font-medium text-gray-700 mb-2">
                Color Scheme
              </label>
              <select
                value={colorScheme}
                onChange={(e) => setColorScheme(e.target.value)}
                className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm"
              >
                {Object.entries(colorSchemes).map(([key, scheme]) => (
                  <option key={key} value={key}>
                    {scheme.name}
                  </option>
                ))}
              </select>
            </div>

            {/* Display Options */}
            <div className="space-y-3">
              <label className="flex items-center">
                <input
                  type="checkbox"
                  checked={showLegend}
                  onChange={(e) => setShowLegend(e.target.checked)}
                  className="mr-2"
                />
                <span className="text-sm">Show Legend</span>
              </label>
            </div>

            {/* Camera Controls - only for basic viewer, Enhanced3DViewer has its own */}
            {!useEnhancedViewer && (
              <div className="mt-4 pt-4 border-t">
                <h4 className="text-sm font-medium mb-2">Camera</h4>
                <div className="space-y-2">
                  <button
                    onClick={focusOnBindingSite}
                    className="w-full px-3 py-2 text-sm bg-blue-600 text-white rounded hover:bg-blue-700"
                  >
                    Focus on Binding Site
                  </button>
                  <button
                    onClick={resetCamera}
                    className="w-full px-3 py-2 text-sm bg-gray-600 text-white rounded hover:bg-gray-700"
                  >
                    Reset View
                  </button>
                </div>
              </div>
            )}

            {/* Docking Summary */}
            {dockingResults?.summary && (
              <div className="mt-4 pt-4 border-t">
                <h4 className="text-sm font-medium mb-2">Summary</h4>
                <div className="text-xs text-gray-600 space-y-1">
                  <div>Best: {dockingResults.summary.best_affinity} kcal/mol</div>
                  <div>Average: {dockingResults.summary.mean_affinity} kcal/mol</div>
                  <div>Time: {dockingResults.summary.calculation_time}s</div>
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  )
}

export default DockingVisualization