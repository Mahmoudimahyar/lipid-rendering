import React from 'react'

const ViewControls = ({ 
  mode, 
  onZoomIn, 
  onZoomOut, 
  onResetView, 
  onExportPNG, 
  onExportSVG, 
  onExportGLTF,
  isLoaded = false 
}) => {
  if (!isLoaded) {
    return null
  }

  return (
    <div className="w-full max-w-4xl mx-auto p-6 bg-white rounded-lg shadow-lg">
      <h2 className="text-2xl font-bold text-gray-800 mb-4">View Controls</h2>
      
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
        {/* Navigation Controls */}
        <div className="space-y-2">
          <h3 className="font-semibold text-gray-700">Navigation</h3>
          <button
            onClick={onZoomIn}
            className="control-button w-full text-sm flex items-center justify-center space-x-2"
            title="Zoom In (or use + key)"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M12 6v6m0 0v6m0-6h6m-6 0H6"/>
            </svg>
            <span>Zoom In</span>
          </button>
          
          <button
            onClick={onZoomOut}
            className="control-button w-full text-sm flex items-center justify-center space-x-2"
            title="Zoom Out (or use - key)"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M18 12H6"/>
            </svg>
            <span>Zoom Out</span>
          </button>
          
          <button
            onClick={onResetView}
            className="control-button w-full text-sm flex items-center justify-center space-x-2"
            title="Reset View (or use Space key)"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15"/>
            </svg>
            <span>Reset View</span>
          </button>
        </div>
        
        {/* Export Controls */}
        <div className="space-y-2">
          <h3 className="font-semibold text-gray-700">Export</h3>
          <button
            onClick={onExportPNG}
            className="control-button w-full text-sm flex items-center justify-center space-x-2"
            title="Export as PNG image"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M4 16l4.586-4.586a2 2 0 012.828 0L16 16m-2-2l1.586-1.586a2 2 0 012.828 0L20 14m-6-6h.01M6 20h12a2 2 0 002-2V6a2 2 0 00-2-2H6a2 2 0 00-2 2v12a2 2 0 002 2z"/>
            </svg>
            <span>Export PNG</span>
          </button>
          
          <button
            onClick={onExportSVG}
            className="control-button w-full text-sm flex items-center justify-center space-x-2"
            title="Export as SVG vector"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z"/>
            </svg>
            <span>Export SVG</span>
          </button>
          
          {mode === '3D' && (
            <button
              onClick={onExportGLTF}
              className="control-button w-full text-sm flex items-center justify-center space-x-2"
              title="Export as GLTF 3D model"
            >
              <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M19 11H5m14 0a2 2 0 012 2v6a2 2 0 01-2 2H5a2 2 0 01-2-2v-6a2 2 0 012-2m14 0V9a2 2 0 00-2-2M5 11V9a2 2 0 012-2m0 0V5a2 2 0 012-2h6a2 2 0 012 2v2M7 7h10"/>
              </svg>
              <span>Export GLTF</span>
            </button>
          )}
        </div>
        
        {/* Additional Information */}
        <div className="md:col-span-2 lg:col-span-2">
          <h3 className="font-semibold text-gray-700 mb-2">Keyboard Shortcuts</h3>
          <div className="text-sm text-gray-600 space-y-1">
            <div className="flex justify-between">
              <span>Zoom In:</span>
              <span className="bg-gray-100 px-2 py-1 rounded">+</span>
            </div>
            <div className="flex justify-between">
              <span>Zoom Out:</span>
              <span className="bg-gray-100 px-2 py-1 rounded">-</span>
            </div>
            <div className="flex justify-between">
              <span>Reset View:</span>
              <span className="bg-gray-100 px-2 py-1 rounded">Space</span>
            </div>
            {mode === '3D' && (
              <>
                <div className="flex justify-between">
                  <span>Rotate:</span>
                  <span className="bg-gray-100 px-2 py-1 rounded">Mouse Drag</span>
                </div>
                <div className="flex justify-between">
                  <span>Pan:</span>
                  <span className="bg-gray-100 px-2 py-1 rounded">Shift + Drag</span>
                </div>
              </>
            )}
          </div>
        </div>
      </div>
      
      {/* Current Mode Info */}
      <div className="mt-4 p-3 bg-blue-50 rounded-lg">
        <p className="text-sm text-blue-800">
          <strong>Current Mode:</strong> {mode === '2D' ? '2D Structure View' : '3D Molecular View'}
          {mode === '3D' && ' - Use mouse to rotate, scroll to zoom'}
        </p>
      </div>
    </div>
  )
}

export default ViewControls 