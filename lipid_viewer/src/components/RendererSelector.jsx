import React from 'react'

const RENDERERS_2D = [
  { id: 'SmilesDrawer', name: 'SmilesDrawer', description: 'Fast client-side 2D rendering' },
  { id: 'RDKit.js', name: 'RDKit.js', description: 'High-quality chemical structures' },
  { id: 'Kekule.js', name: 'Kekule.js', description: 'Comprehensive cheminformatics' },
]

const RENDERERS_3D = [
  { id: '3Dmol.js', name: '3Dmol.js', description: 'WebGL molecular visualization' },
  { id: 'Mol*', name: 'Mol*', description: 'High-performance 3D viewer' },
  { id: 'NGL', name: 'NGL', description: 'Interactive protein structures' },
]

const RendererSelector = ({ mode, renderer, onModeChange, onRendererChange }) => {
  const currentRenderers = mode === '2D' ? RENDERERS_2D : RENDERERS_3D

  return (
    <div className="w-full max-w-4xl mx-auto p-6 bg-white rounded-lg shadow-lg">
      <h2 className="text-2xl font-bold text-gray-800 mb-4">Visualization Mode</h2>
      
      {/* Mode Toggle */}
      <div className="flex space-x-2 mb-6">
        <button
          onClick={() => onModeChange('2D')}
          className={`px-6 py-3 rounded-lg font-semibold transition-colors ${
            mode === '2D'
              ? 'bg-blue-600 text-white'
              : 'bg-gray-200 text-gray-700 hover:bg-gray-300'
          }`}
        >
          2D View
        </button>
        <button
          onClick={() => onModeChange('3D')}
          className={`px-6 py-3 rounded-lg font-semibold transition-colors ${
            mode === '3D'
              ? 'bg-blue-600 text-white'
              : 'bg-gray-200 text-gray-700 hover:bg-gray-300'
          }`}
        >
          3D View
        </button>
      </div>
      
      {/* Renderer Selection */}
      <div>
        <h3 className="text-lg font-semibold text-gray-700 mb-3">
          {mode === '2D' ? '2D Renderer' : '3D Renderer'}
        </h3>
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
          {currentRenderers.map((rendererOption) => (
            <button
              key={rendererOption.id}
              onClick={() => onRendererChange(rendererOption.id)}
              className={`p-4 rounded-lg border-2 text-left transition-colors ${
                renderer === rendererOption.id
                  ? 'border-blue-500 bg-blue-100'
                  : 'border-gray-200 hover:border-gray-300 hover:bg-gray-50'
              }`}
            >
              <div className="font-semibold text-gray-800 mb-1">
                {rendererOption.name}
              </div>
              <div className="text-sm text-gray-600">
                {rendererOption.description}
              </div>
            </button>
          ))}
        </div>
      </div>
      
      {/* Mode-specific information */}
      <div className="mt-6 p-4 bg-gray-50 rounded-lg">
        <h4 className="font-semibold text-gray-700 mb-2">
          {mode === '2D' ? '2D Visualization Features' : '3D Visualization Features'}
        </h4>
        <ul className="text-sm text-gray-600 space-y-1">
          {mode === '2D' ? (
            <>
              <li>• Fast rendering and interaction</li>
              <li>• Export to PNG and SVG formats</li>
              <li>• Chemical structure validation</li>
              <li>• Customizable styling options</li>
            </>
          ) : (
            <>
              <li>• Interactive 3D molecule rotation</li>
              <li>• Multiple rendering styles (stick, ball & stick, surface)</li>
              <li>• Export to GLB format for 3D printing</li>
              <li>• Real-time conformer generation</li>
            </>
          )}
        </ul>
      </div>
    </div>
  )
}

export default RendererSelector 