import React, { useState } from 'react'

/**
 * MetadataPanel - Displays comprehensive docking job metadata
 * Shows engine information, parameters, performance metrics, and accessibility features
 */
function MetadataPanel({ 
  jobData,
  className = '',
  showExportButton = true,
  onExport = null
}) {
  const [isExpanded, setIsExpanded] = useState(false)
  const [exportLoading, setExportLoading] = useState(false)

  if (!jobData) {
    return (
      <div className={`bg-gray-50 rounded-lg p-4 ${className}`}>
        <p className="text-sm text-gray-500">No metadata available</p>
      </div>
    )
  }

  const {
    engine_metadata = {},
    docking_parameters = {},
    performance_metrics = {},
    duration,
    status,
    started_at,
    completed_at,
  } = jobData

  // Determine engine badge style and accessibility
  const getEngineBadgeStyle = () => {
    if (engine_metadata.is_mock) {
      return {
        className: 'bg-red-100 text-red-800 border border-red-200',
        text: 'MOCK (Demo Only)',
        ariaLabel: 'Mock engine - demonstration mode only, not real docking',
        icon: '⚠️'
      }
    } else if (engine_metadata.engine === 'vina') {
      return {
        className: 'bg-green-100 text-green-800 border border-green-200',
        text: `AutoDock Vina ${engine_metadata.version || 'v1.2.5'}`,
        ariaLabel: `Real AutoDock Vina engine version ${engine_metadata.version || '1.2.5'}`,
        icon: '✅'
      }
    } else {
      return {
        className: 'bg-gray-100 text-gray-800 border border-gray-200',
        text: engine_metadata.engine || 'Unknown',
        ariaLabel: `Unknown engine: ${engine_metadata.engine || 'not specified'}`,
        icon: '❓'
      }
    }
  }

  const engineBadge = getEngineBadgeStyle()

  const handleExport = async () => {
    if (!onExport) return
    
    setExportLoading(true)
    try {
      await onExport()
    } catch (error) {
      console.error('Export failed:', error)
    } finally {
      setExportLoading(false)
    }
  }

  const formatDuration = (seconds) => {
    if (!seconds) return 'N/A'
    if (seconds < 60) return `${seconds.toFixed(1)}s`
    const minutes = Math.floor(seconds / 60)
    const remainingSeconds = seconds % 60
    return `${minutes}m ${remainingSeconds.toFixed(1)}s`
  }

  const formatDate = (isoString) => {
    if (!isoString) return 'N/A'
    try {
      return new Date(isoString).toLocaleString()
    } catch {
      return 'Invalid date'
    }
  }

  return (
    <div className={`bg-white rounded-lg shadow-sm border ${className}`}>
      {/* Header with Engine Badge */}
      <div className="flex items-center justify-between p-4 border-b">
        <div className="flex items-center space-x-3">
          <h3 className="text-lg font-medium text-gray-900">Job Metadata</h3>
          <div
            className={`inline-flex items-center px-3 py-1 rounded-full text-sm font-medium ${engineBadge.className}`}
            aria-label={engineBadge.ariaLabel}
            role="status"
          >
            <span className="mr-1" aria-hidden="true">{engineBadge.icon}</span>
            {engineBadge.text}
          </div>
        </div>
        
        <div className="flex items-center space-x-2">
          {showExportButton && onExport && (
            <button
              onClick={handleExport}
              disabled={exportLoading || status !== 'completed'}
              className="inline-flex items-center px-3 py-1.5 text-sm font-medium text-blue-600 bg-blue-50 border border-blue-200 rounded-md hover:bg-blue-100 disabled:opacity-50 disabled:cursor-not-allowed focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2"
              aria-label="Export docking results with metadata"
            >
              {exportLoading ? (
                <>
                  <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-600 mr-2"></div>
                  Exporting...
                </>
              ) : (
                <>
                  <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
                  </svg>
                  Export
                </>
              )}
            </button>
          )}
          
          <button
            onClick={() => setIsExpanded(!isExpanded)}
            className="text-gray-400 hover:text-gray-600 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 rounded"
            aria-label={isExpanded ? 'Collapse metadata details' : 'Expand metadata details'}
            aria-expanded={isExpanded}
          >
            <svg 
              className={`w-5 h-5 transition-transform ${isExpanded ? 'rotate-180' : ''}`} 
              fill="none" 
              stroke="currentColor" 
              viewBox="0 0 24 24"
            >
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
            </svg>
          </button>
        </div>
      </div>

      {/* Quick Summary */}
      <div className="p-4 space-y-3">
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
          <div>
            <dt className="font-medium text-gray-500">Status</dt>
            <dd className={`mt-1 ${
              status === 'completed' ? 'text-green-600' :
              status === 'failed' ? 'text-red-600' :
              status === 'running' ? 'text-blue-600' :
              'text-gray-600'
            }`}>
              {status || 'Unknown'}
            </dd>
          </div>
          
          <div>
            <dt className="font-medium text-gray-500">Runtime</dt>
            <dd className="mt-1 text-gray-900">{formatDuration(duration)}</dd>
          </div>
          
          <div>
            <dt className="font-medium text-gray-500">Seed</dt>
            <dd className="mt-1 text-gray-900">
              {docking_parameters.seed !== null ? docking_parameters.seed : 'Random'}
            </dd>
          </div>
          
          <div>
            <dt className="font-medium text-gray-500">Poses</dt>
            <dd className="mt-1 text-gray-900">
              {performance_metrics.num_poses_generated || 0}
            </dd>
          </div>
        </div>

        {performance_metrics.best_affinity !== undefined && (
          <div className="p-3 bg-blue-50 rounded-lg">
            <div className="text-sm">
              <span className="font-medium text-blue-900">Best Affinity:</span>
              <span className="ml-2 text-blue-700">{performance_metrics.best_affinity} kcal/mol</span>
            </div>
          </div>
        )}
      </div>

      {/* Detailed Information (Expandable) */}
      {isExpanded && (
        <div className="border-t bg-gray-50">
          <div className="p-4 space-y-6">
            
            {/* Engine Information */}
            <div>
              <h4 className="text-sm font-semibold text-gray-900 mb-3">Engine Information</h4>
              <dl className="grid grid-cols-1 md:grid-cols-2 gap-4 text-sm">
                <div>
                  <dt className="font-medium text-gray-500">Method</dt>
                  <dd className="mt-1 text-gray-900">{engine_metadata.method || 'Unknown'}</dd>
                </div>
                <div>
                  <dt className="font-medium text-gray-500">Version</dt>
                  <dd className="mt-1 text-gray-900">{engine_metadata.version || 'Unknown'}</dd>
                </div>
                <div>
                  <dt className="font-medium text-gray-500">Calculation Time</dt>
                  <dd className="mt-1 text-gray-900">{formatDuration(engine_metadata.calculation_time)}</dd>
                </div>
                <div>
                  <dt className="font-medium text-gray-500">Mock Mode</dt>
                  <dd className="mt-1">
                    <span className={`inline-flex items-center px-2 py-0.5 rounded text-xs font-medium ${
                      engine_metadata.is_mock 
                        ? 'bg-red-100 text-red-800' 
                        : 'bg-green-100 text-green-800'
                    }`}>
                      {engine_metadata.is_mock ? 'Yes' : 'No'}
                    </span>
                  </dd>
                </div>
              </dl>
            </div>

            {/* Docking Parameters */}
            <div>
              <h4 className="text-sm font-semibold text-gray-900 mb-3">Docking Parameters</h4>
              <dl className="grid grid-cols-1 md:grid-cols-3 gap-4 text-sm">
                <div>
                  <dt className="font-medium text-gray-500">Exhaustiveness</dt>
                  <dd className="mt-1 text-gray-900">{docking_parameters.exhaustiveness || 'N/A'}</dd>
                </div>
                <div>
                  <dt className="font-medium text-gray-500">Number of Modes</dt>
                  <dd className="mt-1 text-gray-900">{docking_parameters.num_modes || 'N/A'}</dd>
                </div>
                <div>
                  <dt className="font-medium text-gray-500">Random Seed</dt>
                  <dd className="mt-1 text-gray-900">
                    {docking_parameters.seed !== null ? docking_parameters.seed : 'Auto-random'}
                  </dd>
                </div>
              </dl>
              
              <div className="mt-4">
                <h5 className="text-xs font-semibold text-gray-700 mb-2">Binding Site</h5>
                <div className="grid grid-cols-2 md:grid-cols-3 gap-4 text-xs">
                  <div>
                    <dt className="font-medium text-gray-500">Center</dt>
                    <dd className="mt-1 text-gray-900">
                      ({docking_parameters.center_x?.toFixed(1) || 0}, {docking_parameters.center_y?.toFixed(1) || 0}, {docking_parameters.center_z?.toFixed(1) || 0})
                    </dd>
                  </div>
                  <div>
                    <dt className="font-medium text-gray-500">Size</dt>
                    <dd className="mt-1 text-gray-900">
                      {docking_parameters.size_x || 0} × {docking_parameters.size_y || 0} × {docking_parameters.size_z || 0} Å
                    </dd>
                  </div>
                </div>
              </div>
            </div>

            {/* Timing Information */}
            <div>
              <h4 className="text-sm font-semibold text-gray-900 mb-3">Timing Information</h4>
              <dl className="grid grid-cols-1 md:grid-cols-3 gap-4 text-sm">
                <div>
                  <dt className="font-medium text-gray-500">Started At</dt>
                  <dd className="mt-1 text-gray-900 text-xs">{formatDate(started_at)}</dd>
                </div>
                <div>
                  <dt className="font-medium text-gray-500">Completed At</dt>
                  <dd className="mt-1 text-gray-900 text-xs">{formatDate(completed_at)}</dd>
                </div>
                <div>
                  <dt className="font-medium text-gray-500">Total Duration</dt>
                  <dd className="mt-1 text-gray-900">{formatDuration(duration)}</dd>
                </div>
              </dl>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}

export default MetadataPanel
