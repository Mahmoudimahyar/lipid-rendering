import React from 'react'

export default function JobManager({ onRunDocking, onExport, progress, isLoading }) {
  const percent = progress?.percent ?? 0
  const message = progress?.message || (isLoading ? 'Running…' : 'Idle')
  const eta = progress?.estimated_seconds_remaining
  return (
    <div className="bg-white rounded-lg shadow p-4">
      <h3 className="text-lg font-semibold mb-2">Job Manager</h3>
      <div className="flex items-center space-x-3 mb-3">
        <button
          onClick={onRunDocking}
          disabled={isLoading}
          className={`px-4 py-2 rounded text-white ${isLoading ? 'bg-gray-400' : 'bg-green-600 hover:bg-green-700'}`}
        >
          {isLoading ? 'Running…' : 'Run Docking'}
        </button>
        <button
          onClick={onExport}
          className="px-4 py-2 rounded border border-gray-300 hover:bg-gray-50"
        >
          Export Results
        </button>
      </div>
      <div>
        <div className="h-2 bg-gray-200 rounded">
          <div className="h-2 bg-blue-600 rounded" style={{ width: `${percent}%`, transition: 'width 300ms ease' }} />
        </div>
        <div className="mt-2 text-sm text-gray-700 flex justify-between">
          <span>{message}</span>
          <span>{percent.toFixed(0)}%</span>
          {eta != null && <span>ETA: ~{Math.ceil(eta)}s</span>}
        </div>
      </div>
    </div>
  )
}




