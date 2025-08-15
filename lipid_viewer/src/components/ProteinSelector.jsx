import React from 'react'

export default function ProteinSelector({ pdbId, setPdbId, onLoad, isLoading }) {
  return (
    <div className="bg-white rounded-lg shadow p-4">
      <h3 className="text-lg font-semibold mb-2">Protein Selector</h3>
      <div className="flex items-end space-x-3">
        <div className="flex-1">
          <label htmlFor="pdb_id" className="block text-sm font-medium text-gray-700">PDB ID</label>
          <input
            id="pdb_id"
            type="text"
            value={pdbId}
            onChange={(e) => setPdbId(e.target.value.toUpperCase())}
            maxLength={4}
            className="mt-1 block w-full rounded-md border-gray-300 shadow-sm focus:border-blue-500 focus:ring-blue-500"
            placeholder="e.g., 1CRN"
          />
        </div>
        <button
          onClick={onLoad}
          disabled={isLoading}
          className={`inline-flex items-center px-4 py-2 border border-transparent text-sm font-medium rounded-md shadow-sm text-white ${isLoading ? 'bg-gray-400' : 'bg-blue-600 hover:bg-blue-700'} focus:outline-none`}
        >
          {isLoading ? 'Loading…' : 'Load Protein'}
        </button>
      </div>
    </div>
  )
}



