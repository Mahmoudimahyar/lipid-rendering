import React from 'react'

export default function DockingParameters({ params, setParams }) {
  const update = (k, v) => setParams({ ...params, [k]: v })
  return (
    <div className="bg-white rounded-lg shadow p-4">
      <h3 className="text-lg font-semibold mb-2">Docking Parameters</h3>
      <div className="grid grid-cols-3 gap-3">
        <div>
          <label htmlFor="exhaustiveness" className="block text-sm text-gray-700">Exhaustiveness</label>
          <input id="exhaustiveness" type="number" value={params.exhaustiveness}
            onChange={(e)=>update('exhaustiveness', parseInt(e.target.value||8))}
            className="mt-1 w-full rounded border-gray-300" />
        </div>
        <div>
          <label htmlFor="num_modes" className="block text-sm text-gray-700">Num Modes</label>
          <input id="num_modes" type="number" value={params.num_modes}
            onChange={(e)=>update('num_modes', parseInt(e.target.value||9))}
            className="mt-1 w-full rounded border-gray-300" />
        </div>
        <div>
          <label htmlFor="seed" className="block text-sm text-gray-700">Seed</label>
          <input id="seed" type="number" value={params.seed ?? ''}
            onChange={(e)=>update('seed', e.target.value===''? null : parseInt(e.target.value))}
            className="mt-1 w-full rounded border-gray-300" />
        </div>
      </div>
    </div>
  )
}



