import React, { useRef, useEffect } from 'react'

export default function BindingSiteConfigurator({ bindingSite, setBindingSite, refs }) {
  const { centerXRef, centerYRef, centerZRef, sizeXRef, sizeYRef, sizeZRef } = refs

  useEffect(() => {
    if (centerXRef.current) centerXRef.current.value = bindingSite.center_x
    if (centerYRef.current) centerYRef.current.value = bindingSite.center_y
    if (centerZRef.current) centerZRef.current.value = bindingSite.center_z
    if (sizeXRef.current) sizeXRef.current.value = bindingSite.size_x
    if (sizeYRef.current) sizeYRef.current.value = bindingSite.size_y
    if (sizeZRef.current) sizeZRef.current.value = bindingSite.size_z
  }, [bindingSite, centerXRef, centerYRef, centerZRef, sizeXRef, sizeYRef, sizeZRef])

  const handleChange = () => {
    setBindingSite({
      center_x: parseFloat(centerXRef.current?.value || 0),
      center_y: parseFloat(centerYRef.current?.value || 0),
      center_z: parseFloat(centerZRef.current?.value || 0),
      size_x: parseFloat(sizeXRef.current?.value || 20),
      size_y: parseFloat(sizeYRef.current?.value || 20),
      size_z: parseFloat(sizeZRef.current?.value || 20),
    })
  }

  return (
    <div className="bg-white rounded-lg shadow p-4">
      <h3 className="text-lg font-semibold mb-2">Binding Site Configurator</h3>
      <div className="grid grid-cols-6 gap-3">
        <div>
          <label htmlFor="center_x" className="block text-sm text-gray-700">Center X</label>
          <input id="center_x" ref={centerXRef} onChange={handleChange} type="number" step="0.1" className="mt-1 w-full rounded border-gray-300" />
        </div>
        <div>
          <label htmlFor="center_y" className="block text-sm text-gray-700">Center Y</label>
          <input id="center_y" ref={centerYRef} onChange={handleChange} type="number" step="0.1" className="mt-1 w-full rounded border-gray-300" />
        </div>
        <div>
          <label htmlFor="center_z" className="block text-sm text-gray-700">Center Z</label>
          <input id="center_z" ref={centerZRef} onChange={handleChange} type="number" step="0.1" className="mt-1 w-full rounded border-gray-300" />
        </div>
        <div>
          <label htmlFor="size_x" className="block text-sm text-gray-700">Size X</label>
          <input id="size_x" ref={sizeXRef} onChange={handleChange} type="number" step="0.1" className="mt-1 w-full rounded border-gray-300" />
        </div>
        <div>
          <label htmlFor="size_y" className="block text-sm text-gray-700">Size Y</label>
          <input id="size_y" ref={sizeYRef} onChange={handleChange} type="number" step="0.1" className="mt-1 w-full rounded border-gray-300" />
        </div>
        <div>
          <label htmlFor="size_z" className="block text-sm text-gray-700">Size Z</label>
          <input id="size_z" ref={sizeZRef} onChange={handleChange} type="number" step="0.1" className="mt-1 w-full rounded border-gray-300" />
        </div>
      </div>
    </div>
  )
}



