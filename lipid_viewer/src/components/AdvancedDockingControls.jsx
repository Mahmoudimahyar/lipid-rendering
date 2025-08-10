import React, { useState, useEffect } from 'react'

/**
 * Advanced docking controls for expert users
 * Includes job templates, advanced parameters, pocket detection, and GNINA rescoring
 */
function AdvancedDockingControls({
  onAdvancedDocking,
  isLoading,
  disabled,
  currentParams,
  receptorPdbId
}) {
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [jobName, setJobName] = useState('')
  const [jobDescription, setJobDescription] = useState('')
  const [tags, setTags] = useState('')
  const [selectedTemplate, setSelectedTemplate] = useState('')
  const [templates, setTemplates] = useState([])
  const [detectedPockets, setDetectedPockets] = useState([])
  const [selectedPocket, setSelectedPocket] = useState('')
  
  // Advanced options
  const [usePocketDetection, setUsePocketDetection] = useState(false)
  const [useGninaRescoring, setUseGninaRescoring] = useState(false)
  const [advancedParams, setAdvancedParams] = useState({
    exhaustiveness: 8,
    num_modes: 9,
    energy_range: 3.0,
    cpu_threads: 1,
    seed: 0,
    local_only: false,
    randomize_only: false
  })

  // Load templates on component mount
  useEffect(() => {
    loadJobTemplates()
  }, [])

  // Load pockets when receptor changes
  useEffect(() => {
    if (receptorPdbId && usePocketDetection) {
      loadDetectedPockets()
    }
  }, [receptorPdbId, usePocketDetection])

  const loadJobTemplates = async () => {
    try {
      const response = await fetch('/api/templates')
      if (response.ok) {
        const data = await response.json()
        setTemplates(data.templates || [])
      }
    } catch (error) {
      console.error('Failed to load templates:', error)
    }
  }

  const loadDetectedPockets = async () => {
    if (!receptorPdbId) return
    
    try {
      const response = await fetch(`/api/pockets/${receptorPdbId}`)
      if (response.ok) {
        const data = await response.json()
        setDetectedPockets(data.pockets || [])
      }
    } catch (error) {
      console.error('Failed to load pockets:', error)
    }
  }

  const handleTemplateChange = (templateId) => {
    setSelectedTemplate(templateId)
    const template = templates.find(t => t.template_id === templateId)
    if (template) {
      setAdvancedParams({
        ...advancedParams,
        ...template.default_params
      })
      setJobName(template.name)
      setJobDescription(template.description)
      
      // Set advanced features based on template
      if (template.advanced_params) {
        setUsePocketDetection(template.advanced_params.use_pocket_detection || false)
        setUseGninaRescoring(template.advanced_params.use_gnina_rescoring || false)
      }
    }
  }

  const handlePocketSelect = (pocketId) => {
    setSelectedPocket(pocketId)
    const pocket = detectedPockets.find(p => p.pocket_id === pocketId)
    if (pocket && onAdvancedDocking) {
      // Update binding site coordinates based on selected pocket
      const updatedParams = {
        ...currentParams,
        center_x: pocket.center_x,
        center_y: pocket.center_y,
        center_z: pocket.center_z,
        size_x: pocket.radius * 2,
        size_y: pocket.radius * 2,
        size_z: pocket.radius * 2
      }
      // Notify parent component of parameter changes
      onAdvancedDocking(updatedParams, true) // preview mode
    }
  }

  const handleAdvancedDocking = () => {
    if (!onAdvancedDocking) return

    const dockingParams = {
      ...currentParams,
      ...advancedParams,
      job_name: jobName,
      job_description: jobDescription,
      tags: tags.split(',').map(t => t.trim()).filter(t => t),
      template_id: selectedTemplate,
      use_pocket_detection: usePocketDetection,
      use_gnina_rescoring: useGninaRescoring
    }

    onAdvancedDocking(dockingParams, false) // actual run
  }

  const detectPockets = async () => {
    if (!receptorPdbId) return

    try {
      const response = await fetch('/api/pockets/detect', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          pdb_id: receptorPdbId
        })
      })

      if (response.ok) {
        const data = await response.json()
        setDetectedPockets(data.pockets || [])
      } else {
        throw new Error('Failed to detect pockets')
      }
    } catch (error) {
      console.error('Pocket detection failed:', error)
    }
  }

  return (
    <div className="bg-white rounded-lg shadow-lg p-6">
      <div className="flex items-center justify-between mb-4">
        <h2 className="text-xl font-semibold">Advanced Docking</h2>
        <button
          onClick={() => setShowAdvanced(!showAdvanced)}
          className="text-sm text-blue-600 hover:text-blue-700"
        >
          {showAdvanced ? 'Hide Advanced' : 'Show Advanced'}
        </button>
      </div>

      {/* Basic Job Info */}
      <div className="space-y-4 mb-6">
        <div>
          <label className="block text-sm font-medium text-gray-700 mb-1">
            Job Name
          </label>
          <input
            type="text"
            value={jobName}
            onChange={(e) => setJobName(e.target.value)}
            className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm"
            placeholder="Optional job name"
            disabled={isLoading}
          />
        </div>

        <div>
          <label className="block text-sm font-medium text-gray-700 mb-1">
            Description
          </label>
          <textarea
            value={jobDescription}
            onChange={(e) => setJobDescription(e.target.value)}
            className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm"
            rows={2}
            placeholder="Optional job description"
            disabled={isLoading}
          />
        </div>

        <div>
          <label className="block text-sm font-medium text-gray-700 mb-1">
            Tags (comma-separated)
          </label>
          <input
            type="text"
            value={tags}
            onChange={(e) => setTags(e.target.value)}
            className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm"
            placeholder="e.g., fragment, high-throughput, kinase"
            disabled={isLoading}
          />
        </div>
      </div>

      {showAdvanced && (
        <div className="space-y-6">
          {/* Job Templates */}
          <div>
            <label className="block text-sm font-medium text-gray-700 mb-2">
              Job Template
            </label>
            <select
              value={selectedTemplate}
              onChange={(e) => handleTemplateChange(e.target.value)}
              className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm"
              disabled={isLoading}
            >
              <option value="">Select a template (optional)</option>
              {templates.map((template) => (
                <option key={template.template_id} value={template.template_id}>
                  {template.name} - {template.description}
                </option>
              ))}
            </select>
            {selectedTemplate && (
              <p className="text-xs text-gray-500 mt-1">
                Used {templates.find(t => t.template_id === selectedTemplate)?.usage_count || 0} times
              </p>
            )}
          </div>

          {/* Advanced Features */}
          <div>
            <h3 className="text-lg font-medium mb-3">Advanced Features</h3>
            <div className="space-y-3">
              <label className="flex items-center">
                <input
                  type="checkbox"
                  checked={usePocketDetection}
                  onChange={(e) => setUsePocketDetection(e.target.checked)}
                  className="mr-2"
                  disabled={isLoading}
                />
                <span className="text-sm">
                  Use automatic pocket detection
                  <span className="text-gray-500 ml-1">(fpocket algorithm)</span>
                </span>
              </label>

              <label className="flex items-center">
                <input
                  type="checkbox"
                  checked={useGninaRescoring}
                  onChange={(e) => setUseGninaRescoring(e.target.checked)}
                  className="mr-2"
                  disabled={isLoading}
                />
                <span className="text-sm">
                  Use GNINA rescoring
                  <span className="text-gray-500 ml-1">(CNN-based scoring)</span>
                </span>
              </label>
            </div>
          </div>

          {/* Pocket Detection Panel */}
          {usePocketDetection && (
            <div>
              <div className="flex items-center justify-between mb-3">
                <h3 className="text-lg font-medium">Detected Pockets</h3>
                <button
                  onClick={detectPockets}
                  className="px-3 py-1 text-sm bg-blue-100 text-blue-700 rounded hover:bg-blue-200"
                  disabled={isLoading || !receptorPdbId}
                >
                  Detect Pockets
                </button>
              </div>

              {detectedPockets.length > 0 ? (
                <div className="max-h-40 overflow-y-auto space-y-2">
                  {detectedPockets.map((pocket, index) => (
                    <div
                      key={pocket.pocket_id}
                      className={`p-3 border rounded cursor-pointer ${
                        selectedPocket === pocket.pocket_id
                          ? 'border-blue-500 bg-blue-50'
                          : 'border-gray-200 hover:border-gray-300'
                      }`}
                      onClick={() => handlePocketSelect(pocket.pocket_id)}
                    >
                      <div className="flex justify-between items-center">
                        <span className="font-medium text-sm">Pocket {index + 1}</span>
                        <span className="text-sm text-green-600">
                          Druggability: {(pocket.druggability_score * 100).toFixed(0)}%
                        </span>
                      </div>
                      <div className="text-xs text-gray-500 mt-1">
                        Center: ({pocket.center_x.toFixed(1)}, {pocket.center_y.toFixed(1)}, {pocket.center_z.toFixed(1)})
                        | Volume: {pocket.volume.toFixed(0)} Ų
                      </div>
                    </div>
                  ))}
                </div>
              ) : (
                <p className="text-sm text-gray-500">
                  {receptorPdbId 
                    ? 'Click "Detect Pockets" to find binding sites'
                    : 'Load a protein first to detect pockets'
                  }
                </p>
              )}
            </div>
          )}

          {/* Advanced Parameters */}
          <div>
            <h3 className="text-lg font-medium mb-3">Docking Parameters</h3>
            <div className="grid grid-cols-2 gap-4">
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-1">
                  Exhaustiveness
                </label>
                <input
                  type="number"
                  value={advancedParams.exhaustiveness}
                  onChange={(e) => setAdvancedParams({
                    ...advancedParams,
                    exhaustiveness: parseInt(e.target.value)
                  })}
                  min="1"
                  max="32"
                  className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm"
                  disabled={isLoading}
                />
                <p className="text-xs text-gray-500 mt-1">1-32, higher = more thorough</p>
              </div>

              <div>
                <label className="block text-sm font-medium text-gray-700 mb-1">
                  Number of Modes
                </label>
                <input
                  type="number"
                  value={advancedParams.num_modes}
                  onChange={(e) => setAdvancedParams({
                    ...advancedParams,
                    num_modes: parseInt(e.target.value)
                  })}
                  min="1"
                  max="50"
                  className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm"
                  disabled={isLoading}
                />
                <p className="text-xs text-gray-500 mt-1">Number of binding modes</p>
              </div>

              <div>
                <label className="block text-sm font-medium text-gray-700 mb-1">
                  Energy Range (kcal/mol)
                </label>
                <input
                  type="number"
                  value={advancedParams.energy_range}
                  onChange={(e) => setAdvancedParams({
                    ...advancedParams,
                    energy_range: parseFloat(e.target.value)
                  })}
                  min="1"
                  max="10"
                  step="0.1"
                  className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm"
                  disabled={isLoading}
                />
                <p className="text-xs text-gray-500 mt-1">Maximum energy difference</p>
              </div>

              <div>
                <label className="block text-sm font-medium text-gray-700 mb-1">
                  Random Seed
                </label>
                <input
                  type="number"
                  value={advancedParams.seed}
                  onChange={(e) => setAdvancedParams({
                    ...advancedParams,
                    seed: parseInt(e.target.value)
                  })}
                  min="0"
                  className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm"
                  disabled={isLoading}
                />
                <p className="text-xs text-gray-500 mt-1">0 = random, &gt;0 = reproducible</p>
              </div>
            </div>

            <div className="mt-4 space-y-3">
              <label className="flex items-center">
                <input
                  type="checkbox"
                  checked={advancedParams.local_only}
                  onChange={(e) => setAdvancedParams({
                    ...advancedParams,
                    local_only: e.target.checked
                  })}
                  className="mr-2"
                  disabled={isLoading}
                />
                <span className="text-sm">Local search only (no global search)</span>
              </label>

              <label className="flex items-center">
                <input
                  type="checkbox"
                  checked={advancedParams.randomize_only}
                  onChange={(e) => setAdvancedParams({
                    ...advancedParams,
                    randomize_only: e.target.checked
                  })}
                  className="mr-2"
                  disabled={isLoading}
                />
                <span className="text-sm">Randomize ligand only (no docking)</span>
              </label>
            </div>
          </div>
        </div>
      )}

      {/* Run Button */}
      <div className="mt-6 pt-6 border-t">
        <button
          onClick={handleAdvancedDocking}
          disabled={isLoading || disabled}
          className="w-full px-4 py-3 bg-purple-600 text-white rounded-lg hover:bg-purple-700 disabled:opacity-50 disabled:cursor-not-allowed font-medium"
        >
          {isLoading ? 'Running Advanced Docking...' : 'Run Advanced Docking'}
        </button>
        
        {(usePocketDetection || useGninaRescoring) && (
          <div className="mt-2 text-xs text-gray-600 text-center">
            Estimated time: {useGninaRescoring ? '5-10' : '3-6'} minutes
            {usePocketDetection && ' (includes pocket detection)'}
          </div>
        )}
      </div>
    </div>
  )
}

export default AdvancedDockingControls