import React, { useState, useRef, useEffect } from 'react'
import { Link } from 'react-router-dom'
import toast from 'react-hot-toast'
import SMILESInput from '../components/SMILESInput'
import RendererSelector from '../components/RendererSelector'
import MoleculeViewer from '../components/MoleculeViewer'
import ViewControls from '../components/ViewControls'
import DockingVisualization from '../components/DockingVisualization'
import AdvancedDockingControls from '../components/AdvancedDockingControls'
import { DockingAPI } from '../utils/dockingApi'
import { ConnectionHealthChecker } from '../utils/connectionTest'

function DockingPage() {
  const [currentSMILES, setCurrentSMILES] = useState('')
  const [mode, setMode] = useState('3D')
  const [renderer, setRenderer] = useState('3Dmol.js')
  const [isValidSMILES, setIsValidSMILES] = useState(false)
  const [isMoleculeLoaded, setIsMoleculeLoaded] = useState(false)
  const [exportFunctions, setExportFunctions] = useState(null)
  const [pdbId, setPdbId] = useState('1CRN')
  
  // Docking state
  const [proteinInfo, setProteinInfo] = useState(null)
  const [proteinLoaded, setProteinLoaded] = useState(false)
  const [bindingSite, setBindingSite] = useState({
    center_x: 0,
    center_y: 0,
    center_z: 0,
    size_x: 20,
    size_y: 20,
    size_z: 20
  })
  const [dockingJob, setDockingJob] = useState(null)
  const [dockingResults, setDockingResults] = useState(null)
  const [isLoading, setIsLoading] = useState(false)
  const [activeTab, setActiveTab] = useState('ligand') // 'ligand', 'docking', 'advanced'
  
  // Connection health state
  const [backendConnected, setBackendConnected] = useState(true)
  const [connectionChecker] = useState(() => new ConnectionHealthChecker())
  
  // Refs for form inputs
  const centerXRef = useRef(null)
  const centerYRef = useRef(null)
  const centerZRef = useRef(null)
  const sizeXRef = useRef(null)
  const sizeYRef = useRef(null)
  const sizeZRef = useRef(null)
  
  // Connection health check on component mount
  useEffect(() => {
    const checkConnection = async () => {
      const healthResult = await connectionChecker.checkBackendHealth()
      setBackendConnected(healthResult.isConnected)
      
      if (!healthResult.isConnected) {
        toast.error('Backend server is not responding. Please ensure the Django server is running on port 8000.', {
          duration: 6000,
          id: 'backend-connection-error'
        })
      } else {
        toast.dismiss('backend-connection-error')
      }
    }
    
    checkConnection()
  }, [connectionChecker])
  
  // Handler functions
  const handleLoadProtein = async () => {
    if (!pdbId || pdbId.length !== 4) {
      toast.error('Please enter a valid 4-character PDB ID')
      return
    }
    
    // Check backend connection before making API calls
    const healthCheck = await connectionChecker.checkBackendHealth()
    if (!healthCheck.isConnected) {
      setBackendConnected(false)
      toast.error('Cannot connect to backend server. Please ensure the Django server is running on port 8000.', {
        duration: 8000,
        id: 'backend-connection-error'
      })
      return
    }
    
    setBackendConnected(true)
    setIsLoading(true)
    try {
      // Get protein information
      const info = await DockingAPI.getProteinInfo(pdbId)
      setProteinInfo(info)
      
      // Estimate binding site
      const site = await DockingAPI.estimateBindingSite(pdbId, currentSMILES || null)
      setBindingSite(site)
      
      // Update form inputs
      if (centerXRef.current) centerXRef.current.value = site.center_x
      if (centerYRef.current) centerYRef.current.value = site.center_y
      if (centerZRef.current) centerZRef.current.value = site.center_z
      if (sizeXRef.current) sizeXRef.current.value = site.size_x
      if (sizeYRef.current) sizeYRef.current.value = site.size_y
      if (sizeZRef.current) sizeZRef.current.value = site.size_z
      
      setProteinLoaded(true)
      toast.success(`Loaded protein ${pdbId}: ${info.title || 'Unknown'}`)
      toast.dismiss('backend-connection-error')
    } catch (error) {
      console.error('Failed to load protein:', error)
      
      // Check if it's a connection error
      if (error.message.includes('Failed to fetch') || error.message.includes('ERR_CONNECTION_REFUSED')) {
        setBackendConnected(false)
        toast.error('Backend server connection lost. Please restart the Django server on port 8000.', {
          duration: 8000,
          id: 'backend-connection-error'
        })
      } else {
        toast.error(`Failed to load protein: ${error.message}`)
      }
    } finally {
      setIsLoading(false)
    }
  }
  
  const handleRunDocking = async () => {
    if (!currentSMILES) {
      toast.error('Please enter a SMILES string for the ligand')
      return
    }
    
    if (!proteinLoaded) {
      toast.error('Please load a protein first')
      return
    }
    
    // Get current binding site parameters
    const dockingParams = {
      ligand_smiles: currentSMILES,
      receptor_pdb_id: pdbId,
      center_x: parseFloat(centerXRef.current?.value || 0),
      center_y: parseFloat(centerYRef.current?.value || 0),
      center_z: parseFloat(centerZRef.current?.value || 0),
      size_x: parseFloat(sizeXRef.current?.value || 20),
      size_y: parseFloat(sizeYRef.current?.value || 20),
      size_z: parseFloat(sizeZRef.current?.value || 20),
      exhaustiveness: 8,
      num_modes: 9
    }
    
    setIsLoading(true)
    setDockingResults(null)
    
    try {
      // Start docking job
      const job = await DockingAPI.runDocking(dockingParams)
      setDockingJob(job)
      toast.success(`Docking job started: ${job.job_id}`)
      
      // Poll for results
      const results = await DockingAPI.pollDockingJob(
        job.job_id,
        (status) => {
          // Update job status
          setDockingJob(status)
          if (status.status === 'running') {
            toast.loading('Docking calculation in progress...', { id: 'docking-progress' })
          }
        }
      )
      
      setDockingResults(results.results)
      setActiveTab('docking') // Switch to docking results tab
      toast.success('Docking completed successfully!', { id: 'docking-progress' })
    } catch (error) {
      console.error('Docking failed:', error)
      toast.error(`Docking failed: ${error.message}`, { id: 'docking-progress' })
    } finally {
      setIsLoading(false)
    }
  }

  const handleAdvancedDocking = async (params, isPreview = false) => {
    if (isPreview) {
      // Update binding site parameters for preview
      if (centerXRef.current) centerXRef.current.value = params.center_x
      if (centerYRef.current) centerYRef.current.value = params.center_y
      if (centerZRef.current) centerZRef.current.value = params.center_z
      if (sizeXRef.current) sizeXRef.current.value = params.size_x
      if (sizeYRef.current) sizeYRef.current.value = params.size_y
      if (sizeZRef.current) sizeZRef.current.value = params.size_z
      return
    }

    if (!currentSMILES) {
      toast.error('Please enter a SMILES string for the ligand')
      return
    }

    if (!proteinLoaded) {
      toast.error('Please load a protein first')
      return
    }

    setIsLoading(true)
    setDockingResults(null)

    try {
      // Use advanced docking endpoint
      const response = await fetch('/api/dock/advanced/run', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          ligand_smiles: currentSMILES,
          receptor_pdb_id: pdbId,
          center_x: parseFloat(centerXRef.current?.value || 0),
          center_y: parseFloat(centerYRef.current?.value || 0),
          center_z: parseFloat(centerZRef.current?.value || 0),
          size_x: parseFloat(sizeXRef.current?.value || 20),
          size_y: parseFloat(sizeYRef.current?.value || 20),
          size_z: parseFloat(sizeZRef.current?.value || 20),
          ...params // Include all advanced parameters
        })
      })

      if (!response.ok) {
        throw new Error(`Advanced docking failed: ${response.statusText}`)
      }

      const job = await response.json()
      setDockingJob(job)
      toast.success(`Advanced docking job started: ${job.job_name || job.job_id}`)

      // Poll for results using existing polling logic
      const results = await DockingAPI.pollDockingJob(
        job.job_id,
        (status) => {
          setDockingJob(status)
          if (status.status === 'running') {
            toast.loading('Advanced docking in progress...', { id: 'advanced-docking-progress' })
          }
        }
      )

      setDockingResults(results.results)
      setActiveTab('docking') // Switch to docking results tab
      toast.success('Advanced docking completed successfully!', { id: 'advanced-docking-progress' })

    } catch (error) {
      console.error('Advanced docking failed:', error)
      toast.error(`Advanced docking failed: ${error.message}`, { id: 'advanced-docking-progress' })
    } finally {
      setIsLoading(false)
    }
  }

  return (
    <div className="min-h-screen bg-gray-50">
      <header className="bg-white shadow-sm border-b">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-6 flex items-center justify-between">
          <div>
            <h1 className="text-3xl font-bold text-gray-900">Docking</h1>
            <p className="mt-1 text-sm text-gray-600">Protein-Ligand docking with all Lipid Viewer capabilities</p>
          </div>
          <div className="flex items-center space-x-4">
            {/* Backend Connection Status */}
            <div className="flex items-center space-x-2">
              <div className={`w-3 h-3 rounded-full ${backendConnected ? 'bg-green-500' : 'bg-red-500'}`}></div>
              <span className={`text-sm ${backendConnected ? 'text-green-700' : 'text-red-700'}`}>
                Backend {backendConnected ? 'Connected' : 'Disconnected'}
              </span>
            </div>
            <Link to="/" className="text-blue-600 hover:underline">Back to Lipid Viewer</Link>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          <div className="lg:col-span-2 space-y-6">
            {/* Tab Navigation */}
            <div className="bg-white rounded-lg shadow-lg">
              <div className="border-b border-gray-200">
                <nav className="flex space-x-8 px-6 py-3">
                  <button
                    onClick={() => setActiveTab('ligand')}
                    className={`py-2 px-1 border-b-2 font-medium text-sm ${
                      activeTab === 'ligand'
                        ? 'border-blue-500 text-blue-600'
                        : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
                    }`}
                  >
                    Ligand Viewer
                  </button>
                  <button
                    onClick={() => setActiveTab('docking')}
                    className={`py-2 px-1 border-b-2 font-medium text-sm ${
                      activeTab === 'docking'
                        ? 'border-blue-500 text-blue-600'
                        : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
                    }`}
                    disabled={!dockingResults}
                  >
                    Docking Results
                    {dockingResults && (
                      <span className="ml-2 inline-flex items-center px-2 py-0.5 rounded-full text-xs font-medium bg-green-100 text-green-800">
                        {dockingResults.poses?.length || 0} poses
                      </span>
                    )}
                  </button>
                  <button
                    onClick={() => setActiveTab('advanced')}
                    className={`py-2 px-1 border-b-2 font-medium text-sm ${
                      activeTab === 'advanced'
                        ? 'border-blue-500 text-blue-600'
                        : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
                    }`}
                  >
                    Advanced
                    <span className="ml-2 inline-flex items-center px-2 py-0.5 rounded-full text-xs font-medium bg-purple-100 text-purple-800">
                      Expert
                    </span>
                  </button>
                </nav>
              </div>

              {/* Tab Content */}
              <div className="p-6">
                {activeTab === 'ligand' && (
                  <div className="space-y-4">
                    <h2 className="text-xl font-semibold">Ligand Visualization</h2>
                    <SMILESInput onSubmit={(s)=>{setCurrentSMILES(s); setIsMoleculeLoaded(true)}} onValidation={(_,v)=>setIsValidSMILES(v)} isValid={isValidSMILES} />
                    <div className="mt-4">
                      <RendererSelector mode={mode} renderer={renderer} onModeChange={setMode} onRendererChange={setRenderer} />
                    </div>
                    {currentSMILES && activeTab === 'ligand' && (
                      <div className="mt-4">
                        <MoleculeViewer 
                          key={`ligand-viewer-${activeTab}`}
                          smiles={currentSMILES} 
                          mode={mode} 
                          renderer={renderer} 
                          onExport={setExportFunctions} 
                        />
                      </div>
                    )}
                    {isMoleculeLoaded && (
                      <div className="mt-4">
                        <ViewControls mode={mode} onZoomIn={()=>{}} onZoomOut={()=>{}} onResetView={()=>{}} onExportPNG={()=>exportFunctions?.png()} onExportSVG={()=>exportFunctions?.svg()} onExportGLTF={()=>exportFunctions?.gltf()} isLoaded={isMoleculeLoaded} />
                      </div>
                    )}
                  </div>
                )}

                {activeTab === 'docking' && (
                  <div className="space-y-4">
                    <div className="flex items-center justify-between">
                      <h2 className="text-xl font-semibold">Docking Visualization</h2>
                      {dockingResults && (
                        <div className="text-sm text-gray-600">
                          Best Score: {dockingResults.summary?.best_affinity} kcal/mol
                        </div>
                      )}
                    </div>
                    
                    {dockingResults && activeTab === 'docking' ? (
                      <DockingVisualization
                        key={`docking-viewer-${activeTab}`}
                        ligandSmiles={currentSMILES}
                        receptorPdbId={pdbId}
                        dockingResults={dockingResults}
                        mode="3D"
                        renderer="3Dmol.js"
                        className="w-full"
                      />
                    ) : (
                      <div className="text-center py-12 text-gray-500">
                        <div className="text-lg font-medium mb-2">No Docking Results</div>
                        <div className="text-sm">Run a docking calculation to see results here</div>
                      </div>
                    )}
                  </div>
                )}

                {activeTab === 'advanced' && (
                  <AdvancedDockingControls
                    onAdvancedDocking={handleAdvancedDocking}
                    isLoading={isLoading}
                    disabled={!currentSMILES || !proteinLoaded}
                    currentParams={{
                      center_x: parseFloat(centerXRef.current?.value || 0),
                      center_y: parseFloat(centerYRef.current?.value || 0),
                      center_z: parseFloat(centerZRef.current?.value || 0),
                      size_x: parseFloat(sizeXRef.current?.value || 20),
                      size_y: parseFloat(sizeYRef.current?.value || 20),
                      size_z: parseFloat(sizeZRef.current?.value || 20)
                    }}
                    receptorPdbId={pdbId}
                  />
                )}
              </div>
            </div>
          </div>
          <div className="space-y-6">
            <div className="bg-white rounded-lg shadow-lg p-6">
              <h2 className="text-xl font-semibold mb-4">Protein Loader</h2>
              <label className="block text-sm font-medium text-gray-700 mb-2">PDB ID</label>
              <div className="flex space-x-2">
                <input 
                  value={pdbId} 
                  onChange={(e)=>setPdbId(e.target.value.toUpperCase())} 
                  className="flex-1 border rounded px-3 py-2" 
                  placeholder="e.g., 1CRN"
                  disabled={isLoading}
                />
                <button 
                  onClick={handleLoadProtein}
                  disabled={isLoading}
                  className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed"
                >
                  {isLoading ? 'Loading...' : 'Load'}
                </button>
              </div>
              {proteinInfo && (
                <div className="mt-3 p-3 bg-green-50 rounded">
                  <p className="text-sm text-green-800">
                    <strong>{proteinInfo.title || 'Unknown protein'}</strong>
                  </p>
                  {proteinInfo.molecular_weight && (
                    <p className="text-xs text-green-600">
                      MW: {proteinInfo.molecular_weight} Da
                    </p>
                  )}
                </div>
              )}
              <p className="text-xs text-gray-500 mt-2">Backend: /api/pdb/{'{pdbId}'}/info</p>
            </div>
            <div className="bg-white rounded-lg shadow-lg p-6">
              <h2 className="text-xl font-semibold mb-4">Docking Controls</h2>
              <div className="grid grid-cols-3 gap-3 text-sm">
                <div>
                  <label className="block text-gray-600 mb-1">Center X</label>
                  <input 
                    ref={centerXRef}
                    className="w-full border rounded px-2 py-1" 
                    defaultValue={bindingSite.center_x}
                    type="number"
                    step="0.1"
                    disabled={isLoading}
                  />
                </div>
                <div>
                  <label className="block text-gray-600 mb-1">Center Y</label>
                  <input 
                    ref={centerYRef}
                    className="w-full border rounded px-2 py-1" 
                    defaultValue={bindingSite.center_y}
                    type="number"
                    step="0.1"
                    disabled={isLoading}
                  />
                </div>
                <div>
                  <label className="block text-gray-600 mb-1">Center Z</label>
                  <input 
                    ref={centerZRef}
                    className="w-full border rounded px-2 py-1" 
                    defaultValue={bindingSite.center_z}
                    type="number"
                    step="0.1"
                    disabled={isLoading}
                  />
                </div>
                <div>
                  <label className="block text-gray-600 mb-1">Size X</label>
                  <input 
                    ref={sizeXRef}
                    className="w-full border rounded px-2 py-1" 
                    defaultValue={bindingSite.size_x}
                    type="number"
                    step="0.1"
                    min="1"
                    disabled={isLoading}
                  />
                </div>
                <div>
                  <label className="block text-gray-600 mb-1">Size Y</label>
                  <input 
                    ref={sizeYRef}
                    className="w-full border rounded px-2 py-1" 
                    defaultValue={bindingSite.size_y}
                    type="number"
                    step="0.1"
                    min="1"
                    disabled={isLoading}
                  />
                </div>
                <div>
                  <label className="block text-gray-600 mb-1">Size Z</label>
                  <input 
                    ref={sizeZRef}
                    className="w-full border rounded px-2 py-1" 
                    defaultValue={bindingSite.size_z}
                    type="number"
                    step="0.1"
                    min="1"
                    disabled={isLoading}
                  />
                </div>
              </div>
              <button 
                onClick={handleRunDocking}
                disabled={isLoading || !proteinLoaded || !currentSMILES}
                className="mt-4 w-full px-4 py-2 bg-emerald-600 text-white rounded hover:bg-emerald-700 disabled:opacity-50 disabled:cursor-not-allowed"
              >
                {isLoading ? 'Running...' : 'Run Docking'}
              </button>
              {dockingJob && (
                <div className="mt-3 p-3 bg-blue-50 rounded">
                  <p className="text-sm text-blue-800">
                    Job: {dockingJob.job_id ? dockingJob.job_id.substring(0, 8) : 'N/A'}...
                  </p>
                  <p className="text-xs text-blue-600">
                    Status: {dockingJob.status}
                  </p>
                </div>
              )}
            </div>
            <div className="bg-white rounded-lg shadow-lg p-6">
              <h2 className="text-xl font-semibold mb-4">Results</h2>
              {!dockingResults ? (
                <p className="text-sm text-gray-500">No results yet.</p>
              ) : (
                <div className="space-y-4">
                  <div className="p-3 bg-gray-50 rounded">
                    <h3 className="text-sm font-medium">Summary</h3>
                    <p className="text-xs text-gray-600">
                      Best Affinity: {dockingResults.summary?.best_affinity} kcal/mol
                    </p>
                    <p className="text-xs text-gray-600">
                      Poses: {dockingResults.summary?.num_poses}
                    </p>
                    <p className="text-xs text-gray-600">
                      Time: {dockingResults.summary?.calculation_time}s
                    </p>
                  </div>
                  
                  <div className="max-h-60 overflow-y-auto">
                    <h3 className="text-sm font-medium mb-2">Top Poses</h3>
                    {dockingResults.poses?.slice(0, 5).map((pose, index) => (
                      <div key={index} className="p-2 border rounded mb-2 text-xs">
                        <div className="flex justify-between items-center">
                          <span className="font-medium">Mode {pose.mode}</span>
                          <span className="text-green-600 font-bold">
                            {pose.affinity} kcal/mol
                          </span>
                        </div>
                        <div className="text-gray-500 mt-1">
                          RMSD: {pose.rmsd_lb} - {pose.rmsd_ub} Å
                        </div>
                        <div className="text-gray-500">
                          Center: ({pose.center_x}, {pose.center_y}, {pose.center_z})
                        </div>
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </div>
          </div>
        </div>
      </main>
    </div>
  )
}

export default DockingPage