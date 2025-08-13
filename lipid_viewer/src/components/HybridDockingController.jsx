import React, { useState, useEffect } from 'react';
import WebinaIntegration from './WebinaIntegration';

/**
 * Hybrid Docking Controller
 * 
 * Intelligently switches between:
 * 1. Client-side docking (Webina) - Fast, immediate results
 * 2. Server-side docking (AWS/GPU) - Powerful, high-accuracy results
 * 3. Fallback server (CPU) - Reliable backup option
 */
const HybridDockingController = ({
  ligandSmiles,
  receptorPdbId,
  dockingParams,
  onResults,
  onError,
  onProgress
}) => {
  const [dockingMode, setDockingMode] = useState('auto'); // 'auto', 'client', 'server-gpu', 'server-cpu'
  const [availableEngines, setAvailableEngines] = useState({
    client: false,
    serverGpu: false,
    serverCpu: false
  });
  const [engineCapabilities, setEngineCapabilities] = useState(null);
  const [currentExecution, setCurrentExecution] = useState(null);
  const [executionHistory, setExecutionHistory] = useState([]);

  // Check available engines on component mount
  useEffect(() => {
    checkEngineAvailability();
  }, []);

  const checkEngineAvailability = async () => {
    onProgress?.({ stage: 'checking', message: 'Checking available docking engines...' });

    const engines = {
      client: checkClientSideSupport(),
      serverGpu: false,
      serverCpu: false
    };

    // Check server capabilities
    try {
      const response = await fetch('/api/dock/capabilities');
      const capabilities = await response.json();
      
      setEngineCapabilities(capabilities);
      engines.serverGpu = capabilities.cuda_available || false;
      engines.serverCpu = capabilities.vina_available || capabilities.mock_allowed || false;
      
      console.log('🔍 Engine capabilities detected:', capabilities);
    } catch (error) {
      console.warn('⚠️ Server capabilities check failed:', error);
      engines.serverCpu = false; // Assume server unavailable
    }

    setAvailableEngines(engines);
    
    // Auto-select best available engine
    if (dockingMode === 'auto') {
      selectOptimalEngine(engines);
    }
  };

  const checkClientSideSupport = () => {
    const hasWebAssembly = 'WebAssembly' in window;
    const hasSharedArrayBuffer = 'SharedArrayBuffer' in window;
    
    // Minimum requirement: WebAssembly support
    return hasWebAssembly;
  };

  const selectOptimalEngine = (engines) => {
    // Priority order: Client (fast) > Server GPU (powerful) > Server CPU (reliable)
    if (engines.client) {
      setDockingMode('client');
      console.log('🎯 Auto-selected: Client-side docking');
    } else if (engines.serverGpu) {
      setDockingMode('server-gpu');
      console.log('🎯 Auto-selected: Server GPU docking');
    } else if (engines.serverCpu) {
      setDockingMode('server-cpu');
      console.log('🎯 Auto-selected: Server CPU docking');
    } else {
      onError?.({
        type: 'no_engines',
        message: 'No docking engines available. Please check your browser compatibility or server status.'
      });
    }
  };

  const runDocking = async () => {
    if (!dockingMode || dockingMode === 'auto') {
      onError?.({
        type: 'no_mode',
        message: 'No docking mode selected'
      });
      return;
    }

    const execution = {
      id: `exec_${Date.now()}`,
      mode: dockingMode,
      startTime: Date.now(),
      params: dockingParams
    };

    setCurrentExecution(execution);
    onProgress?.({ stage: 'starting', message: `Starting ${dockingMode} docking...` });

    try {
      let results;
      
      switch (dockingMode) {
        case 'client':
          results = await runClientSideDocking();
          break;
        case 'server-gpu':
          results = await runServerDocking({ useGpu: true });
          break;
        case 'server-cpu':
          results = await runServerDocking({ useGpu: false });
          break;
        default:
          throw new Error(`Unknown docking mode: ${dockingMode}`);
      }

      // Record successful execution
      execution.endTime = Date.now();
      execution.duration = execution.endTime - execution.startTime;
      execution.success = true;
      execution.resultCount = results.poses?.length || 0;

      setExecutionHistory(prev => [execution, ...prev.slice(0, 9)]); // Keep last 10
      setCurrentExecution(null);

      onResults?.(results);

    } catch (error) {
      console.error(`❌ ${dockingMode} docking failed:`, error);
      
      execution.endTime = Date.now();
      execution.duration = execution.endTime - execution.startTime;
      execution.success = false;
      execution.error = error.message;

      setExecutionHistory(prev => [execution, ...prev.slice(0, 9)]);
      setCurrentExecution(null);

      // Try fallback if auto mode
      if (dockingMode === 'auto') {
        await tryFallbackEngine(error);
      } else {
        onError?.(error);
      }
    }
  };

  const runClientSideDocking = () => {
    return new Promise((resolve, reject) => {
      // This will be handled by WebinaIntegration component
      // For now, return a promise that resolves when client-side completes
      reject(new Error('Client-side docking integration in progress'));
    });
  };

  const runServerDocking = async ({ useGpu }) => {
    const endpoint = '/api/dock/run';
    const params = {
      ...dockingParams,
      ligand_smiles: ligandSmiles,
      receptor_pdb_id: receptorPdbId,
      runtime_preference: useGpu ? 'gpu' : 'cpu'
    };

    onProgress?.({ 
      stage: 'submitting', 
      message: `Submitting to server (${useGpu ? 'GPU' : 'CPU'})...` 
    });

    const response = await fetch(endpoint, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(params)
    });

    if (!response.ok) {
      const errorData = await response.json().catch(() => ({}));
      throw new Error(errorData.error || `Server error: ${response.status}`);
    }

    const jobData = await response.json();
    
    if (jobData.error) {
      throw new Error(jobData.error);
    }

    // Poll for results
    return await pollForResults(jobData.job_id);
  };

  const pollForResults = async (jobId) => {
    const maxAttempts = 120; // 10 minutes max
    const pollInterval = 5000; // 5 seconds

    for (let attempt = 0; attempt < maxAttempts; attempt++) {
      onProgress?.({ 
        stage: 'processing', 
        message: `Processing... (${attempt + 1}/${maxAttempts})`,
        percentage: (attempt / maxAttempts) * 100
      });

      await new Promise(resolve => setTimeout(resolve, pollInterval));

      try {
        const response = await fetch(`/api/dock/status/${jobId}`);
        const statusData = await response.json();

        if (statusData.status === 'completed') {
          onProgress?.({ stage: 'completed', message: 'Docking completed successfully!' });
          return statusData;
        } else if (statusData.status === 'failed') {
          throw new Error(statusData.error || 'Docking job failed');
        }
        // Continue polling if status is 'pending' or 'running'
        
      } catch (error) {
        console.warn(`Polling attempt ${attempt + 1} failed:`, error);
        if (attempt === maxAttempts - 1) {
          throw error;
        }
      }
    }

    throw new Error('Docking job timed out');
  };

  const tryFallbackEngine = async (originalError) => {
    // Try fallback engines in order
    const fallbackOrder = ['client', 'server-cpu', 'server-gpu'];
    const currentIndex = fallbackOrder.indexOf(dockingMode);
    
    for (let i = currentIndex + 1; i < fallbackOrder.length; i++) {
      const fallbackMode = fallbackOrder[i];
      if (availableEngines[fallbackMode.replace('-', '')]) {
        console.log(`🔄 Trying fallback: ${fallbackMode}`);
        setDockingMode(fallbackMode);
        await runDocking();
        return;
      }
    }

    // No fallbacks available
    onError?.({
      type: 'all_failed',
      message: `All docking engines failed. Original error: ${originalError.message}`,
      originalError
    });
  };

  const getEngineStatusIcon = (engineKey) => {
    if (!availableEngines[engineKey]) return '❌';
    if (dockingMode === engineKey || (dockingMode.includes(engineKey))) return '🟢';
    return '🟡';
  };

  const formatExecutionTime = (duration) => {
    if (duration < 1000) return `${duration}ms`;
    if (duration < 60000) return `${(duration / 1000).toFixed(1)}s`;
    return `${(duration / 60000).toFixed(1)}m`;
  };

  return (
    <div className="hybrid-docking-controller space-y-4">
      
      {/* Engine Selection */}
      <div className="bg-gray-50 p-4 rounded-lg">
        <h3 className="font-semibold mb-3">Docking Engine Selection</h3>
        
        <div className="grid grid-cols-1 md:grid-cols-4 gap-2 mb-4">
          <button
            onClick={() => setDockingMode('auto')}
            className={`p-2 rounded text-sm ${
              dockingMode === 'auto' ? 'bg-blue-600 text-white' : 'bg-white border hover:bg-blue-50'
            }`}
          >
            🤖 Auto Select
          </button>
          
          <button
            onClick={() => setDockingMode('client')}
            disabled={!availableEngines.client}
            className={`p-2 rounded text-sm ${
              dockingMode === 'client' ? 'bg-green-600 text-white' : 
              availableEngines.client ? 'bg-white border hover:bg-green-50' : 'bg-gray-200 text-gray-500'
            }`}
          >
            {getEngineStatusIcon('client')} Client-Side
          </button>
          
          <button
            onClick={() => setDockingMode('server-gpu')}
            disabled={!availableEngines.serverGpu}
            className={`p-2 rounded text-sm ${
              dockingMode === 'server-gpu' ? 'bg-purple-600 text-white' : 
              availableEngines.serverGpu ? 'bg-white border hover:bg-purple-50' : 'bg-gray-200 text-gray-500'
            }`}
          >
            {getEngineStatusIcon('serverGpu')} Server GPU
          </button>
          
          <button
            onClick={() => setDockingMode('server-cpu')}
            disabled={!availableEngines.serverCpu}
            className={`p-2 rounded text-sm ${
              dockingMode === 'server-cpu' ? 'bg-orange-600 text-white' : 
              availableEngines.serverCpu ? 'bg-white border hover:bg-orange-50' : 'bg-gray-200 text-gray-500'
            }`}
          >
            {getEngineStatusIcon('serverCpu')} Server CPU
          </button>
        </div>

        {/* Engine Details */}
        {engineCapabilities && (
          <div className="text-xs text-gray-600 bg-white p-2 rounded">
            <strong>Server Capabilities:</strong> 
            {engineCapabilities.vina_available && ` AutoDock Vina ${engineCapabilities.vina_version || 'unknown'}`}
            {engineCapabilities.cuda_available && ` • CUDA GPU`}
            {engineCapabilities.mock_allowed && ` • Mock Mode`}
            {!engineCapabilities.vina_available && !engineCapabilities.mock_allowed && ` Unavailable`}
          </div>
        )}
      </div>

      {/* Webina Integration (for client-side) */}
      {dockingMode === 'client' && (
        <WebinaIntegration
          ligandSmiles={ligandSmiles}
          receptorPdbId={receptorPdbId}
          dockingParams={dockingParams}
          onResults={onResults}
          onError={onError}
          onProgress={onProgress}
        />
      )}

      {/* Main Docking Button */}
      <div className="flex space-x-2">
        <button
          onClick={runDocking}
          disabled={currentExecution || !dockingMode || dockingMode === 'auto'}
          className={`px-6 py-3 rounded-lg font-medium ${
            currentExecution ? 'bg-gray-400 text-gray-700 cursor-not-allowed' :
            dockingMode && dockingMode !== 'auto' ? 'bg-blue-600 hover:bg-blue-700 text-white' :
            'bg-gray-300 text-gray-500 cursor-not-allowed'
          }`}
        >
          {currentExecution ? '⏳ Running...' : `🚀 Run ${dockingMode} Docking`}
        </button>
      </div>

      {/* Execution History */}
      {executionHistory.length > 0 && (
        <div className="bg-gray-50 p-3 rounded">
          <h4 className="font-medium mb-2">Recent Executions</h4>
          <div className="space-y-1 text-sm">
            {executionHistory.slice(0, 3).map((exec, index) => (
              <div key={exec.id} className="flex justify-between items-center">
                <span className={exec.success ? 'text-green-600' : 'text-red-600'}>
                  {exec.success ? '✅' : '❌'} {exec.mode} 
                  {exec.success && ` (${exec.resultCount} poses)`}
                </span>
                <span className="text-gray-500">
                  {formatExecutionTime(exec.duration)}
                </span>
              </div>
            ))}
          </div>
        </div>
      )}

    </div>
  );
};

export default HybridDockingController;
