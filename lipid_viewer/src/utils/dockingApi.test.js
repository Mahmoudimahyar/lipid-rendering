import { DockingAPI } from './dockingApi';

// Mock fetch globally
global.fetch = jest.fn();

describe('DockingAPI', () => {
  beforeEach(() => {
    fetch.mockClear();
  });

  describe('estimateBindingSite', () => {
    test('should estimate binding site successfully', async () => {
      const mockResponse = {
        center_x: 5.0,
        center_y: 5.0,
        center_z: 5.0,
        size_x: 20.0,
        size_y: 20.0,
        size_z: 20.0,
        confidence: 'medium'
      };

      fetch.mockResolvedValueOnce({
        ok: true,
        json: () => Promise.resolve(mockResponse)
      });

      const result = await DockingAPI.estimateBindingSite('1CRN', 'CCO');

      expect(fetch).toHaveBeenCalledWith(
        '/api/binding-site/estimate',
        {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            pdb_id: '1CRN',
            ligand_smiles: 'CCO'
          })
        }
      );

      expect(result).toEqual(mockResponse);
    });

    test('should handle binding site estimation error', async () => {
      fetch.mockResolvedValueOnce({
        ok: false,
        text: () => Promise.resolve('PDB not found')
      });

      await expect(DockingAPI.estimateBindingSite('INVALID'))
        .rejects.toThrow('Failed to estimate binding site: PDB not found');
    });
  });

  describe('getProteinInfo', () => {
    test('should get protein info successfully', async () => {
      const mockInfo = {
        title: 'Test Protein',
        molecular_weight: 10000
      };

      fetch.mockResolvedValueOnce({
        ok: true,
        json: () => Promise.resolve(mockInfo)
      });

      const result = await DockingAPI.getProteinInfo('1CRN');

      expect(fetch).toHaveBeenCalledWith('/api/pdb/1CRN/info');
      expect(result).toEqual(mockInfo);
    });

    test('should handle protein info error', async () => {
      fetch.mockResolvedValueOnce({
        ok: false,
        text: () => Promise.resolve('Not found')
      });

      await expect(DockingAPI.getProteinInfo('INVALID'))
        .rejects.toThrow('Failed to get protein info: Not found');
    });
  });

  describe('prepareLigand', () => {
    test('should prepare ligand successfully', async () => {
      const mockResult = {
        smiles: 'CCO',
        molecular_formula: 'C2H6O',
        molecular_weight: 46.07
      };

      fetch.mockResolvedValueOnce({
        ok: true,
        json: () => Promise.resolve(mockResult)
      });

      const result = await DockingAPI.prepareLigand('CCO');

      expect(fetch).toHaveBeenCalledWith(
        '/api/ligand/prepare',
        {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ smiles: 'CCO' })
        }
      );

      expect(result).toEqual(mockResult);
    });
  });

  describe('prepareReceptor', () => {
    test('should prepare receptor successfully', async () => {
      const mockResult = {
        pdb_id: '1CRN',
        protein_info: { title: 'Test Protein' }
      };

      fetch.mockResolvedValueOnce({
        ok: true,
        json: () => Promise.resolve(mockResult)
      });

      const result = await DockingAPI.prepareReceptor('1CRN');

      expect(fetch).toHaveBeenCalledWith(
        '/api/receptor/prepare',
        {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ pdb_id: '1CRN' })
        }
      );

      expect(result).toEqual(mockResult);
    });
  });

  describe('runDocking', () => {
    test('should start docking job successfully', async () => {
      const mockJob = {
        job_id: 'test-job-id',
        status: 'pending',
        message: 'Docking job started'
      };

      const dockingParams = {
        ligand_smiles: 'CCO',
        receptor_pdb_id: '1CRN',
        center_x: 5.0,
        center_y: 5.0,
        center_z: 5.0,
        size_x: 20.0,
        size_y: 20.0,
        size_z: 20.0
      };

      fetch.mockResolvedValueOnce({
        ok: true,
        json: () => Promise.resolve(mockJob)
      });

      const result = await DockingAPI.runDocking(dockingParams);

      expect(fetch).toHaveBeenCalledWith(
        '/api/dock/run',
        {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(dockingParams)
        }
      );

      expect(result).toEqual(mockJob);
    });

    test('should handle docking start error', async () => {
      const dockingParams = {
        ligand_smiles: 'INVALID',
        receptor_pdb_id: '1CRN'
      };

      fetch.mockResolvedValueOnce({
        ok: false,
        text: () => Promise.resolve('Invalid parameters')
      });

      await expect(DockingAPI.runDocking(dockingParams))
        .rejects.toThrow('Failed to start docking: Invalid parameters');
    });
  });

  describe('getDockingStatus', () => {
    test('should get docking status successfully', async () => {
      const mockStatus = {
        job_id: 'test-job-id',
        status: 'completed',
        results: { poses: [] }
      };

      fetch.mockResolvedValueOnce({
        ok: true,
        json: () => Promise.resolve(mockStatus)
      });

      const result = await DockingAPI.getDockingStatus('test-job-id');

      expect(fetch).toHaveBeenCalledWith('/api/dock/status/test-job-id');
      expect(result).toEqual(mockStatus);
    });
  });

  describe('listDockingJobs', () => {
    test('should list docking jobs successfully', async () => {
      const mockJobs = {
        jobs: [
          { job_id: 'job1', status: 'completed' },
          { job_id: 'job2', status: 'running' }
        ],
        total_count: 2
      };

      fetch.mockResolvedValueOnce({
        ok: true,
        json: () => Promise.resolve(mockJobs)
      });

      const result = await DockingAPI.listDockingJobs();

      expect(fetch).toHaveBeenCalledWith('/api/dock/jobs');
      expect(result).toEqual(mockJobs);
    });
  });

  describe('pollDockingJob', () => {
    beforeEach(() => {
      jest.useFakeTimers();
    });

    afterEach(() => {
      jest.useRealTimers();
    });

    test('should poll until completion', async () => {
      const mockProgressStatus = {
        job_id: 'test-job-id',
        status: 'running'
      };

      const mockCompletedStatus = {
        job_id: 'test-job-id',
        status: 'completed',
        results: { poses: [] }
      };

      // Mock getDockingStatus to return running, then completed
      jest.spyOn(DockingAPI, 'getDockingStatus')
        .mockResolvedValueOnce(mockProgressStatus)
        .mockResolvedValueOnce(mockCompletedStatus);

      const onProgress = jest.fn();
      
      // Start polling
      const pollPromise = DockingAPI.pollDockingJob('test-job-id', onProgress, 1000);

      // Advance timer to trigger first poll
      jest.advanceTimersByTime(1000);
      await Promise.resolve(); // Let the first poll complete

      // Advance timer to trigger second poll
      jest.advanceTimersByTime(1000);
      await Promise.resolve(); // Let the second poll complete

      const result = await pollPromise;

      expect(onProgress).toHaveBeenCalledWith(mockProgressStatus);
      expect(result).toEqual(mockCompletedStatus);
    });

    test('should handle failed job', async () => {
      const mockFailedStatus = {
        job_id: 'test-job-id',
        status: 'failed',
        error_message: 'Docking failed'
      };

      jest.spyOn(DockingAPI, 'getDockingStatus')
        .mockResolvedValueOnce(mockFailedStatus);

      const pollPromise = DockingAPI.pollDockingJob('test-job-id', null, 1000);

      // Advance timer to trigger poll
      jest.advanceTimersByTime(1000);
      await Promise.resolve();

      await expect(pollPromise).rejects.toThrow('Docking failed');
    });
  });
});
