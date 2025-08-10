/**
 * API utilities for molecular docking backend integration
 */

// Use relative URLs since frontend and backend are served from the same domain
const API_BASE_URL = '/api';

export class DockingAPI {
  /**
   * Estimate binding site for a given PDB ID
   */
  static async estimateBindingSite(pdbId, ligandSmiles = null) {
    const response = await fetch(`${API_BASE_URL}/binding-site/estimate`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        pdb_id: pdbId,
        ligand_smiles: ligandSmiles
      })
    });

    if (!response.ok) {
      const error = await response.text();
      throw new Error(`Failed to estimate binding site: ${error}`);
    }

    return response.json();
  }

  /**
   * Get protein information for a PDB ID
   */
  static async getProteinInfo(pdbId) {
    const response = await fetch(`${API_BASE_URL}/pdb/${pdbId}/info`);

    if (!response.ok) {
      const error = await response.text();
      throw new Error(`Failed to get protein info: ${error}`);
    }

    return response.json();
  }

  /**
   * Prepare ligand for docking
   */
  static async prepareLigand(smiles) {
    const response = await fetch(`${API_BASE_URL}/ligand/prepare`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        smiles: smiles
      })
    });

    if (!response.ok) {
      const error = await response.text();
      throw new Error(`Failed to prepare ligand: ${error}`);
    }

    return response.json();
  }

  /**
   * Prepare receptor for docking
   */
  static async prepareReceptor(pdbId) {
    const response = await fetch(`${API_BASE_URL}/receptor/prepare`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        pdb_id: pdbId
      })
    });

    if (!response.ok) {
      const error = await response.text();
      throw new Error(`Failed to prepare receptor: ${error}`);
    }

    return response.json();
  }

  /**
   * Start docking job
   */
  static async runDocking(params) {
    const response = await fetch(`${API_BASE_URL}/dock/run`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(params)
    });

    if (!response.ok) {
      const error = await response.text();
      throw new Error(`Failed to start docking: ${error}`);
    }

    return response.json();
  }

  /**
   * Get docking job status
   */
  static async getDockingStatus(jobId) {
    const response = await fetch(`${API_BASE_URL}/dock/status/${jobId}`);

    if (!response.ok) {
      const error = await response.text();
      throw new Error(`Failed to get docking status: ${error}`);
    }

    return response.json();
  }

  /**
   * List recent docking jobs
   */
  static async listDockingJobs() {
    const response = await fetch(`${API_BASE_URL}/dock/jobs`);

    if (!response.ok) {
      const error = await response.text();
      throw new Error(`Failed to list docking jobs: ${error}`);
    }

    return response.json();
  }

  /**
   * Poll docking job until completion
   */
  static async pollDockingJob(jobId, onProgress = null, pollInterval = 2000) {
    return new Promise((resolve, reject) => {
      const pollTimer = setInterval(async () => {
        try {
          const status = await this.getDockingStatus(jobId);
          
          if (onProgress) {
            onProgress(status);
          }

          if (status.status === 'completed') {
            clearInterval(pollTimer);
            resolve(status);
          } else if (status.status === 'failed' || status.status === 'cancelled') {
            clearInterval(pollTimer);
            reject(new Error(status.error_message || `Docking job ${status.status}`));
          }
          // Continue polling if status is 'pending' or 'running'
        } catch (error) {
          clearInterval(pollTimer);
          reject(error);
        }
      }, pollInterval);

      // Set a timeout to avoid infinite polling
      setTimeout(() => {
        clearInterval(pollTimer);
        reject(new Error('Docking job timeout'));
      }, 300000); // 5 minutes timeout
    });
  }
}
