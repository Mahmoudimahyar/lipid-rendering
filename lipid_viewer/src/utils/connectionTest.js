/**
 * Utility class for checking backend connection health
 */
export class ConnectionHealthChecker {
  constructor() {
    this.lastHealthCheck = null
    this.checkInterval = 30000 // 30 seconds
  }

  /**
   * Check if the backend server is responding
   * @returns {Promise<{isConnected: boolean, timestamp: number, error?: string}>}
   */
  async checkBackendHealth() {
    try {
      const response = await fetch('/api/healthz', {
        method: 'GET',
        headers: {
          'Content-Type': 'application/json',
        },
        timeout: 5000 // 5 second timeout
      })

      const isConnected = response.ok
      const result = {
        isConnected,
        timestamp: Date.now()
      }

      if (!isConnected) {
        result.error = `Server responded with status ${response.status}: ${response.statusText}`
      }

      this.lastHealthCheck = result
      return result

    } catch (error) {
      const result = {
        isConnected: false,
        timestamp: Date.now(),
        error: error.message || 'Connection failed'
      }

      this.lastHealthCheck = result
      return result
    }
  }

  /**
   * Get the last health check result
   * @returns {object|null}
   */
  getLastHealthCheck() {
    return this.lastHealthCheck
  }

  /**
   * Check if the last health check is still valid (not expired)
   * @returns {boolean}
   */
  isLastCheckValid() {
    if (!this.lastHealthCheck) return false
    const now = Date.now()
    return (now - this.lastHealthCheck.timestamp) < this.checkInterval
  }

  /**
   * Start periodic health checks
   * @param {function} callback - Callback function to handle health check results
   * @param {number} interval - Check interval in milliseconds (default: 30000)
   */
  startPeriodicChecks(callback, interval = this.checkInterval) {
    this.stopPeriodicChecks() // Stop any existing checks
    
    this.periodicCheckTimer = setInterval(async () => {
      const result = await this.checkBackendHealth()
      if (callback) callback(result)
    }, interval)
  }

  /**
   * Stop periodic health checks
   */
  stopPeriodicChecks() {
    if (this.periodicCheckTimer) {
      clearInterval(this.periodicCheckTimer)
      this.periodicCheckTimer = null
    }
  }
}

export default ConnectionHealthChecker
