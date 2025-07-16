// AWS Configuration for Lipid Rendering App
export const awsConfig = {
  // Environment detection
  isProduction: import.meta.env.PROD,
  isDevelopment: import.meta.env.DEV,
  
  // App configuration
  appName: import.meta.env.VITE_APP_NAME || 'Lipid Rendering',
  appVersion: import.meta.env.VITE_APP_VERSION || '1.0.0',
  
  // API configuration
  apiBaseUrl: import.meta.env.VITE_API_BASE_URL || 'http://localhost:3001',
  
  // AWS specific settings
  region: import.meta.env.VITE_AWS_REGION || 'us-east-1',
  
  // Performance settings for molecular rendering
  performance: {
    enableWebGL: true,
    maxMoleculeSize: 10000, // atoms
    enableLazyLoading: true,
    chunkSizeLimit: 2048, // KB
  },
  
  // Security settings
  security: {
    enableCSP: true,
    allowedDomains: [
      'cdnjs.cloudflare.com',
      'cdn.jsdelivr.net',
      'unpkg.com'
    ]
  }
};

// Export for use in components
export default awsConfig; 