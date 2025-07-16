import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  server: {
    port: 3000,
    open: true
  },
  build: {
    outDir: 'dist',
    assetsDir: 'assets',
    // Optimize for AWS Amplify
    sourcemap: false,
    rollupOptions: {
      external: [
        // These libraries are loaded via CDN, don't bundle them
        '3dmol',
        'ngl',
        'molstar',
        'kekule'
      ],
      output: {
        manualChunks: {
          vendor: ['react', 'react-dom', 'react-hot-toast'],
          utils: ['smiles-drawer']
        },
        globals: {
          // Map CDN libraries to their global variables
          '3dmol': '$3Dmol',
          'ngl': 'NGL',
          'molstar': 'molstar',
          'kekule': 'Kekule'
        }
      }
    }
  },
  optimizeDeps: {
    include: ['smiles-drawer'],
    exclude: [
      // Exclude CDN libraries from optimization
      '3dmol',
      'ngl', 
      'molstar',
      'kekule',
      'rdkit-js'
    ]
  },
  // For SPA routing
  preview: {
    port: 3000
  }
}) 