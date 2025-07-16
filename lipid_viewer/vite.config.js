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
      output: {
        manualChunks: {
          vendor: ['react', 'react-dom'],
          molecular: ['3dmol', 'ngl', 'molstar', 'kekule', 'smiles-drawer']
        }
      }
    }
  },
  optimizeDeps: {
    include: ['rdkit-js']
  },
  // For SPA routing
  preview: {
    port: 3000
  }
}) 