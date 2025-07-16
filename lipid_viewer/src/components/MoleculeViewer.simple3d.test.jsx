import React from 'react'
import { render, screen, waitFor } from '@testing-library/react'
import '@testing-library/jest-dom'
import MoleculeViewer from './MoleculeViewer'

// Mock all 3D libraries before any imports
const mock3DmolViewer = {
  clear: jest.fn(),
  addModel: jest.fn(),
  setStyle: jest.fn(),
  zoomTo: jest.fn(),
  center: jest.fn(),
  rotate: jest.fn(),
  render: jest.fn()
}

// Set up global mocks before tests run
beforeAll(() => {
  // Mock 3Dmol.js
  Object.defineProperty(window, '$3Dmol', {
    writable: true,
    value: {
      createViewer: jest.fn().mockReturnValue(mock3DmolViewer),
      rasmolElementColors: {}
    }
  })

  // Mock NGL
  Object.defineProperty(window, 'NGL', {
    writable: true,
    value: {
      Stage: jest.fn().mockImplementation(() => ({
        loadFile: jest.fn().mockResolvedValue({
          addRepresentation: jest.fn().mockReturnThis(),
          autoView: jest.fn()
        }),
        dispose: jest.fn()
      }))
    }
  })

  // Mock Molstar
  Object.defineProperty(window, 'molstar', {
    writable: true,
    value: {
      createPlugin: jest.fn().mockResolvedValue({
        builders: {
          data: { download: jest.fn().mockResolvedValue({}) },
          structure: {
            parseTrajectory: jest.fn().mockResolvedValue({}),
            createModel: jest.fn().mockResolvedValue({}),
            createStructure: jest.fn().mockResolvedValue({}),
            representation: { addRepresentation: jest.fn().mockResolvedValue({}) }
          }
        },
        dataTransaction: jest.fn((callback) => callback()),
        managers: {
          camera: { focusLoci: jest.fn() },
          structure: { hierarchy: { current: { structures: [] } } }
        },
        dispose: jest.fn()
      })
    }
  })

  // Mock DOM APIs
  global.Blob = jest.fn().mockImplementation((parts, properties) => ({
    size: parts[0]?.length || 0,
    type: properties?.type || 'application/octet-stream'
  }))

  global.File = jest.fn().mockImplementation((parts, filename, properties) => ({
    name: filename,
    size: parts[0]?.length || 0,
    type: properties?.type || 'application/octet-stream'
  }))

  Object.defineProperty(window.URL, 'createObjectURL', {
    writable: true,
    value: jest.fn().mockReturnValue('blob:test-url')
  })

  Object.defineProperty(window.URL, 'revokeObjectURL', {
    writable: true,
    value: jest.fn()
  })
})

describe('MoleculeViewer 3D Basic Tests', () => {
  beforeEach(() => {
    jest.clearAllMocks()
    // Reduce console noise
    jest.spyOn(console, 'log').mockImplementation(() => {})
    jest.spyOn(console, 'warn').mockImplementation(() => {})
    jest.spyOn(console, 'error').mockImplementation(() => {})
  })

  afterEach(() => {
    jest.restoreAllMocks()
  })

  test('should render 3Dmol.js viewer when libraries are available', async () => {
    render(
      <MoleculeViewer 
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C" 
        mode="3D" 
        renderer="3Dmol.js" 
      />
    )

    // Wait for component to render
    await waitFor(() => {
      expect(screen.getByTestId('molecule-viewer')).toBeInTheDocument()
    }, { timeout: 5000 })

    // Check that 3Dmol was called
    await waitFor(() => {
      expect(window.$3Dmol.createViewer).toHaveBeenCalled()
    }, { timeout: 3000 })

    expect(mock3DmolViewer.addModel).toHaveBeenCalled()
    expect(mock3DmolViewer.setStyle).toHaveBeenCalled()
    expect(mock3DmolViewer.render).toHaveBeenCalled()
  })

  test('should render NGL viewer when libraries are available', async () => {
    render(
      <MoleculeViewer 
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C" 
        mode="3D" 
        renderer="NGL" 
      />
    )

    // Wait for component to render
    await waitFor(() => {
      expect(screen.getByTestId('molecule-viewer')).toBeInTheDocument()
    }, { timeout: 5000 })

    // Check that NGL was called
    await waitFor(() => {
      expect(window.NGL.Stage).toHaveBeenCalled()
    }, { timeout: 3000 })
  })

  test('should render Molstar viewer when libraries are available', async () => {
    render(
      <MoleculeViewer 
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C" 
        mode="3D" 
        renderer="Mol*" 
      />
    )

    // Wait for component to render
    await waitFor(() => {
      expect(screen.getByTestId('molecule-viewer')).toBeInTheDocument()
    }, { timeout: 5000 })

    // Check that Molstar was called
    await waitFor(() => {
      expect(window.molstar.createPlugin).toHaveBeenCalled()
    }, { timeout: 3000 })
  })

  test('should handle empty SMILES gracefully', async () => {
    render(
      <MoleculeViewer 
        smiles="" 
        mode="3D" 
        renderer="3Dmol.js" 
      />
    )

    // Should show error message
    await waitFor(() => {
      expect(screen.getByText(/No SMILES string provided/)).toBeInTheDocument()
    })

    // Should not call 3D libraries
    expect(window.$3Dmol.createViewer).not.toHaveBeenCalled()
  })

  test('should cleanup 3D viewers on unmount', async () => {
    const { unmount } = render(
      <MoleculeViewer 
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C" 
        mode="3D" 
        renderer="3Dmol.js" 
      />
    )

    await waitFor(() => {
      expect(window.$3Dmol.createViewer).toHaveBeenCalled()
    })

    unmount()

    // Should call cleanup methods
    expect(mock3DmolViewer.clear).toHaveBeenCalled()
  })

  test('should switch between 3D renderers correctly', async () => {
    const { rerender } = render(
      <MoleculeViewer 
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C" 
        mode="3D" 
        renderer="3Dmol.js" 
      />
    )

    // Wait for 3Dmol to load
    await waitFor(() => {
      expect(window.$3Dmol.createViewer).toHaveBeenCalled()
    })

    // Switch to NGL
    rerender(
      <MoleculeViewer 
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C" 
        mode="3D" 
        renderer="NGL" 
      />
    )

    // Should cleanup 3Dmol and start NGL
    await waitFor(() => {
      expect(mock3DmolViewer.clear).toHaveBeenCalled()
      expect(window.NGL.Stage).toHaveBeenCalled()
    })
  })
}) 