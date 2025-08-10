import React from 'react'
import { render, screen, fireEvent, waitFor } from '@testing-library/react'
import '@testing-library/jest-dom'
import Enhanced3DViewer from './Enhanced3DViewer'

// Mock 3Dmol.js
const mockViewer = {
  addModel: jest.fn(),
  setStyle: jest.fn(),
  addSurface: jest.fn(),
  addLabels: jest.fn(),
  zoomTo: jest.fn(),
  center: jest.fn(),
  render: jest.fn(),
  clear: jest.fn(),
  setBackgroundColor: jest.fn(),
  pngURI: jest.fn(() => 'data:image/png;base64,mock-image'),
  getModel: jest.fn(() => ({
    getCentroid: jest.fn(() => ({ x: 0, y: 0, z: 0 }))
  }))
}

const mockCreateViewer = jest.fn(() => mockViewer)

// Mock window.$3Dmol
Object.defineProperty(window, '$3Dmol', {
  value: {
    createViewer: mockCreateViewer,
    SurfaceType: {
      VDW: 'VDW'
    }
  },
  writable: true
})

// Mock fetch for API calls
global.fetch = jest.fn()

describe('Enhanced3DViewer', () => {
  const defaultProps = {
    ligandSmiles: 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C',
    receptorPdbId: '1CRN'
  }

  beforeEach(() => {
    jest.clearAllMocks()
    
    // Mock successful protein fetch
    fetch.mockImplementation((url) => {
      if (url.includes('rcsb.org')) {
        return Promise.resolve({
          ok: true,
          text: () => Promise.resolve('MOCK PDB DATA')
        })
      }
      
      // Mock ligand preparation API
      if (url.includes('/api/ligand/prepare')) {
        return Promise.resolve({
          ok: true,
          json: () => Promise.resolve({
            sdf_content: 'MOCK SDF DATA'
          })
        })
      }
      
      return Promise.reject(new Error('Unknown URL'))
    })
  })

  it('renders enhanced viewer with control panel', async () => {
    render(<Enhanced3DViewer {...defaultProps} />)
    
    // Check for control panel elements
    expect(screen.getByText('Visualization Controls')).toBeInTheDocument()
    expect(screen.getByText('Protein Style')).toBeInTheDocument()
    expect(screen.getByText('Ligand Style')).toBeInTheDocument()
    expect(screen.getByText('Color Scheme')).toBeInTheDocument()
    expect(screen.getByText('Background')).toBeInTheDocument()
    
    // Check for action buttons
    expect(screen.getByText('Reset View')).toBeInTheDocument()
    expect(screen.getByText('Save Image')).toBeInTheDocument()
  })

  it('initializes 3D viewer and loads molecules', async () => {
    render(<Enhanced3DViewer {...defaultProps} />)
    
    await waitFor(() => {
      expect(mockCreateViewer).toHaveBeenCalled()
    })
    
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalledTimes(2) // Protein + Ligand
    })
    
    expect(mockViewer.setStyle).toHaveBeenCalled()
    expect(mockViewer.render).toHaveBeenCalled()
  })

  it('changes protein visualization style', async () => {
    render(<Enhanced3DViewer {...defaultProps} />)
    
    // Wait for initial load
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalled()
    })
    
    // Change protein style
    const proteinStyleSelect = screen.getByDisplayValue('Cartoon')
    fireEvent.change(proteinStyleSelect, { target: { value: 'surface' } })
    
    await waitFor(() => {
      expect(mockViewer.setStyle).toHaveBeenCalledWith(
        { model: 0 },
        expect.objectContaining({
          surface: expect.any(Object)
        })
      )
    })
  })

  it('changes ligand visualization style', async () => {
    render(<Enhanced3DViewer {...defaultProps} />)
    
    // Wait for initial load
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalled()
    })
    
    // Change ligand style
    const ligandStyleSelect = screen.getByDisplayValue('Stick')
    fireEvent.change(ligandStyleSelect, { target: { value: 'sphere' } })
    
    await waitFor(() => {
      expect(mockViewer.setStyle).toHaveBeenCalledWith(
        { model: 1 },
        expect.objectContaining({
          sphere: expect.any(Object)
        })
      )
    })
  })

  it('changes background color', async () => {
    render(<Enhanced3DViewer {...defaultProps} />)
    
    // Wait for initial load
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalled()
    })
    
    // Find and click a background color button (white)
    const whiteBackgroundButton = screen.getByTitle('White')
    fireEvent.click(whiteBackgroundButton)
    
    await waitFor(() => {
      expect(mockViewer.setBackgroundColor).toHaveBeenCalledWith('#ffffff')
    })
  })

  it('toggles surface display', async () => {
    render(<Enhanced3DViewer {...defaultProps} />)
    
    // Wait for initial load
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalled()
    })
    
    // Toggle surface checkbox
    const surfaceCheckbox = screen.getByLabelText('Show Surface')
    fireEvent.click(surfaceCheckbox)
    
    await waitFor(() => {
      expect(mockViewer.addSurface).toHaveBeenCalled()
    })
  })

  it('toggles labels display', async () => {
    render(<Enhanced3DViewer {...defaultProps} />)
    
    // Wait for initial load
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalled()
    })
    
    // Toggle labels checkbox
    const labelsCheckbox = screen.getByLabelText('Show Labels')
    fireEvent.click(labelsCheckbox)
    
    await waitFor(() => {
      expect(mockViewer.addLabels).toHaveBeenCalled()
    })
  })

  it('resets view when reset button clicked', async () => {
    render(<Enhanced3DViewer {...defaultProps} />)
    
    // Wait for initial load
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalled()
    })
    
    // Click reset view button
    const resetButton = screen.getByText('Reset View')
    fireEvent.click(resetButton)
    
    expect(mockViewer.zoomTo).toHaveBeenCalled()
    expect(mockViewer.center).toHaveBeenCalled()
    expect(mockViewer.render).toHaveBeenCalled()
  })

  it('saves image when save button clicked', async () => {
    // Mock document.createElement and click
    const mockLink = {
      click: jest.fn(),
      download: '',
      href: ''
    }
    jest.spyOn(document, 'createElement').mockReturnValue(mockLink)
    
    render(<Enhanced3DViewer {...defaultProps} />)
    
    // Wait for initial load
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalled()
    })
    
    // Click save image button
    const saveButton = screen.getByText('Save Image')
    fireEvent.click(saveButton)
    
    expect(mockViewer.pngURI).toHaveBeenCalled()
    expect(mockLink.click).toHaveBeenCalled()
    expect(mockLink.download).toMatch(/molecular-docking-\d+\.png/)
  })

  it('handles errors gracefully', async () => {
    // Make fetch fail
    fetch.mockRejectedValue(new Error('Network error'))
    
    render(<Enhanced3DViewer {...defaultProps} />)
    
    await waitFor(() => {
      expect(screen.getByText('Visualization Error')).toBeInTheDocument()
    })
    
    expect(screen.getByText('Reload Page')).toBeInTheDocument()
  })

  it('displays loading state initially', () => {
    render(<Enhanced3DViewer {...defaultProps} />)
    
    expect(screen.getByText('Loading molecular structures...')).toBeInTheDocument()
  })

  it('displays status information', async () => {
    render(<Enhanced3DViewer {...defaultProps} />)
    
    // Check for status display
    expect(screen.getByText(/Protein: 1CRN/)).toBeInTheDocument()
    expect(screen.getByText(/Ligand: CC\(C\)CCCC/)).toBeInTheDocument()
  })

  it('handles missing 3Dmol library', () => {
    // Temporarily remove 3Dmol
    const original3Dmol = window.$3Dmol
    delete window.$3Dmol
    
    render(<Enhanced3DViewer {...defaultProps} />)
    
    expect(screen.getByText('3Dmol.js library not loaded. Please ensure the script is included.')).toBeInTheDocument()
    
    // Restore 3Dmol
    window.$3Dmol = original3Dmol
  })

  it('cleans up on unmount', async () => {
    const { unmount } = render(<Enhanced3DViewer {...defaultProps} />)
    
    // Wait for initial load
    await waitFor(() => {
      expect(mockViewer.addModel).toHaveBeenCalled()
    })
    
    // Unmount component
    unmount()
    
    expect(mockViewer.clear).toHaveBeenCalled()
  })
})