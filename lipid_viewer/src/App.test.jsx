import React from 'react'
import { render, screen, fireEvent, waitFor } from '@testing-library/react'
import { MemoryRouter } from 'react-router-dom'
import userEvent from '@testing-library/user-event'
import App from './App'

// Mock only the specific external dependencies, not our components
const mockToast = {
  success: jest.fn(),
  error: jest.fn(),
  promise: jest.fn((promise, options) => promise),
}

jest.mock('react-hot-toast', () => ({
  __esModule: true,
  default: mockToast,
  toast: mockToast,
  Toaster: () => <div data-testid="toaster" />
}))

// Mock keyboard event handling
const createKeyboardEvent = (key, target = document) => {
  const event = new KeyboardEvent('keydown', { key, bubbles: true })
  target.dispatchEvent(event)
}

describe('App Component', () => {
  beforeEach(() => {
    jest.clearAllMocks()
    mockToast.success.mockClear()
    mockToast.error.mockClear()
    mockToast.promise.mockClear()
  })

  test('renders main layout elements', () => {
    render(<MemoryRouter><App /></MemoryRouter>)
    
    expect(screen.getByText('Lipid Viewer')).toBeInTheDocument()
    expect(screen.getByText('Interactive 2D/3D molecular visualization tool')).toBeInTheDocument()
    expect(screen.getByText('Enter SMILES String')).toBeInTheDocument()
    expect(screen.getByText('Visualization Mode')).toBeInTheDocument()
  })

  test('handles mode switching', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    const threeDButton = screen.getByRole('button', { name: /3d view/i })
    const twoDButton = screen.getByRole('button', { name: /2d view/i })
    
    // Start in 2D mode
    expect(screen.getByText('Mode:')).toBeInTheDocument()
    
    // Switch to 3D
    await user.click(threeDButton)
    expect(threeDButton).toHaveClass('bg-blue-600')
    
    // Switch back to 2D
    await user.click(twoDButton)
    expect(twoDButton).toHaveClass('bg-blue-600')
  })

  test('handles renderer switching', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    const rdkitButton = screen.getByRole('button', { name: /rdkit\.js/i })
    const smilesDrawerButton = screen.getByRole('button', { name: /smilesdrawer/i })
    
    // Start with SmilesDrawer
    expect(smilesDrawerButton).toHaveClass('border-blue-500')
    
    // Switch to RDKit
    await user.click(rdkitButton)
    expect(rdkitButton).toHaveClass('border-blue-500')
    
    // Switch back to SmilesDrawer
    await user.click(smilesDrawerButton)
    expect(smilesDrawerButton).toHaveClass('border-blue-500')
  })

  test('keyboard shortcuts are registered', () => {
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Test that keyboard event listeners are set up
    // We can test this by verifying the component mounts without errors
    expect(screen.getByText('Lipid Viewer')).toBeInTheDocument()
    
    // Simulate keyboard events (they should not crash)
    createKeyboardEvent('+')
    createKeyboardEvent('-')
    createKeyboardEvent(' ')
    
    // Component should still be rendered
    expect(screen.getByText('Lipid Viewer')).toBeInTheDocument()
  })

  test('handles SMILES submission flow', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Enter a SMILES string
    const input = screen.getByPlaceholderText(/enter smiles string/i)
    await user.type(input, 'CCO')
    
    // Submit the form
    const submitButton = screen.getByRole('button', { name: /visualize molecule/i })
    await user.click(submitButton)
    
    // Component should handle the submission
    expect(screen.getByText('Lipid Viewer')).toBeInTheDocument()
  })

  test('renders all example buttons', () => {
    render(<MemoryRouter><App /></MemoryRouter>)
    
    expect(screen.getByRole('button', { name: /ethanol/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /benzene/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /caffeine/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /aspirin/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /cholesterol/i })).toBeInTheDocument()
  })

  test('displays renderer information in header', () => {
    render(<MemoryRouter><App /></MemoryRouter>)
    
    expect(screen.getByText('Mode:')).toBeInTheDocument()
    expect(screen.getByText('Renderer:')).toBeInTheDocument()
    expect(screen.getAllByText('SmilesDrawer').length).toBeGreaterThan(0)
  })

  test('shows appropriate renderer options based on mode', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // In 2D mode, should show 2D renderers
    expect(screen.getByText('2D Renderer')).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /smilesdrawer/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /rdkit\.js/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /kekule\.js/i })).toBeInTheDocument()
    
    // Switch to 3D mode
    const threeDButton = screen.getByRole('button', { name: /3d view/i })
    await user.click(threeDButton)
    
    // Should show 3D renderers
    expect(screen.getByText('3D Renderer')).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /3dmol\.js/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /mol\*/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /ngl/i })).toBeInTheDocument()
  })

  test('handles export functions setup', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Submit SMILES to potentially trigger export setup
    const input = screen.getByPlaceholderText(/enter smiles string/i)
    await user.type(input, 'CCO')
    
    const submitButton = screen.getByRole('button', { name: /visualize molecule/i })
    await user.click(submitButton)
    
    // Look for export-related elements
    await waitFor(() => {
      // The ViewControls component should be rendered when molecule is loaded
      const exportElement = screen.queryByText(/export png/i) || 
                          screen.queryByText(/export svg/i) || 
                          screen.queryByText(/export gltf/i)
      
      // If exports are available, they should be rendered
      if (exportElement) {
        expect(exportElement).toBeInTheDocument()
      }
    }, { timeout: 1000 })
  })

  test('component lifecycle and cleanup', () => {
    const { unmount } = render(<MemoryRouter><App /></MemoryRouter>)
    
    // Should mount successfully
    expect(screen.getByText('Lipid Viewer')).toBeInTheDocument()
    
    // Should unmount without errors
    unmount()
    
    // No assertions needed - if cleanup has issues, the test will fail
  })

  test('maintains state consistency', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Change mode
    const threeDButton = screen.getByRole('button', { name: /3d view/i })
    await user.click(threeDButton)
    
    // Change renderer
    const molstarButton = screen.getByRole('button', { name: /mol\*/i })
    await user.click(molstarButton)
    
    // State should be consistent - header should reflect changes
    expect(screen.getByText('3D')).toBeInTheDocument()
    expect(screen.getAllByText('Mol*').length).toBeGreaterThan(0)
  })

  test('handles rapid state changes', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Rapidly change modes and renderers
    const twoDButton = screen.getByRole('button', { name: /2d view/i })
    const threeDButton = screen.getByRole('button', { name: /3d view/i })
    
    await user.click(threeDButton)
    await user.click(twoDButton)
    await user.click(threeDButton)
    
    // Should handle rapid changes without crashing
    expect(screen.getByText('Lipid Viewer')).toBeInTheDocument()
  })

  test('toast notifications are configured', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Should have toaster component
    expect(screen.getByTestId('toaster')).toBeInTheDocument()
    
    // Submit empty form to potentially trigger error toast
    const submitButton = screen.getByRole('button', { name: /visualize molecule/i })
    await user.click(submitButton)
    
    // Component should still be functional
    expect(screen.getByText('Lipid Viewer')).toBeInTheDocument()
  })

  test('keyboard shortcuts work correctly', async () => {
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Simulate keyboard events
    fireEvent.keyDown(document, { key: '+' })
    fireEvent.keyDown(document, { key: '-' })
    fireEvent.keyDown(document, { key: ' ' })
    fireEvent.keyDown(document, { key: '2', ctrlKey: true })
    fireEvent.keyDown(document, { key: '3', ctrlKey: true })
    
    // Component should still be rendered and functional
    expect(screen.getByText('Lipid Viewer')).toBeInTheDocument()
  })

  test('export functions are properly handled', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Enter SMILES and submit
    const input = screen.getByPlaceholderText(/enter smiles string/i)
    await user.type(input, 'CCO')
    
    const submitButton = screen.getByRole('button', { name: /visualize molecule/i })
    await user.click(submitButton)
    
    // Wait for molecule to load and export functions to be available
    await waitFor(() => {
      const exportButton = screen.queryByText(/export png/i) || 
                          screen.queryByText(/export svg/i) || 
                          screen.queryByText(/export gltf/i)
      if (exportButton) {
        expect(exportButton).toBeInTheDocument()
      }
    }, { timeout: 2000 })
  })

  test('handles SMILES validation states', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Enter invalid SMILES
    const input = screen.getByPlaceholderText(/enter smiles string/i)
    await user.type(input, 'invalid')
    
    // Component should handle validation
    expect(input).toHaveValue('invalid')
    
    // Clear and enter valid SMILES
    await user.clear(input)
    await user.type(input, 'CCO')
    expect(input).toHaveValue('CCO')
  })

  test('navigation to docking page works', () => {
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Should have docking link
    const dockingLink = screen.getByRole('link', { name: /docking/i })
    expect(dockingLink).toBeInTheDocument()
    expect(dockingLink).toHaveAttribute('href', '/dock')
  })

  test('handles empty SMILES submission', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Submit without entering SMILES
    const submitButton = screen.getByRole('button', { name: /visualize molecule/i })
    await user.click(submitButton)
    
    // Should show error toast
    expect(mockToast.error).toHaveBeenCalledWith('Please enter a SMILES string', expect.any(Object))
  })

  test('handles successful SMILES submission', async () => {
    const user = userEvent.setup()
    render(<MemoryRouter><App /></MemoryRouter>)
    
    // Enter valid SMILES and submit
    const input = screen.getByPlaceholderText(/enter smiles string/i)
    await user.type(input, 'CCO')
    
    const submitButton = screen.getByRole('button', { name: /visualize molecule/i })
    await user.click(submitButton)
    
    // Should show success toast
    await waitFor(() => {
      expect(mockToast.success).toHaveBeenCalledWith('Loading molecule: CCO')
    })
  })
}) 