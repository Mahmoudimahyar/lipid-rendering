import React from 'react'
import { render, screen, fireEvent } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import ViewControls from './ViewControls'

describe('ViewControls Component', () => {
  const defaultProps = {
    mode: '2D',
    onZoomIn: jest.fn(),
    onZoomOut: jest.fn(),
    onResetView: jest.fn(),
    onExportPNG: jest.fn(),
    onExportSVG: jest.fn(),
    onExportGLTF: jest.fn(),
    isLoaded: true,
  }

  beforeEach(() => {
    jest.clearAllMocks()
  })

  test('renders all control buttons when molecule is loaded', () => {
    render(<ViewControls {...defaultProps} />)
    
    expect(screen.getByRole('button', { name: /Reset View/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /Zoom In/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /Zoom Out/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /Export PNG/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /Export SVG/i })).toBeInTheDocument()
  })

  test('does not render controls when molecule is not loaded', () => {
    render(<ViewControls {...defaultProps} isLoaded={false} />)
    
    expect(screen.queryByRole('button', { name: /Reset View/i })).not.toBeInTheDocument()
  })

  test('calls zoom in function when zoom in button is clicked', async () => {
    const user = userEvent.setup()
    render(<ViewControls {...defaultProps} />)
    
    const zoomInButton = screen.getByRole('button', { name: /Zoom In/i })
    await user.click(zoomInButton)
    
    expect(defaultProps.onZoomIn).toHaveBeenCalledTimes(1)
  })

  test('calls zoom out function when zoom out button is clicked', async () => {
    const user = userEvent.setup()
    render(<ViewControls {...defaultProps} />)
    
    const zoomOutButton = screen.getByRole('button', { name: /Zoom Out/i })
    await user.click(zoomOutButton)
    
    expect(defaultProps.onZoomOut).toHaveBeenCalledTimes(1)
  })

  test('calls reset view function when reset button is clicked', async () => {
    const user = userEvent.setup()
    render(<ViewControls {...defaultProps} />)
    
    const resetButton = screen.getByRole('button', { name: /Reset View/i })
    await user.click(resetButton)
    
    expect(defaultProps.onResetView).toHaveBeenCalledTimes(1)
  })

  test('calls PNG export function when PNG button is clicked', async () => {
    const user = userEvent.setup()
    render(<ViewControls {...defaultProps} />)
    
    const pngButton = screen.getByRole('button', { name: /Export PNG/i })
    await user.click(pngButton)
    
    expect(defaultProps.onExportPNG).toHaveBeenCalledTimes(1)
  })

  test('calls SVG export function when SVG button is clicked', async () => {
    const user = userEvent.setup()
    render(<ViewControls {...defaultProps} />)
    
    const svgButton = screen.getByRole('button', { name: /Export SVG/i })
    await user.click(svgButton)
    
    expect(defaultProps.onExportSVG).toHaveBeenCalledTimes(1)
  })

  test('shows GLTF export button only in 3D mode', () => {
    const { rerender } = render(<ViewControls {...defaultProps} mode="2D" />)
    
    expect(screen.queryByRole('button', { name: /Export GLTF/i })).not.toBeInTheDocument()
    
    rerender(<ViewControls {...defaultProps} mode="3D" />)
    
    expect(screen.getByRole('button', { name: /Export GLTF/i })).toBeInTheDocument()
  })

  test('calls GLTF export function when GLTF button is clicked in 3D mode', async () => {
    const user = userEvent.setup()
    render(<ViewControls {...defaultProps} mode="3D" />)
    
    const gltfButton = screen.getByRole('button', { name: /Export GLTF/i })
    await user.click(gltfButton)
    
    expect(defaultProps.onExportGLTF).toHaveBeenCalledTimes(1)
  })

  test('displays keyboard shortcuts hint', () => {
    render(<ViewControls {...defaultProps} />)
    
    expect(screen.getByText(/keyboard shortcuts/i)).toBeInTheDocument()
  })

  test('handles multiple rapid clicks gracefully', async () => {
    const user = userEvent.setup()
    render(<ViewControls {...defaultProps} />)
    
    const zoomInButton = screen.getByRole('button', { name: /Zoom In/i })
    
    // Rapid clicks
    await user.click(zoomInButton)
    await user.click(zoomInButton)
    await user.click(zoomInButton)
    
    expect(defaultProps.onZoomIn).toHaveBeenCalledTimes(3)
  })

  test('has proper accessibility with titles', () => {
    render(<ViewControls {...defaultProps} />)
    
    const zoomInButton = screen.getByRole('button', { name: /Zoom In/i })
    expect(zoomInButton).toHaveAttribute('title')
  })

  test('displays correct title and sections', () => {
    render(<ViewControls {...defaultProps} />)
    
    expect(screen.getByText('View Controls')).toBeInTheDocument()
    expect(screen.getByText('Navigation')).toBeInTheDocument()
    expect(screen.getByText('Export')).toBeInTheDocument()
    expect(screen.getByText('Keyboard Shortcuts')).toBeInTheDocument()
  })

  test('shows consistent export options for different modes', () => {
    const { rerender } = render(<ViewControls {...defaultProps} mode="2D" />)
    
    // Check initial 2D mode exports
    expect(screen.getByRole('button', { name: /Export PNG/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /Export SVG/i })).toBeInTheDocument()
    
    // Switch to 3D mode
    rerender(<ViewControls {...defaultProps} mode="3D" />)
    expect(screen.getByRole('button', { name: /Export PNG/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /Export GLTF/i })).toBeInTheDocument()
  })
}) 