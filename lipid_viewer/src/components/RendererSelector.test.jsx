import React from 'react'
import { render, screen } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import RendererSelector from './RendererSelector'

describe('RendererSelector Component', () => {
  const mockOnModeChange = jest.fn()
  const mockOnRendererChange = jest.fn()

  beforeEach(() => {
    jest.clearAllMocks()
  })

  test('renders 2D/3D mode toggle buttons', () => {
    render(
      <RendererSelector 
        mode="2D" 
        renderer="SmilesDrawer"
        onModeChange={mockOnModeChange}
        onRendererChange={mockOnRendererChange}
      />
    )
    
    expect(screen.getByRole('button', { name: /2D View/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /3D View/i })).toBeInTheDocument()
  })

  test('shows 2D mode as active when in 2D mode', () => {
    render(
      <RendererSelector 
        mode="2D" 
        renderer="SmilesDrawer"
        onModeChange={mockOnModeChange}
        onRendererChange={mockOnRendererChange}
      />
    )
    
    const twoDButton = screen.getByRole('button', { name: /2D View/i })
    expect(twoDButton).toHaveClass('bg-blue-600')
  })

  test('shows 3D mode as active when in 3D mode', () => {
    render(
      <RendererSelector 
        mode="3D" 
        renderer="3Dmol.js"
        onModeChange={mockOnModeChange}
        onRendererChange={mockOnRendererChange}
      />
    )
    
    // Get the 3D View toggle button specifically (not renderer buttons)
    const threeDButton = screen.getByRole('button', { name: '3D View' })
    expect(threeDButton).toHaveClass('bg-blue-600')
  })

  test('calls onModeChange when switching to 3D', async () => {
    const user = userEvent.setup()
    render(
      <RendererSelector 
        mode="2D" 
        renderer="SmilesDrawer"
        onModeChange={mockOnModeChange}
        onRendererChange={mockOnRendererChange}
      />
    )
    
    const threeDButton = screen.getByRole('button', { name: /3D View/i })
    await user.click(threeDButton)
    
    expect(mockOnModeChange).toHaveBeenCalledWith('3D')
  })

  test('displays 2D renderer options when in 2D mode', () => {
    render(
      <RendererSelector 
        mode="2D" 
        renderer="SmilesDrawer"
        onModeChange={mockOnModeChange}
        onRendererChange={mockOnRendererChange}
      />
    )
    
    expect(screen.getByText('SmilesDrawer')).toBeInTheDocument()
    expect(screen.getByText('RDKit.js')).toBeInTheDocument()
    expect(screen.getByText('Kekule.js')).toBeInTheDocument()
  })

  test('displays 3D renderer options when in 3D mode', () => {
    render(
      <RendererSelector 
        mode="3D" 
        renderer="3Dmol.js"
        onModeChange={mockOnModeChange}
        onRendererChange={mockOnRendererChange}
      />
    )
    
    expect(screen.getByText('3Dmol.js')).toBeInTheDocument()
    expect(screen.getByText('Mol*')).toBeInTheDocument()
    expect(screen.getByText('NGL')).toBeInTheDocument()
  })

  test('calls onRendererChange when selecting different 2D renderer', async () => {
    const user = userEvent.setup()
    render(
      <RendererSelector 
        mode="2D" 
        renderer="SmilesDrawer"
        onModeChange={mockOnModeChange}
        onRendererChange={mockOnRendererChange}
      />
    )
    
    const rdkitButton = screen.getByText('RDKit.js')
    await user.click(rdkitButton)
    
    expect(mockOnRendererChange).toHaveBeenCalledWith('RDKit.js')
  })

  test('calls onRendererChange when selecting different 3D renderer', async () => {
    const user = userEvent.setup()
    render(
      <RendererSelector 
        mode="3D" 
        renderer="3Dmol.js"
        onModeChange={mockOnModeChange}
        onRendererChange={mockOnRendererChange}
      />
    )
    
    const molstarButton = screen.getByText('Mol*')
    await user.click(molstarButton)
    
    expect(mockOnRendererChange).toHaveBeenCalledWith('Mol*')
  })

  test('highlights active renderer', () => {
    render(
      <RendererSelector 
        mode="2D" 
        renderer="RDKit.js"
        onModeChange={mockOnModeChange}
        onRendererChange={mockOnRendererChange}
      />
    )
    
    // Find the button containing RDKit.js text
    const activeRenderer = screen.getByText('RDKit.js').closest('button')
    expect(activeRenderer).toHaveClass('bg-blue-100')
  })
}) 