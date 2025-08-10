import React from 'react'
import { render, screen } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { MemoryRouter } from 'react-router-dom'
import DockingPage from './DockingPage'

describe('DockingPage', () => {
  test('renders docking layout and navigation', () => {
    render(
      <MemoryRouter>
        <DockingPage />
      </MemoryRouter>
    )

    expect(screen.getByText('Docking')).toBeInTheDocument()
    expect(screen.getByText(/Protein-Ligand docking/i)).toBeInTheDocument()
    expect(screen.getByRole('link', { name: /back to lipid viewer/i })).toBeInTheDocument()
    expect(screen.getByText(/Protein Loader/i)).toBeInTheDocument()
    expect(screen.getByText(/Docking Controls/i)).toBeInTheDocument()
    expect(screen.getByText(/Results/i)).toBeInTheDocument()
  })

  test('handles PDB ID input and loading', async () => {
    const user = userEvent.setup()
    render(
      <MemoryRouter>
        <DockingPage />
      </MemoryRouter>
    )

    const pdbInput = screen.getByPlaceholderText(/e.g., 1CRN/i)
    const loadButton = screen.getByRole('button', { name: /load/i })

    // Should have default value
    expect(pdbInput).toHaveValue('1CRN')

    // Should be able to change PDB ID
    await user.clear(pdbInput)
    await user.type(pdbInput, '2ABC')
    expect(pdbInput).toHaveValue('2ABC')

    // Load button should be clickable
    await user.click(loadButton)
    expect(loadButton).toBeInTheDocument()
  })

  test('handles docking controls input', async () => {
    const user = userEvent.setup()
    render(
      <MemoryRouter>
        <DockingPage />
      </MemoryRouter>
    )

    // Check all docking control inputs
    const centerX = screen.getByDisplayValue('0')
    const centerY = screen.getAllByDisplayValue('0')[1] 
    const centerZ = screen.getAllByDisplayValue('0')[2]
    const sizeX = screen.getByDisplayValue('20')
    const sizeY = screen.getAllByDisplayValue('20')[1]
    const sizeZ = screen.getAllByDisplayValue('20')[2]

    // Should be able to change values
    await user.clear(centerX)
    await user.type(centerX, '5')
    expect(centerX).toHaveValue('5')

    await user.clear(sizeX)
    await user.type(sizeX, '25')
    expect(sizeX).toHaveValue('25')
  })

  test('handles ligand input and visualization', async () => {
    const user = userEvent.setup()
    render(
      <MemoryRouter>
        <DockingPage />
      </MemoryRouter>
    )

    // Should have SMILES input
    const smilesInput = screen.getByPlaceholderText(/enter smiles string/i)
    expect(smilesInput).toBeInTheDocument()

    // Should be able to enter SMILES
    await user.type(smilesInput, 'CCO')
    expect(smilesInput).toHaveValue('CCO')

    // Should have visualization button
    const visualizeButton = screen.getByRole('button', { name: /visualize molecule/i })
    await user.click(visualizeButton)
  })

  test('handles mode and renderer selection', async () => {
    const user = userEvent.setup()
    render(
      <MemoryRouter>
        <DockingPage />
      </MemoryRouter>
    )

    // Should default to 3D mode
    expect(screen.getByText('3D Renderer')).toBeInTheDocument()

    // Should have 3D renderers
    expect(screen.getByRole('button', { name: /3dmol\.js/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /mol\*/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /ngl/i })).toBeInTheDocument()

    // Should be able to switch to 2D
    const twoDButton = screen.getByRole('button', { name: /2d view/i })
    await user.click(twoDButton)
    expect(screen.getByText('2D Renderer')).toBeInTheDocument()
  })

  test('run docking button is functional', async () => {
    const user = userEvent.setup()
    render(
      <MemoryRouter>
        <DockingPage />
      </MemoryRouter>
    )

    const runDockingButton = screen.getByRole('button', { name: /run docking/i })
    expect(runDockingButton).toBeInTheDocument()

    await user.click(runDockingButton)
    // Button should still be present after click
    expect(runDockingButton).toBeInTheDocument()
  })

  test('results section shows no results initially', () => {
    render(
      <MemoryRouter>
        <DockingPage />
      </MemoryRouter>
    )

    expect(screen.getByText('No results yet.')).toBeInTheDocument()
  })

  test('displays backend API endpoint information', () => {
    render(
      <MemoryRouter>
        <DockingPage />
      </MemoryRouter>
    )

    expect(screen.getByText(/Backend: \/api\/pdb\//)).toBeInTheDocument()
  })
})


