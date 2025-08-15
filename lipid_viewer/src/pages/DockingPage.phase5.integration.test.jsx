import React from 'react'
import { render, screen } from '@testing-library/react'
import { MemoryRouter } from 'react-router-dom'
import DockingPage from './DockingPage'

describe('DockingPage Phase 5 Integration', () => {
  test('renders new modular components', () => {
    render(
      <MemoryRouter>
        <DockingPage />
      </MemoryRouter>
    )

    // ProteinSelector
    expect(screen.getByText(/Protein Selector/i)).toBeInTheDocument()
    expect(screen.getByPlaceholderText(/e\.g\., 1CRN/i)).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /Load Protein/i })).toBeInTheDocument()

    // BindingSiteConfigurator
    expect(screen.getByText(/Binding Site Configurator/i)).toBeInTheDocument()
    expect(screen.getByLabelText(/Center X/i)).toBeInTheDocument()
    expect(screen.getByLabelText(/Size Z/i)).toBeInTheDocument()

    // DockingParameters
    expect(screen.getByText(/Docking Parameters/i)).toBeInTheDocument()
    expect(screen.getByLabelText(/Exhaustiveness/i)).toBeInTheDocument()

    // JobManager
    expect(screen.getByRole('button', { name: /Run Docking/i })).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /Export Results/i })).toBeInTheDocument()
  })
})



