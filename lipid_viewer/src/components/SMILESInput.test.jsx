import React from 'react'
import { render, screen, fireEvent, waitFor } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import SMILESInput from './SMILESInput'

// Mock the validation utility with more comprehensive behavior
jest.mock('../utils/smilesValidator', () => {
  const validateSMILES = jest.fn()
  return {
    __esModule: true,
    validateSMILES,
    getExampleMolecules: jest.fn(() => [
      { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' },
      { name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O' }
    ])
  }
})

import { validateSMILES as mockValidateSMILES } from '../utils/smilesValidator'

describe('SMILESInput Component', () => {
  const mockOnSubmit = jest.fn()
  const mockOnValidation = jest.fn()

  beforeEach(() => {
    jest.clearAllMocks()
    // Default mock behavior
    mockValidateSMILES.mockResolvedValue({ isValid: true, validationType: 'basic' })
  })

  test('renders input field and visualize button', () => {
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    expect(screen.getByPlaceholderText(/Enter SMILES string/i)).toBeInTheDocument()
    expect(screen.getByRole('button', { name: /visualize/i })).toBeInTheDocument()
  })

  test('renders example molecule buttons', () => {
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    expect(screen.getByText('Ethanol')).toBeInTheDocument()
    expect(screen.getByText('Benzene')).toBeInTheDocument()
    expect(screen.getByText('Caffeine')).toBeInTheDocument()
  })

  test('allows user to type in input field', async () => {
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    await user.type(input, 'CCO')
    
    expect(input).toHaveValue('CCO')
  })

  test('button is disabled when no input', () => {
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const button = screen.getByRole('button', { name: /visualize/i })
    expect(button).toBeDisabled()
  })

  test('example molecule buttons populate input', async () => {
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    const ethanolButton = screen.getByText('Ethanol')
    
    await user.click(ethanolButton)
    expect(input).toHaveValue('CCO')
  })

  test('form submission triggers onSubmit callback', async () => {
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    const button = screen.getByRole('button', { name: /visualize/i })
    
    await user.type(input, 'CCO')

    // wait for debounce + validation to enable button
    await waitFor(() => expect(button).not.toBeDisabled(), { timeout: 1500 })

    await user.click(button)
    
    expect(mockOnSubmit).toHaveBeenCalledWith('CCO')
  })

  test('displays SMILES information text', () => {
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    expect(screen.getAllByText(/SMILES/i).length).toBeGreaterThan(0)
    expect(screen.getByText(/Simplified Molecular Input Line Entry System/i)).toBeInTheDocument()
  })

  test('input field has proper styling classes', () => {
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    expect(input).toHaveClass('w-full', 'px-4', 'py-3')
  })

  test('button has proper styling when disabled', () => {
    render(<SMILESInput onSubmit={mockOnValidation} onValidation={mockOnValidation} />)
    
    const button = screen.getByRole('button', { name: /visualize/i })
    expect(button).toHaveClass('control-button')
    expect(button).toBeDisabled()
  })

  test('clears input when example button clicked multiple times', async () => {
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    const ethanolButton = screen.getByText('Ethanol')
    const benzeneButton = screen.getByText('Benzene')
    
    await user.click(ethanolButton)
    expect(input).toHaveValue('CCO')
    
    await user.click(benzeneButton)
    expect(input).toHaveValue('c1ccccc1')
  })

  test('input accepts complex SMILES strings', async () => {
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    const complexSmiles = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
    
    await user.type(input, complexSmiles)
    expect(input).toHaveValue(complexSmiles)
  })

  test('form has proper structure and classes', () => {
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const button = screen.getByRole('button', { name: /visualize/i })
    const form = button.closest('form')
    expect(form).toBeInTheDocument()
    expect(form).toHaveClass('space-y-4')
  })

  test('example section has proper heading', () => {
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    expect(screen.getByText('Try examples:')).toBeInTheDocument()
  })

  test('all example molecule buttons are present', () => {
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const expectedButtons = ['Ethanol', 'Benzene', 'Caffeine', 'Aspirin', 'Cholesterol']
    expectedButtons.forEach(buttonText => {
      expect(screen.getByText(buttonText)).toBeInTheDocument()
    })
  })

  test('input has correct placeholder text', () => {
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText('Enter SMILES string (e.g., CCO for ethanol)')
    expect(input).toBeInTheDocument()
  })

  test('component renders without crashing with minimal props', () => {
    render(<SMILESInput onSubmit={jest.fn()} />)
    expect(screen.getByPlaceholderText(/Enter SMILES string/i)).toBeInTheDocument()
  })

  test('handles validation with invalid SMILES', async () => {
    mockValidateSMILES.mockResolvedValue({ isValid: false, error: 'Invalid SMILES', validationType: 'basic' })
    
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    await user.type(input, 'invalid')
    
    // Allow time for validation
    await waitFor(() => {
      expect(mockValidateSMILES).toHaveBeenCalledWith('invalid')
    }, { timeout: 1500 })
  })

  test('handles validation errors gracefully', async () => {
    mockValidateSMILES.mockRejectedValue(new Error('Validation failed'))
    
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    await user.type(input, 'CCO')
    
    // Should not crash when validation fails
    expect(input).toBeInTheDocument()
  })

  test('calls onValidation callback when provided', async () => {
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    await user.type(input, 'CCO')
    
    // Allow time for validation debounce
    await waitFor(() => {
      expect(mockOnValidation).toHaveBeenCalled()
    }, { timeout: 1000 })
  })

  test('works without onValidation callback', async () => {
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    await user.type(input, 'CCO')
    
    expect(input).toHaveValue('CCO')
  })

  test('shows loading state during validation', async () => {
    // Make validation take longer
    mockValidateSMILES.mockImplementation(() => 
      new Promise(resolve => setTimeout(() => resolve({ isValid: true, validationType: 'basic' }), 200))
    )
    
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    await user.type(input, 'CCO')
    
    // Component should handle the loading state
    expect(input).toBeInTheDocument()
  })

  test('prevents form submission with empty input', () => {
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const button = screen.getByRole('button', { name: /visualize/i })
    expect(button).toBeDisabled()
    
    // Clicking disabled button should not trigger submission
    fireEvent.click(button)
    expect(mockOnSubmit).not.toHaveBeenCalled()
  })

  test('handles Enter key submission', async () => {
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    const button = screen.getByRole('button', { name: /visualize/i })

    await user.type(input, 'CCO')
    await waitFor(() => expect(button).not.toBeDisabled(), { timeout: 1500 })

    await user.keyboard('{Enter}')
    
    expect(mockOnSubmit).toHaveBeenCalledWith('CCO')
  })

  test('button becomes enabled with valid input', async () => {
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const input = screen.getByPlaceholderText(/Enter SMILES string/i)
    const button = screen.getByRole('button', { name: /visualize/i })
    
    await user.type(input, 'CCO')

    await waitFor(() => expect(button).not.toBeDisabled(), { timeout: 1500 })
  })

  test('example buttons trigger validation and submission', async () => {
    const user = userEvent.setup()
    render(<SMILESInput onSubmit={mockOnSubmit} onValidation={mockOnValidation} />)
    
    const ethanolButton = screen.getByText('Ethanol')
    await user.click(ethanolButton)

    // validation called and then submit after enabling
    await waitFor(() => expect(mockValidateSMILES).toHaveBeenCalled(), { timeout: 1500 })

    const button = screen.getByRole('button', { name: /visualize/i })
    await waitFor(() => expect(button).not.toBeDisabled(), { timeout: 1500 })
    await user.click(button)
    
    expect(mockOnSubmit).toHaveBeenCalledWith('CCO')
  })
})
