import { test, expect } from '@playwright/test'

test.describe('Lipid Viewer E2E Tests', () => {
  test.beforeEach(async ({ page }) => {
    await page.goto('http://localhost:3000')
  })

  test('loads application with correct title and interface', async ({ page }) => {
    await expect(page).toHaveTitle(/Lipid Viewer/)
    await expect(page.getByText('Lipid Viewer')).toBeVisible()
    await expect(page.getByPlaceholder(/Enter SMILES string/i)).toBeVisible()
    await expect(page.getByRole('button', { name: /visualize/i })).toBeVisible()
  })

  test('validates SMILES input and shows error for invalid input', async ({ page }) => {
    await page.getByPlaceholder(/Enter SMILES string/i).fill('invalid_smiles')
    await page.getByRole('button', { name: /visualize/i }).click()
    
    await expect(page.getByText(/Invalid SMILES string/i)).toBeVisible()
  })

  test('successfully visualizes valid SMILES molecule', async ({ page }) => {
    await page.getByPlaceholder(/Enter SMILES string/i).fill('CCO')
    await page.getByRole('button', { name: /visualize/i }).click()
    
    await expect(page.getByTestId('molecule-viewer')).toBeVisible()
    await expect(page.getByRole('button', { name: /Export PNG/i })).toBeVisible()
  })

  test('switches between 2D and 3D visualization modes', async ({ page }) => {
    // Start with 2D mode
    await expect(page.getByRole('button', { name: /2D View/i })).toHaveClass(/bg-blue-600/)
    
    // Switch to 3D mode
    await page.getByRole('button', { name: /3D View/i }).click()
    await expect(page.getByRole('button', { name: /3D View/i })).toHaveClass(/bg-blue-600/)
    
    // Check 3D renderer options appear
    await expect(page.getByText('3Dmol.js')).toBeVisible()
    await expect(page.getByText('Mol*')).toBeVisible()
    await expect(page.getByText('NGL')).toBeVisible()
  })

  test('changes 2D renderers successfully', async ({ page }) => {
    // Ensure we're in 2D mode
    await page.getByRole('button', { name: /2D View/i }).click()
    
    // Load a molecule first
    await page.getByPlaceholder(/Enter SMILES string/i).fill('CCO')
    await page.getByRole('button', { name: /visualize/i }).click()
    await expect(page.getByTestId('molecule-viewer')).toBeVisible()
    
    // Switch between 2D renderers
    await page.getByText('RDKit.js').click()
    await expect(page.getByText('RDKit.js')).toHaveClass(/bg-blue-100/)
    
    await page.getByText('Kekule.js').click()
    await expect(page.getByText('Kekule.js')).toHaveClass(/bg-blue-100/)
  })

  test('changes 3D renderers successfully', async ({ page }) => {
    // Switch to 3D mode
    await page.getByRole('button', { name: /3D View/i }).click()
    
    // Load a molecule
    await page.getByPlaceholder(/Enter SMILES string/i).fill('CCO')
    await page.getByRole('button', { name: /visualize/i }).click()
    await expect(page.getByTestId('molecule-viewer')).toBeVisible()
    
    // Switch between 3D renderers
    await page.getByText('Mol*').click()
    await expect(page.getByText('Mol*')).toHaveClass(/bg-blue-100/)
    
    await page.getByText('NGL').click()
    await expect(page.getByText('NGL')).toHaveClass(/bg-blue-100/)
  })

  test('uses view controls for zoom and reset', async ({ page }) => {
    // Load a molecule
    await page.getByPlaceholder(/Enter SMILES string/i).fill('CCO')
    await page.getByRole('button', { name: /visualize/i }).click()
    await expect(page.getByTestId('molecule-viewer')).toBeVisible()
    
    // Test view controls
    await expect(page.getByRole('button', { name: /Zoom In/i })).toBeVisible()
    await expect(page.getByRole('button', { name: /Zoom Out/i })).toBeVisible()
    await expect(page.getByRole('button', { name: /Reset View/i })).toBeVisible()
    
    // Click controls (functionality testing)
    await page.getByRole('button', { name: /Zoom In/i }).click()
    await page.getByRole('button', { name: /Zoom Out/i }).click()
    await page.getByRole('button', { name: /Reset View/i }).click()
  })

  test('exports molecule in different formats', async ({ page }) => {
    // Load a molecule
    await page.getByPlaceholder(/Enter SMILES string/i).fill('CCO')
    await page.getByRole('button', { name: /visualize/i }).click()
    await expect(page.getByTestId('molecule-viewer')).toBeVisible()
    
    // Test 2D export buttons
    await expect(page.getByRole('button', { name: /Export PNG/i })).toBeVisible()
    await expect(page.getByRole('button', { name: /Export SVG/i })).toBeVisible()
    
    // Switch to 3D and test 3D export
    await page.getByRole('button', { name: /3D View/i }).click()
    await expect(page.getByRole('button', { name: /Export GLB/i })).toBeVisible()
  })

  test('loads example SMILES molecules', async ({ page }) => {
    // Check examples are available
    await expect(page.getByText(/Try examples:/i)).toBeVisible()
    await expect(page.getByText('Ethanol')).toBeVisible()
    await expect(page.getByText('Benzene')).toBeVisible()
    
    // Click on example
    await page.getByText('Ethanol').click()
    
    // Check input is populated
    const input = page.getByPlaceholder(/Enter SMILES string/i)
    await expect(input).toHaveValue('CCO')
  })

  test('handles keyboard shortcuts', async ({ page }) => {
    // Load a molecule first
    await page.getByPlaceholder(/Enter SMILES string/i).fill('CCO')
    await page.getByRole('button', { name: /visualize/i }).click()
    await expect(page.getByTestId('molecule-viewer')).toBeVisible()
    
    // Test keyboard shortcuts (if implemented)
    await page.keyboard.press('Space') // Reset view
    await page.keyboard.press('Plus')  // Zoom in
    await page.keyboard.press('Minus') // Zoom out
  })

  test('persists state when switching modes', async ({ page }) => {
    // Load molecule in 2D
    await page.getByPlaceholder(/Enter SMILES string/i).fill('c1ccccc1')
    await page.getByRole('button', { name: /visualize/i }).click()
    await expect(page.getByTestId('molecule-viewer')).toBeVisible()
    
    // Switch to 3D
    await page.getByRole('button', { name: /3D View/i }).click()
    await expect(page.getByTestId('molecule-viewer')).toBeVisible()
    
    // Switch back to 2D
    await page.getByRole('button', { name: /2D View/i }).click()
    await expect(page.getByTestId('molecule-viewer')).toBeVisible()
    
    // Input should still contain the SMILES
    const input = page.getByPlaceholder(/Enter SMILES string/i)
    await expect(input).toHaveValue('c1ccccc1')
  })
}) 