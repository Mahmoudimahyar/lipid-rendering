import React from 'react'
import { render, screen, fireEvent, waitFor } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import DockingVisualization from './DockingVisualization'

// Mock data for testing
const mockDockingResults = {
  poses: [
    {
      mode: 1,
      affinity: -8.5,
      rmsd_lb: 0.0,
      rmsd_ub: 0.0,
      center_x: 5.123,
      center_y: 10.456,
      center_z: -2.789,
      coordinates: 'MOCK_COORDINATES_1',
      sdf: 'MOCK_SDF_1'
    },
    {
      mode: 2,
      affinity: -7.2,
      rmsd_lb: 1.5,
      rmsd_ub: 2.0,
      center_x: 4.888,
      center_y: 11.234,
      center_z: -3.456,
      coordinates: 'MOCK_COORDINATES_2',
      sdf: 'MOCK_SDF_2'
    },
    {
      mode: 3,
      affinity: -6.8,
      rmsd_lb: 2.1,
      rmsd_ub: 2.5,
      center_x: 6.234,
      center_y: 9.876,
      center_z: -1.234,
      coordinates: 'MOCK_COORDINATES_3',
      sdf: 'MOCK_SDF_3'
    }
  ],
  summary: {
    num_poses: 3,
    best_affinity: -8.5,
    mean_affinity: -7.5,
    calculation_time: 3.2,
    software: 'AutoDock Vina (Mock)',
    version: '1.2.0 (simulated)'
  }
}

const defaultProps = {
  ligandSmiles: 'CCO',
  receptorPdbId: '1CRN',
  dockingResults: mockDockingResults,
  mode: '3D',
  renderer: '3Dmol.js'
}

describe('DockingVisualization', () => {
  test('renders docking visualization with pose information', () => {
    render(<DockingVisualization {...defaultProps} />)
    
    // Check main visualization area - now shows MoleculeViewer
    expect(screen.getByTestId('molecule-viewer')).toBeInTheDocument()
    expect(screen.getByText('Ligand: CCO')).toBeInTheDocument()
    expect(screen.getByText('Receptor: 1CRN')).toBeInTheDocument()
    
    // Check poses panel
    expect(screen.getByText('Poses')).toBeInTheDocument()
    expect(screen.getByText('Mode 1')).toBeInTheDocument()
    expect(screen.getByText('-8.5 kcal/mol')).toBeInTheDocument()
    
    // Check visualization controls
    expect(screen.getByText('Visualization')).toBeInTheDocument()
    expect(screen.getByText('Color Scheme')).toBeInTheDocument()
  })

  test('displays all poses with correct information', () => {
    render(<DockingVisualization {...defaultProps} />)
    
    // Should show all 3 poses
    expect(screen.getByText('Mode 1')).toBeInTheDocument()
    expect(screen.getByText('Mode 2')).toBeInTheDocument()
    expect(screen.getByText('Mode 3')).toBeInTheDocument()
    
    // Check affinity scores
    expect(screen.getByText('-8.5 kcal/mol')).toBeInTheDocument()
    expect(screen.getByText('-7.2 kcal/mol')).toBeInTheDocument()
    expect(screen.getByText('-6.8 kcal/mol')).toBeInTheDocument()
    
    // Check RMSD values
    expect(screen.getByText('RMSD: 0 - 0 Å')).toBeInTheDocument()
    expect(screen.getByText('RMSD: 1.5 - 2 Å')).toBeInTheDocument()
    expect(screen.getByText('RMSD: 2.1 - 2.5 Å')).toBeInTheDocument()
  })

  test('handles pose selection', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Click on Mode 2 pose
    const mode2Pose = screen.getByText('Mode 2').closest('div')
    await user.click(mode2Pose)
    
    // Should update selected pose (visual indication)
    await waitFor(() => {
      expect(mode2Pose.closest('.p-2')).toHaveClass('border-blue-500', 'bg-blue-50')
    })
  })

  test('handles pose visibility toggle', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Find visibility toggle buttons (eye icons)
    const visibilityButtons = screen.getAllByText('👁️')
    expect(visibilityButtons).toHaveLength(1) // Only first pose initially visible
    
    // Toggle visibility of first pose
    await user.click(visibilityButtons[0])
    
    // Should change the icon (simplified check)
    await waitFor(() => {
      expect(screen.getAllByText('👁️‍🗨️')).toHaveLength(1)
    })
  })

  test('handles show all and hide all pose controls', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Click Hide All
    const hideAllButton = screen.getByText('Hide All')
    await user.click(hideAllButton)
    
    // All poses should be hidden
    await waitFor(() => {
      expect(screen.getAllByText('👁️‍🗨️')).toHaveLength(3)
    })
    
    // Click Show All
    const showAllButton = screen.getByText('Show All')
    await user.click(showAllButton)
    
    // All poses should be visible again
    await waitFor(() => {
      expect(screen.getAllByText('👁️')).toHaveLength(3)
    })
  })

  test('handles color scheme changes', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Find color scheme selector
    const colorSchemeSelect = screen.getByDisplayValue('Binding Affinity')
    
    // Change to RMSD color scheme
    await user.selectOptions(colorSchemeSelect, 'rmsd')
    
    // Should update the legend
    await waitFor(() => {
      expect(screen.getByText('Color Scheme: RMSD')).toBeInTheDocument()
      expect(screen.getByText('Blue = Low RMSD, Yellow = High RMSD')).toBeInTheDocument()
    })
    
    // Change to Rainbow color scheme
    await user.selectOptions(colorSchemeSelect, 'rainbow')
    
    await waitFor(() => {
      expect(screen.getByText('Color Scheme: Pose Order')).toBeInTheDocument()
      expect(screen.getByText('Different colors for each pose')).toBeInTheDocument()
    })
  })

  test('handles camera controls', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Test focus on binding site
    const focusButton = screen.getByText('Focus Site')
    await user.click(focusButton)
    
    // Should update button state
    await waitFor(() => {
      expect(focusButton).toHaveClass('bg-blue-600', 'text-white')
    })
    
    // Test reset view (there are multiple Reset View buttons)
    const resetButtons = screen.getAllByText('Reset View')
    await user.click(resetButtons[0])
    
    // Focus button should return to normal state
    await waitFor(() => {
      expect(focusButton).toHaveClass('bg-white', 'text-gray-700')
    })
  })

  test('handles legend toggle', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Legend should be visible initially
    expect(screen.getByText('Color Scheme: Binding Affinity')).toBeInTheDocument()
    
    // Toggle legend off via checkbox
    const legendCheckbox = screen.getByRole('checkbox', { name: /show legend/i })
    await user.click(legendCheckbox)
    
    // Legend should be hidden
    await waitFor(() => {
      expect(screen.queryByText('Color Scheme: Binding Affinity')).not.toBeInTheDocument()
    })
    
    // Toggle legend back on
    await user.click(legendCheckbox)
    
    // Legend should be visible again
    await waitFor(() => {
      expect(screen.getByText('Color Scheme: Binding Affinity')).toBeInTheDocument()
    })
  })

  test('displays docking summary information', () => {
    render(<DockingVisualization {...defaultProps} />)
    
    // Check summary section
    expect(screen.getByText('Summary')).toBeInTheDocument()
    expect(screen.getByText('Best: -8.5 kcal/mol')).toBeInTheDocument()
    expect(screen.getByText('Average: -7.5 kcal/mol')).toBeInTheDocument()
    expect(screen.getByText('Time: 3.2s')).toBeInTheDocument()
  })

  test('handles camera control buttons in control panel', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Test binding site focus button in controls panel
    const focusButtons = screen.getAllByText('Focus on Binding Site')
    expect(focusButtons.length).toBeGreaterThan(0)
    
    await user.click(focusButtons[0])
    
    // Test reset view button in controls panel  
    const resetButtons = screen.getAllByText('Reset View')
    expect(resetButtons.length).toBeGreaterThan(0)
    
    await user.click(resetButtons[0])
  })

  test('renders without docking results', () => {
    const propsWithoutResults = {
      ...defaultProps,
      dockingResults: null
    }
    
    render(<DockingVisualization {...propsWithoutResults} />)
    
    // Should still render basic structure - now shows MoleculeViewer
    expect(screen.getByTestId('molecule-viewer')).toBeInTheDocument()
    expect(screen.getByText(/Visible poses:/)).toBeInTheDocument()
    
    // Controls should be present but might be disabled
    expect(screen.getByText('Poses')).toBeInTheDocument()
    expect(screen.getByText('Visualization')).toBeInTheDocument()
  })

  test('renders with empty poses array', () => {
    const propsWithEmptyPoses = {
      ...defaultProps,
      dockingResults: {
        poses: [],
        summary: {
          num_poses: 0,
          best_affinity: 0,
          mean_affinity: 0,
          calculation_time: 0
        }
      }
    }
    
    render(<DockingVisualization {...propsWithEmptyPoses} />)
    
    expect(screen.getByText(/Visible poses:/)).toBeInTheDocument()
    expect(screen.getByText('Best: 0 kcal/mol')).toBeInTheDocument()
  })

  test('handles missing summary data gracefully', () => {
    const propsWithoutSummary = {
      ...defaultProps,
      dockingResults: {
        poses: mockDockingResults.poses,
        summary: null
      }
    }
    
    render(<DockingVisualization {...propsWithoutSummary} />)
    
    // Should still render poses
    expect(screen.getByText('Mode 1')).toBeInTheDocument()
    expect(screen.getByText('-8.5 kcal/mol')).toBeInTheDocument()
    
    // Summary section may not be present without summary data
    expect(screen.getByText('Visualization')).toBeInTheDocument()
  })
})



import userEvent from '@testing-library/user-event'
import DockingVisualization from './DockingVisualization'

// Mock data for testing
const mockDockingResults = {
  poses: [
    {
      mode: 1,
      affinity: -8.5,
      rmsd_lb: 0.0,
      rmsd_ub: 0.0,
      center_x: 5.123,
      center_y: 10.456,
      center_z: -2.789,
      coordinates: 'MOCK_COORDINATES_1',
      sdf: 'MOCK_SDF_1'
    },
    {
      mode: 2,
      affinity: -7.2,
      rmsd_lb: 1.5,
      rmsd_ub: 2.0,
      center_x: 4.888,
      center_y: 11.234,
      center_z: -3.456,
      coordinates: 'MOCK_COORDINATES_2',
      sdf: 'MOCK_SDF_2'
    },
    {
      mode: 3,
      affinity: -6.8,
      rmsd_lb: 2.1,
      rmsd_ub: 2.5,
      center_x: 6.234,
      center_y: 9.876,
      center_z: -1.234,
      coordinates: 'MOCK_COORDINATES_3',
      sdf: 'MOCK_SDF_3'
    }
  ],
  summary: {
    num_poses: 3,
    best_affinity: -8.5,
    mean_affinity: -7.5,
    calculation_time: 3.2,
    software: 'AutoDock Vina (Mock)',
    version: '1.2.0 (simulated)'
  }
}

const defaultProps = {
  ligandSmiles: 'CCO',
  receptorPdbId: '1CRN',
  dockingResults: mockDockingResults,
  mode: '3D',
  renderer: '3Dmol.js'
}

describe('DockingVisualization', () => {
  test('renders docking visualization with pose information', () => {
    render(<DockingVisualization {...defaultProps} />)
    
    // Check main visualization area - now shows MoleculeViewer
    expect(screen.getByTestId('molecule-viewer')).toBeInTheDocument()
    expect(screen.getByText('Ligand: CCO')).toBeInTheDocument()
    expect(screen.getByText('Receptor: 1CRN')).toBeInTheDocument()
    
    // Check poses panel
    expect(screen.getByText('Poses')).toBeInTheDocument()
    expect(screen.getByText('Mode 1')).toBeInTheDocument()
    expect(screen.getByText('-8.5 kcal/mol')).toBeInTheDocument()
    
    // Check visualization controls
    expect(screen.getByText('Visualization')).toBeInTheDocument()
    expect(screen.getByText('Color Scheme')).toBeInTheDocument()
  })

  test('displays all poses with correct information', () => {
    render(<DockingVisualization {...defaultProps} />)
    
    // Should show all 3 poses
    expect(screen.getByText('Mode 1')).toBeInTheDocument()
    expect(screen.getByText('Mode 2')).toBeInTheDocument()
    expect(screen.getByText('Mode 3')).toBeInTheDocument()
    
    // Check affinity scores
    expect(screen.getByText('-8.5 kcal/mol')).toBeInTheDocument()
    expect(screen.getByText('-7.2 kcal/mol')).toBeInTheDocument()
    expect(screen.getByText('-6.8 kcal/mol')).toBeInTheDocument()
    
    // Check RMSD values
    expect(screen.getByText('RMSD: 0 - 0 Å')).toBeInTheDocument()
    expect(screen.getByText('RMSD: 1.5 - 2 Å')).toBeInTheDocument()
    expect(screen.getByText('RMSD: 2.1 - 2.5 Å')).toBeInTheDocument()
  })

  test('handles pose selection', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Click on Mode 2 pose
    const mode2Pose = screen.getByText('Mode 2').closest('div')
    await user.click(mode2Pose)
    
    // Should update selected pose (visual indication)
    await waitFor(() => {
      expect(mode2Pose.closest('.p-2')).toHaveClass('border-blue-500', 'bg-blue-50')
    })
  })

  test('handles pose visibility toggle', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Find visibility toggle buttons (eye icons)
    const visibilityButtons = screen.getAllByText('👁️')
    expect(visibilityButtons).toHaveLength(1) // Only first pose initially visible
    
    // Toggle visibility of first pose
    await user.click(visibilityButtons[0])
    
    // Should change the icon (simplified check)
    await waitFor(() => {
      expect(screen.getAllByText('👁️‍🗨️')).toHaveLength(1)
    })
  })

  test('handles show all and hide all pose controls', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Click Hide All
    const hideAllButton = screen.getByText('Hide All')
    await user.click(hideAllButton)
    
    // All poses should be hidden
    await waitFor(() => {
      expect(screen.getAllByText('👁️‍🗨️')).toHaveLength(3)
    })
    
    // Click Show All
    const showAllButton = screen.getByText('Show All')
    await user.click(showAllButton)
    
    // All poses should be visible again
    await waitFor(() => {
      expect(screen.getAllByText('👁️')).toHaveLength(3)
    })
  })

  test('handles color scheme changes', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Find color scheme selector
    const colorSchemeSelect = screen.getByDisplayValue('Binding Affinity')
    
    // Change to RMSD color scheme
    await user.selectOptions(colorSchemeSelect, 'rmsd')
    
    // Should update the legend
    await waitFor(() => {
      expect(screen.getByText('Color Scheme: RMSD')).toBeInTheDocument()
      expect(screen.getByText('Blue = Low RMSD, Yellow = High RMSD')).toBeInTheDocument()
    })
    
    // Change to Rainbow color scheme
    await user.selectOptions(colorSchemeSelect, 'rainbow')
    
    await waitFor(() => {
      expect(screen.getByText('Color Scheme: Pose Order')).toBeInTheDocument()
      expect(screen.getByText('Different colors for each pose')).toBeInTheDocument()
    })
  })

  test('handles camera controls', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Test focus on binding site
    const focusButton = screen.getByText('Focus Site')
    await user.click(focusButton)
    
    // Should update button state
    await waitFor(() => {
      expect(focusButton).toHaveClass('bg-blue-600', 'text-white')
    })
    
    // Test reset view (there are multiple Reset View buttons)
    const resetButtons = screen.getAllByText('Reset View')
    await user.click(resetButtons[0])
    
    // Focus button should return to normal state
    await waitFor(() => {
      expect(focusButton).toHaveClass('bg-white', 'text-gray-700')
    })
  })

  test('handles legend toggle', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Legend should be visible initially
    expect(screen.getByText('Color Scheme: Binding Affinity')).toBeInTheDocument()
    
    // Toggle legend off via checkbox
    const legendCheckbox = screen.getByRole('checkbox', { name: /show legend/i })
    await user.click(legendCheckbox)
    
    // Legend should be hidden
    await waitFor(() => {
      expect(screen.queryByText('Color Scheme: Binding Affinity')).not.toBeInTheDocument()
    })
    
    // Toggle legend back on
    await user.click(legendCheckbox)
    
    // Legend should be visible again
    await waitFor(() => {
      expect(screen.getByText('Color Scheme: Binding Affinity')).toBeInTheDocument()
    })
  })

  test('displays docking summary information', () => {
    render(<DockingVisualization {...defaultProps} />)
    
    // Check summary section
    expect(screen.getByText('Summary')).toBeInTheDocument()
    expect(screen.getByText('Best: -8.5 kcal/mol')).toBeInTheDocument()
    expect(screen.getByText('Average: -7.5 kcal/mol')).toBeInTheDocument()
    expect(screen.getByText('Time: 3.2s')).toBeInTheDocument()
  })

  test('handles camera control buttons in control panel', async () => {
    const user = userEvent.setup()
    render(<DockingVisualization {...defaultProps} />)
    
    // Test binding site focus button in controls panel
    const focusButtons = screen.getAllByText('Focus on Binding Site')
    expect(focusButtons.length).toBeGreaterThan(0)
    
    await user.click(focusButtons[0])
    
    // Test reset view button in controls panel  
    const resetButtons = screen.getAllByText('Reset View')
    expect(resetButtons.length).toBeGreaterThan(0)
    
    await user.click(resetButtons[0])
  })

  test('renders without docking results', () => {
    const propsWithoutResults = {
      ...defaultProps,
      dockingResults: null
    }
    
    render(<DockingVisualization {...propsWithoutResults} />)
    
    // Should still render basic structure - now shows MoleculeViewer
    expect(screen.getByTestId('molecule-viewer')).toBeInTheDocument()
    expect(screen.getByText(/Visible poses:/)).toBeInTheDocument()
    
    // Controls should be present but might be disabled
    expect(screen.getByText('Poses')).toBeInTheDocument()
    expect(screen.getByText('Visualization')).toBeInTheDocument()
  })

  test('renders with empty poses array', () => {
    const propsWithEmptyPoses = {
      ...defaultProps,
      dockingResults: {
        poses: [],
        summary: {
          num_poses: 0,
          best_affinity: 0,
          mean_affinity: 0,
          calculation_time: 0
        }
      }
    }
    
    render(<DockingVisualization {...propsWithEmptyPoses} />)
    
    expect(screen.getByText(/Visible poses:/)).toBeInTheDocument()
    expect(screen.getByText('Best: 0 kcal/mol')).toBeInTheDocument()
  })

  test('handles missing summary data gracefully', () => {
    const propsWithoutSummary = {
      ...defaultProps,
      dockingResults: {
        poses: mockDockingResults.poses,
        summary: null
      }
    }
    
    render(<DockingVisualization {...propsWithoutSummary} />)
    
    // Should still render poses
    expect(screen.getByText('Mode 1')).toBeInTheDocument()
    expect(screen.getByText('-8.5 kcal/mol')).toBeInTheDocument()
    
    // Summary section may not be present without summary data
    expect(screen.getByText('Visualization')).toBeInTheDocument()
  })
})


