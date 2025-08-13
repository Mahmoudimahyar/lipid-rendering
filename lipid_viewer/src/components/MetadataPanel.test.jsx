import React from 'react'
import { render, screen, fireEvent, waitFor } from '@testing-library/react'
import { vi } from 'vitest'
import MetadataPanel from './MetadataPanel'

// Mock job data for testing
const mockJobDataReal = {
  job_id: 'test-job-123',
  status: 'completed',
  started_at: '2024-01-15T10:00:00Z',
  completed_at: '2024-01-15T10:05:30Z',
  duration: 330,
  docking_parameters: {
    center_x: 1.5,
    center_y: 2.5,
    center_z: 3.5,
    size_x: 20.0,
    size_y: 20.0,
    size_z: 20.0,
    exhaustiveness: 8,
    num_modes: 5,
    seed: 42
  },
  engine_metadata: {
    engine: 'vina',
    is_mock: false,
    method: 'AutoDock Vina',
    version: '1.2.5',
    calculation_time: 280.5
  },
  performance_metrics: {
    best_affinity: -8.5,
    num_poses_generated: 5,
    calculation_time: 280.5
  }
}

const mockJobDataMock = {
  ...mockJobDataReal,
  engine_metadata: {
    engine: 'mock',
    is_mock: true,
    method: 'Mock Docking Engine',
    version: 'mock-1.0',
    calculation_time: 15.2
  },
  performance_metrics: {
    best_affinity: -7.2,
    num_poses_generated: 3,
    calculation_time: 15.2
  }
}

describe('MetadataPanel', () => {
  
  describe('Basic Rendering', () => {
    
    test('renders without crashing with no data', () => {
      render(<MetadataPanel jobData={null} />)
      expect(screen.getByText('No metadata available')).toBeInTheDocument()
    })
    
    test('renders job metadata header', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      expect(screen.getByText('Job Metadata')).toBeInTheDocument()
    })
    
    test('shows real engine badge for real Vina', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      const badge = screen.getByText('AutoDock Vina 1.2.5')
      expect(badge).toBeInTheDocument()
      expect(badge).toHaveClass('bg-green-100', 'text-green-800')
      
      const realBadge = screen.getByLabelText(/Real AutoDock Vina engine version 1.2.5/)
      expect(realBadge).toBeInTheDocument()
    })
    
    test('shows mock engine badge for mock engine', () => {
      render(<MetadataPanel jobData={mockJobDataMock} />)
      
      const badge = screen.getByText('MOCK (Demo Only)')
      expect(badge).toBeInTheDocument()
      expect(badge).toHaveClass('bg-red-100', 'text-red-800')
      
      const mockBadge = screen.getByLabelText(/Mock engine - demonstration mode only, not real docking/)
      expect(mockBadge).toBeInTheDocument()
    })
    
  })
  
  describe('Quick Summary Display', () => {
    
    test('displays status information', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      expect(screen.getByText('Status')).toBeInTheDocument()
      expect(screen.getByText('completed')).toBeInTheDocument()
    })
    
    test('displays runtime information', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      expect(screen.getByText('Runtime')).toBeInTheDocument()
      expect(screen.getByText('5m 30.0s')).toBeInTheDocument()
    })
    
    test('displays seed information', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      expect(screen.getByText('Seed')).toBeInTheDocument()
      expect(screen.getByText('42')).toBeInTheDocument()
    })
    
    test('displays poses count', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      expect(screen.getByText('Poses')).toBeInTheDocument()
      expect(screen.getByText('5')).toBeInTheDocument()
    })
    
    test('displays best affinity when available', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      expect(screen.getByText('Best Affinity:')).toBeInTheDocument()
      expect(screen.getByText('-8.5 kcal/mol')).toBeInTheDocument()
    })
    
    test('handles missing seed (shows Random)', () => {
      const jobDataNoSeed = {
        ...mockJobDataReal,
        docking_parameters: {
          ...mockJobDataReal.docking_parameters,
          seed: null
        }
      }
      
      render(<MetadataPanel jobData={jobDataNoSeed} />)
      expect(screen.getByText('Random')).toBeInTheDocument()
    })
    
  })
  
  describe('Expandable Details', () => {
    
    test('expands and collapses details on button click', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      const expandButton = screen.getByLabelText('Expand metadata details')
      expect(expandButton).toBeInTheDocument()
      
      // Initially collapsed
      expect(screen.queryByText('Engine Information')).not.toBeInTheDocument()
      
      // Expand
      fireEvent.click(expandButton)
      expect(screen.getByText('Engine Information')).toBeInTheDocument()
      expect(screen.getByText('Docking Parameters')).toBeInTheDocument()
      expect(screen.getByText('Timing Information')).toBeInTheDocument()
      
      // Check aria-expanded attribute
      expect(expandButton).toHaveAttribute('aria-expanded', 'true')
      
      // Collapse
      const collapseButton = screen.getByLabelText('Collapse metadata details')
      fireEvent.click(collapseButton)
      expect(screen.queryByText('Engine Information')).not.toBeInTheDocument()
    })
    
    test('shows engine details when expanded', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      const expandButton = screen.getByLabelText('Expand metadata details')
      fireEvent.click(expandButton)
      
      // Check engine information
      expect(screen.getByText('Method')).toBeInTheDocument()
      expect(screen.getByText('AutoDock Vina')).toBeInTheDocument()
      expect(screen.getByText('Version')).toBeInTheDocument()
      expect(screen.getByText('1.2.5')).toBeInTheDocument()
      expect(screen.getByText('Mock Mode')).toBeInTheDocument()
      expect(screen.getByText('No')).toBeInTheDocument()
    })
    
    test('shows docking parameters when expanded', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      const expandButton = screen.getByLabelText('Expand metadata details')
      fireEvent.click(expandButton)
      
      // Check docking parameters
      expect(screen.getByText('Exhaustiveness')).toBeInTheDocument()
      expect(screen.getByText('8')).toBeInTheDocument()
      expect(screen.getByText('Number of Modes')).toBeInTheDocument()
      expect(screen.getByText('5')).toBeInTheDocument()
      
      // Check binding site information
      expect(screen.getByText('Center')).toBeInTheDocument()
      expect(screen.getByText('(1.5, 2.5, 3.5)')).toBeInTheDocument()
      expect(screen.getByText('Size')).toBeInTheDocument()
      expect(screen.getByText('20 × 20 × 20 Å')).toBeInTheDocument()
    })
    
    test('shows timing information when expanded', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      const expandButton = screen.getByLabelText('Expand metadata details')
      fireEvent.click(expandButton)
      
      // Check timing information
      expect(screen.getByText('Started At')).toBeInTheDocument()
      expect(screen.getByText('Completed At')).toBeInTheDocument()
      expect(screen.getByText('Total Duration')).toBeInTheDocument()
    })
    
  })
  
  describe('Export Functionality', () => {
    
    test('shows export button when enabled', () => {
      const mockExport = vi.fn()
      render(
        <MetadataPanel 
          jobData={mockJobDataReal} 
          showExportButton={true}
          onExport={mockExport}
        />
      )
      
      const exportButton = screen.getByLabelText('Export docking results with metadata')
      expect(exportButton).toBeInTheDocument()
      expect(exportButton).toHaveTextContent('Export')
    })
    
    test('hides export button when showExportButton is false', () => {
      const mockExport = vi.fn()
      render(
        <MetadataPanel 
          jobData={mockJobDataReal} 
          showExportButton={false}
          onExport={mockExport}
        />
      )
      
      expect(screen.queryByLabelText('Export docking results with metadata')).not.toBeInTheDocument()
    })
    
    test('disables export button for incomplete jobs', () => {
      const mockExport = vi.fn()
      const incompleteJob = {
        ...mockJobDataReal,
        status: 'running'
      }
      
      render(
        <MetadataPanel 
          jobData={incompleteJob} 
          onExport={mockExport}
        />
      )
      
      const exportButton = screen.getByLabelText('Export docking results with metadata')
      expect(exportButton).toBeDisabled()
    })
    
    test('calls export function when button clicked', async () => {
      const mockExport = vi.fn().mockResolvedValue()
      render(
        <MetadataPanel 
          jobData={mockJobDataReal} 
          onExport={mockExport}
        />
      )
      
      const exportButton = screen.getByLabelText('Export docking results with metadata')
      fireEvent.click(exportButton)
      
      expect(mockExport).toHaveBeenCalledTimes(1)
    })
    
    test('shows loading state during export', async () => {
      let resolveExport
      const mockExport = vi.fn(() => new Promise(resolve => {
        resolveExport = resolve
      }))
      
      render(
        <MetadataPanel 
          jobData={mockJobDataReal} 
          onExport={mockExport}
        />
      )
      
      const exportButton = screen.getByLabelText('Export docking results with metadata')
      fireEvent.click(exportButton)
      
      // Should show loading state
      expect(screen.getByText('Exporting...')).toBeInTheDocument()
      expect(exportButton).toBeDisabled()
      
      // Complete export
      resolveExport()
      await waitFor(() => {
        expect(screen.getByText('Export')).toBeInTheDocument()
      })
    })
    
  })
  
  describe('Accessibility Features', () => {
    
    test('has proper ARIA labels for engine badges', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      const badge = screen.getByRole('status')
      expect(badge).toHaveAttribute('aria-label', 'Real AutoDock Vina engine version 1.2.5')
    })
    
    test('has proper ARIA labels for mock engine', () => {
      render(<MetadataPanel jobData={mockJobDataMock} />)
      
      const badge = screen.getByRole('status')
      expect(badge).toHaveAttribute('aria-label', 'Mock engine - demonstration mode only, not real docking')
    })
    
    test('has keyboard navigation for expand/collapse', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      const expandButton = screen.getByLabelText('Expand metadata details')
      expect(expandButton).toHaveAttribute('aria-expanded', 'false')
      
      // Test keyboard activation
      expandButton.focus()
      fireEvent.keyDown(expandButton, { key: 'Enter' })
      
      expect(expandButton).toHaveAttribute('aria-expanded', 'true')
    })
    
    test('has proper focus management for export button', () => {
      const mockExport = vi.fn()
      render(
        <MetadataPanel 
          jobData={mockJobDataReal} 
          onExport={mockExport}
        />
      )
      
      const exportButton = screen.getByLabelText('Export docking results with metadata')
      expect(exportButton).toHaveClass('focus:outline-none', 'focus:ring-2', 'focus:ring-blue-500', 'focus:ring-offset-2')
    })
    
  })
  
  describe('Color-blind Friendly Indicators', () => {
    
    test('uses both color and text for mock/real distinction', () => {
      // Real engine
      render(<MetadataPanel jobData={mockJobDataReal} />)
      const realIcon = screen.getByText('✅')
      expect(realIcon).toBeInTheDocument()
      
      // Mock engine  
      render(<MetadataPanel jobData={mockJobDataMock} />)
      const mockIcon = screen.getByText('⚠️')
      expect(mockIcon).toBeInTheDocument()
    })
    
    test('provides text alternatives for all color-coded information', () => {
      render(<MetadataPanel jobData={mockJobDataReal} />)
      
      // Status uses both color and text
      expect(screen.getByText('completed')).toBeInTheDocument()
      
      // Engine type uses both badge color and text content
      expect(screen.getByText('AutoDock Vina 1.2.5')).toBeInTheDocument()
      
      // Mock mode indicator uses text
      const expandButton = screen.getByLabelText('Expand metadata details')
      fireEvent.click(expandButton)
      expect(screen.getByText('No')).toBeInTheDocument() // Mock Mode: No
    })
    
  })
  
  describe('Data Formatting', () => {
    
    test('formats duration correctly', () => {
      const testCases = [
        { duration: 30, expected: '30.0s' },
        { duration: 90, expected: '1m 30.0s' },
        { duration: 3661, expected: '61m 1.0s' }
      ]
      
      testCases.forEach(({ duration, expected }) => {
        const jobData = { ...mockJobDataReal, duration }
        render(<MetadataPanel jobData={jobData} />)
        expect(screen.getByText(expected)).toBeInTheDocument()
      })
    })
    
    test('handles null/undefined values gracefully', () => {
      const incompleteJobData = {
        status: 'running',
        docking_parameters: {
          seed: null,
          exhaustiveness: undefined
        },
        engine_metadata: {
          is_mock: false
        },
        performance_metrics: {}
      }
      
      render(<MetadataPanel jobData={incompleteJobData} />)
      
      // Should not crash and should show appropriate defaults
      expect(screen.getByText('Random')).toBeInTheDocument() // for null seed
      expect(screen.getByText('running')).toBeInTheDocument()
    })
    
  })
  
})
