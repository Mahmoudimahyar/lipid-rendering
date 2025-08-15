import React from 'react'
import { render, screen } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import JobManager from './JobManager'

describe('JobManager', () => {
  test('renders progress and buttons', async () => {
    const user = userEvent.setup()
    const onRunDocking = jest.fn()
    const onExport = jest.fn()
    const progress = { percent: 42, message: 'Running', estimated_seconds_remaining: 12 }

    render(
      <JobManager onRunDocking={onRunDocking} onExport={onExport} progress={progress} isLoading={false} />
    )

    expect(screen.getByText('Running')).toBeInTheDocument()
    expect(screen.getByText(/42%/)).toBeInTheDocument()
    expect(screen.getByText(/ETA: ~12s/)).toBeInTheDocument()

    await user.click(screen.getByRole('button', { name: /Run Docking/i }))
    expect(onRunDocking).toHaveBeenCalled()

    await user.click(screen.getByRole('button', { name: /Export Results/i }))
    expect(onExport).toHaveBeenCalled()
  })
})



