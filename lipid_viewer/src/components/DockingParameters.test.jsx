import React from 'react'
import { render, screen } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import DockingParameters from './DockingParameters'

describe('DockingParameters', () => {
  test('renders and updates parameters', async () => {
    const user = userEvent.setup()
    const setParams = jest.fn()
    const params = { exhaustiveness: 8, num_modes: 9, seed: null }

    render(<DockingParameters params={params} setParams={setParams} />)

    const ex = screen.getByLabelText(/Exhaustiveness/i)
    await user.clear(ex)
    await user.type(ex, '16')
    expect(setParams).toHaveBeenCalled()
  })
})



