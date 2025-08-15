import React, { createRef } from 'react'
import { render, screen } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import BindingSiteConfigurator from './BindingSiteConfigurator'

function makeRefs() {
  return {
    centerXRef: createRef(),
    centerYRef: createRef(),
    centerZRef: createRef(),
    sizeXRef: createRef(),
    sizeYRef: createRef(),
    sizeZRef: createRef(),
  }
}

describe('BindingSiteConfigurator', () => {
  test('renders inputs and updates bindingSite', async () => {
    const user = userEvent.setup()
    const refs = makeRefs()
    const setBindingSite = jest.fn()
    const site = { center_x: 0, center_y: 0, center_z: 0, size_x: 20, size_y: 20, size_z: 20 }

    render(
      <BindingSiteConfigurator bindingSite={site} setBindingSite={setBindingSite} refs={refs} />
    )

    const cx = screen.getByLabelText(/Center X/i)
    await user.clear(cx)
    await user.type(cx, '5')
    expect(setBindingSite).toHaveBeenCalled()
  })
})



