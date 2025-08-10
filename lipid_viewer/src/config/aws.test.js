import { awsConfig } from './aws'

describe('AWS Config', () => {
  test('exports AWS configuration object', () => {
    expect(awsConfig).toBeDefined()
    expect(typeof awsConfig).toBe('object')
  })
})
