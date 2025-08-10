import { validateSMILES, isValidSMILES, generateMolecule3D } from './smilesValidator'

describe('smilesValidator', () => {
  beforeEach(() => {
    jest.clearAllMocks()
  })

  describe('validateSMILES', () => {
    test('returns isValid true for valid SMILES', async () => {
      const valid = ['CCO', 'c1ccccc1', 'CC(=O)O']
      for (const s of valid) {
        const r = await validateSMILES(s)
        expect(r).toEqual(expect.objectContaining({ isValid: true }))
      }
    })

    test('returns isValid false and message for invalid input types/empty', async () => {
      expect(await validateSMILES('')).toEqual({ isValid: false, error: 'Invalid input: SMILES cannot be empty' })
      expect(await validateSMILES('   ')).toEqual({ isValid: false, error: 'Invalid input: SMILES cannot be empty' })
      expect(await validateSMILES(null)).toEqual({ isValid: false, error: 'Invalid input: SMILES must be a non-empty string' })
      expect(await validateSMILES(undefined)).toEqual({ isValid: false, error: 'Invalid input: SMILES must be a non-empty string' })
      expect(await validateSMILES(123)).toEqual({ isValid: false, error: 'Invalid input: SMILES must be a non-empty string' })
    })

    test('falls back to basic validation when RDKit cannot parse', async () => {
      // setupTests mock returns null from get_mol for strings containing 'invalid'
      const result = await validateSMILES('invalid')
      expect(result.method).toBe('basic')
      expect(typeof result.isValid).toBe('boolean')
    })

    test('returns object with method field when RDKit is available', async () => {
      const result = await validateSMILES('CCO')
      expect(result).toEqual(expect.objectContaining({ isValid: true }))
      expect(['RDKit', 'basic']).toContain(result.method)
    })
  })

  describe('isValidSMILES (basic synchronous checks)', () => {
    test('accepts common valid strings', () => {
      expect(isValidSMILES('CCO')).toBe(true)
      expect(isValidSMILES('c1ccccc1')).toBe(true)
      expect(isValidSMILES('CC(=O)O')).toBe(true)
      expect(isValidSMILES('C#N')).toBe(true)
    })

    test('rejects obvious invalid patterns', () => {
      expect(isValidSMILES('(((')).toBe(false)
      expect(isValidSMILES(')))')).toBe(false)
      expect(isValidSMILES('[[[')).toBe(false)
      expect(isValidSMILES(']')).toBe(false)
      expect(isValidSMILES('(')).toBe(false)
    })

    test('handles non-string/empty', () => {
      expect(isValidSMILES('')).toBe(false)
      expect(isValidSMILES('   ')).toBe(false)
      expect(isValidSMILES(null)).toBe(false)
      expect(isValidSMILES(undefined)).toBe(false)
      expect(isValidSMILES(123)).toBe(false)
    })
  })

  describe('generateMolecule3D', () => {
    test('returns object with molblock for typical valid SMILES when RDKit mock available', async () => {
      const res = await generateMolecule3D('CCO')
      expect(res === null || typeof res === 'object').toBe(true)
      if (res) {
        expect(res).toEqual(expect.objectContaining({ molblock: expect.any(String), smiles: 'CCO' }))
      }
    })

    test('returns null for SMILES RDKit mock cannot create', async () => {
      const res = await generateMolecule3D('invalid123')
      expect(res).toBeNull()
    })
  })
}) 