import React, { useState, useEffect } from 'react'
import { isValidSMILES, validateSMILES } from '../utils/smilesValidator'
import toast from 'react-hot-toast'

const EXAMPLE_MOLECULES = [
  { name: 'Ethanol', smiles: 'CCO' },
  { name: 'Benzene', smiles: 'c1ccccc1' },
  { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' },
  { name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O' },
  { name: 'Cholesterol', smiles: 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C' },
]

const SMILESInput = ({ onSubmit, onValidation, isValid = null }) => {
  const [smiles, setSmiles] = useState('')
  const [isValidating, setIsValidating] = useState(false)
  const [validationResult, setValidationResult] = useState(null)
  const [isValidSMILES, setIsValidSMILES] = useState(false)
  const [validationMessage, setValidationMessage] = useState('')

  // Real-time validation
  useEffect(() => {
    const performValidation = async () => {
      if (smiles.trim()) {
        setIsValidating(true)
        setValidationMessage('')
        
        try {
          const validation = await validateSMILES(smiles.trim())
          setIsValidSMILES(validation.isValid)
          setValidationResult(validation.isValid)
          
          if (validation.isValid) {
            setValidationMessage(`Valid SMILES (${validation.method || 'basic'} validation)`)
            toast.success(`Molecule loaded: ${smiles.trim()}`)
          } else {
            setValidationMessage(validation.error || 'Invalid SMILES string')
            toast.error('Invalid SMILES string')
          }
          
          // Call onValidation callback if provided
          if (onValidation) {
            onValidation(smiles.trim(), validation.isValid)
          }
        } catch (error) {
          console.error('Validation error:', error)
          setIsValidSMILES(false)
          setValidationResult(false)
          setValidationMessage('Validation failed')
          toast.error('Validation failed')
          
          if (onValidation) {
            onValidation(smiles.trim(), false)
          }
        } finally {
          setIsValidating(false)
        }
      } else {
        setIsValidSMILES(false)
        setValidationResult(null)
        setValidationMessage('')
        
        if (onValidation) {
          onValidation('', false)
        }
      }
    }

    const timeoutId = setTimeout(performValidation, 500)
    return () => clearTimeout(timeoutId)
  }, [smiles, onValidation])

  const handleSubmit = (e) => {
    e.preventDefault()
    if (smiles.trim() && isValidSMILES) {
      onSubmit(smiles.trim())
    }
  }

  const handleExampleClick = (exampleSmiles) => {
    setSmiles(exampleSmiles)
  }

  const currentValidState = isValid !== null ? isValid : isValidSMILES

  return (
    <div className="w-full max-w-4xl mx-auto p-6 bg-white rounded-lg shadow-lg">
      <h2 className="text-2xl font-bold text-gray-800 mb-4">Enter SMILES String</h2>
      
      <form onSubmit={handleSubmit} className="space-y-4">
        <div className="relative">
          <input
            type="text"
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            placeholder="Enter SMILES string (e.g., CCO for ethanol)"
            className={`w-full px-4 py-3 text-lg border-2 rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500 transition-colors ${
              smiles && currentValidState === false
                ? 'border-red-400 bg-red-50'
                : smiles && currentValidState === true
                ? 'border-green-400 bg-green-50'
                : 'border-gray-300'
            }`}
            autoComplete="off"
            spellCheck="false"
          />
          
          {isValidating && (
            <div className="absolute right-3 top-1/2 transform -translate-y-1/2">
              <div className="animate-spin rounded-full h-5 w-5 border-b-2 border-blue-600"></div>
            </div>
          )}
          
          {smiles && !isValidating && currentValidState !== null && (
            <div className="absolute right-3 top-1/2 transform -translate-y-1/2">
              {currentValidState ? (
                <svg className="h-5 w-5 text-green-500" fill="none" strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" viewBox="0 0 24 24" stroke="currentColor">
                  <path d="M5 13l4 4L19 7"></path>
                </svg>
              ) : (
                <svg className="h-5 w-5 text-red-500" fill="none" strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" viewBox="0 0 24 24" stroke="currentColor">
                  <path d="M6 18L18 6M6 6l12 12"></path>
                </svg>
              )}
            </div>
          )}
        </div>
        
        {validationMessage && (
          <p className={`text-sm ${isValidSMILES ? 'text-green-600' : 'text-red-600'}`}>
            {validationMessage}
          </p>
        )}
        
        <button
          type="submit"
          disabled={!smiles.trim() || !isValidSMILES || isValidating}
          className="control-button w-full py-3 text-lg font-semibold"
        >
          {isValidating ? 'Validating...' : 'Visualize Molecule'}
        </button>
      </form>
      
      <div className="mt-6">
        <h3 className="text-lg font-semibold text-gray-700 mb-3">Try examples:</h3>
        <div className="flex flex-wrap gap-2">
          {EXAMPLE_MOLECULES.map((molecule) => (
            <button
              key={molecule.name}
              onClick={() => handleExampleClick(molecule.smiles)}
              className="px-3 py-2 bg-gray-100 hover:bg-blue-100 text-gray-700 rounded-md transition-colors text-sm font-medium"
            >
              {molecule.name}
            </button>
          ))}
        </div>
      </div>
      
      <div className="mt-4 text-sm text-gray-600">
        <p>
          <strong>SMILES</strong> (Simplified Molecular Input Line Entry System) is a notation for describing chemical structures.
          Examples: <code className="bg-gray-100 px-1 rounded">CCO</code> (ethanol), 
          <code className="bg-gray-100 px-1 rounded">c1ccccc1</code> (benzene).
        </p>
      </div>
    </div>
  )
}

export default SMILESInput 