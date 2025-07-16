// Initialize RDKit from CDN
let rdkitInstance = null
let rdkitPromise = null

const initRDKit = async () => {
  // Return cached instance if available
  if (rdkitInstance) {
    return rdkitInstance
  }
  
  // Return existing promise if initialization is in progress
  if (rdkitPromise) {
    return rdkitPromise
  }

  rdkitPromise = (async () => {
    try {
      // Check if RDKit is already loaded
      if (window.RDKit && window.RDKit.get_mol) {
        rdkitInstance = window.RDKit
        return rdkitInstance
      }

      if (!window.initRDKitModule) {
        // Try multiple CDN sources for better reliability
        const cdnUrls = [
          'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js',
          'https://cdn.jsdelivr.net/npm/@rdkit/rdkit/dist/RDKit_minimal.js',
          'https://unpkg.com/@rdkit/rdkit@2024.3.1/dist/RDKit_minimal.js'
        ]

        let loadSuccess = false
        for (const url of cdnUrls) {
          try {
            console.log(`Attempting to load RDKit from: ${url}`)
            const script = document.createElement('script')
            script.src = url
            script.crossOrigin = 'anonymous'
            document.head.appendChild(script)
            
            await new Promise((resolve, reject) => {
              script.onload = () => {
                console.log(`Successfully loaded RDKit from: ${url}`)
                resolve()
              }
              script.onerror = (error) => {
                console.warn(`Failed to load RDKit from ${url}:`, error)
                reject(error)
              }
              // Add timeout to prevent hanging
              setTimeout(() => reject(new Error('RDKit loading timeout')), 5000)
            })
            
            loadSuccess = true
            break
          } catch (error) {
            console.warn(`Failed to load from ${url}, trying next...`)
            continue
          }
        }
        
        if (!loadSuccess) {
          throw new Error('All RDKit CDN sources failed to load')
        }
      }
      
      if (window.initRDKitModule) {
        console.log('Initializing RDKit module...')
        rdkitInstance = await window.initRDKitModule()
        console.log('RDKit initialized successfully')
        return rdkitInstance
      }
      
      throw new Error('RDKit module not available after script load')
    } catch (error) {
      console.warn('RDKit initialization failed:', error)
      rdkitPromise = null // Reset promise so we can try again later
      return null
    }
  })()

  return rdkitPromise
}

/**
 * Validates a SMILES string using RDKit.js
 * @param {string} smiles - The SMILES string to validate
 * @returns {Promise<boolean>} - True if valid, false otherwise
 */
export const validateSMILES = async (smiles) => {
  if (!smiles || typeof smiles !== 'string') {
    return { isValid: false, error: 'Invalid input: SMILES must be a non-empty string' }
  }

  const trimmedSmiles = smiles.trim()
  if (!trimmedSmiles) {
    return { isValid: false, error: 'Invalid input: SMILES cannot be empty' }
  }

  // First try basic validation (always available)
  const basicValidation = basicValidateSMILES(trimmedSmiles)
  
  try {
    // Try to use RDKit for advanced validation
    console.log('Attempting RDKit validation for:', trimmedSmiles)
    const RDKit = await initRDKit()
    
    if (RDKit) {
      console.log('RDKit available, performing advanced validation')
      const mol = RDKit.get_mol(trimmedSmiles)
      if (mol) {
        const isValid = mol.is_valid()
        mol.delete() // Clean up
        console.log('RDKit validation result:', isValid)
        return { 
          isValid, 
          error: isValid ? null : 'Invalid SMILES structure (RDKit validation failed)',
          method: 'RDKit'
        }
      } else {
        console.warn('RDKit could not parse molecule, falling back to basic validation')
        return { 
          isValid: basicValidation.isValid, 
          error: basicValidation.isValid ? null : 'Invalid SMILES structure (RDKit could not parse)',
          method: 'basic'
        }
      }
    } else {
      console.warn('RDKit not available, using basic validation')
      return { 
        isValid: basicValidation.isValid, 
        error: basicValidation.isValid ? null : basicValidation.error,
        method: 'basic'
      }
    }
  } catch (error) {
    console.warn('RDKit validation failed, falling back to basic validation:', error)
    return { 
      isValid: basicValidation.isValid, 
      error: basicValidation.isValid ? null : `Validation error: ${error.message}`,
      method: 'basic'
    }
  }
}

/**
 * Synchronous basic SMILES validation for quick feedback
 * @param {string} smiles - The SMILES string to validate
 * @returns {boolean} - True if passes basic validation
 */
export const isValidSMILES = (smiles) => {
  if (!smiles || typeof smiles !== 'string') {
    return false
  }
  
  const trimmedSmiles = smiles.trim()
  if (!trimmedSmiles) {
    return false
  }
  
  // Check for balanced parentheses and brackets
  const parenthesesCount = (trimmedSmiles.match(/\(/g) || []).length - (trimmedSmiles.match(/\)/g) || []).length
  const bracketsCount = (trimmedSmiles.match(/\[/g) || []).length - (trimmedSmiles.match(/\]/g) || []).length
  
  if (parenthesesCount !== 0 || bracketsCount !== 0) {
    return false
  }
  
  // Check for valid atom symbols (basic check)
  const invalidAtomPattern = /[^A-Za-z0-9()[\]=#@+-:/.\\]/
  if (invalidAtomPattern.test(trimmedSmiles)) {
    return false
  }
  
  // Check for some obviously invalid patterns
  const invalidPatterns = [
    /^\(/, // Starts with (
    /\)$/, // Ends with )
    /^\[/, // Starts with [
    /\]$/, // Ends with ]
    /[XZ]/, // Invalid atom symbols
  ]
  
  return !invalidPatterns.some(pattern => pattern.test(trimmedSmiles))
}

/**
 * Generate 3D coordinates for a SMILES string
 * @param {string} smiles - The SMILES string
 * @returns {Promise<Object|null>} - 3D molecule data or null if failed
 */
export const generateMolecule3D = async (smiles) => {
  try {
    const rdkit = await initRDKit()
    if (!rdkit) {
      return null
    }
    
    console.log('generateMolecule3D: Starting 3D generation for SMILES:', smiles.substring(0, 50) + '...')
    
    const mol = rdkit.get_mol(smiles)
    if (!mol) {
      console.warn('generateMolecule3D: Failed to create molecule from SMILES')
      return null
    }
    
    console.log('generateMolecule3D: Molecule created successfully')
    
    // Generate 3D coordinates using RDKit
    let mol3D
    try {
      // First add hydrogens to the molecule (important for 3D)
      const molWithH = rdkit.get_mol(smiles, JSON.stringify({ 
        "removeHs": false,
        "sanitize": true 
      }))
      
      if (molWithH) {
        console.log('generateMolecule3D: Added hydrogens to molecule')
        
        // Generate 3D conformer - this is the key step for real 3D coordinates
        try {
          // Method 1: Try addHs() if available
          if (typeof molWithH.addHs === 'function') {
            molWithH.addHs()
            console.log('generateMolecule3D: Added hydrogens using addHs()')
          }
          
          // Method 2: Try generate_aligned_coords() for 3D coordinates
          if (typeof molWithH.generate_aligned_coords === 'function') {
            molWithH.generate_aligned_coords(mol, true, true, false)
            console.log('generateMolecule3D: Generated 3D coordinates using generate_aligned_coords')
          }
          // Method 3: Try get_new_coords() for conformer generation  
          else if (typeof molWithH.get_new_coords === 'function') {
            const coords = molWithH.get_new_coords(true) // true for 3D
            if (coords) {
              console.log('generateMolecule3D: Generated 3D coordinates using get_new_coords')
            }
          }
          // Method 4: Try conformer-based approach
          else if (typeof rdkit.get_mol_copy === 'function') {
            const mol3DCopy = rdkit.get_mol_copy(molWithH)
            if (mol3DCopy) {
              console.log('generateMolecule3D: Created 3D molecule copy')
              molWithH.delete?.()
              mol3D = mol3DCopy
            }
          }
          
          if (!mol3D) {
            mol3D = molWithH
            console.log('generateMolecule3D: Using molecule with hydrogens (attempting 3D coordinates)')
          }
        } catch (coord3DError) {
          console.warn('generateMolecule3D: 3D coordinate generation failed, using molecule with hydrogens:', coord3DError)
          mol3D = molWithH
        }
      } else {
        console.warn('generateMolecule3D: Failed to add hydrogens, using original molecule')
        mol3D = mol
      }
    } catch (error) {
      console.warn('generateMolecule3D: 3D molecule processing failed, using original molecule:', error)
      mol3D = mol
    }
    
    if (!mol3D) {
      mol.delete?.()
      return null
    }
    
    // Generate the molblock (SDF format)
    let molblock
    try {
      molblock = mol3D.get_molblock()
      
      if (molblock) {
        // Update the header to indicate 3D coordinates if we processed them
        if (mol3D !== mol) {
          molblock = molblock.replace('     RDKit          2D', '     RDKit          3D')
          console.log('generateMolecule3D: Updated molblock header to indicate 3D coordinates')
        }
        
        const result = {
          molblock,
          smiles,
          mol: mol3D,
          has3D: mol3D !== mol,
          atomCount: molblock.split('\n')[3]?.trim().split(/\s+/)[0] || 'unknown'
        }
        
        console.log('generateMolecule3D: Successfully generated molecule data')
        return result
      }
    } catch (molblockError) {
      console.error('generateMolecule3D: Failed to generate molblock:', molblockError)
    }
    
    // Clean up
    mol.delete?.()
    mol3D.delete?.()
    
    return null
  } catch (error) {
    console.warn('3D molecule generation error:', error)
    return null
  }
} 

// Basic SMILES validation without RDKit
const basicValidateSMILES = (smiles) => {
  if (!smiles || typeof smiles !== 'string') {
    return { isValid: false, error: 'SMILES must be a non-empty string' }
  }
  
  const trimmed = smiles.trim()
  if (trimmed.length === 0) {
    return { isValid: false, error: 'SMILES cannot be empty' }
  }
  
  // Basic SMILES character validation
  const validChars = /^[A-Za-z0-9()\[\]@+\-=#$:/\\%*.]+$/
  if (!validChars.test(trimmed)) {
    return { isValid: false, error: 'SMILES contains invalid characters' }
  }
  
  // Check for basic structural validity
  if (trimmed.length < 1) {
    return { isValid: false, error: 'SMILES too short' }
  }
  
  // Very basic checks - could be expanded
  const openParens = (trimmed.match(/\(/g) || []).length
  const closeParens = (trimmed.match(/\)/g) || []).length
  const openBrackets = (trimmed.match(/\[/g) || []).length
  const closeBrackets = (trimmed.match(/\]/g) || []).length
  
  if (openParens !== closeParens) {
    return { isValid: false, error: 'Unmatched parentheses in SMILES' }
  }
  
  if (openBrackets !== closeBrackets) {
    return { isValid: false, error: 'Unmatched brackets in SMILES' }
  }
  
  // If it passes basic checks, consider it potentially valid
  return { isValid: true, error: null }
} 