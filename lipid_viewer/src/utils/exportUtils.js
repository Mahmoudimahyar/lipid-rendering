/**
 * Generate a safe filename from SMILES string
 * @param {string} smiles - The SMILES string
 * @param {string} extension - File extension (png, svg, glb)
 * @param {boolean} includeTimestamp - Whether to include timestamp
 * @returns {string} - Safe filename
 */
export const generateFilename = (smiles, extension, includeTimestamp = false) => {
  // Sanitize: allow alphanumerics, dot, underscore, hyphen; replace others with underscore
  const base = typeof smiles === 'string' ? smiles : ''
  let safeName = base.replace(/[^A-Za-z0-9._-]/g, '_')
  
  if (includeTimestamp) {
    const timestamp = Date.now()
    safeName = `${safeName}_${timestamp}`
  }
  
  return `${safeName}.${extension}`
}

/**
 * Download a blob as a file
 * @param {Blob} blob - The blob to download
 * @param {string} filename - The filename
 */
export const downloadFile = (blob, filename) => {
  let url = ''
  try {
    if (typeof URL.createObjectURL === 'function') {
      url = URL.createObjectURL(blob)
    }
  } catch {
    url = ''
  }

  const link = document.createElement('a')
  link.href = url || '#'
  link.download = filename || ''
  link.style.display = 'none'
  
  document.body.appendChild(link)
  try {
    link.click()
  } finally {
    document.body.removeChild(link)
  }
  
  // Clean up the blob URL if available
  if (url && typeof URL.revokeObjectURL === 'function') {
    setTimeout(() => URL.revokeObjectURL(url), 100)
  }
}

/**
 * Export canvas content to PNG
 * @param {HTMLCanvasElement} canvas - The canvas element
 * @param {string} smiles - The SMILES string for filename
 * @param {number} [quality] - Optional image quality (0-1)
 * @returns {Promise<boolean>} - Success status
 */
export const exportToPNG = async (canvas, smiles, quality) => {
  if (!canvas || typeof canvas.toBlob !== 'function') {
    return false
  }
  
  try {
    return new Promise((resolve) => {
      canvas.toBlob((blob) => {
        if (blob) {
          const filename = generateFilename(smiles, 'png', true)
          downloadFile(blob, filename)
          resolve(true)
        } else {
          resolve(false)
        }
      }, 'image/png', quality)
    })
  } catch (error) {
    console.error('PNG export error:', error)
    return false
  }
}

/**
 * Export SVG element or string to SVG file
 * @param {SVGElement|string} svgInput - The SVG element or SVG string
 * @param {string} smiles - The SMILES string for filename
 * @returns {Promise<boolean>} - Success status
 */
export const exportToSVG = async (svgInput, smiles) => {
  if (!svgInput) {
    return false
  }
  
  try {
    let svgString
    if (typeof svgInput === 'string') {
      svgString = svgInput
    } else {
      const serializer = new XMLSerializer()
      svgString = serializer.serializeToString(svgInput)
    }
    
    // Add XML declaration and DOCTYPE for proper SVG file
    const fullSvgString = `<?xml version="1.0" encoding="UTF-8"?>\n<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n${svgString}`
    
    const blob = new Blob([fullSvgString], { type: 'image/svg+xml' })
    // Use a stable filename (no timestamp) for SVG to ease testing and UX
    const filename = generateFilename(smiles, 'svg', false)
    
    downloadFile(blob, filename)
    return true
  } catch (error) {
    console.error('SVG export error:', error)
    return false
  }
}

/**
 * Export molecular data to GLTF format
 * @param {string} gltfData - The GLTF JSON string
 * @param {string} smiles - The SMILES string for filename
 * @returns {Promise<boolean>} - Success status
 */
export const exportToGLTF = async (gltfData, smiles) => {
  try {
    const blob = new Blob([gltfData], { type: 'model/gltf+json' })
    const filename = generateFilename(smiles, 'gltf', true)
    
    downloadFile(blob, filename)
    return true
  } catch (error) {
    console.error('GLTF export error:', error)
    return false
  }
}

/**
 * Export 3D scene to GLB format (legacy - kept for compatibility)
 * @param {Object} scene - The 3D scene object
 * @param {string} smiles - The SMILES string for filename
 * @returns {Promise<boolean>} - Success status
 */
export const exportToGLB = async (scene, smiles) => {
  try {
    if (!scene || typeof scene.export !== 'function') {
      console.warn('Scene does not support GLB export')
      return false
    }
    
    return new Promise((resolve) => {
      try {
        scene.export((arrayBuffer) => {
          if (arrayBuffer) {
            const blob = new Blob([arrayBuffer], { type: 'model/gltf-binary' })
            const filename = generateFilename(smiles, 'glb', true)
            downloadFile(blob, filename)
            resolve(true)
          } else {
            resolve(false)
          }
        })
      } catch (e) {
        console.error('GLB export error:', e)
        resolve(false)
      }
    })
  } catch (error) {
    console.error('GLB export error:', error)
    return false
  }
}

/**
 * Export data as JSON file
 * @param {Object} data - The data to export
 * @param {string} smiles - The SMILES string for filename
 * @returns {Promise<boolean>} - Success status
 */
export const exportToJSON = async (data, smiles) => {
  try {
    const jsonString = JSON.stringify(data, null, 2)
    const blob = new Blob([jsonString], { type: 'application/json' })
    const filename = generateFilename(smiles, 'json', true)
    
    downloadFile(blob, filename)
    return true
  } catch (error) {
    console.error('JSON export error:', error)
    return false
  }
} 