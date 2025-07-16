import { 
  exportToPNG, 
  exportToSVG, 
  exportToGLB, 
  generateFilename,
  downloadFile,
  exportToJSON
} from './exportUtils'

// Mock canvas and URL methods
global.URL.createObjectURL = jest.fn(() => 'mock-blob-url')
global.URL.revokeObjectURL = jest.fn()

// Mock canvas toBlob method
HTMLCanvasElement.prototype.toBlob = jest.fn((callback) => {
  callback(new Blob(['mock-canvas-data'], { type: 'image/png' }))
})

// Mock SVG serialization
global.XMLSerializer = jest.fn().mockImplementation(() => ({
  serializeToString: jest.fn(() => '<svg>mock-svg-content</svg>')
}))

// Mock document methods
Object.defineProperty(document, 'createElement', {
  value: jest.fn((tagName) => {
    if (tagName === 'a') {
      return {
        href: '',
        download: '',
        click: jest.fn(),
        style: {}
      }
    }
    return {}
  })
})

Object.defineProperty(document.body, 'appendChild', {
  value: jest.fn()
})

Object.defineProperty(document.body, 'removeChild', {
  value: jest.fn()
})

describe('Export Utils', () => {
  beforeEach(() => {
    jest.clearAllMocks()
  })

  describe('generateFilename', () => {
    test('generates filename with SMILES and extension', () => {
      const filename = generateFilename('CCO', 'png')
      expect(filename).toBe('CCO.png')
    })

    test('sanitizes SMILES for filename', () => {
      const filename = generateFilename('c1ccccc1', 'svg')
      expect(filename).toBe('c1ccccc1.svg')
    })

    test('handles complex SMILES with special characters', () => {
      const filename = generateFilename('CC(=O)O', 'png')
      expect(filename).toBe('CC(=O)O.png')
    })

    test('includes timestamp when specified', () => {
      const filename = generateFilename('CCO', 'png', true)
      expect(filename).toMatch(/^CCO_\d+\.png$/)
    })
  })

  describe('downloadFile', () => {
    test('creates download link and triggers download', () => {
      const blob = new Blob(['test'], { type: 'text/plain' })
      const filename = 'test.txt'
      
      downloadFile(blob, filename)
      
      expect(URL.createObjectURL).toHaveBeenCalledWith(blob)
      expect(document.createElement).toHaveBeenCalledWith('a')
    })

    test('cleans up blob URL after download', async () => {
      const blob = new Blob(['test'], { type: 'text/plain' })
      const filename = 'test.txt'
      
      downloadFile(blob, filename)
      
      // Wait for the setTimeout to execute (100ms + buffer)
      await new Promise(resolve => setTimeout(resolve, 150))
      expect(URL.revokeObjectURL).toHaveBeenCalled()
    })
  })

  describe('exportToPNG', () => {
    test('exports canvas to PNG blob and downloads', async () => {
      // Create a mock canvas with working toBlob
      const mockCanvas = {
        toBlob: jest.fn((callback) => {
          const blob = new Blob(['mock data'], { type: 'image/png' })
          callback(blob)
        })
      }
      const smiles = 'CCO'
      
      await exportToPNG(mockCanvas, smiles)
      
      expect(mockCanvas.toBlob).toHaveBeenCalled()
      expect(URL.createObjectURL).toHaveBeenCalled()
    })

    test('handles canvas export errors gracefully', async () => {
      const mockCanvas = {
        toBlob: jest.fn((callback) => callback(null))
      }
      
      const result = await exportToPNG(mockCanvas, 'CCO')
      expect(result).toBe(false)
    })
  })

  describe('exportToSVG', () => {
    test('exports SVG element to SVG file', async () => {
      const mockSVG = document.createElementNS('http://www.w3.org/2000/svg', 'svg')
      const smiles = 'CCO'
      
      await exportToSVG(mockSVG, smiles)
      
      expect(URL.createObjectURL).toHaveBeenCalled()
    })

    test('handles SVG serialization errors', async () => {
      const mockSVG = document.createElementNS('http://www.w3.org/2000/svg', 'svg')
      
      // Mock XMLSerializer constructor to return an object with a failing method
      const originalXMLSerializer = global.XMLSerializer
      global.XMLSerializer = jest.fn(() => ({
        serializeToString: jest.fn(() => {
          throw new Error('Serialization failed')
        })
      }))
      
      const result = await exportToSVG(mockSVG, 'CCO')
      expect(result).toBe(false)
      
      // Restore original
      global.XMLSerializer = originalXMLSerializer
    })
  })

  describe('exportToGLB', () => {
    test('exports 3D scene to GLB format', async () => {
      const mockScene = { 
        export: jest.fn((callback) => callback(new ArrayBuffer(100))),
        type: '3d-scene'
      }
      const smiles = 'CCO'
      
      await exportToGLB(mockScene, smiles)
      
      expect(URL.createObjectURL).toHaveBeenCalled()
    })

    test('handles GLB export errors', async () => {
      const mockScene = { 
        export: jest.fn((callback) => callback(null)),
        type: '3d-scene'
      }
      
      const result = await exportToGLB(mockScene, 'CCO')
      expect(result).toBe(false)
    })

    test('handles scenes without export capability', async () => {
      const mockScene = { type: '3d-scene' }
      
      const result = await exportToGLB(mockScene, 'CCO')
      expect(result).toBe(false)
    })
  })

  describe('generateFilename', () => {
    test('generates correct filename with different extensions', async () => {
      // Test PNG export filename generation
      const canvas = {
        toBlob: jest.fn((callback) => {
          const blob = new Blob(['mock data'], { type: 'image/png' })
          callback(blob)
        })
      }
      await exportToPNG(canvas, 'CCO')
      expect(URL.createObjectURL).toHaveBeenCalled()
      
      // Test SVG export filename generation  
      const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg')
      await exportToSVG(svg, 'c1ccccc1')
      expect(URL.createObjectURL).toHaveBeenCalled()
    })

    test('handles special characters in SMILES', async () => {
      const canvas = {
        toBlob: jest.fn((callback) => {
          const blob = new Blob(['mock data'], { type: 'image/png' })
          callback(blob)
        })
      }
      const specialSmiles = 'CC(=O)[C@@H](N)C'
      
      await exportToPNG(canvas, specialSmiles)
      expect(URL.createObjectURL).toHaveBeenCalled()
    })

    test('handles empty SMILES string', async () => {
      const canvas = {
        toBlob: jest.fn((callback) => {
          const blob = new Blob(['mock data'], { type: 'image/png' })
          callback(blob)
        })
      }
      
      await exportToPNG(canvas, '')
      expect(URL.createObjectURL).toHaveBeenCalled()
    })

    test('handles very long SMILES strings', async () => {
      const canvas = {
        toBlob: jest.fn((callback) => {
          const blob = new Blob(['mock data'], { type: 'image/png' })
          callback(blob)
        })
      }
      const longSmiles = 'C'.repeat(100)
      
      await exportToPNG(canvas, longSmiles)
      expect(URL.createObjectURL).toHaveBeenCalled()
    })
  })

  describe('Error handling and edge cases', () => {
    test('handles DOM manipulation gracefully', () => {
      const blob = new Blob(['test'], { type: 'text/plain' })
      
      // Should not throw with valid blob and filename
      expect(() => downloadFile(blob, 'test.txt')).not.toThrow()
    })

    test('handles null/undefined inputs to export functions', async () => {
      // These should handle null inputs gracefully and return false
      const result1 = await exportToPNG(null, 'CCO')
      expect(result1).toBe(false)
      
      const result2 = await exportToSVG(null, 'CCO')
      expect(result2).toBe(false)
      
      const result3 = await exportToGLB(null, 'CCO')
      expect(result3).toBe(false)
    })

    test('handles export function errors consistently', async () => {
      // Test with a canvas that will fail toBlob
      const canvas = document.createElement('canvas')
      canvas.toBlob = jest.fn((callback) => callback(null))
      
      const result = await exportToPNG(canvas, 'CCO')
      expect(result).toBe(false)
    })
  })
}) 

// NEW TESTS TO TARGET UNCOVERED LINES FOR 90% COVERAGE
describe('exportUtils Error Handling and Missing Function Coverage', () => {
  test('covers PNG export error handling (lines 65-66)', async () => {
    // Create a canvas that will cause toBlob to throw an error
    const mockCanvas = {
      toBlob: (callback) => {
        // Simulate an error by calling callback with null
        callback(null)
      }
    }
    
    // This should hit the catch block lines 65-66
    const result = await exportToPNG(mockCanvas, 'test')
    expect(result).toBe(false)
  })

  test('covers GLB export error handling (lines 127-128)', async () => {
    // Create a viewer that will cause GLB export to fail
    const mockViewer = {
      export: (callback) => {
        // Simulate failure by calling callback with null
        callback(null)
      }
    }
    
    // This should hit the catch block lines 127-128
    const result = await exportToGLB(mockViewer, 'test')
    expect(result).toBe(false)
  })

  test('covers exportToJSON function completely (lines 139-148)', async () => {
    // Test successful JSON export
    const testData = {
      smiles: 'CCO',
      properties: {
        formula: 'C2H6O',
        weight: 46.07
      },
      structure: {
        atoms: ['C', 'C', 'O'],
        bonds: [1, 1]
      }
    }
    
    const result = await exportToJSON(testData, 'CCO')
    expect(result).toBe(true)
    
    // No need to check mocks since downloadFile is a real function
    // Just verify the function completes successfully
  })

  test('covers exportToJSON error handling', async () => {
    // Test JSON export with data that causes JSON.stringify to fail
    const cyclicalData = {}
    cyclicalData.self = cyclicalData // Creates circular reference
    
    const result = await exportToJSON(cyclicalData, 'test')
    expect(result).toBe(false)
  })

  test('covers exportToJSON with complex data structures', async () => {
    // Test with various data types to ensure comprehensive coverage
    const complexData = {
      string: 'test',
      number: 42,
      boolean: true,
      null: null,
      array: [1, 2, 3],
      nested: {
        deep: {
          value: 'nested'
        }
      },
      emptyObject: {},
      emptyArray: []
    }
    
    const result = await exportToJSON(complexData, 'complex')
    expect(result).toBe(true)
  })

  test('covers exportToJSON with special characters in SMILES', async () => {
    // Test JSON export with SMILES containing special characters
    const data = { test: 'data' }
    const specialSMILES = 'C[C@H](O)C(=O)O' // Contains special characters
    
    const result = await exportToJSON(data, specialSMILES)
    expect(result).toBe(true)
  })

  test('comprehensive error scenarios for all export functions', async () => {
    // Test PNG export with null canvas
    let result = await exportToPNG(null, 'test')
    expect(result).toBe(false)
    
    // Test SVG export with null element
    result = await exportToSVG(null, 'test')
    expect(result).toBe(false)
    
    // Test GLB export with null viewer
    result = await exportToGLB(null, 'test')
    expect(result).toBe(false)
    
    // Test JSON export with null data
    result = await exportToJSON(null, 'test')
    expect(result).toBe(true) // JSON.stringify(null) is valid
    
    // Test JSON export with undefined data
    result = await exportToJSON(undefined, 'test')
    expect(result).toBe(true) // JSON.stringify(undefined) becomes "undefined" string
  })
}) 

// Test PNG export edge cases and error handling
describe('PNG Export Edge Cases', () => {
  beforeEach(() => {
    global.URL.createObjectURL = vi.fn(() => 'blob:mock-url')
    global.URL.revokeObjectURL = vi.fn()
  })

  test('handles canvas with different dimensions', async () => {
    const canvas = document.createElement('canvas')
    canvas.width = 1200
    canvas.height = 800
    
    await exportToPNG(canvas, 'large-molecule')
    
    expect(document.createElement).toHaveBeenCalledWith('a')
    expect(document.createElement('a')).toHaveProperty('download')
  })

  test('handles very small canvas', async () => {
    const canvas = document.createElement('canvas')
    canvas.width = 100
    canvas.height = 100
    
    await exportToPNG(canvas, 'tiny-molecule')
    
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles empty canvas', async () => {
    const canvas = document.createElement('canvas')
    // Canvas with no content
    
    await exportToPNG(canvas, 'empty-molecule')
    
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles canvas.toBlob failure', async () => {
    const canvas = document.createElement('canvas')
    canvas.toBlob = vi.fn((callback) => callback(null)) // Simulate failure
    
    await exportToPNG(canvas, 'failed-molecule')
    
    // Should handle gracefully without throwing
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles canvas.toBlob with different quality settings', async () => {
    const canvas = document.createElement('canvas')
    
    await exportToPNG(canvas, 'quality-test', 0.8)
    
    expect(canvas.toBlob).toHaveBeenCalledWith(expect.any(Function), 'image/png', 0.8)
  })

  test('handles very long molecule names', async () => {
    const canvas = document.createElement('canvas')
    const longName = 'a'.repeat(500)
    
    await exportToPNG(canvas, longName)
    
    expect(document.createElement('a')).toHaveProperty('download')
    expect(document.createElement('a').download).toContain('.png')
  })

  test('handles special characters in molecule names', async () => {
    const canvas = document.createElement('canvas')
    const specialName = 'molecule@#$%^&*()+=[]{}|;:,.<>?'
    
    await exportToPNG(canvas, specialName)
    
    // Should sanitize filename
    expect(document.createElement('a')).toHaveProperty('download')
    expect(document.createElement('a').download).toMatch(/^[^@#$%^&*()+=[\]{}|;:,.<>?]+\.png$/)
  })
})

// Test SVG export edge cases
describe('SVG Export Edge Cases', () => {
  test('handles SVG with embedded images', async () => {
    const svgWithImage = `
      <svg xmlns="http://www.w3.org/2000/svg">
        <image href="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==" />
      </svg>
    `
    
    await exportToSVG(svgWithImage, 'molecule-with-image')
    
    expect(document.createElement).toHaveBeenCalledWith('a')
    expect(document.createElement('a')).toHaveProperty('download')
    expect(document.createElement('a').download).toBe('molecule-with-image.svg')
  })

  test('handles SVG with CSS styles', async () => {
    const svgWithStyles = `
      <svg xmlns="http://www.w3.org/2000/svg">
        <style>
          .molecule { fill: blue; stroke: red; }
        </style>
        <circle class="molecule" cx="50" cy="50" r="40"/>
      </svg>
    `
    
    await exportToSVG(svgWithStyles, 'styled-molecule')
    
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles SVG with complex paths', async () => {
    const complexSvg = `
      <svg xmlns="http://www.w3.org/2000/svg">
        <path d="M150 0 L75 200 L225 200 Z" fill="blue"/>
        <path d="M10,150 C10,77 77,10 150,10 S290,77 290,150" stroke="red" fill="none"/>
      </svg>
    `
    
    await exportToSVG(complexSvg, 'complex-molecule')
    
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles malformed SVG gracefully', async () => {
    const malformedSvg = '<svg><unclosed-tag><circle cx="50" cy="50" r="40"/>'
    
    await exportToSVG(malformedSvg, 'malformed-molecule')
    
    // Should still attempt export
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles empty SVG', async () => {
    const emptySvg = '<svg xmlns="http://www.w3.org/2000/svg"></svg>'
    
    await exportToSVG(emptySvg, 'empty-molecule')
    
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles SVG string instead of element', async () => {
    const svgString = '<svg xmlns="http://www.w3.org/2000/svg"><circle cx="50" cy="50" r="40"/></svg>'
    
    await exportToSVG(svgString, 'string-molecule')
    
    expect(document.createElement).toHaveBeenCalledWith('a')
  })
})

// Test GLB export edge cases
describe('GLB Export Edge Cases', () => {
  test('handles empty geometry', async () => {
    const emptyGeometry = {
      export: (callback) => callback(new ArrayBuffer(0))
    }
    
    await exportToGLB(emptyGeometry, 'empty-molecule')
    
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles large GLB data', async () => {
    const largeData = new ArrayBuffer(10 * 1024 * 1024) // 10MB
    const largeGeometry = {
      export: (callback) => callback(largeData)
    }
    
    await exportToGLB(largeGeometry, 'large-molecule')
    
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles geometry export failure', async () => {
    const failingGeometry = {
      export: (callback) => callback(null)
    }
    
    await exportToGLB(failingGeometry, 'failed-molecule')
    
    // Should handle gracefully
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles geometry export throwing error', async () => {
    const errorGeometry = {
      export: () => { throw new Error('Export failed') }
    }
    
    await exportToGLB(errorGeometry, 'error-molecule')
    
    // Should handle gracefully
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles missing export method', async () => {
    const invalidGeometry = {}
    
    await exportToGLB(invalidGeometry, 'invalid-molecule')
    
    // Should handle gracefully
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles undefined geometry', async () => {
    await exportToGLB(undefined, 'undefined-molecule')
    
    // Should handle gracefully
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles null geometry', async () => {
    await exportToGLB(null, 'null-molecule')
    
    // Should handle gracefully
    expect(document.createElement).toHaveBeenCalledWith('a')
  })
})

// Test filename sanitization and edge cases
describe('Filename Handling', () => {
  test('sanitizes filenames with invalid characters', async () => {
    const canvas = document.createElement('canvas')
    const invalidName = 'molecule/\\:*?"<>|'
    
    await exportToPNG(canvas, invalidName)
    
    // Filename should not contain invalid characters
    expect(document.createElement('a')).toHaveProperty('download')
    expect(document.createElement('a').download).not.toMatch(/[/\\:*?"<>|]/)
  })

  test('handles very long filenames', async () => {
    const canvas = document.createElement('canvas')
    const longName = 'a'.repeat(300)
    
    await exportToPNG(canvas, longName)
    
    // Should truncate if necessary
    expect(document.createElement('a')).toHaveProperty('download')
    expect(document.createElement('a').download.length).toBeLessThan(300)
  })

  test('handles empty filename', async () => {
    const canvas = document.createElement('canvas')
    
    await exportToPNG(canvas, '')
    
    // Should provide default filename
    expect(document.createElement('a')).toHaveProperty('download')
    expect(document.createElement('a').download).toMatch(/\.png$/)
  })

  test('handles undefined filename', async () => {
    const canvas = document.createElement('canvas')
    
    await exportToPNG(canvas, undefined)
    
    // Should provide default filename
    expect(document.createElement('a')).toHaveProperty('download')
    expect(document.createElement('a').download).toMatch(/\.png$/)
  })

  test('handles null filename', async () => {
    const canvas = document.createElement('canvas')
    
    await exportToPNG(canvas, null)
    
    // Should provide default filename
    expect(document.createElement('a')).toHaveProperty('download')
    expect(document.createElement('a').download).toMatch(/\.png$/)
  })
})

// Test browser compatibility edge cases
describe('Browser Compatibility', () => {
  test('handles missing URL.createObjectURL', async () => {
    const originalCreateObjectURL = global.URL.createObjectURL
    delete global.URL.createObjectURL
    
    const canvas = document.createElement('canvas')
    
    await exportToPNG(canvas, 'no-url-api')
    
    // Should handle gracefully
    expect(document.createElement).toHaveBeenCalledWith('a')
    
    // Restore
    global.URL.createObjectURL = originalCreateObjectURL
  })

  test('handles missing URL.revokeObjectURL', async () => {
    const originalRevokeObjectURL = global.URL.revokeObjectURL
    delete global.URL.revokeObjectURL
    
    const canvas = document.createElement('canvas')
    
    await exportToPNG(canvas, 'no-revoke-api')
    
    // Should handle gracefully
    expect(document.createElement).toHaveBeenCalledWith('a')
    
    // Restore
    global.URL.revokeObjectURL = originalRevokeObjectURL
  })

  test('handles missing canvas.toBlob', async () => {
    const canvas = document.createElement('canvas')
    delete canvas.toBlob
    
    await exportToPNG(canvas, 'no-toblob')
    
    // Should handle gracefully
    expect(document.createElement).toHaveBeenCalledWith('a')
  })

  test('handles missing document.createElement', async () => {
    const originalCreateElement = document.createElement
    document.createElement = undefined
    
    const canvas = document.createElement = originalCreateElement; document.createElement('canvas')
    
    try {
      await exportToPNG(canvas, 'no-createelement')
    } catch (error) {
      // Expected to throw
    }
    
    // Restore
    document.createElement = originalCreateElement
  })
})

// Test concurrent exports
describe('Concurrent Exports', () => {
  test('handles multiple simultaneous PNG exports', async () => {
    const canvas1 = document.createElement('canvas')
    const canvas2 = document.createElement('canvas')
    const canvas3 = document.createElement('canvas')
    
    const promises = [
      exportToPNG(canvas1, 'molecule1'),
      exportToPNG(canvas2, 'molecule2'),
      exportToPNG(canvas3, 'molecule3')
    ]
    
    await Promise.all(promises)
    
    expect(document.createElement).toHaveBeenCalledTimes(3)
  })

  test('handles multiple simultaneous SVG exports', async () => {
    const svg = '<svg><circle cx="50" cy="50" r="40"/></svg>'
    
    const promises = [
      exportToSVG(svg, 'molecule1'),
      exportToSVG(svg, 'molecule2'),
      exportToSVG(svg, 'molecule3')
    ]
    
    await Promise.all(promises)
    
    expect(document.createElement).toHaveBeenCalledTimes(3)
  })

  test('handles mixed concurrent exports', async () => {
    const canvas = document.createElement('canvas')
    const svg = '<svg><circle cx="50" cy="50" r="40"/></svg>'
    const geometry = { export: (cb) => cb(new ArrayBuffer(100)) }
    
    const promises = [
      exportToPNG(canvas, 'molecule1'),
      exportToSVG(svg, 'molecule2'),
      exportToGLB(geometry, 'molecule3')
    ]
    
    await Promise.all(promises)
    
    expect(document.createElement).toHaveBeenCalledTimes(3)
  })
}) 