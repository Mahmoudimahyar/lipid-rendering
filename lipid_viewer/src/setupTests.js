import '@testing-library/jest-dom'

// Mock for matchMedia (required by react-hot-toast)
Object.defineProperty(window, 'matchMedia', {
  writable: true,
  value: jest.fn().mockImplementation(query => ({
    matches: false,
    media: query,
    onchange: null,
    addListener: jest.fn(), // deprecated
    removeListener: jest.fn(), // deprecated
    addEventListener: jest.fn(),
    removeEventListener: jest.fn(),
    dispatchEvent: jest.fn(),
  })),
})

// Mock for molecular visualization libraries
global.ResizeObserver = jest.fn().mockImplementation(() => ({
  observe: jest.fn(),
  unobserve: jest.fn(),
  disconnect: jest.fn(),
}))

// Mock WebGL context
HTMLCanvasElement.prototype.getContext = jest.fn()

// Create canvas mock functions  
const mockCanvasToBlob = jest.fn(function(callback) {
  const blob = new Blob(['mock canvas data'], { type: 'image/png' })
  setTimeout(() => callback(blob), 0)
})

const mockCanvasToDataURL = jest.fn(() => 
  'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg=='
)

const mockCanvasGetContext = jest.fn(() => ({ fillRect: jest.fn(), drawImage: jest.fn(), fillStyle: '' }))

// Mock HTMLCanvasElement constructor
global.HTMLCanvasElement = class HTMLCanvasElement extends global.HTMLCanvasElement {
  constructor() {
    super()
    this.toBlob = mockCanvasToBlob
    this.toDataURL = mockCanvasToDataURL
    this.getContext = mockCanvasGetContext
  }
}

// Also mock the prototype for existing elements
HTMLCanvasElement.prototype.toBlob = mockCanvasToBlob
HTMLCanvasElement.prototype.toDataURL = mockCanvasToDataURL  
HTMLCanvasElement.prototype.getContext = mockCanvasGetContext

// Mock URL methods for file download testing
global.URL.createObjectURL = jest.fn(() => 'blob:http://localhost:3000/mock-url')
global.URL.revokeObjectURL = jest.fn()

// Mock document.createElement for canvas  
const originalCreateElement = document.createElement
document.createElement = jest.fn((tagName) => {
  if (tagName === 'canvas') {
    const canvas = originalCreateElement.call(document, tagName)
    // Override the methods on the actual canvas element
    canvas.toBlob = mockCanvasToBlob
    canvas.toDataURL = mockCanvasToDataURL
    canvas.getContext = mockCanvasGetContext
    return canvas
  }
  if (tagName === 'a') {
    const link = originalCreateElement.call(document, tagName)
    link.click = jest.fn()
    return link
  }
  return originalCreateElement.call(document, tagName)
})

// Mock for 3D libraries that require WebGL
Object.defineProperty(window, 'WebGLRenderingContext', {
  value: jest.fn()
})

// Mock for RDKit WASM (loaded via CDN)
Object.defineProperty(window, 'initRDKitModule', {
  value: jest.fn(() => Promise.resolve({
    get_mol: jest.fn((smiles) => {
      // Return null for invalid SMILES patterns
      const invalidPatterns = ['invalid', '((', '))', '[[[', 'xyz', 'invalid123', '~~~']
      if (!smiles || invalidPatterns.some(pattern => String(smiles).includes(pattern))) {
        return null
      }
      // Return mock molecule object for valid SMILES
      return {
        is_valid: jest.fn(() => true),
        get_svg: jest.fn().mockReturnValue('<svg viewBox="0 0 800 600"><circle r="10"/></svg>'),
        get_molblock: jest.fn().mockReturnValue('     RDKit          2D\n\n 14 15  0  0  0  0  0  0  0  0999 V2000\n'),
        get_json: jest.fn().mockReturnValue('{"atoms":[]}'),
        addHs: jest.fn(),
        generate_aligned_coords: jest.fn(),
        get_new_coords: jest.fn(() => 'coords'),
        delete: jest.fn()
      }
    }),
  }))
})

// Provide global shims for external libraries to avoid module resolution issues
// 3Dmol.js
window.$3Dmol = {
  rasmolElementColors: {},
  createViewer: jest.fn(() => ({
    setBackgroundColor: jest.fn(),
    addModel: jest.fn(),
    setStyle: jest.fn(),
    zoomTo: jest.fn(),
    center: jest.fn(),
    rotate: jest.fn(),
    render: jest.fn(),
    resize: jest.fn(),
    clear: jest.fn(),
    getModel: jest.fn(() => ({ selectedAtoms: jest.fn(() => new Array(10)) }))
  }))
}

// SmilesDrawer
window.SmilesDrawer = {
  Drawer: function() { this.draw = jest.fn() },
  parse: (smiles, callback) => callback && callback({})
}

// Kekule
window.Kekule = {
  ChemWidget: {
    SmilesParser: function() { this.parse = jest.fn(() => ({ success: true })) }
  }
}

// NGL
window.NGL = {
  Stage: jest.fn(() => ({
    loadFile: jest.fn(() => Promise.resolve({ addRepresentation: jest.fn() })),
    dispose: jest.fn()
  }))
}

// Mol* (molstar)
window.molstar = {
  createPlugin: jest.fn(() => Promise.resolve({
    clear: jest.fn(),
    dispose: jest.fn()
  }))
}

// Suppress console warnings during tests
global.console = {
  ...console,
  warn: jest.fn(),
  error: jest.fn(),
} 