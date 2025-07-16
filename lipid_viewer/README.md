# Lipid Viewer

A single-page React 18 + Vite + Tailwind molecular visualization application for interactive 2D/3D rendering of lipids and other molecules from SMILES strings.

## 🧬 Features

### Core Functionality
- **SMILES Input Interface**: Real-time validation with toast error notifications
- **2D/3D Visualization Toggle**: Seamless switching between visualization modes
- **Multiple Renderers**: 
  - **2D**: SmilesDrawer @2.1.1, RDKit.js @2024.09, Kekule.js @0.9.6
  - **3D**: 3Dmol.js @2.4.0, Mol* @3.0.0, NGL @2.0.0
- **In-browser SMILES→3D Conversion**: Via RDKit.js WASM
- **Interactive Controls**: Rotate, zoom, pan with mouse and keyboard
- **Export Functionality**: PNG, SVG, GLB formats
- **Framework Preparation**: Ready for electrostatic surface visualization

### Technical Excellence
- **≥90% Test Coverage**: Enforced across all metrics (branches, functions, lines, statements)
- **TDD Methodology**: Tests written first, comprehensive coverage
- **CI/CD Pipeline**: GitHub Actions with lint → test → build → zip workflow
- **Modern Stack**: React 18, Vite, Tailwind CSS, Jest, Playwright

## 🚀 Quick Start

### Prerequisites
- Node.js 18+ 
- npm 9+

### Installation

```bash
# Clone the repository
git clone <repository-url>
cd lipid-viewer

# Install dependencies
npm install

# Start development server
npm run dev
```

Open [http://localhost:3000](http://localhost:3000) in your browser.

### Example Usage

1. Enter a SMILES string (e.g., `CCO` for ethanol)
2. Choose 2D or 3D visualization mode
3. Select your preferred renderer
4. Click "Visualize Molecule"
5. Use controls to interact and export

## 🧪 Testing

The project maintains ≥90% test coverage across all metrics:

```bash
# Run unit tests with coverage
npm run test

# Run E2E tests
npm run test:e2e

# Run all tests
npm run test:all
```

### Test Coverage Breakdown
- **Unit Tests**: Jest + React Testing Library
- **E2E Tests**: Playwright across Chrome, Firefox, Safari
- **Coverage Threshold**: 90% enforced for branches, functions, lines, statements

## 🏗️ Architecture

### Component Structure
```
src/
├── components/
│   ├── SMILESInput.jsx          # Input validation & examples
│   ├── RendererSelector.jsx     # 2D/3D mode & renderer selection
│   ├── MoleculeViewer.jsx       # Core visualization engine
│   └── ViewControls.jsx         # Zoom, rotate, export controls
├── utils/
│   ├── smilesValidator.js       # RDKit.js SMILES validation
│   └── exportUtils.js           # PNG/SVG/GLB export functionality
└── App.jsx                      # Main orchestrator component
```

### Renderer Integration
- **SmilesDrawer**: Client-side 2D rendering with SVG output
- **RDKit.js**: WASM-based validation and high-quality 2D structures
- **Kekule.js**: Comprehensive cheminformatics with canvas rendering
- **3Dmol.js**: WebGL molecular visualization with interactive 3D
- **Mol***: High-performance structural biology viewer
- **NGL**: Advanced protein and molecular visualization

## 🛠️ Build & Deployment

### Development
```bash
npm run dev          # Start dev server
npm run build        # Production build
npm run preview      # Preview production build
npm run lint         # ESLint checking
```

### Production Deployment
The build creates a static `dist/` folder that can be deployed to any web server:

```bash
npm run build
# Deploy dist/ folder to your hosting platform
```

### GitHub Actions CI/CD
Automated pipeline on push to main:
1. **Lint**: ESLint code quality checks
2. **Test**: Jest unit tests + Playwright E2E tests
3. **Build**: Production Vite build
4. **Package**: Create ZIP with interactive `index.html`

## 📋 Requirements Compliance

### Functional Requirements ✅
- [x] SMILES string input interface
- [x] 2D/3D visualization toggle  
- [x] Multiple renderer support (6 libraries)
- [x] In-browser SMILES→3D conversion via RDKit.js WASM
- [x] Rotate/zoom controls + export PNG/SVG/GLB functionality
- [x] Inline SMILES validation with toast error notifications
- [x] Framework preparation for electrostatic-surface toggle

### Non-functional Requirements ✅
- [x] ≥90% test coverage (Jest + Playwright) with enforced coverageThreshold
- [x] GitHub Actions pipeline: lint → test → build → zip
- [x] Interactive rendering via index.html at project root  
- [x] TDD methodology with tests written first

### Library Versions ✅
All dependencies pinned to exact specified versions:
- SmilesDrawer @2.1.1
- RDKit.js @2024.09
- Kekule.js @0.9.6  
- 3Dmol.js @2.4.0
- Mol* @3.0.0
- NGL @2.0.0

## 🎮 Usage Examples

### Basic SMILES Visualization
```javascript
// Example molecules included:
"CCO"                    // Ethanol
"c1ccccc1"              // Benzene  
"CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  // Caffeine
"CC(=O)OC1=CC=CC=C1C(=O)O"       // Aspirin
```

### Keyboard Shortcuts
- `+` / `-`: Zoom in/out
- `Space`: Reset view
- `Ctrl+2`: Switch to 2D mode
- `Ctrl+3`: Switch to 3D mode

### Export Formats
- **PNG**: Raster image for presentations
- **SVG**: Vector graphics for publications  
- **GLB**: 3D model for printing/AR/VR

## 🔬 Framework Preparation

The application architecture supports future electrostatic surface visualization:

- Modular renderer system allows easy addition of surface rendering
- RDKit.js integration provides molecular property calculation
- 3D libraries (3Dmol.js, Mol*, NGL) support surface representations
- Export system ready for surface data formats

## 📄 License

[MIT License](LICENSE) - see LICENSE file for details.

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/amazing-feature`
3. Follow TDD: Write tests first
4. Ensure ≥90% coverage: `npm run test`
5. Commit changes: `git commit -m 'Add amazing feature'`
6. Push to branch: `git push origin feature/amazing-feature`
7. Open a Pull Request

## 📞 Support

For issues and questions:
- Create an issue in the GitHub repository
- Check existing documentation and examples
- Review test files for usage patterns

---

Built with ❤️ for the molecular visualization community. 