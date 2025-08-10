# Lipid Rendering

A modern React application for visualizing lipid molecules using SMILES notation. This interactive web tool allows users to input molecular structures and view them in various 2D and 3D representations.

## Quick Runbook (avoid common pitfalls)

- Always run commands inside `lipid_viewer` (not the repo root)
- Dev server: runs on http://localhost:3000
- Preview server (production build): runs on http://localhost:4173

Steps:

```bash
# 1) From repository root, cd into the app folder
cd lipid_viewer

# 2) Install dependencies
npm install

# 3) Start the dev server (Vite)
npm run dev
# open http://localhost:3000

# If you prefer serving the production build locally:
npm run build            # creates dist/
npm run preview          # serves dist/ at http://localhost:4173
```

Troubleshooting:
- If you see “Could not read package.json” errors, you are in the wrong folder. Change directory to `lipid_viewer` and retry.
- If `npm run preview` says “dist does not exist,” run `npm run build` first.
- If a port is busy, use: `npm run dev -- --port 3001` or `npm run preview -- --port 4174`.

## Features

- 🧪 **SMILES Input**: Enter molecular structures using SMILES notation
- 🎨 **Multiple Renderers**: Support for different visualization engines
- 🔍 **Interactive Views**: Zoom, rotate, and explore molecular structures
- 📱 **Responsive Design**: Works on desktop and mobile devices
- ⚡ **Real-time Rendering**: Instant visualization updates
- 📊 **Export Options**: Save molecular visualizations

## Demo

The application is accessible at: [https://github.com/Mahmoudimahyar/lipid-rendering](https://github.com/Mahmoudimahyar/lipid-rendering)

## Tech Stack

- **Frontend**: React 18 with Vite
- **Styling**: Tailwind CSS
- **Testing**: Jest & React Testing Library
- **E2E Testing**: Playwright
- **Molecular Visualization**: RDKit integration
- **Build Tool**: Vite

## Getting Started

### Prerequisites

- Node.js (version 18+)
- npm

### Installation

```bash
git clone https://github.com/Mahmoudimahyar/lipid-rendering.git
cd lipid-rendering
cd lipid_viewer
npm install
```

### Development

```bash
npm run dev
# open http://localhost:3000
```

### Production Preview

```bash
npm run build
npm run preview
# open http://localhost:4173
```

## Available Scripts

In the `lipid_viewer` directory, you can run:

- `npm run dev` - Starts the development server on port 3000 (host exposed)
- `npm run build` - Builds the app for production
- `npm run preview` - Preview the production build on port 4173 (host exposed)
- `npm test` - Runs the test suite
- `npm run test:watch` - Runs tests in watch mode
- `npm run test:coverage` - Runs tests with coverage report
- `npm run lint` - Runs ESLint
- `npm run e2e` - Runs end-to-end tests with Playwright

## Project Structure

```
lipid_viewer/
├── src/
│   ├── components/          # React components
│   │   ├── MoleculeViewer.jsx
│   │   ├── SMILESInput.jsx
│   │   ├── RendererSelector.jsx
│   │   └── ViewControls.jsx
│   ├── utils/              # Utility functions
│   │   ├── smilesValidator.js
│   │   └── exportUtils.js
│   ├── App.jsx             # Main application component
│   └── main.jsx           # Application entry point
├── public/                # Static assets
├── dist/                  # Production build output (after npm run build)
└── package.json
```

## Testing

```bash
npm test
```

---

For detailed feature docs, see `lipid_viewer/README.md`.

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- RDKit for molecular visualization capabilities
- React community for excellent tooling and libraries
- Contributors and users who help improve this project

## Contact

Mahyar Mahmoudi - mahmoudimahyar@gmail.com

Project Link: [https://github.com/Mahmoudimahyar/lipid-rendering](https://github.com/Mahmoudimahyar/lipid-rendering) 