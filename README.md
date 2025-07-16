# Lipid Rendering

A modern React application for visualizing lipid molecules using SMILES notation. This interactive web tool allows users to input molecular structures and view them in various 2D and 3D representations.

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

- Node.js (version 16 or higher)
- npm or yarn

### Installation

1. Clone the repository:
```bash
git clone https://github.com/Mahmoudimahyar/lipid-rendering.git
cd lipid-rendering
```

2. Navigate to the application directory:
```bash
cd lipid_viewer
```

3. Install dependencies:
```bash
npm install
```

4. Start the development server:
```bash
npm run dev
```

5. Open your browser and visit `http://localhost:3001`

## Available Scripts

In the `lipid_viewer` directory, you can run:

- `npm run dev` - Starts the development server
- `npm run build` - Builds the app for production
- `npm run preview` - Preview the production build locally
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
├── tests/                  # E2E tests
├── public/                # Static assets
└── package.json
```

## Usage

1. **Enter SMILES**: Input a valid SMILES notation in the text field
2. **Select Renderer**: Choose your preferred visualization engine
3. **Adjust Controls**: Use the view controls to customize the display
4. **Export**: Save your molecular visualization

### Example SMILES

Try these example SMILES notations:
- `CCO` (Ethanol)
- `CC(=O)OC1=CC=CC=C1C(=O)O` (Aspirin)
- `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` (Caffeine)

## Testing

The project includes comprehensive testing:

- **Unit Tests**: Component and utility function tests
- **Integration Tests**: Component interaction tests
- **E2E Tests**: Full user workflow tests

Run all tests:
```bash
npm test
```

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