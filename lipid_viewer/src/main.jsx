import React from 'react'
import ReactDOM from 'react-dom/client'
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom'
import App from './App.jsx'
import DockingPage from './pages/DockingPage.jsx'
import './index.css'

console.log('=== main.jsx: Starting React application ===')

ReactDOM.createRoot(document.getElementById('root')).render(
  <React.StrictMode>
    <Router>
      <Routes>
        <Route path="/" element={<App />} />
        <Route path="/dock" element={<DockingPage />} />
      </Routes>
    </Router>
  </React.StrictMode>,
)

console.log('=== main.jsx: React application rendered ===')