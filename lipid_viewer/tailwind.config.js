/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        'molecule-bg': '#1e293b',
        'molecule-surface': '#334155',
      }
    },
  },
  plugins: [],
} 