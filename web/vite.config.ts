import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  base: '/Hackthon_tahoe_exp/',
  server: {
    proxy: {
      '/api': 'http://localhost:8000',
    },
  },
})
