import { defineConfig } from 'vitest/config'
import react from '@vitejs/plugin-react'
import { loadEnv } from 'vite'

// Permit only the exact tunnel hostname supplied by the researcher. Never use
// allowedHosts: true (which would expose dev source through DNS rebinding).
export default defineConfig(({ mode }) => {
  const tunnelHost = loadEnv(mode, '.', 'GENE2BRAIN_').GENE2BRAIN_TUNNEL_HOST
  const allowedHosts = tunnelHost ? [tunnelHost] : []
  return {
    plugins: [react()],
    base: '/',
    server: { host: '127.0.0.1', port: 8000, strictPort: true, allowedHosts },
    preview: { allowedHosts },
    test: {
      environment: 'jsdom',
      setupFiles: './src/test/setup.ts',
      include: ['src/**/*.test.ts', 'src/**/*.test.tsx'],
    },
  }
})
