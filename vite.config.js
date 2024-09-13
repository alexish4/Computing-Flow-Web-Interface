import { defineConfig } from 'vite';

export default defineConfig({
  root: './static',
  build: {
    outDir: '../static/dist',
    emptyOutDir: true,
  },
});
