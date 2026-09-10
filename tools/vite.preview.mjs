import { defineConfig } from 'vite';
export default defineConfig({
  root: '_site',
  server: { host: '0.0.0.0', allowedHosts: ['terminal.local'] },
  plugins: [{
    // Review-only route: a real 390px iframe viewport, never copied into _site.
    name: 'narrow-layout-review',
    configureServer(server) {
      server.middlewares.use('/__review/mobile', (_req, res) => {
        res.setHeader('Content-Type', 'text/html; charset=utf-8');
        res.end('<!doctype html><html lang="en"><title>DrFARM narrow layout review</title><body style="margin:0;background:#e9eeee"><p style="margin:12px 24px;font:14px system-ui">390 px viewport · responsive layout review</p><iframe title="Narrow DrFARM preview" src="/index.html" style="display:block;margin:0 24px;width:390px;height:780px;border:0;background:white"></iframe></body></html>');
      });
    },
  }],
});
