/**
 * Minimal static file server, so IonTrace can be opened locally.
 *
 * ES modules are blocked over file:// by the browser's origin rules, so the
 * site needs to come over HTTP even when everything is on disk. This uses only
 * Node's standard library - the project has no dependencies and this file does
 * not add any.
 *
 *     node serve.js [port]
 */

import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { extname, join, normalize, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = fileURLToPath(new URL('.', import.meta.url));
const PORT = Number(process.argv[2]) || 8080;

const TYPES = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.md': 'text/markdown; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.png': 'image/png',
  '.ico': 'image/x-icon',
};

const server = createServer(async (req, res) => {
  try {
    const url = new URL(req.url, `http://localhost:${PORT}`);
    let rel = decodeURIComponent(url.pathname);
    if (rel.endsWith('/')) rel += 'index.html';

    // Resolve inside ROOT and refuse anything that escapes it, so a crafted
    // path cannot read files from elsewhere on the machine.
    const path = normalize(join(ROOT, rel));
    if (!path.startsWith(ROOT.endsWith(sep) ? ROOT : ROOT + sep)) {
      res.writeHead(403).end('Forbidden');
      return;
    }

    const body = await readFile(path);
    res.writeHead(200, {
      'Content-Type': TYPES[extname(path).toLowerCase()] ?? 'application/octet-stream',
      'Cache-Control': 'no-cache',
    });
    res.end(body);
  } catch (err) {
    const code = err.code === 'ENOENT' || err.code === 'EISDIR' ? 404 : 500;
    res.writeHead(code).end(code === 404 ? 'Not found' : 'Server error');
  }
});

server.listen(PORT, () => {
  console.log(`IonTrace serving ${ROOT}`);
  console.log(`  http://localhost:${PORT}/`);
  console.log(`  http://localhost:${PORT}/tests/   (physics tests)`);
});
