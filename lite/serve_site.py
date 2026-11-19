# Preview server for the built JupyterLite site, for use while developing.
#
# It serves the current folder on port 8000, like `python -m http.server`,
# but tells the browser not to keep copies of anything. Without that, a
# browser can go on using an old pyodide-lock.json (the list of packages the
# kernel may load) for up to an hour after the site is rebuilt.
#
# Run from the repository root, after building the site into ./_output:
#   docker run --rm -p 8000:8000 -v $(pwd)/_output:/srv -v $(pwd)/lite/serve_site.py:/serve_site.py:ro -w /srv python:3.13-slim python /serve_site.py
# then open http://localhost:8000/
import http.server


class NoCacheHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # Added to every response: do not store this file for later reuse.
        self.send_header('Cache-Control', 'no-store')
        super().end_headers()


server = http.server.ThreadingHTTPServer(('0.0.0.0', 8000), NoCacheHandler)
print('Serving on http://localhost:8000/ with browser caching switched off')
server.serve_forever()
