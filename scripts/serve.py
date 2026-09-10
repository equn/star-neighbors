"""Dev server for web/ that never serves a stale file, and only to this machine.

`python3 -m http.server` has two problems for this project. It sends
Last-Modified with no cache directive, so a browser may reuse a stale index.html
or app.js after an edit without asking - which looks exactly like a broken
feature. And it binds every interface, so the directory is offered to whatever
network you happen to be on.

This serves the same directory with:
  * `Cache-Control: no-cache`, which permits storage but requires revalidation,
    so edits always land and an unchanged 1.9 MB cheb.bin still costs a 304;
  * a bind to 127.0.0.1 unless you ask otherwise;
  * directory listings off;
  * the response headers a hosted deployment should also set.

    python3 scripts/serve.py [port] [--public]
"""

from __future__ import annotations

import functools
import http.server
import pathlib
import socketserver
import sys

ROOT = pathlib.Path(__file__).parent.parent / "web"

# The page loads nothing from anywhere else: no CDN, no fonts, no analytics, no
# outbound requests at all. Saying so in a header means a script injected into
# the bundle could not phone home either. 'unsafe-inline' for styles is needed
# only because the stylesheet lives in a <style> block in index.html; scripts
# have no such exemption.
CSP = ("default-src 'none'; "
       "script-src 'self'; "
       "style-src 'self' 'unsafe-inline'; "
       "img-src 'self' data:; "
       "connect-src 'self'; "
       "base-uri 'none'; "
       "form-action 'none'; "
       "frame-ancestors 'none'")


class Handler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # `no-cache` rather than `no-store`: the browser may keep the response
        # but must revalidate before reusing it, so an edited file is always
        # picked up while an unchanged 1.9 MB cheb.bin costs a 304 instead of a
        # full refetch. `no-store` forbids keeping it at all, which makes every
        # reload re-download the trajectories for no correctness gain -
        # SimpleHTTPRequestHandler already answers If-Modified-Since.
        self.send_header("Cache-Control", "no-cache, must-revalidate")
        self.send_header("Pragma", "no-cache")
        self.send_header("Content-Security-Policy", CSP)
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        super().end_headers()

    def list_directory(self, path):
        """No directory listings - serve files, not an index of the build."""
        self.send_error(404, "No permission to list directory")
        return None

    def log_message(self, fmt, *args):        # keep the console readable
        if "304" not in fmt % args:
            super().log_message(fmt, *args)


def main() -> None:
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    port = int(args[0]) if args else 8777
    # Loopback by default. A dev server has no authentication and no rate limit,
    # and "" would offer web/ to every machine on the café Wi-Fi.
    host = "0.0.0.0" if "--public" in sys.argv else "127.0.0.1"
    handler = functools.partial(Handler, directory=str(ROOT))
    socketserver.TCPServer.allow_reuse_address = True
    with socketserver.TCPServer((host, port), handler) as httpd:
        where = "all interfaces" if host == "0.0.0.0" else "localhost only"
        print(f"serving {ROOT} at http://localhost:{port}  ({where}, no-cache)")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print()


if __name__ == "__main__":
    main()
