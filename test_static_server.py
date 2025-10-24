#!/usr/bin/env python3
"""
Simple HTTP server for testing static version locally
"""
import http.server
import socketserver
import os
from pathlib import Path

# Configuration
PORT = 8000
DIRECTORY = Path(__file__).parent.absolute()


class CustomHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=DIRECTORY, **kwargs)

    def end_headers(self):
        # Add CORS headers for local testing
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        super().end_headers()


def main():
    print("="*60)
    print("ClaimSpy UI - Static Version Test Server")
    print("="*60)
    print(f"Serving directory: {DIRECTORY}")
    print(f"Server: http://localhost:{PORT}")
    print(f"Interface: http://localhost:{PORT}/evaluation-interface/")
    print("="*60)
    print("Press Ctrl+C to stop the server")

    with socketserver.TCPServer(("", PORT), CustomHTTPRequestHandler) as httpd:
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\n\n🛑 Server stopped")


if __name__ == "__main__":
    main()
