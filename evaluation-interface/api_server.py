#!/usr/bin/env python3
"""
Scientific Claim Assessment Interface API Server

Provides REST API endpoints for:
- Loading assessment result folders
- Serving JSON data files (assessment, trajectory, mapping, predictions)
- Serving evidence files (PDF, images, text files)
- Serving static files (HTML, CSS, JS) with modular architecture support
- CORS handling for cross-origin requests

Supports the modular frontend architecture with separated HTML/CSS/JS files.
"""

import json
import os
import sys
from pathlib import Path
from http.server import HTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs
import mimetypes
import logging
from datetime import datetime

# Configuration
# Use relative paths based on the script location
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_BASE_PATH = PROJECT_ROOT / 'data' / \
    'sonnet45-full-run-codex-gpt5-1010' / 'output'
SANDBOX_BASE_PATH = PROJECT_ROOT / 'data' / \
    'sonnet45-full-run-codex-gpt5-1010'
STATIC_BASE_PATH = SCRIPT_DIR
PORT = 8080

# Allowed data files for security
ALLOWED_DATA_FILES = {
    'assessment.json',
    'trajectory.json',
    'mapping.json',
    'likert_score_prediction.json',
    'evidence_source.json'
}

# Allowed simulation file extensions
ALLOWED_SIMULATION_EXTENSIONS = {'.csv', '.json', '.py', '.png'}

# Allowed evidence file extensions
ALLOWED_EVIDENCE_EXTENSIONS = {'.pdf', '.png', '.jpg', '.jpeg', '.txt', '.csv'}

# PDF.js viewer path
PDFJS_VIEWER_PATH = STATIC_BASE_PATH / 'pdfjs' / 'web' / 'viewer.html'

# Sandbox directory mapping
SANDBOX_MAPPING = {
    'dry-run': 'sandbox-dry-run',
    'drop1': 'sandbox-drop1',
    'drop2': 'sandbox-drop2',
    'drop4': 'sandbox-drop4',
    'continuous-release': 'sandbox-continuous-release'
}

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


class APIHandler(BaseHTTPRequestHandler):
    def log_message(self, format, *args):
        """Override to use logging instead of stderr"""
        logging.info(f"{self.address_string()} - {format % args}")

    def send_cors_headers(self):
        """Send common CORS headers"""
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers',
                         'Content-Type, Authorization')
        self.send_header('Cache-Control', 'no-cache')

    def do_GET(self):
        parsed_url = urlparse(self.path)
        path = parsed_url.path
        query = parse_qs(parsed_url.query)

        try:
            if parsed_url.path == '/api/folders':
                self.handle_get_folders()
            elif parsed_url.path.startswith('/api/data/'):
                # Check if this is an evidence file request
                if '/evidences/' in parsed_url.path:
                    self.handle_get_evidence_file(parsed_url.path)
                else:
                    self.handle_get_data(parsed_url.path)
            elif parsed_url.path.startswith('/api/simulation-files/'):
                self.handle_get_simulation_files(parsed_url.path)
            elif parsed_url.path.startswith('/api/simulation-file/'):
                self.handle_get_simulation_file(parsed_url.path)
            else:
                # Serve static files
                self.serve_static_file(parsed_url.path)
        except Exception as e:
            logging.error(f"Error handling request {path}: {e}")
            self.send_error(500, f"Internal server error: {e}")

    def do_OPTIONS(self):
        """Handle preflight CORS requests"""
        self.send_response(200)
        self.send_cors_headers()
        self.end_headers()

    def handle_health_check(self):
        """Simple health check endpoint"""
        self.send_response(200)
        self.send_header('Content-Type', 'application/json')
        self.send_cors_headers()
        self.end_headers()

        response = {
            'status': 'healthy',
            'timestamp': datetime.now().isoformat(),
            'data_path_exists': Path(DATA_BASE_PATH).exists(),
            'static_path_exists': Path(STATIC_BASE_PATH).exists()
        }
        self.wfile.write(json.dumps(response).encode('utf-8'))

    def handle_get_folders(self):
        """Return list of available result folders with metadata from all subdirectories"""
        try:
            folders = []
            data_path = Path(DATA_BASE_PATH)

            if not data_path.exists():
                logging.warning(f"Data path does not exist: {DATA_BASE_PATH}")
                self.send_response(200)
                self.send_header('Content-Type', 'application/json')
                self.send_cors_headers()
                self.end_headers()
                response = {'folders': [],
                            'warning': 'Data directory not found'}
                self.wfile.write(json.dumps(response).encode('utf-8'))
                return

            # Scan all subdirectories for folders with assessment.json
            for run_type_dir in data_path.iterdir():
                if run_type_dir.is_dir():
                    # Scan each run type directory (dry-run, drop1, drop2, drop4, etc.)
                    for item in run_type_dir.iterdir():
                        if item.is_dir():
                            # Check if folder contains required files
                            assessment_file = item / 'assessment.json'
                            if assessment_file.exists():
                                folder_info = {
                                    'name': item.name,
                                    'run_type': run_type_dir.name,
                                    'has_trajectory': (item / 'trajectory.json').exists(),
                                    'has_mapping': (item / 'mapping.json').exists(),
                                    'has_prediction': (item / 'likert_score_prediction.json').exists()
                                }
                                folders.append(folder_info)

            # Sort by folder name
            folders.sort(key=lambda x: x['name'])

            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.send_cors_headers()
            self.end_headers()

            response = {'folders': folders}
            self.wfile.write(json.dumps(response).encode('utf-8'))
            logging.info(f"Served {len(folders)} folders")

        except Exception as e:
            logging.error(f"Error getting folders: {e}")
            self.send_error(500, f"Error reading folders: {e}")

    def handle_get_data(self, path):
        """Handle data file requests like /api/data/computational_tools_0001/assessment.json"""
        try:
            # Extract path components
            path_parts = path.strip('/').split('/')
            if len(path_parts) != 4 or path_parts[0] != 'api' or path_parts[1] != 'data':
                logging.warning(f"Invalid data path format: {path}")
                self.send_error(
                    400, "Invalid data path format. Expected: /api/data/{folder}/{filename}")
                return

            folder = path_parts[2]
            filename = path_parts[3]

            # Check if data file is allowed
            if filename not in ALLOWED_DATA_FILES:
                logging.warning(
                    f"Attempt to access disallowed file: {filename}")
                self.send_error(
                    400, f"File {filename} not allowed. Allowed files: {list(ALLOWED_DATA_FILES)}")
                return

            # Validate folder name (basic security check)
            if not folder.replace('_', '').replace('-', '').isalnum():
                logging.warning(f"Invalid folder name: {folder}")
                self.send_error(400, "Invalid folder name")
                return

            # Search for the folder in all run type directories
            file_path = None
            for run_type_dir in Path(DATA_BASE_PATH).iterdir():
                if run_type_dir.is_dir():
                    potential_path = run_type_dir / folder / filename
                    if potential_path.exists():
                        file_path = potential_path
                        break

            if file_path is None:
                logging.info(
                    f"File not found: {folder}/{filename} in any run type directory")
                self.send_error(404, f"File not found: {filename}")
                return

            # Handle JSON files
            with open(file_path, 'r', encoding='utf-8') as f:
                data = json.load(f)

            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.send_cors_headers()
            self.end_headers()

            self.wfile.write(json.dumps(
                data, ensure_ascii=False).encode('utf-8'))
            logging.info(f"Served data file: {file_path}")

        except json.JSONDecodeError as e:
            logging.error(f"Invalid JSON in file {file_path}: {e}")
            self.send_error(500, f"Invalid JSON in file: {e}")
        except UnicodeDecodeError as e:
            logging.error(f"Encoding error in file {file_path}: {e}")
            self.send_error(500, f"File encoding error: {e}")
        except Exception as e:
            logging.error(f"Error reading data file {path}: {e}")
            self.send_error(500, f"Error reading file: {e}")

    def handle_get_evidence_file(self, path):
        """Handle evidence file requests like /api/data/alloys_0003/evidences/ev_baseline_kic.pdf"""
        try:
            # Extract path components
            path_parts = path.strip('/').split('/')
            if len(path_parts) != 5 or path_parts[0] != 'api' or path_parts[1] != 'data' or path_parts[3] != 'evidences':
                logging.warning(f"Invalid evidence path format: {path}")
                self.send_error(
                    400, "Invalid evidence path format. Expected: /api/data/{folder}/evidences/{filename}")
                return

            folder = path_parts[2]
            filename = path_parts[4]

            # Check if evidence file extension is allowed
            file_extension = Path(filename).suffix.lower()
            if file_extension not in ALLOWED_EVIDENCE_EXTENSIONS:
                logging.warning(
                    f"Attempt to access disallowed evidence file: {filename}")
                self.send_error(
                    400, f"File type {file_extension} not allowed. Allowed extensions: {list(ALLOWED_EVIDENCE_EXTENSIONS)}")
                return

            # Validate folder name (basic security check)
            if not folder.replace('_', '').replace('-', '').isalnum():
                logging.warning(f"Invalid folder name: {folder}")
                self.send_error(400, "Invalid folder name")
                return

            # Search for the folder in all run type directories
            file_path = None
            for run_type_dir in Path(DATA_BASE_PATH).iterdir():
                if run_type_dir.is_dir():
                    potential_path = run_type_dir / folder / 'evidences' / filename
                    if potential_path.exists():
                        file_path = potential_path
                        break

            if file_path is None:
                logging.info(
                    f"Evidence file not found: {folder}/evidences/{filename} in any run type directory")
                self.send_error(404, f"Evidence file not found: {filename}")
                return

            # Determine content type
            content_type, _ = mimetypes.guess_type(str(file_path))
            if content_type is None:
                content_type = 'application/octet-stream'

            # Serve the file
            self.send_response(200)
            self.send_header('Content-Type', content_type)
            self.send_header('Content-Length', str(file_path.stat().st_size))
            self.send_cors_headers()
            self.end_headers()

            with open(file_path, 'rb') as f:
                self.wfile.write(f.read())

            logging.info(f"Served evidence file: {file_path}")

        except Exception as e:
            logging.error(f"Error serving evidence file {path}: {e}")
            self.send_error(500, f"Error serving evidence file: {e}")

    def handle_get_simulation_files(self, path):
        """Handle simulation files list requests like /api/simulation-files/dry-run/alloys_0003"""
        try:
            # Extract run_type and folder from path
            # Path format: /api/simulation-files/{run_type}/{folder}
            path_parts = path.strip('/').split('/')
            if len(path_parts) != 4 or path_parts[0] != 'api' or path_parts[1] != 'simulation-files':
                logging.warning(
                    f"Invalid simulation files path format: {path}")
                self.send_error(
                    400, "Invalid path format. Expected: /api/simulation-files/{run_type}/{folder}")
                return

            run_type = path_parts[2]
            folder = path_parts[3]

            # Validate run_type
            if run_type not in SANDBOX_MAPPING:
                logging.warning(f"Invalid run type: {run_type}")
                self.send_error(
                    400, f"Invalid run type. Allowed: {list(SANDBOX_MAPPING.keys())}")
                return

            # Validate folder name
            if not folder.replace('_', '').replace('-', '').isalnum():
                logging.warning(f"Invalid folder name: {folder}")
                self.send_error(400, "Invalid folder name")
                return

            # Construct sandbox path
            sandbox_dir = SANDBOX_MAPPING[run_type]
            sandbox_path = Path(SANDBOX_BASE_PATH) / sandbox_dir / folder

            if not sandbox_path.exists():
                logging.info(f"Sandbox folder not found: {sandbox_path}")
                self.send_error(
                    404, f"Sandbox folder not found: {sandbox_dir}/{folder}")
                return

            # Get list of allowed simulation files
            files = []
            for file_path in sandbox_path.iterdir():
                if file_path.is_file() and file_path.suffix in ALLOWED_SIMULATION_EXTENSIONS:
                    files.append({
                        'name': file_path.name,
                        'size': file_path.stat().st_size,
                        'extension': file_path.suffix
                    })

            # Sort by name
            files.sort(key=lambda x: x['name'])

            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.send_cors_headers()
            self.end_headers()

            response = {'files': files}
            self.wfile.write(json.dumps(response).encode('utf-8'))
            logging.info(
                f"Served {len(files)} simulation files for {run_type}/{folder}")

        except Exception as e:
            logging.error(f"Error getting simulation files {path}: {e}")
            self.send_error(500, f"Error reading simulation files: {e}")

    def handle_get_simulation_file(self, path):
        """Handle individual simulation file requests like /api/simulation-file/dry-run/alloys_0003/file.csv"""
        try:
            # Extract run_type, folder, and filename from path
            # Path format: /api/simulation-file/{run_type}/{folder}/{filename}
            path_parts = path.strip('/').split('/')
            if len(path_parts) != 5 or path_parts[0] != 'api' or path_parts[1] != 'simulation-file':
                logging.warning(f"Invalid simulation file path format: {path}")
                self.send_error(
                    400, "Invalid path format. Expected: /api/simulation-file/{run_type}/{folder}/{filename}")
                return

            run_type = path_parts[2]
            folder = path_parts[3]
            filename = path_parts[4]

            # Validate run_type
            if run_type not in SANDBOX_MAPPING:
                logging.warning(f"Invalid run type: {run_type}")
                self.send_error(
                    400, f"Invalid run type. Allowed: {list(SANDBOX_MAPPING.keys())}")
                return

            # Validate folder and filename
            if not folder.replace('_', '').replace('-', '').isalnum():
                logging.warning(f"Invalid folder name: {folder}")
                self.send_error(400, "Invalid folder name")
                return

            # Check file extension
            file_ext = Path(filename).suffix
            if file_ext not in ALLOWED_SIMULATION_EXTENSIONS:
                logging.warning(f"Disallowed file extension: {file_ext}")
                self.send_error(
                    400, f"File extension not allowed. Allowed: {list(ALLOWED_SIMULATION_EXTENSIONS)}")
                return

            # Construct file path
            sandbox_dir = SANDBOX_MAPPING[run_type]
            file_path = Path(SANDBOX_BASE_PATH) / \
                sandbox_dir / folder / filename

            # Security check: ensure path is within sandbox
            try:
                file_path.resolve().relative_to(Path(SANDBOX_BASE_PATH).resolve())
            except ValueError:
                logging.error(f"Path traversal attempt: {file_path}")
                self.send_error(403, "Access denied")
                return

            if not file_path.exists():
                logging.info(f"Simulation file not found: {file_path}")
                self.send_error(404, f"File not found: {filename}")
                return

            # Determine content type and read file
            if file_ext == '.png':
                content_type = 'image/png'
                # Read binary file
                with open(file_path, 'rb') as f:
                    content = f.read()

                self.send_response(200)
                self.send_header('Content-Type', content_type)
                self.send_header('Content-Length', str(len(content)))
                self.send_cors_headers()
                self.end_headers()

                self.wfile.write(content)
            else:
                # Handle text files
                content_type = 'text/plain; charset=utf-8'
                if file_ext == '.json':
                    content_type = 'application/json; charset=utf-8'
                elif file_ext == '.csv':
                    content_type = 'text/csv; charset=utf-8'
                elif file_ext == '.py':
                    content_type = 'text/x-python; charset=utf-8'

                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.read()

                self.send_response(200)
                self.send_header('Content-Type', content_type)
                self.send_cors_headers()
                self.end_headers()

                self.wfile.write(content.encode('utf-8'))
            logging.info(
                f"Served simulation file: {run_type}/{folder}/{filename}")

        except UnicodeDecodeError:
            logging.error(f"Encoding error in file {file_path}")
            self.send_error(500, "File encoding error")
        except Exception as e:
            logging.error(f"Error reading simulation file {path}: {e}")
            self.send_error(500, f"Error reading file: {e}")

    def serve_static_file(self, path):
        """Serve static files (HTML, CSS, JS) with security checks"""
        try:
            # Default to index.html for root path
            if path == '/' or path == '':
                path = '/index.html'

            # Security: remove any path traversal attempts
            path = path.lstrip('/')
            if '..' in path or path.startswith('/'):
                logging.warning(f"Path traversal attempt blocked: {path}")
                self.send_error(403, "Access denied")
                return

            file_path = Path(STATIC_BASE_PATH) / path

            # Security check: ensure path is within STATIC_BASE_PATH
            try:
                file_path.resolve().relative_to(Path(STATIC_BASE_PATH).resolve())
            except ValueError:
                logging.error(f"Path traversal attempt: {file_path}")
                self.send_error(403, "Access denied")
                return

            if not file_path.exists():
                logging.info(f"Static file not found: {file_path}")
                self.send_error(404, "File not found")
                return

            # Determine content type
            content_type, encoding = mimetypes.guess_type(str(file_path))
            if content_type is None:
                content_type = 'application/octet-stream'

            # Set appropriate headers for different file types
            self.send_response(200)
            self.send_header('Content-Type', content_type)

            # Add encoding for text files
            if content_type.startswith('text/') or content_type == 'application/javascript':
                self.send_header(
                    'Content-Type', f'{content_type}; charset=utf-8')

            self.send_cors_headers()

            # Cache control for static assets
            if path.endswith(('.css', '.js')):
                self.send_header(
                    'Cache-Control', 'public, max-age=3600')  # 1 hour
            else:
                self.send_header('Cache-Control', 'no-cache')

            self.end_headers()

            # Read and send file
            with open(file_path, 'rb') as f:
                self.wfile.write(f.read())

            logging.info(f"Served static file: {path}")

        except Exception as e:
            logging.error(f"Error serving static file {path}: {e}")
            self.send_error(500, f"Error serving file: {e}")


def main():
    """Start the API server with proper initialization and error handling"""
    print("="*60)
    print("Scientific Claim Assessment Interface API Server")
    print("="*60)
    print(f"Server: http://localhost:{PORT}")
    print(f"Data directory: {DATA_BASE_PATH}")
    print(f"Static files: {STATIC_BASE_PATH}")
    print("Endpoints:")
    print(f"  GET  /                     - Main interface")
    print(f"  GET  /api/health           - Health check")
    print(f"  GET  /api/folders          - Available folders")
    print(f"  GET  /api/data/{{folder}}/{{file}} - Data files")
    print(f"  GET  /api/data/{{folder}}/evidences/{{file}} - Evidence files")
    print("="*60)

    # Validate configuration
    data_path = Path(DATA_BASE_PATH)
    static_path = Path(STATIC_BASE_PATH)

    if not data_path.exists():
        logging.warning(f"Data directory does not exist: {DATA_BASE_PATH}")
        print(f"⚠️  WARNING: Data directory not found")
    else:
        folder_count = len([d for d in data_path.iterdir() if d.is_dir()])
        logging.info(f"Data directory contains {folder_count} folders")
        print(f"✅ Data directory found ({folder_count} folders)")

    if not static_path.exists():
        logging.error(
            f"Static files directory does not exist: {STATIC_BASE_PATH}")
        print(f"❌ ERROR: Static files directory not found")
        return 1
    else:
        print(f"✅ Static files directory found")

    # Check for required static files
    required_files = ['index.html', 'styles.css', 'script.js']
    missing_files = []
    for file in required_files:
        if not (static_path / file).exists():
            missing_files.append(file)

    if missing_files:
        logging.error(f"Missing required static files: {missing_files}")
        print(f"❌ ERROR: Missing files: {', '.join(missing_files)}")
        return 1
    else:
        print(f"✅ All required static files found")

    print("\n🚀 Starting server...")
    server = HTTPServer(('localhost', PORT), APIHandler)
    logging.info(f"Server started on http://localhost:{PORT}")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n\n🛑 Shutting down server...")
        logging.info("Server shutdown requested")
        server.server_close()
        print("✅ Server stopped successfully")
        return 0
    except Exception as e:
        logging.error(f"Server error: {e}")
        print(f"❌ Server error: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
