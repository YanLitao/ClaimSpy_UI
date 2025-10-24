#!/usr/bin/env python3
"""
Export Dynamic Data to Static Structure for GitHub Pages

This script converts the current dynamic API-based data structure to a static structure
that can be deployed on GitHub Pages without requiring a Python server.

The script will:
1. Scan all available assessment folders and their data files
2. Copy all JSON data files to a flat static structure  
3. Generate a folders.json configuration file
4. Create a static data directory structure that can be served directly

Usage:
    python export_static_data.py
"""

import json
import os
import shutil
from pathlib import Path
from datetime import datetime

# Configuration - same as API server
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_ROOT = SCRIPT_DIR
DATA_BASE_PATH = PROJECT_ROOT / 'data' / \
    'sonnet45-full-run-codex-gpt5-1010' / 'output'
SANDBOX_BASE_PATH = PROJECT_ROOT / 'data' / 'sonnet45-full-run-codex-gpt5-1010'
EVALUATION_INTERFACE_PATH = PROJECT_ROOT / 'evaluation-interface'

# Output configuration for static deployment
STATIC_OUTPUT_PATH = EVALUATION_INTERFACE_PATH / 'static-data'
STATIC_FOLDERS_FILE = STATIC_OUTPUT_PATH / 'folders.json'

# Allowed data files (same as API server)
ALLOWED_DATA_FILES = {
    'assessment.json',
    'trajectory.json',
    'mapping.json',
    'likert_score_prediction.json',
    'evidence_source.json'
}

# Sandbox directory mapping (same as API server)
SANDBOX_MAPPING = {
    'dry-run': 'sandbox-dry-run',
    'drop1': 'sandbox-drop1',
    'drop2': 'sandbox-drop2',
    'drop4': 'sandbox-drop4',
    'continuous-release': 'sandbox-continuous-release'
}


def scan_available_folders():
    """Scan all available result folders and return metadata"""
    folders = []
    data_path = Path(DATA_BASE_PATH)

    if not data_path.exists():
        print(f"❌ Data path does not exist: {DATA_BASE_PATH}")
        return folders

    print(f"📁 Scanning data directory: {DATA_BASE_PATH}")

    # Scan all subdirectories for folders with assessment.json
    for run_type_dir in data_path.iterdir():
        if run_type_dir.is_dir():
            print(f"  📂 Scanning run type: {run_type_dir.name}")
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
                            'has_prediction': (item / 'likert_score_prediction.json').exists(),
                            'has_evidence_source': (item / 'evidence_source.json').exists(),
                            'source_path': str(item)
                        }
                        folders.append(folder_info)
                        print(f"    ✅ Found: {item.name}")

    # Sort by folder name
    folders.sort(key=lambda x: x['name'])
    print(f" Total folders found: {len(folders)}")
    return folders


def copy_data_files(folders):
    """Copy all data files to static structure"""
    print(f"\n📋 Copying data files to: {STATIC_OUTPUT_PATH}")

    # Create output directory
    STATIC_OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

    # Create folders for each assessment
    for folder_info in folders:
        folder_name = folder_info['name']
        run_type = folder_info['run_type']
        source_path = Path(folder_info['source_path'])

        # Create destination folder
        dest_folder = STATIC_OUTPUT_PATH / folder_name
        dest_folder.mkdir(parents=True, exist_ok=True)

        print(f"  📁 Processing {folder_name} ({run_type})")

        # Copy each allowed data file
        for filename in ALLOWED_DATA_FILES:
            source_file = source_path / filename
            if source_file.exists():
                dest_file = dest_folder / filename
                shutil.copy2(source_file, dest_file)
                print(f"    ✅ Copied {filename}")
            else:
                print(f"    ⚠️  Missing {filename}")

        # Copy evidences directory if it exists
        evidences_source = source_path / 'evidences'
        if evidences_source.exists() and evidences_source.is_dir():
            evidences_dest = dest_folder / 'evidences'
            if evidences_dest.exists():
                shutil.rmtree(evidences_dest)
            shutil.copytree(evidences_source, evidences_dest)
            evidence_count = len(list(evidences_dest.iterdir()))
            print(f"    📄 Copied evidences directory ({evidence_count} files)")


def copy_simulation_files(folders):
    """Copy simulation files from sandbox directories"""
    print(f"\n🔬 Copying simulation files...")

    for folder_info in folders:
        folder_name = folder_info['name']
        run_type = folder_info['run_type']

        # Check if there's a corresponding sandbox directory
        if run_type in SANDBOX_MAPPING:
            sandbox_dir = SANDBOX_MAPPING[run_type]
            sandbox_path = Path(SANDBOX_BASE_PATH) / sandbox_dir / folder_name

            if sandbox_path.exists():
                # Create destination
                dest_folder = STATIC_OUTPUT_PATH / folder_name / 'simulation-files'
                dest_folder.mkdir(parents=True, exist_ok=True)

                print(f"  📁 Processing simulation files for {folder_name}")

                # Copy simulation files and create manifest
                allowed_extensions = {'.csv', '.json', '.py', '.png'}
                copied_files = []

                for file_path in sandbox_path.iterdir():
                    if file_path.is_file() and file_path.suffix in allowed_extensions:
                        dest_file = dest_folder / file_path.name
                        shutil.copy2(file_path, dest_file)

                        # Add to manifest
                        copied_files.append({
                            'name': file_path.name,
                            'size': file_path.stat().st_size,
                            'extension': file_path.suffix
                        })

                # Create manifest file for static loading
                if copied_files:
                    manifest = {
                        'files': sorted(copied_files, key=lambda x: x['name']),
                        'total_count': len(copied_files),
                        'folder': folder_name,
                        'run_type': run_type
                    }

                    manifest_file = dest_folder / 'manifest.json'
                    with open(manifest_file, 'w', encoding='utf-8') as f:
                        json.dump(manifest, f, indent=2, ensure_ascii=False)

                    print(f"    ✅ Copied {len(copied_files)} simulation files")
                    print(f"    📋 Created manifest.json")
                else:
                    print(f"    ⚠️  No simulation files found")


def generate_folders_config(folders):
    """Generate the static folders.json configuration file"""
    print(f"\n⚙️ Generating folders configuration...")

    # Remove source_path from folder info for client use
    client_folders = []
    for folder_info in folders:
        client_folder = {
            'name': folder_info['name'],
            'run_type': folder_info['run_type'],
            'has_trajectory': folder_info['has_trajectory'],
            'has_mapping': folder_info['has_mapping'],
            'has_prediction': folder_info['has_prediction'],
            'has_evidence_source': folder_info.get('has_evidence_source', False)
        }
        client_folders.append(client_folder)

    # Create folders configuration
    folders_config = {
        'folders': client_folders,
        'generated_at': datetime.now().isoformat(),
        'total_count': len(client_folders),
        'version': '1.0.0',
        'note': 'Static export for GitHub Pages deployment'
    }

    # Write configuration file
    with open(STATIC_FOLDERS_FILE, 'w', encoding='utf-8') as f:
        json.dump(folders_config, f, indent=2, ensure_ascii=False)

    print(f"     Generated {STATIC_FOLDERS_FILE}")
    print(f"     {len(client_folders)} folders configured")


def create_index_redirect():
    """Create a simple index.html redirect in the root for GitHub Pages"""
    index_content = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ClaimSpy UI - Redirecting...</title>
    <meta http-equiv="refresh" content="0; url=./evaluation-interface/index.html">
</head>
<body>
    <p>Redirecting to the evaluation interface...</p>
    <p>If you are not redirected automatically, <a href="./evaluation-interface/index.html">click here</a>.</p>
</body>
</html>'''

    root_index = PROJECT_ROOT / 'index.html'
    with open(root_index, 'w', encoding='utf-8') as f:
        f.write(index_content)

    print(f"✅ Created root index.html redirect")


def main():
    """Main export process"""
    print("="*60)
    print("ClaimSpy UI - Static Data Export for GitHub Pages")
    print("="*60)
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Data source: {DATA_BASE_PATH}")
    print(f"Static output: {STATIC_OUTPUT_PATH}")
    print("="*60)

    # Step 1: Scan available folders
    print("🔍 Step 1: Scanning available folders...")
    folders = scan_available_folders()

    if not folders:
        print("❌ No folders found! Exiting.")
        return 1

    # Step 2: Copy data files
    print("📋 Step 2: Copying data files...")
    copy_data_files(folders)

    # Step 3: Copy simulation files
    print("🔬 Step 3: Copying simulation files...")
    copy_simulation_files(folders)

    # Step 4: Generate configuration
    print("⚙️ Step 4: Generating configuration...")
    generate_folders_config(folders)

    # Step 5: Create root index redirect
    print("🔗 Step 5: Creating root index redirect...")
    create_index_redirect()

    print("\n" + "="*60)
    print("✅ Static export completed successfully!")
    print("="*60)
    print("Next steps:")
    print("1. Modify evaluation-interface/script.js to use static data")
    print("2. Test the static version locally")
    print("3. Deploy to GitHub Pages")
    print("="*60)

    return 0


if __name__ == '__main__':
    exit(main())
