#!/usr/bin/env python3
"""
Validate Static Export for ClaimSpy UI

This script validates that the static export was successful and all
necessary files are present for GitHub Pages deployment.
"""

import json
import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.absolute()
STATIC_DATA_PATH = SCRIPT_DIR / 'evaluation-interface' / 'static-data'
EVALUATION_INTERFACE_PATH = SCRIPT_DIR / 'evaluation-interface'


def validate_static_export():
    """Validate the static export"""
    print("="*60)
    print("ClaimSpy UI - Static Export Validation")
    print("="*60)

    errors = []
    warnings = []

    # 1. Check folders.json exists
    folders_file = STATIC_DATA_PATH / 'folders.json'
    if not folders_file.exists():
        errors.append("❌ Missing folders.json configuration file")
    else:
        print("✅ folders.json configuration file found")

        # Validate folders.json content
        try:
            with open(folders_file, 'r') as f:
                folders_config = json.load(f)

            if 'folders' not in folders_config:
                errors.append("❌ Invalid folders.json - missing 'folders' key")
            else:
                folder_count = len(folders_config['folders'])
                print(f"✅ {folder_count} folders configured")

                # Check each folder
                for folder_info in folders_config['folders']:
                    folder_name = folder_info['name']
                    folder_path = STATIC_DATA_PATH / folder_name

                    if not folder_path.exists():
                        errors.append(
                            f"❌ Missing folder directory: {folder_name}")
                        continue

                    print(f"  📁 Validating {folder_name}...")

                    # Check required JSON files
                    required_files = ['assessment.json', 'trajectory.json',
                                      'mapping.json', 'likert_score_prediction.json']
                    for filename in required_files:
                        file_path = folder_path / filename
                        if not file_path.exists():
                            if folder_info.get(f'has_{filename.split(".")[0]}', True):
                                errors.append(
                                    f"❌ Missing {filename} in {folder_name}")
                            else:
                                warnings.append(
                                    f"⚠️  Optional {filename} missing in {folder_name}")
                        else:
                            # Try to parse JSON
                            try:
                                with open(file_path, 'r') as f:
                                    json.load(f)
                                print(f"    ✅ {filename}")
                            except json.JSONDecodeError as e:
                                errors.append(
                                    f"❌ Invalid JSON in {folder_name}/{filename}: {e}")

                    # Check evidences directory
                    evidences_dir = folder_path / 'evidences'
                    if evidences_dir.exists():
                        evidence_count = len(list(evidences_dir.iterdir()))
                        print(
                            f"    ✅ evidences/ directory ({evidence_count} files)")
                    else:
                        if folder_info.get('has_evidence_source', False):
                            warnings.append(
                                f"⚠️  No evidences directory in {folder_name}")
                        else:
                            print(f"    ℹ️  No evidences directory (expected)")

                    # Check simulation files
                    sim_files_dir = folder_path / 'simulation-files'
                    if sim_files_dir.exists():
                        manifest_file = sim_files_dir / 'manifest.json'
                        if not manifest_file.exists():
                            errors.append(
                                f"❌ Missing manifest.json in {folder_name}/simulation-files/")
                        else:
                            try:
                                with open(manifest_file, 'r') as f:
                                    manifest = json.load(f)
                                sim_count = len(manifest.get('files', []))
                                print(
                                    f"    ✅ simulation-files/ directory ({sim_count} files)")
                            except json.JSONDecodeError as e:
                                errors.append(
                                    f"❌ Invalid manifest.json in {folder_name}: {e}")
                    else:
                        warnings.append(
                            f"⚠️  No simulation-files directory in {folder_name}")

        except json.JSONDecodeError as e:
            errors.append(f"❌ Invalid folders.json: {e}")

    # 2. Check static-config.js
    static_config_file = EVALUATION_INTERFACE_PATH / 'static-config.js'
    if not static_config_file.exists():
        errors.append("❌ Missing static-config.js")
    else:
        print("✅ static-config.js found")

        # Check if static mode is enabled
        with open(static_config_file, 'r') as f:
            config_content = f.read()

        if 'isStaticMode: true' in config_content:
            print("✅ Static mode is enabled")
        elif 'isStaticMode: false' in config_content:
            warnings.append(
                "⚠️  Static mode is disabled - set isStaticMode: true for GitHub Pages")
        else:
            errors.append("❌ Cannot determine static mode setting")

    # 3. Check HTML file includes static-config.js
    html_file = EVALUATION_INTERFACE_PATH / 'index.html'
    if not html_file.exists():
        errors.append("❌ Missing index.html")
    else:
        with open(html_file, 'r') as f:
            html_content = f.read()

        if 'static-config.js' in html_content:
            print("✅ index.html includes static-config.js")
        else:
            errors.append("❌ index.html does not include static-config.js")

    # 4. Check root index.html redirect
    root_index = SCRIPT_DIR / 'index.html'
    if not root_index.exists():
        warnings.append("⚠️  No root index.html redirect for GitHub Pages")
    else:
        print("✅ Root index.html redirect found")

    # 5. Check GitHub Actions workflow
    workflow_file = SCRIPT_DIR / '.github' / 'workflows' / 'deploy.yml'
    if not workflow_file.exists():
        warnings.append("⚠️  No GitHub Actions workflow found")
    else:
        print("✅ GitHub Actions workflow found")

    # Summary
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)

    if not errors:
        print("🎉 All validation checks passed!")
        if warnings:
            print(f"\n📋 {len(warnings)} warnings:")
            for warning in warnings:
                print(f"  {warning}")
        print("\n✅ Ready for GitHub Pages deployment!")
        return 0
    else:
        print(f"❌ {len(errors)} errors found:")
        for error in errors:
            print(f"  {error}")

        if warnings:
            print(f"\n📋 {len(warnings)} warnings:")
            for warning in warnings:
                print(f"  {warning}")

        print("\n🔧 Please fix the errors before deploying to GitHub Pages.")
        return 1


if __name__ == '__main__':
    exit(validate_static_export())
