// Simulation Manager Module - Handles all simulation-related functionality

// Load simulation files from sandbox folder for a given evidence item
async function loadSandboxFiles(evidenceId) {
    // The API server automatically handles sandbox folder mapping
    // So we can use the regular loadSimulationFiles function
    return loadSimulationFiles(evidenceId);
}

// Load simulation files for a given evidence item
async function loadSimulationFiles(evidenceId) {
    // Access these functions from the main script
    const runType = (typeof getRunTypeFromCurrentFolder === 'function') ?
        getRunTypeFromCurrentFolder() : null;
    const problemFolder = (typeof getProblemFolderFromCurrentFolder === 'function') ?
        getProblemFolderFromCurrentFolder() : null;

    if (!runType || !problemFolder) {
        console.error('Cannot determine run type or problem folder');
        return [];
    }

    try {
        const response = await fetch(`/api/simulation-files/${runType}/${problemFolder}`);

        if (!response.ok) {
            if (response.status === 404) {
                console.log(`🚫 No simulation files found for ${runType}/${problemFolder}`);
                return [];
            }
            throw new Error(`HTTP ${response.status}: ${response.statusText}`);
        }

        const data = await response.json();
        const allFiles = data.files || [];

        // Only return data files (CSV, JSON, PY), PNG files are handled as visualizations
        // Also filter out explanation files from display
        const dataFiles = allFiles.filter(file => {
            const isDataFile = ['.csv', '.json', '.py'].includes(file.extension);
            // Check if this is an explanation file - use window.codeExplanation if available
            const isExplanation = (typeof window !== 'undefined' && window.codeExplanation && typeof window.codeExplanation.isExplanationFile === 'function') ?
                window.codeExplanation.isExplanationFile(file.name) : file.name.endsWith('_explanation.json');

            return isDataFile && !isExplanation;
        });

        const pngFiles = allFiles.filter(file => file.extension === '.png');

        // Associate PNG files with their corresponding data files
        const filesWithVisualizations = dataFiles.map(file => {
            const baseName = file.name.replace(/\.[^/.]+$/, ''); // Remove extension
            const correspondingPng = pngFiles.find(png => {
                const pngBaseName = png.name.replace(/\.[^/.]+$/, '');
                return pngBaseName === baseName;
            });

            return {
                ...file,
                hasVisualization: !!correspondingPng,
                visualizationFile: correspondingPng?.name
            };
        });

        return filesWithVisualizations;
    } catch (error) {
        console.error('Error loading simulation files:', error);
        return [];
    }
}

// Display simulation file content in right panel
async function showSimulationFile(filename) {
    // Use the new enhanced function with explanation support if available
    if (typeof window !== 'undefined' && window.codeExplanation && typeof window.codeExplanation.showSimulationFileWithExplanation === 'function') {
        return window.codeExplanation.showSimulationFileWithExplanation(filename);
    }

    // Fallback to original implementation
    const runType = (typeof getRunTypeFromCurrentFolder === 'function') ?
        getRunTypeFromCurrentFolder() : null;
    const problemFolder = (typeof getProblemFolderFromCurrentFolder === 'function') ?
        getProblemFolderFromCurrentFolder() : null;

    if (!runType || !problemFolder) {
        console.error('Cannot determine run type or problem folder');
        return;
    }

    // Show right panel
    const rightPanel = document.getElementById('rightPanel');
    if (!rightPanel.classList.contains('expanded')) {
        rightPanel.classList.add('expanded');
    }

    // Show loading state
    const loadingOverlay = document.getElementById('loadingOverlay');
    const readingModeViewer = document.getElementById('readingModeViewer');
    const webpageViewer = document.getElementById('webpageViewer');
    const welcomeMessage = document.getElementById('welcomeMessage');

    loadingOverlay.style.display = 'flex';
    readingModeViewer.style.display = 'none';
    webpageViewer.style.display = 'none';
    welcomeMessage.style.display = 'none';

    try {
        // API server will automatically check sandbox folders
        const response = await fetch(`/api/simulation-file/${runType}/${problemFolder}/${filename}`);

        if (!response.ok) {
            throw new Error(`HTTP ${response.status}: ${response.statusText}`);
        }

        const content = await response.text();

        // Visualization display has been removed from right panel

        let contentHtml = '';

        if (filename.endsWith('.csv')) {
            // Parse and display CSV as a table
            contentHtml = formatCsvAsTable(content, filename);
        } else {
            // Handle other file types
            let displayContent = content;
            let language = 'text';

            if (filename.endsWith('.json')) {
                try {
                    const jsonData = JSON.parse(content);
                    displayContent = JSON.stringify(jsonData, null, 2);
                    language = 'json';
                } catch (e) {
                    language = 'text';
                }
            } else if (filename.endsWith('.py')) {
                language = 'python';
            }

            contentHtml = `
                <div class="content-display">
                    <pre><code class="language-${language}">${displayContent}</code></pre>
                </div>
            `;
        }

        // Display in reading mode viewer
        readingModeViewer.innerHTML = `
            <div class="source-info">
                <div class="source-title">Simulation File: ${filename}</div>
                <div style="margin-top: 0.5rem; font-style: italic; color: #666;">
                    Location: ${runType}/${problemFolder}/${filename}
                </div>
            </div>
            <div style="margin-top: 1rem;">
                ${contentHtml}
            </div>
        `;

        readingModeViewer.style.display = 'block';

    } catch (error) {
        console.error('Error loading simulation file:', error);
        readingModeViewer.innerHTML = `
            <div class="source-info">
                <div class="source-title">Error Loading File</div>
            </div>
            <div class="error-content" style="margin-top: 1rem; padding: 1rem; background: #fee; border: 1px solid #fcc; border-radius: 4px; color: #c33;">
                Failed to load ${filename}: ${error.message}
            </div>
        `;
        readingModeViewer.style.display = 'block';
    } finally {
        loadingOverlay.style.display = 'none';
    }
}

// Format CSV content as an HTML table
function formatCsvAsTable(csvContent, filename) {
    try {
        const lines = csvContent.trim().split('\n');
        if (lines.length === 0) {
            return '<p>Empty CSV file</p>';
        }

        // Parse CSV headers
        const headers = parseCsvLine(lines[0]);

        let tableHtml = `
            <div class="csv-table-container">
                <h4 style="color: #2c3e50; margin-bottom: 1rem;">📋 Data Table</h4>
                <div style="overflow-x: auto; background: white; border-radius: 8px; border: 1px solid #e9ecef;">
                    <table style="width: 100%; border-collapse: collapse; font-size: 0.9rem;">
                        <thead>
                            <tr style="background: #f8f9fa;">
        `;

        // Add headers
        headers.forEach(header => {
            tableHtml += `<th style="padding: 0.75rem; text-align: left; border-bottom: 2px solid #dee2e6; font-weight: 600; color: #495057;">${header}</th>`;
        });

        tableHtml += '</tr></thead><tbody>';

        // Add data rows
        for (let i = 1; i < lines.length; i++) {
            if (lines[i].trim()) {
                const cells = parseCsvLine(lines[i]);
                const rowStyle = i % 2 === 0 ? 'background: #f8f9fa;' : 'background: white;';
                tableHtml += `<tr style="${rowStyle}">`;

                cells.forEach((cell, index) => {
                    const cellContent = cell || '-';
                    const isNumeric = !isNaN(cell) && cell !== '';
                    const textAlign = isNumeric ? 'text-align: right;' : 'text-align: left;';
                    tableHtml += `<td style="padding: 0.75rem; border-bottom: 1px solid #dee2e6; ${textAlign}">${cellContent}</td>`;
                });

                tableHtml += '</tr>';
            }
        }

        tableHtml += `
                    </tbody>
                </table>
            </div>
            <div style="margin-top: 0.5rem; font-size: 0.85rem; color: #666;">
                Rows: ${lines.length - 1} | Columns: ${headers.length}
            </div>
        </div>
        `;

        return tableHtml;

    } catch (error) {
        console.error('Error parsing CSV:', error);
        return `
            <div class="content-display">
                <h4 style="color: #2c3e50; margin-bottom: 1rem;">📋 Raw CSV Content</h4>
                <pre style="background: #f8f9fa; border: 1px solid #e9ecef; border-radius: 4px; padding: 1rem; overflow-x: auto; white-space: pre-wrap;">${csvContent}</pre>
            </div>
        `;
    }
}

// Simple CSV line parser (handles basic cases)
function parseCsvLine(line) {
    const result = [];
    let current = '';
    let isInQuotes = false;

    for (let i = 0; i < line.length; i++) {
        const char = line[i];

        if (char === '"') {
            isInQuotes = !isInQuotes;
        } else if (char === ',' && !isInQuotes) {
            result.push(current.trim());
            current = '';
        } else {
            current += char;
        }
    }

    result.push(current.trim());
    return result;
}

// Format file size for display
function formatFileSize(bytes) {
    if (bytes === 0) return '0 B';
    const k = 1024;
    const sizes = ['B', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(1)) + ' ' + sizes[i];
}

// Show error in reading mode
function showErrorInReadingMode(message) {
    const rightPanel = document.getElementById('rightPanel');
    const readingModeViewer = document.getElementById('readingModeViewer');
    const loadingOverlay = document.getElementById('loadingOverlay');

    if (!rightPanel.classList.contains('expanded')) {
        rightPanel.classList.add('expanded');
    }

    readingModeViewer.innerHTML = `
        <div class="source-info">
            <div class="source-title">Error</div>
        </div>
        <div class="error-content" style="margin-top: 1rem; padding: 1rem; background: #fee; border: 1px solid #fcc; border-radius: 4px; color: #c33;">
            ${message}
        </div>
    `;
    readingModeViewer.style.display = 'block';
    loadingOverlay.style.display = 'none';
}

// Export functions for global access
if (typeof window !== 'undefined') {
    window.simulationManager = {
        loadSandboxFiles,
        loadSimulationFiles,
        showSimulationFile,
        formatCsvAsTable,
        parseCsvLine,
        formatFileSize,
        showErrorInReadingMode
    };
}