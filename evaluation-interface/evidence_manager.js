// Evidence Manager Module - Handles all evidence-related functionality

// Global variables for evidence management
let evidenceData = null;
let evidenceCache = new Map(); // Cache for scraped content and AI analysis

// Helper function to get appropriate CSS class for evidence type
function getEvidenceTypeClass(type) {
    if (!type) return 'unknown';

    const typeMap = {
        'web search': 'web-search',
        'simulation': 'simulation',
        'simulation result': 'simulation',
        'reasoning': 'reasoning',
        'experiment': 'experiment',
        'experimental': 'experiment',
        'literature': 'literature',
        'literature review': 'literature',
        'database': 'database',
        'data': 'database'
    };

    const lowerType = type.toLowerCase();
    return typeMap[lowerType] || 'web-search'; // default to web-search style
}

// Helper function to validate if a string is a valid URL
function isValidUrl(string) {
    if (!string || typeof string !== 'string') return false;
    try {
        const url = new URL(string);
        return url.protocol === 'http:' || url.protocol === 'https:';
    } catch (_) {
        return false;
    }
}

// Wrapper functions for simulation functionality - delegates to simulation manager
function loadSandboxFiles(evidenceId) {
    if (window.simulationManager && window.simulationManager.loadSandboxFiles) {
        return window.simulationManager.loadSandboxFiles(evidenceId);
    } else {
        console.error('Simulation manager not loaded');
        return [];
    }
}

function loadSimulationFiles(evidenceId) {
    if (window.simulationManager && window.simulationManager.loadSimulationFiles) {
        return window.simulationManager.loadSimulationFiles(evidenceId);
    } else {
        console.error('Simulation manager not loaded');
        return [];
    }
}

// Wrapper function for showSimulationFile - delegates to simulation manager
function showSimulationFile(filename) {
    if (window.simulationManager && window.simulationManager.showSimulationFile) {
        return window.simulationManager.showSimulationFile(filename);
    } else {
        console.error('Simulation manager not loaded');
    }
}

// Wrapper functions for CSV and file utility functions - delegates to simulation manager
function formatCsvAsTable(csvContent, filename) {
    if (window.simulationManager && window.simulationManager.formatCsvAsTable) {
        return window.simulationManager.formatCsvAsTable(csvContent, filename);
    } else {
        console.error('Simulation manager not loaded');
        return '<p>Error: Simulation manager not available</p>';
    }
}

function parseCsvLine(line) {
    if (window.simulationManager && window.simulationManager.parseCsvLine) {
        return window.simulationManager.parseCsvLine(line);
    } else {
        console.error('Simulation manager not loaded');
        return [];
    }
}

function formatFileSize(bytes) {
    if (window.simulationManager && window.simulationManager.formatFileSize) {
        return window.simulationManager.formatFileSize(bytes);
    } else {
        console.error('Simulation manager not loaded');
        return '0 B';
    }
}

// Highlight active evidence item
function highlightActiveEvidence(evidenceId) {
    // Remove previous active states
    document.querySelectorAll('.evidence-item.active').forEach(item => {
        item.classList.remove('active');
    });

    // Find and highlight the current evidence item by looking for the jump button
    const jumpButton = document.getElementById(`jumpBtn${evidenceId}`);
    if (jumpButton) {
        const evidenceElement = jumpButton.closest('.evidence-item');
        if (evidenceElement) {
            evidenceElement.classList.add('active');
        }
    }
}

// View local evidence file
async function viewLocalEvidence(evidenceId) {
    try {
        // Load evidence source mapping
        const runType = (typeof getRunTypeFromCurrentFolder === 'function') ?
            getRunTypeFromCurrentFolder() : null;
        const problemFolder = (typeof getProblemFolderFromCurrentFolder === 'function') ?
            getProblemFolderFromCurrentFolder() : null;

        if (!runType || !problemFolder) {
            showErrorInReadingMode('Cannot determine current folder context');
            return;
        }

        // Load evidence_source.json from the same folder as assessment.json
        const evidenceSourceData = (typeof loadJSONFromAPI === 'function') ?
            await loadJSONFromAPI(problemFolder, 'evidence_source.json') : null;

        if (!evidenceSourceData || !evidenceSourceData.mappings) {
            showErrorInReadingMode('No evidence source mapping found');
            return;
        }

        const evidenceMapping = evidenceSourceData.mappings[evidenceId];
        if (!evidenceMapping) {
            showErrorInReadingMode(`No local file mapping found for evidence: ${evidenceId}`);
            return;
        }

        // Show the PDF in the right panel
        await displayLocalPDF(evidenceMapping, evidenceId);

        // Highlight the active evidence
        highlightActiveEvidence(evidenceId);

    } catch (error) {
        console.error('Error loading local evidence:', error);
        showErrorInReadingMode('Failed to load local evidence file');
    }
}

// Show local PDF in the right panel (using pdf_controls.js)
async function displayLocalPDF(evidenceMapping, evidenceId) {
    // Call the enhanced PDF viewer from pdf_controls.js
    // The function in pdf_controls.js handles PDF.js integration and highlighting
    if (typeof window.showLocalPDF === 'function') {
        return await window.showLocalPDF(evidenceMapping, evidenceId);
    } else {
        console.warn('PDF controls not loaded, falling back to basic PDF viewer');
        // Fallback implementation
        const rightPanel = document.getElementById('rightPanel');
        if (!rightPanel.classList.contains('expanded')) {
            rightPanel.classList.add('expanded');
        }

        const readingModeViewer = document.getElementById('readingModeViewer');
        const loadingOverlay = document.getElementById('loadingOverlay');

        const runType = (typeof getRunTypeFromCurrentFolder === 'function') ?
            getRunTypeFromCurrentFolder() : null;
        const problemFolder = (typeof getProblemFolderFromCurrentFolder === 'function') ?
            getProblemFolderFromCurrentFolder() : null;
        const pdfPath = `/api/data/${runType}/${problemFolder}/evidences/${evidenceMapping.filename}`;

        readingModeViewer.innerHTML = `
            <div class="pdf-viewer-container">
                <h4>📄 ${evidenceMapping.filename}</h4>
                <embed src="${pdfPath}" type="application/pdf" width="100%" height="600px" />
            </div>
        `;
        readingModeViewer.style.display = 'block';
        loadingOverlay.style.display = 'none';
    }
}

// Open PDF in new tab
function openPDFInNewTab(pdfPath) {
    window.open(pdfPath, '_blank');
}

// Update local file button texts
function updateLocalFileButtonTexts() {
    // Update button titles for local file viewing
    document.querySelectorAll('.local-file-btn').forEach(btn => {
        btn.title = 'View local evidence file';
    });
}

// Set evidence data (called from main script)
function setEvidenceData(data) {
    evidenceData = data;
}

// Get evidence data
function getEvidenceData() {
    return evidenceData;
}

// Toggle evidence section
function toggleEvidenceSection(index) {
    const evidenceSection = document.getElementById(`evidenceSection${index}`);
    const isExpanding = !evidenceSection.classList.contains('expanded');
    evidenceSection.classList.toggle('expanded');

    // If trajectory flow is visible, handle evidence connections
    if (window.trajectoryFlow && window.trajectoryFlow.isVisible()) {
        if (isExpanding) {
            // Evidence section is being expanded - draw connections
            setTimeout(() => {
                window.trajectoryFlow.drawEvidenceConnections(index);
            }, 100); // Wait for expansion animation
        } else {
            // Evidence section is being collapsed - clear connections
            window.trajectoryFlow.clearEvidenceConnections();
        }
    }
}

// Wrapper function for showErrorInReadingMode - delegates to simulation manager
function showErrorInReadingMode(message) {
    if (window.simulationManager && window.simulationManager.showErrorInReadingMode) {
        window.simulationManager.showErrorInReadingMode(message);
    } else {
        console.error('Simulation manager not loaded');
    }
}

// Export functions for global access
if (typeof window !== 'undefined') {
    window.evidenceManager = {
        setEvidenceData,
        getEvidenceData,
        getEvidenceTypeClass,
        isValidUrl,
        loadSandboxFiles,
        loadSimulationFiles,
        showSimulationFile,
        formatCsvAsTable,
        parseCsvLine,
        formatFileSize,
        highlightActiveEvidence,
        viewLocalEvidence,
        displayLocalPDF,
        openPDFInNewTab,
        updateLocalFileButtonTexts,
        toggleEvidenceSection,
        showErrorInReadingMode
    };
}