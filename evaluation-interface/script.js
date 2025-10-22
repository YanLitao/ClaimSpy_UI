// Scientific Claim Assessment Interface - Main JavaScript

// Global variables
let currentData = null;
let trajectoryData = null;
let mappingData = null;
let userApiKey = null;
let apiKeyRemembered = false;
let folderMetadata = new Map(); // Store folder metadata including run_type

// Data directory configuration
// Note: This path is now handled by the API server, so this constant is kept for reference only
let availableFolders = [];

// Initialize
document.addEventListener('DOMContentLoaded', function () {
    // Initialize trajectory flow system
    if (window.trajectoryFlow && typeof window.trajectoryFlow.init === 'function') {
        window.trajectoryFlow.init();
    }

    // Load available folders on page load
    loadAvailableFolders();
});

// Folder selection handling
async function loadAvailableFolders() {
    try {
        const response = await fetch('/api/folders');
        if (response.ok) {
            const data = await response.json();

            // Handle new API format with folder metadata
            if (Array.isArray(data.folders) && data.folders.length > 0) {
                if (typeof data.folders[0] === 'string') {
                    // Old format: array of strings
                    availableFolders = data.folders;
                    // Assume dry-run for old format
                    data.folders.forEach(folder => {
                        folderMetadata.set(folder, { run_type: 'dry-run' });
                    });
                } else {
                    // New format: array of objects with metadata
                    availableFolders = data.folders.map(folder => folder.name || folder);
                    // Store metadata
                    data.folders.forEach(folder => {
                        folderMetadata.set(folder.name || folder, {
                            run_type: folder.run_type || 'dry-run',
                            has_trajectory: folder.has_trajectory,
                            has_mapping: folder.has_mapping,
                            has_prediction: folder.has_prediction
                        });
                    });
                }
            } else {
                availableFolders = [];
            }

            populateFolderDropdown();

            if (data.warning) {
                showError(`Warning: ${data.warning}`);
            }
        } else {
            throw new Error(`Failed to load folders: ${response.statusText}`);
        }
    } catch (error) {
        console.error('Error loading available folders:', error);
        // Fallback to simulated data for demo
        const simulatedFolders = [
            'alloys_0001', 'alloys_0002', 'alloys_0003', 'alloys_0004', 'alloys_0005', 'alloys_0006',
            'alloys_0007', 'alloys_0008', 'alloys_0009', 'alloys_0010', 'alloys_0011', 'alloys_0012',
            'batteries_0001', 'batteries_0002', 'batteries_0003', 'batteries_0004',
            'computational_tools_0001', 'computational_tools_0002', 'computational_tools_0003'
        ];
        availableFolders = simulatedFolders;
        populateFolderDropdown();
        showError('Using demo data - start API server for real data.');
    }
}

function populateFolderDropdown() {
    const dropdown = document.getElementById('folderDropdown');
    dropdown.innerHTML = '<option value="">-- Select a result folder --</option>';

    availableFolders.forEach(folder => {
        const option = document.createElement('option');
        option.value = folder;
        option.textContent = folder;
        dropdown.appendChild(option);
    });
}

async function loadSelectedFolder() {
    const dropdown = document.getElementById('folderDropdown');
    const selectedFolder = dropdown.value;

    if (!selectedFolder) return;

    const loadingStatus = document.getElementById('loadingStatus');
    loadingStatus.style.display = 'block';
    loadingStatus.textContent = `Loading data from ${selectedFolder}...`;

    try {
        // Load all files using API endpoints
        const [assessmentData, trajectoryDataResult, mappingDataResult] = await Promise.all([
            loadJSONFromAPI(selectedFolder, 'assessment.json'),
            loadJSONFromAPI(selectedFolder, 'trajectory.json'),
            loadJSONFromAPI(selectedFolder, 'mapping.json')
        ]);

        // Set global variables
        if (trajectoryDataResult) {
            trajectoryData = trajectoryDataResult;
        }

        if (mappingDataResult) {
            mappingData = mappingDataResult;
        }

        // Load data into trajectory flow module
        if (window.trajectoryFlow && trajectoryDataResult && mappingDataResult) {
            window.trajectoryFlow.loadData(trajectoryDataResult, mappingDataResult);
        }

        // Load assessment data (this will trigger the UI update)
        if (assessmentData) {
            // Add problem/claim information to assessment data for display
            const problemClaim = getProblemClaimForFolder(selectedFolder);
            const dataWithProblem = {
                ...assessmentData,
                problem: problemClaim
            };

            await loadAssessmentData(dataWithProblem);

            // Update folder selection UI
            const folderSelection = document.getElementById('folderSelection');
            folderSelection.classList.add('has-selection', 'collapsed');
        } else {
            throw new Error('Assessment data not found');
        }

    } catch (error) {
        console.error('Error loading folder data:', error);
        showError(`Error loading data from ${selectedFolder}: ${error.message}`);
    } finally {
        loadingStatus.style.display = 'none';
    }
}

function getProblemClaimForFolder(folderName) {
    // Define claims for different problem folders
    const problemClaims = {
        'alloys_0003': 'A novel variation of the WE43 lightweight magnesium alloy has been demonstrated for bioimplants in additively manufactured parts using laser powder bed fusion under high partial pressure of oxygen, attaining concentrations of 10000 ppm of oxygen for higher hardness while maintaining a K_Ic fracture toughness of 16 MPa m^(1/2).',
        'alloys_0001': 'A new high-strength aluminum alloy demonstrates superior tensile strength of 600 MPa while maintaining excellent corrosion resistance.',
        'alloys_0002': 'Novel titanium-vanadium alloy exhibits biocompatibility and enhanced mechanical properties for medical implants.',
        'batteries_0001': 'Lithium-ion battery with silicon nanowire anodes achieves 4000 mAh/g capacity with stable cycling performance.',
        'batteries_0002': 'Solid-state electrolyte enables safer battery operation with energy density exceeding 400 Wh/kg.',
        'computational_tools_0001': 'Al20Zn80 at 870K is a solid at equilibrium',
        'computational_tools_0002': 'Quantum Monte Carlo simulations accurately predict electronic band gaps in semiconductors.'
    };

    return {
        claim: problemClaims[folderName] || `Assessment for ${folderName}`,
        problem_id: folderName
    };
}

// Load JSON data from API endpoint
async function loadJSONFromAPI(folder, filename) {
    try {
        const apiUrl = `/api/data/${folder}/${filename}`;

        const response = await fetch(apiUrl);
        if (response.ok) {
            const data = await response.json();
            return data;
        } else if (response.status === 404) {
            console.log(`File not found: ${filename} in ${folder}`);
            return null;
        } else {
            throw new Error(`API error: ${response.statusText}`);
        }

    } catch (error) {
        console.error(`Error loading ${filename} from ${folder}:`, error);

        // Fallback to simulated data for demo purposes (only for alloys_0003)
        if (folder === 'alloys_0003') {
            try {
                if (filename === 'assessment.json') {
                    return getSampleAssessmentData();
                } else if (filename === 'trajectory.json') {
                    return getSampleTrajectoryData();
                } else if (filename === 'mapping.json') {
                    return getSampleMappingData();
                }
            } catch (fallbackError) {
                console.error(`Fallback data error for ${folder}/${filename}:`, fallbackError);
            }
        }

        return null;
    }
}

// Legacy function - kept for compatibility
async function loadJSONFromPath(filePath) {
    console.warn('loadJSONFromPath is deprecated, use loadJSONFromAPI instead');
    // Extract folder and filename from path for compatibility
    const pathParts = filePath.split('/');
    const folder = pathParts[pathParts.length - 2];
    const filename = pathParts[pathParts.length - 1];
    return loadJSONFromAPI(folder, filename);
}

// Load assessment data
async function loadAssessmentData(data) {
    currentData = data;

    // Update folder selection indicator
    const folderSelection = document.getElementById('folderSelection');
    if (folderSelection) {
        folderSelection.classList.add('has-selection', 'collapsed');
    }

    // Parse assessment data - handle multiple formats
    let assessment;
    try {
        // Format 1: Direct assessment object (assessment.json)
        if (data.type === 'assessment' && data.explanation && data.evidence) {
            assessment = data;
            if (window.evidenceManager) {
                window.evidenceManager.setEvidenceData(assessment.evidence || {});
            }
        }
        // Format 2: Nested solution.assessment format (alloys_0006.json, assessment_1.0_*.json)
        else if (data.solution && data.solution.assessment) {
            assessment = data.solution.assessment;
            if (window.evidenceManager) {
                window.evidenceManager.setEvidenceData(assessment.evidence || {});
            }
        }
        // Format 3: Old solution.json_output format (legacy PoVE-multiagent.json)
        else if (data.solution && data.solution.json_output) {
            assessment = JSON.parse(data.solution.json_output);
            if (window.evidenceManager) {
                window.evidenceManager.setEvidenceData(assessment.evidence || {});
            }
        }
        else {
            throw new Error('Unsupported JSON format. Please ensure your file contains:\n' +
                '• Direct assessment format: {type: "assessment", explanation: [...], evidence: {...}}\n' +
                '• Nested format: {solution: {assessment: {...}}}\n' +
                '• Legacy format: {solution: {json_output: "..."}}');
        }
    } catch (error) {
        showError('Error parsing assessment data: ' + error.message);
        throw error; // Re-throw the error so it can be caught by the caller
    }

    // Parse problem data - handle multiple formats
    let problem = null;
    try {
        // Check if problem exists and parse accordingly
        if (typeof data.problem === 'string') {
            problem = JSON.parse(data.problem);
        } else if (typeof data.problem === 'object' && data.problem !== null) {
            problem = data.problem;
        } else if (data.problem_id) {
            // For direct assessment format, create a minimal problem object
            problem = {
                claim: `Assessment for problem: ${data.problem_id}`,
                problem_id: data.problem_id
            };
        } else {
            // If no problem data, create a default
            problem = {
                claim: 'Assessment data (no specific claim provided)',
                problem_id: 'unknown'
            };
        }
    } catch (error) {
        showError('Error parsing problem data: ' + error.message + '. Using default.');
        problem = {
            claim: 'Assessment data (parsing error)',
            problem_id: 'error'
        };
    }

    // Display claim
    displayClaim(problem.claim);

    // Display claim card interpretation (if available)
    displayClaimCard(assessment);

    // Display system scores
    displaySystemScores(assessment);

    // Display explanations - now with full file loading
    if (window.explanationDisplay) {
        await window.explanationDisplay.displayExplanations(assessment.explanation || []);
    }

    // Update JSON context panel
    updateJSONContext(data);

    // Show content area
    document.getElementById('contentArea').classList.remove('hidden');
}

// Helper function to get the correct assessment path based on data format
function getAssessmentPath() {
    if (!currentData) return [];

    // Format 1: Direct assessment object (assessment.json)
    if (currentData.type === 'assessment' && currentData.explanation && currentData.evidence) {
        return []; // Root level
    }
    // Format 2: Nested solution.assessment format
    else if (currentData.solution && currentData.solution.assessment) {
        return ['solution', 'assessment'];
    }
    // Format 3: Old solution.json_output format
    else if (currentData.solution && currentData.solution.json_output) {
        return ['solution', 'json_output'];
    }

    return []; // default fallback for direct format
}

// Get run type from current folder selection
function getRunTypeFromCurrentFolder() {
    const dropdown = document.getElementById('folderDropdown');
    const selectedFolder = dropdown.value;
    if (!selectedFolder) return null;

    // Get run type from stored metadata
    const metadata = folderMetadata.get(selectedFolder);
    return metadata ? metadata.run_type : 'dry-run';
}

// Get problem folder name from current selection
function getProblemFolderFromCurrentFolder() {
    const dropdown = document.getElementById('folderDropdown');
    const selectedFolder = dropdown.value;
    return selectedFolder; // The dropdown value is already the problem folder name
}

// Wrapper function for showSimulationFile - delegates to simulation manager
function showSimulationFile(filename) {
    if (window.simulationManager && window.simulationManager.showSimulationFile) {
        return window.simulationManager.showSimulationFile(filename);
    } else {
        console.error('Simulation manager not loaded');
    }
}

// Display claim
function displayClaim(claim) {
    const claimElement = document.getElementById('claimText');
    claimElement.innerHTML = `${claim}`;
    claimElement.classList.add('clickable-element');
    claimElement.onclick = () => showInJsonPanel(['problem'], 'claim');
}

// Display claim card interpretation
function displayClaimCard(assessment) {
    const claimCardSection = document.getElementById('claimCardSection');
    const claimCardContent = document.getElementById('claimCardContent');

    // Look for claim card evidence
    let claimCardEvidence = null;
    let claimCardId = null;

    if (assessment.evidence) {
        for (const [evidenceId, evidence] of Object.entries(assessment.evidence)) {
            if (evidence.source === 'claim_card') {
                claimCardEvidence = evidence;
                claimCardId = evidenceId;
                break;
            }
        }
    }

    // If no claim card found, hide the section
    if (!claimCardEvidence) {
        claimCardSection.style.display = 'none';
        return;
    }

    // Show the section
    claimCardSection.style.display = 'block';

    // Parse the claim card JSON from citation
    let parsedClaimCard = null;
    try {
        parsedClaimCard = JSON.parse(claimCardEvidence.citation);
    } catch (error) {
        console.warn('Failed to parse claim card JSON:', error);
        // Show raw citation if JSON parsing fails
        claimCardContent.innerHTML = `
            <div class="claim-card-raw">
                <div class="claim-card-raw-label">Raw Claim Card Data:</div>
                <div class="claim-card-raw-content">${claimCardEvidence.citation}</div>
            </div>
        `;
        return;
    }

    // Create structured display of claim card fields
    const assessmentPath = getAssessmentPath();
    const evidencePath = assessmentPath.length > 0 ? [...assessmentPath, 'evidence', claimCardId] : ['evidence', claimCardId];
    const pathStr = evidencePath.length > 0 ? `'${evidencePath.join("', '")}'` : '';

    // Define the order and labels for claim card fields
    const fieldLabels = {
        'SystemEntity': 'System Entity',
        'StatePhase': 'State/Phase',
        'Conditions': 'Conditions',
        'Manipulation': 'Manipulation',
        'ObservableFoM': 'Observable/FoM',
        'Units': 'Units',
        'Quantifier': 'Quantifier',
        'Threshold': 'Threshold'
    };

    // Build table rows
    let tableRows = [];

    // Display structured fields
    for (const [field, label] of Object.entries(fieldLabels)) {
        if (parsedClaimCard[field] !== undefined) {
            tableRows.push(`
                <tr class="clickable-element" onclick="showInJsonPanel([${pathStr}], 'citation')">
                    <td class="claim-card-field">${label}</td>
                    <td class="claim-card-value">${parsedClaimCard[field]}</td>
                </tr>
            `);
        }
    }

    // Add any additional fields not in the predefined list
    for (const [field, value] of Object.entries(parsedClaimCard)) {
        if (!fieldLabels[field]) {
            tableRows.push(`
                <tr class="clickable-element" onclick="showInJsonPanel([${pathStr}], 'citation')">
                    <td class="claim-card-field">${field}</td>
                    <td class="claim-card-value">${value}</td>
                </tr>
            `);
        }
    }

    // Create the table
    const tableHtml = `
        <table class="claim-card-table">
            <thead>
                <tr>
                    <th>Field</th>
                    <th>Value</th>
                </tr>
            </thead>
            <tbody>
                ${tableRows.join('')}
            </tbody>
        </table>
    `;

    claimCardContent.innerHTML = tableHtml;
}

// Display system scores
function displaySystemScores(assessment) {
    const container = document.getElementById('systemScores');
    const assessmentPath = getAssessmentPath();
    const pathStr = assessmentPath.length > 0 ? `'${assessmentPath.join("', '")}'` : '';

    container.innerHTML = `
        <div class="score-item clickable-element" onclick="showInJsonPanel([${pathStr}], 'likert_score')">
            <div class="score-value">
                ${assessment.likert_score !== undefined ? assessment.likert_score : 'N/A'}
            </div>
            <div class="score-label">Feasibility Score</div>
        </div>
        <div class="score-item clickable-element" onclick="showInJsonPanel([${pathStr}], 'confidence')">
            <div class="score-value">
                ${assessment.confidence !== undefined && assessment.confidence !== null ? assessment.confidence.toFixed(2) : 'N/A'}
            </div>
            <div class="score-label">Confidence</div>
        </div>
    `;
}

// Wrapper functions for evidence and explanation functionality

// Toggle explanation content - delegates to explanation display module
function toggleExplanationContent(index) {
    if (window.explanationDisplay && window.explanationDisplay.toggleExplanationContent) {
        window.explanationDisplay.toggleExplanationContent(index);
    } else {
        // Fallback implementation
        const content = document.getElementById(`explanationContent${index}`);
        content.classList.toggle('expanded');
    }
}

// Mark explanation - delegates to explanation display module
function markExplanation(index, status) {
    if (window.explanationDisplay && window.explanationDisplay.markExplanation) {
        window.explanationDisplay.markExplanation(index, status);
    } else {
        // Fallback implementation
        const item = document.querySelectorAll('.explanation-item')[index];
        const buttons = item.querySelectorAll('.btn');

        // Reset all buttons
        buttons.forEach(btn => {
            btn.classList.remove('checked', 'issue');
        });

        // Mark the clicked button
        event.target.classList.add(status);
    }
}

// Toggle evidence section - delegates to evidence manager
function toggleEvidenceSection(index) {
    if (window.evidenceManager && window.evidenceManager.toggleEvidenceSection) {
        window.evidenceManager.toggleEvidenceSection(index);
    } else {
        // Fallback implementation
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
}

// View local evidence - delegates to evidence manager
function viewLocalEvidence(evidenceId) {
    if (window.evidenceManager && window.evidenceManager.viewLocalEvidence) {
        window.evidenceManager.viewLocalEvidence(evidenceId);
    } else {
        console.error('Evidence manager not loaded');
    }
}

// Show trajectory flow in left panel
function showTrajectoryFlow(explanationIndex = null) {
    // Check if trajectory is currently visible
    const isCurrentlyVisible = window.trajectoryFlow && window.trajectoryFlow.isVisible();

    if (window.trajectoryFlow && typeof window.trajectoryFlow.show === 'function') {
        window.trajectoryFlow.show(explanationIndex);

        // If trajectory was just shown (not hidden), expand evidence section
        if (!isCurrentlyVisible && explanationIndex !== null) {
            const evidenceSection = document.getElementById(`evidenceSection${explanationIndex}`);
            if (evidenceSection && !evidenceSection.classList.contains('expanded')) {
                toggleEvidenceSection(explanationIndex);
            }
        }
    } else {
        console.error('Trajectory flow module not loaded');
    }
}

// Wrapper functions for formatting - delegates to explanation display module
function formatTrajectoryText(text) {
    if (window.explanationDisplay && window.explanationDisplay.formatTrajectoryText) {
        return window.explanationDisplay.formatTrajectoryText(text);
    }
    // Fallback - simple text conversion
    return String(text || '');
}

function formatSearchResults(results) {
    if (window.explanationDisplay && window.explanationDisplay.formatSearchResults) {
        return window.explanationDisplay.formatSearchResults(results);
    }
    // Fallback - simple JSON representation
    return Array.isArray(results) ? results.map(r => JSON.stringify(r)).join('<br>') : '';
}



// Update JSON context
function updateJSONContext(data) {
    const container = document.getElementById('jsonContent');
    // Store the full data for reference
    window.fullJsonData = data;
    window.currentJsonData = data;
    window.currentJsonType = 'assessment';

    container.innerHTML = `
        <div class="json-viewer" id="fullJsonViewer">${formatJsonForHighlight(data, [])}</div>
    `;
}

// Show trajectory data in JSON panel
function showTrajectoryInJsonPanel(trajectoryKey) {
    if (!trajectoryData) {
        console.error('Trajectory data not available');
        return;
    }

    // Open left panel if not already open
    const leftPanel = document.getElementById('leftPanel');
    if (!leftPanel.classList.contains('expanded')) {
        toggleLeftPanel();
    }

    // Update JSON context to show trajectory data
    const container = document.getElementById('jsonContent');
    window.currentJsonData = trajectoryData;
    window.currentJsonType = 'trajectory';

    container.innerHTML = `
        <div class="json-viewer" id="fullJsonViewer">${formatJsonForHighlight(trajectoryData, [])}</div>
    `;

    // Clear previous highlights
    setTimeout(() => {
        document.querySelectorAll('.json-highlight').forEach(el => {
            el.classList.remove('json-highlight');
        });

        // Find and highlight the target trajectory element
        const targetElement = document.querySelector(`[data-path*="trajectory.${trajectoryKey}"]`);
        if (targetElement) {
            targetElement.classList.add('json-highlight');
            targetElement.scrollIntoView({ behavior: 'smooth', block: 'center' });

            // Remove highlight after 3 seconds
            setTimeout(() => {
                targetElement.classList.remove('json-highlight');
            }, 3000);
        } else {
            // Try alternative path formats
            const altElement = document.querySelector(`[data-path="trajectory.${trajectoryKey}"]`) ||
                document.querySelector(`[data-path*="${trajectoryKey}"]`);
            if (altElement) {
                altElement.classList.add('json-highlight');
                altElement.scrollIntoView({ behavior: 'smooth', block: 'center' });
                setTimeout(() => {
                    altElement.classList.remove('json-highlight');
                }, 3000);
            }
        }
    }, 300);
}

// Format JSON with hierarchical structure for highlighting
function formatJsonForHighlight(obj, path = [], depth = 0) {
    const indent = '  '.repeat(depth);
    const pathStr = path.join('.');

    if (typeof obj !== 'object' || obj === null) {
        if (typeof obj === 'string' && obj.length > 100) {
            // Truncate long strings but keep them searchable
            const truncated = obj.substring(0, 100) + '...';
            return `<span class="json-value" data-path="${pathStr}" title="${obj.replace(/"/g, '&quot;')}">${JSON.stringify(truncated)}</span>`;
        }
        return `<span class="json-value" data-path="${pathStr}">${JSON.stringify(obj)}</span>`;
    }

    if (Array.isArray(obj)) {
        if (obj.length === 0) return `<span data-path="${pathStr}">[]</span>`;
        let result = `<span data-path="${pathStr}">[\n`;
        obj.forEach((item, index) => {
            const itemPath = [...path, index];
            result += `${indent}  <span class="json-item" data-path="${itemPath.join('.')}">${formatJsonForHighlight(item, itemPath, depth + 1)}</span>`;
            if (index < obj.length - 1) result += ',';
            result += '\n';
        });
        result += `${indent}]</span>`;
        return result;
    }

    const keys = Object.keys(obj);
    if (keys.length === 0) return `<span data-path="${pathStr}">{}</span>`;

    let result = `<span data-path="${pathStr}">{\n`;
    keys.forEach((key, index) => {
        const keyPath = [...path, key];
        const keyPathStr = keyPath.join('.');
        const value = obj[key];

        // Handle nested JSON strings
        if (typeof value === 'string' && (key === 'problem' || key === 'json_output')) {
            try {
                const parsed = JSON.parse(value);
                result += `${indent}  <span class="json-key" data-path="${keyPathStr}">"${key}": ${formatJsonForHighlight(parsed, keyPath, depth + 1)}</span>`;
            } catch {
                result += `${indent}  <span class="json-key" data-path="${keyPathStr}">"${key}": ${formatJsonForHighlight(value, keyPath, depth + 1)}</span>`;
            }
        } else {
            result += `${indent}  <span class="json-key" data-path="${keyPathStr}">"${key}": ${formatJsonForHighlight(value, keyPath, depth + 1)}</span>`;
        }

        if (index < keys.length - 1) result += ',';
        result += '\n';
    });
    result += `${indent}}</span>`;
    return result;
}

// Show element in JSON panel with highlighting
function showInJsonPanel(basePath, targetKey, index = null) {
    // Open left panel if not already open
    const leftPanel = document.getElementById('leftPanel');
    if (!leftPanel.classList.contains('expanded')) {
        toggleLeftPanel();
    }

    // Switch back to assessment data if we're currently showing trajectory data
    if (window.currentJsonType === 'trajectory') {
        const container = document.getElementById('jsonContent');
        window.currentJsonData = window.fullJsonData;
        window.currentJsonType = 'assessment';
        container.innerHTML = `
            <div class="json-viewer" id="fullJsonViewer">${formatJsonForHighlight(window.fullJsonData, [])}</div>
        `;
    }

    // Clear previous highlights
    document.querySelectorAll('.json-highlight').forEach(el => {
        el.classList.remove('json-highlight');
    });

    // Ensure JSON viewer is populated
    const jsonViewer = document.getElementById('fullJsonViewer');
    if (!jsonViewer || jsonViewer.innerHTML.trim().length === 0) {
        const container = document.getElementById('jsonContent');
        if (container) {
            container.innerHTML = `
                <div class="json-viewer" id="fullJsonViewer">${formatJsonForHighlight(window.currentJsonData || window.fullJsonData, [])}</div>
            `;
        }
    }

    // Find and highlight the target element with longer timeout
    setTimeout(() => {
        let targetElement = null;
        let searchPaths = [];

        // Build search paths based on the parameters
        if (targetKey === 'evidence' && typeof index === 'string') {
            // For evidence items, construct proper paths
            searchPaths = [
                `evidence.${index}`,
                `solution.assessment.evidence.${index}`,
                `assessment.evidence.${index}`,
                index // fallback to just the evidence ID
            ];
        } else if (targetKey === 'explanation' && typeof index === 'number') {
            // For explanation items
            searchPaths = [
                `explanation.${index}`,
                `solution.assessment.explanation.${index}`,
                `assessment.explanation.${index}`
            ];
        } else if (basePath && basePath.length > 0) {
            // For paths with base path
            const fullPath = [...basePath, targetKey].join('.');
            searchPaths = [
                fullPath,
                targetKey // fallback to just the key
            ];
        } else {
            // For simple keys
            searchPaths = [
                targetKey,
                `solution.assessment.${targetKey}`,
                `assessment.${targetKey}`
            ];
        }

        // Try each search path
        for (const searchPath of searchPaths) {
            // Try exact match first
            targetElement = document.querySelector(`[data-path="${searchPath}"]`);
            if (targetElement) break;

            // Try partial match
            targetElement = document.querySelector(`[data-path*="${searchPath}"]`);
            if (targetElement) break;

            // For evidence IDs, also try searching in the content
            if (targetKey === 'evidence' && typeof index === 'string') {
                const allElements = document.querySelectorAll('[data-path*="evidence"]');
                for (const el of allElements) {
                    if (el.textContent && el.textContent.includes(index)) {
                        targetElement = el;
                        break;
                    }
                }
                if (targetElement) break;
            }
        }

        if (targetElement) {
            targetElement.classList.add('json-highlight');
            targetElement.scrollIntoView({ behavior: 'smooth', block: 'center' });

            // Remove highlight after 3 seconds
            setTimeout(() => {
                targetElement.classList.remove('json-highlight');
            }, 3000);
        } else {

            // Comprehensive debug info
            const allElements = document.querySelectorAll('[data-path]');
            const allPaths = Array.from(allElements).map(el => el.getAttribute('data-path'));
            const jsonViewer = document.getElementById('fullJsonViewer');
            // Show filtered paths
            const filteredPaths = allPaths.filter(path => path && (path.includes(targetKey) || (index && path.includes(index.toString()))));

            // If no exact matches, show paths that might be related
            if (filteredPaths.length === 0) {
                const relatedPaths = allPaths.filter(path => {
                    if (!path) return false;
                    const pathLower = path.toLowerCase();
                    const targetLower = targetKey.toLowerCase();
                    const indexStr = index ? index.toString().toLowerCase() : '';
                    return pathLower.includes(targetLower) || (indexStr && pathLower.includes(indexStr));
                });
            }
        }
    }, 500); // Wait longer for panel to expand and DOM to be ready
}

// Close left panel
function closeLeftPanel() {
    const leftPanel = document.getElementById('leftPanel');
    if (leftPanel) {
        leftPanel.classList.remove('expanded');
        // Clear any inline width style that may have been set during resizing
        leftPanel.style.width = '';
    }
}

// Close right panel
function closeRightPanel() {
    const rightPanel = document.getElementById('rightPanel');
    if (rightPanel) {
        rightPanel.classList.remove('expanded');
        // Clear any inline width style that may have been set during resizing
        rightPanel.style.width = '';
        // Clear evidence highlighting when panel is closed
        document.querySelectorAll('.evidence-item.active').forEach(item => {
            item.classList.remove('active');
        });

        // Hide all viewers
        const readingModeViewer = document.getElementById('readingModeViewer');
        const webpageViewer = document.getElementById('webpageViewer');
        const welcomeMessage = document.getElementById('welcomeMessage');

        if (readingModeViewer) readingModeViewer.style.display = 'none';
        if (webpageViewer) webpageViewer.style.display = 'none';
        if (welcomeMessage) welcomeMessage.style.display = 'block';
    }
}

// Make close functions globally accessible
window.closeLeftPanel = closeLeftPanel;
window.closeRightPanel = closeRightPanel;

// Panel toggles
function toggleLeftPanel() {
    const leftPanel = document.getElementById('leftPanel');
    const wasExpanded = leftPanel.classList.contains('expanded');
    leftPanel.classList.toggle('expanded');

    // Clear any inline width style when closing
    if (wasExpanded) {
        leftPanel.style.width = '';
    }
}

function toggleRightPanel() {
    const rightPanel = document.getElementById('rightPanel');
    const wasExpanded = rightPanel.classList.contains('expanded');
    rightPanel.classList.toggle('expanded');

    // Clear evidence highlighting when panel is closed
    if (wasExpanded) {
        // Clear any inline width style when closing
        rightPanel.style.width = '';
        document.querySelectorAll('.evidence-item.active').forEach(item => {
            item.classList.remove('active');
        });
        // Hide all viewers and show welcome message
        document.getElementById('readingModeViewer').style.display = 'none';
        document.getElementById('webpageViewer').style.display = 'none';
        document.getElementById('welcomeMessage').style.display = 'block';
    }
}

// Utility functions
function showError(message) {
    const errorDiv = document.createElement('div');
    errorDiv.className = 'error';
    errorDiv.textContent = message;
    document.querySelector('.main-content').insertBefore(errorDiv, document.querySelector('.main-content').firstChild);
    setTimeout(() => errorDiv.remove(), 5000);
}

