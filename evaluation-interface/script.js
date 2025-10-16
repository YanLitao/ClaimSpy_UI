// Scientific Claim Assessment Interface - Main JavaScript

// Global variables
let currentData = null;
let evidenceData = null;
let trajectoryData = null;
let mappingData = null;
let userApiKey = null;
let apiKeyRemembered = false;
let evidenceCache = new Map(); // Cache for scraped content and AI analysis

// Data directory configuration
// Note: This path is now handled by the API server, so this constant is kept for reference only
let availableFolders = [];

// Initialize
document.addEventListener('DOMContentLoaded', function () {
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
                } else {
                    // New format: array of objects with metadata
                    availableFolders = data.folders.map(folder => folder.name || folder);
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
            console.log('Trajectory data loaded successfully');
        }

        if (mappingDataResult) {
            mappingData = mappingDataResult;
            console.log('Mapping data loaded successfully');
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
        'computational_tools_0001': 'Machine learning model predicts material properties with 95% accuracy using only composition data.',
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
        console.log(`Loading ${filename} from ${folder}...`);

        const response = await fetch(apiUrl);
        if (response.ok) {
            const data = await response.json();
            console.log(`✅ Successfully loaded ${filename} from ${folder}`);
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
                    console.log('🎭 Using demo assessment data');
                    return getSampleAssessmentData();
                } else if (filename === 'trajectory.json') {
                    console.log('🎭 Using demo trajectory data');
                    return getSampleTrajectoryData();
                } else if (filename === 'mapping.json') {
                    console.log('🎭 Using demo mapping data');
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
            evidenceData = assessment.evidence || {};
        }
        // Format 2: Nested solution.assessment format (alloys_0006.json, assessment_1.0_*.json)
        else if (data.solution && data.solution.assessment) {
            assessment = data.solution.assessment;
            evidenceData = assessment.evidence || {};
        }
        // Format 3: Old solution.json_output format (legacy PoVE-multiagent.json)
        else if (data.solution && data.solution.json_output) {
            assessment = JSON.parse(data.solution.json_output);
            evidenceData = assessment.evidence || {};
        }
        else {
            console.log('Data structure:', Object.keys(data));
            throw new Error('Unsupported JSON format. Please ensure your file contains:\n' +
                '• Direct assessment format: {type: "assessment", explanation: [...], evidence: {...}}\n' +
                '• Nested format: {solution: {assessment: {...}}}\n' +
                '• Legacy format: {solution: {json_output: "..."}}');
        }
    } catch (error) {
        showError('Error parsing assessment data: ' + error.message);
        return;
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

    // Display system scores
    displaySystemScores(assessment);

    // Display explanations - now with full file loading
    await displayExplanations(assessment.explanation || []);

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

// Get run type from current folder selection
function getRunTypeFromCurrentFolder() {
    const dropdown = document.getElementById('folderDropdown');
    const selectedFolder = dropdown.value;
    if (!selectedFolder) return null;

    // For now, all folders in the current data directory are from dry-run
    // In the future, this could be made more dynamic by including run type info in the API response
    return 'dry-run';
}

// Get problem folder name from current selection
function getProblemFolderFromCurrentFolder() {
    const dropdown = document.getElementById('folderDropdown');
    const selectedFolder = dropdown.value;
    return selectedFolder; // The dropdown value is already the problem folder name
}

// Load simulation files for a given evidence item
async function loadSimulationFiles(evidenceId) {
    const runType = getRunTypeFromCurrentFolder();
    const problemFolder = getProblemFolderFromCurrentFolder();

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
            const isExplanation = typeof isExplanationFile === 'function' ?
                isExplanationFile(file.name) : file.name.endsWith('_explanation.json');



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
}// Display simulation file content in right panel
async function showSimulationFile(filename) {
    // Use the new enhanced function with explanation support
    if (typeof showSimulationFileWithExplanation === 'function') {
        return showSimulationFileWithExplanation(filename);
    }

    // Fallback to original implementation
    const runType = getRunTypeFromCurrentFolder();
    const problemFolder = getProblemFolderFromCurrentFolder();

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
        const response = await fetch(`/api/simulation-file/${runType}/${problemFolder}/${filename}`);
        if (!response.ok) {
            throw new Error(`HTTP ${response.status}: ${response.statusText}`);
        }

        const content = await response.text();

        // Check if there's a corresponding visualization (PNG file)
        const baseName = filename.replace(/\.[^/.]+$/, ''); // Remove extension
        const visualizationFilename = `${baseName}.png`;

        console.log(`Looking for visualization: ${visualizationFilename} for file: ${filename}`);

        let visualizationHtml = '';
        try {
            const vizResponse = await fetch(`/api/simulation-file/${runType}/${problemFolder}/${visualizationFilename}`);
            console.log(`Visualization response status: ${vizResponse.status} for ${visualizationFilename}`);

            if (vizResponse.ok) {
                console.log(`✅ Found visualization: ${visualizationFilename}`);
                // PNG file exists, create image element
                visualizationHtml = `
                    <div class="visualization-section" style="margin-bottom: 2rem;">
                        <h4 style="color: #2c3e50; margin-bottom: 1rem;">📊 Visualization</h4>
                        <div style="text-align: center; background: #f8f9fa; padding: 1rem; border-radius: 8px; border: 1px solid #e9ecef;">
                            <img src="/api/simulation-file/${runType}/${problemFolder}/${visualizationFilename}" 
                                 alt="Visualization for ${filename}" 
                                 style="max-width: 100%; height: auto; border-radius: 4px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);" />
                        </div>
                    </div>
                `;
            } else {
                console.log(`❌ No visualization found: ${vizResponse.status} for ${visualizationFilename}`);
            }
        } catch (vizError) {
            // No visualization available, continue without it
            console.log(`❌ Error fetching visualization for ${filename}:`, vizError);
        }

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
                ${visualizationHtml}
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

// Display claim
function displayClaim(claim) {
    const claimElement = document.getElementById('claimText');
    claimElement.textContent = claim;
    claimElement.classList.add('clickable-element');
    claimElement.onclick = () => showInJsonPanel(['problem'], 'claim');
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
            <div class="score-label">Likert Score</div>
        </div>
        <div class="score-item clickable-element" onclick="showInJsonPanel([${pathStr}], 'continuous_score')">
            <div class="score-value">
                ${assessment.continuous_score !== undefined && assessment.continuous_score !== null ? assessment.continuous_score.toFixed(2) : 'N/A'}
            </div>
            <div class="score-label">Continuous Score</div>
        </div>
        <div class="score-item clickable-element" onclick="showInJsonPanel([${pathStr}], 'confidence')">
            <div class="score-value">
                ${assessment.confidence !== undefined && assessment.confidence !== null ? assessment.confidence.toFixed(2) : 'N/A'}
            </div>
            <div class="score-label">Confidence</div>
        </div>
    `;
}

// Display explanations
async function displayExplanations(explanations) {
    const container = document.getElementById('explanationsList');
    container.innerHTML = '';
    const assessmentPath = getAssessmentPath();

    // Handle case where explanations might be undefined or empty
    if (!explanations || !Array.isArray(explanations)) {
        container.innerHTML = '<p style="color: #666; padding: 1rem;">No explanations available.</p>';
        return;
    }

    // Pre-load all simulation files for all simulation evidence
    const simulationFilesCache = {};
    if (evidenceData) {
        const simulationEvidence = Object.keys(evidenceData).filter(evidenceId => {
            const evidence = evidenceData[evidenceId];
            return evidence && evidence.type && evidence.type.toLowerCase().includes('simulation');
        });

        // Load all simulation files at once
        for (const evidenceId of simulationEvidence) {
            try {
                const files = await loadSimulationFiles(evidenceId);
                simulationFilesCache[evidenceId] = files;
            } catch (error) {
                console.error(`Error loading simulation files for ${evidenceId}:`, error);
                simulationFilesCache[evidenceId] = [];
            }
        }
    }

    explanations.forEach((explanation, index) => {
        const explanationEl = document.createElement('div');
        explanationEl.className = 'explanation-item';
        const pathStr = assessmentPath.length > 0 ? `'${assessmentPath.join("', '")}'` : '';

        // Check if explanation has valid evidence
        const validEvidence = (explanation.evidence || []).filter(id => evidenceData[id] && evidenceData[id].citation);
        const hasEvidence = validEvidence.length > 0;

        explanationEl.innerHTML = `
            <div class="trajectory-section" id="trajectorySection${index}" style="display: none;">
                <div class="evidence-header-toggle" onclick="toggleTrajectoryDetails(${index})">
                    <h4>🔍 Trajectory Details</h4>
                    <span class="collapse-icon">▼</span>
                </div>
                <div class="trajectory-content" id="trajectoryContent${index}">
                    ${mappingData ? generateTrajectoryContent(index + 1) : '<p>Trajectory data not available</p>'}
                </div>
            </div>
            <div class="explanation-header">
                <div class="explanation-text clickable-element" onclick="showInJsonPanel([${pathStr}], 'explanation', ${index})">${explanation.text}</div>
                <div class="explanation-actions">
                    ${hasEvidence ? `<button class="btn" onclick="toggleEvidenceSection(${index})">Evidence (${validEvidence.length})</button>` : ''}
                    ${mappingData ? `<button class="btn" onclick="toggleTrajectoryContent(${index})">Show Trajectory</button>` : ''}
                </div>
            </div>
            <div class="explanation-content" id="explanationContent${index}">
                <textarea class="comment-box" placeholder="Add your comments about this explanation..."></textarea>
            </div>
            ${hasEvidence ? `
            <div class="evidence-section" id="evidenceSection${index}">
                <div class="evidence-header-toggle" onclick="toggleEvidenceSection(${index})">
                    <h4>Evidence (${validEvidence.length})</h4>
                    <span class="collapse-icon">▼</span>
                </div>
                <div id="evidenceContent${index}">
                    ${(explanation.evidence || []).map(evidenceId => {
            const evidence = evidenceData[evidenceId];
            if (!evidence || !evidence.citation) return '';
            const hasValidUrl = isValidUrl(evidence.source);

            // Generate simulation files HTML if this is simulation evidence
            let simulationFilesHtml = '';
            if (evidence.type && evidence.type.toLowerCase().includes('simulation') && simulationFilesCache[evidenceId]) {
                const files = simulationFilesCache[evidenceId];
                if (files.length > 0) {
                    const fileBoxesHtml = files.map(file => {
                        const visualizationIndicator = file.hasVisualization ? ' 📊' : '';
                        const tooltipText = file.hasVisualization ?
                            `View ${file.name} with visualization (${formatFileSize(file.size)})` :
                            `View ${file.name} (${formatFileSize(file.size)})`;

                        return `
                            <div class="simulation-file-box" onclick="event.stopPropagation(); showSimulationFile('${file.name}')" 
                                 style="display: inline-block; margin: 2px; padding: 6px 10px; background: #f8f9fa; 
                                        border: 1px solid #dee2e6; border-radius: 4px; cursor: pointer; font-size: 0.85rem;
                                        transition: background-color 0.2s;"
                                 onmouseover="this.style.backgroundColor='#e9ecef'" 
                                 onmouseout="this.style.backgroundColor='#f8f9fa'"
                                 title="${tooltipText}">
                                ${file.name}${visualizationIndicator}
                            </div>
                        `;
                    }).join('');

                    simulationFilesHtml = `
                        <div class="simulation-files-section" id="simFiles${evidenceId}" style="margin-top: 0.5rem;">
                            <div class="simulation-files-header" style="font-size: 0.9rem; color: #666; margin-bottom: 0.5rem;">Simulation Files:</div>
                            <div class="simulation-files-content" id="simFilesContent${evidenceId}">${fileBoxesHtml}</div>
                        </div>
                    `;
                } else {
                    simulationFilesHtml = `
                        <div class="simulation-files-section" id="simFiles${evidenceId}" style="margin-top: 0.5rem;">
                            <div class="simulation-files-header" style="font-size: 0.9rem; color: #666; margin-bottom: 0.5rem;">Simulation Files:</div>
                            <div class="simulation-files-content" id="simFilesContent${evidenceId}" style="color: #999;">No simulation files found</div>
                        </div>
                    `;
                }
            }

            return `
                            <div class="evidence-item" onclick="showInJsonPanel([${pathStr}], 'evidence', '${evidenceId}')">
                                <div class="evidence-header">
                                    <span class="evidence-type ${getEvidenceTypeClass(evidence.type)}">${evidence.type || 'Unknown'}</span>
                                </div>
                                <div class="evidence-citation-row">
                                    <div class="evidence-citation">${evidence.citation}</div>
                                    <button class="local-file-btn" onclick="event.stopPropagation(); viewLocalEvidence('${evidenceId}')" id="localBtn${evidenceId}" title="View local evidence file">
                                        <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                                            <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path>
                                            <polyline points="14,2 14,8 20,8"></polyline>
                                            <line x1="16" y1="13" x2="8" y2="13"></line>
                                            <line x1="16" y1="17" x2="8" y2="17"></line>
                                            <polyline points="10,9 9,9 8,9"></polyline>
                                        </svg>
                                        📄
                                    </button>
                                </div>
                            </div>
                            ${simulationFilesHtml}
                        `;
        }).join('')}
                </div>
            </div>
            ` : ''}
        `;
        container.appendChild(explanationEl);
    });

    // Update local file button texts
    updateLocalFileButtonTexts();
}

// Reset API key function (kept for compatibility)
function resetApiKey() {
    alert('API functionality has been replaced with local file viewing.');
}

// Update local file button texts
function updateLocalFileButtonTexts() {
    // Update button titles for local file viewing
    document.querySelectorAll('.local-file-btn').forEach(btn => {
        btn.title = 'View local evidence file';
    });
}

// Toggle explanation content
function toggleExplanationContent(index) {
    const content = document.getElementById(`explanationContent${index}`);
    content.classList.toggle('expanded');
}

// Mark explanation
function markExplanation(index, status) {
    const item = document.querySelectorAll('.explanation-item')[index];
    const buttons = item.querySelectorAll('.btn');

    // Reset all buttons
    buttons.forEach(btn => {
        btn.classList.remove('checked', 'issue');
    });

    // Mark the clicked button
    event.target.classList.add(status);
}

// Toggle evidence section
function toggleEvidenceSection(index) {
    const evidenceSection = document.getElementById(`evidenceSection${index}`);
    evidenceSection.classList.toggle('expanded');
}

// Toggle trajectory content visibility
function toggleTrajectoryContent(index) {
    const trajectorySection = document.getElementById(`trajectorySection${index}`);

    // Toggle the trajectory section display
    if (trajectorySection.style.display === 'none') {
        trajectorySection.style.display = 'block';
        trajectorySection.classList.add('expanded');
    } else {
        trajectorySection.style.display = 'none';
        trajectorySection.classList.remove('expanded');
    }
}

// Toggle trajectory details section
function toggleTrajectoryDetails(index) {
    const section = document.getElementById(`trajectorySection${index}`);
    section.classList.toggle('expanded');
}

// Format trajectory text with markdown-like formatting and handle objects
function formatTrajectoryText(text) {
    if (!text) return '';

    // Handle different data types
    if (typeof text === 'object') {
        if (Array.isArray(text)) {
            // Handle arrays of objects (like search results)
            return formatSearchResults(text);
        } else {
            // Handle single objects - show as formatted JSON
            return `<pre>${JSON.stringify(text, null, 2)}</pre>`;
        }
    }

    // Convert to string if not already
    let formattedText = String(text);

    // Apply markdown-like formatting
    formattedText = formattedText
        // Bold text: **text** -> <strong>text</strong>
        .replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>')
        // Bullet points: - text -> • text
        .replace(/^- (.+)$/gm, '• $1')
        // Preserve line breaks
        .replace(/\n/g, '<br>')
        // Handle numbered lists better
        .replace(/^(\d+\.\s)/gm, '<strong>$1</strong>');

    return formattedText;
}

// Format search results as structured cards
function formatSearchResults(results) {
    if (!Array.isArray(results)) return '';

    return results.map((result, index) => {
        if (typeof result === 'object' && result.source && result.citation) {
            return `
                <div class="search-result-card">
                    <div class="search-result-header">
                        <strong>Result ${index + 1}</strong>
                        <a href="${result.source}" target="_blank" class="source-link">🔗 Source</a>
                    </div>
                    <div class="search-result-citation">${result.citation}</div>
                    ${result.why ? `<div class="search-result-why"><strong>Why relevant:</strong> ${result.why}</div>` : ''}
                </div>
            `;
        } else {
            return `<div class="search-result-simple">${typeof result === 'object' ? JSON.stringify(result, null, 2) : String(result)}</div>`;
        }
    }).join('');
}

// Generate trajectory content for a given explanation
function generateTrajectoryContent(explanationNumber) {
    if (!trajectoryData || !mappingData) {
        return '<p>Trajectory data not available</p>';
    }

    const stepIndex = mappingData.explanation_to_trajectory[explanationNumber.toString()];
    if (stepIndex === undefined || stepIndex === null) {
        return '<p>No trajectory mapping found for this explanation</p>';
    }

    const trajectory = trajectoryData.trajectory;
    let html = '';

    const thoughtKey = `thought_${stepIndex}`;
    const obsKey = `observation_${stepIndex}`;
    const toolNameKey = `tool_name_${stepIndex}`;
    const toolArgsKey = `tool_args_${stepIndex}`;

    // Show Thought
    if (trajectory[thoughtKey]) {
        html += `
            <div class="trajectory-item">
                <div class="trajectory-step">
                    <span class="trajectory-type thought">Thought</span>
                    Step ${stepIndex}
                </div>
                <div class="trajectory-text clickable-element" onclick="showTrajectoryInJsonPanel('${thoughtKey}')">${formatTrajectoryText(trajectory[thoughtKey])}</div>
            </div>
        `;
    }

    // Show Tool (if exists for this step)
    if (trajectory[toolNameKey]) {
        const toolName = trajectory[toolNameKey];
        const toolArgs = trajectory[toolArgsKey];

        // Format tool args in a user-friendly way
        let formattedArgs = '';
        if (toolArgs && typeof toolArgs === 'object') {
            if (toolArgs.query) {
                formattedArgs += `<strong>Query:</strong> ${toolArgs.query}`;
                if (toolArgs.max_results) {
                    formattedArgs += `<br><strong>Max Results:</strong> ${toolArgs.max_results}`;
                }
                // Add other common parameters if they exist
                Object.keys(toolArgs).forEach(key => {
                    if (key !== 'query' && key !== 'max_results') {
                        formattedArgs += `<br><strong>${key}:</strong> ${toolArgs[key]}`;
                    }
                });
            } else {
                // Fallback to formatted JSON for other tool types
                formattedArgs = `<pre>${JSON.stringify(toolArgs, null, 2)}</pre>`;
            }
        } else {
            formattedArgs = String(toolArgs || 'No arguments');
        }

        html += `
            <div class="trajectory-item">
                <div class="trajectory-step">
                    <span class="trajectory-type tool">Tool</span>
                    ${toolName} (Step ${stepIndex})
                </div>
                <div class="trajectory-text clickable-element" onclick="showTrajectoryInJsonPanel('${toolArgsKey}')">${formattedArgs}</div>
            </div>
        `;
    }

    // Show Observation
    if (trajectory[obsKey]) {
        html += `
            <div class="trajectory-item">
                <div class="trajectory-step">
                    <span class="trajectory-type observation">Observation</span>
                    Step ${stepIndex}
                </div>
                <div class="trajectory-text clickable-element" onclick="showTrajectoryInJsonPanel('${obsKey}')">${formatTrajectoryText(trajectory[obsKey])}</div>
            </div>
        `;
    }

    return html || '<p>No trajectory data found for this explanation</p>';
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

    if (typeof obj !== 'object' || obj === null) {
        if (typeof obj === 'string' && obj.length > 100) {
            // Truncate long strings but keep them searchable
            const truncated = obj.substring(0, 100) + '...';
            return `<span class="json-value" data-path="${path.join('.')}" title="${obj.replace(/"/g, '&quot;')}">${JSON.stringify(truncated)}</span>`;
        }
        return `<span class="json-value" data-path="${path.join('.')}">${JSON.stringify(obj)}</span>`;
    }

    if (Array.isArray(obj)) {
        if (obj.length === 0) return '[]';
        let result = '[\n';
        obj.forEach((item, index) => {
            const itemPath = [...path, index];
            result += `${indent}  <span class="json-item" data-path="${itemPath.join('.')}">${formatJsonForHighlight(item, itemPath, depth + 1)}</span>`;
            if (index < obj.length - 1) result += ',';
            result += '\n';
        });
        result += `${indent}]`;
        return result;
    }

    const keys = Object.keys(obj);
    if (keys.length === 0) return '{}';

    let result = '{\n';
    keys.forEach((key, index) => {
        const keyPath = [...path, key];
        const value = obj[key];

        // Handle nested JSON strings
        if (typeof value === 'string' && (key === 'problem' || key === 'json_output')) {
            try {
                const parsed = JSON.parse(value);
                result += `${indent}  <span class="json-key" data-path="${keyPath.join('.')}">"${key}": ${formatJsonForHighlight(parsed, keyPath, depth + 1)}</span>`;
            } catch {
                result += `${indent}  <span class="json-key" data-path="${keyPath.join('.')}">"${key}": ${formatJsonForHighlight(value, keyPath, depth + 1)}</span>`;
            }
        } else {
            result += `${indent}  <span class="json-key" data-path="${keyPath.join('.')}">"${key}": ${formatJsonForHighlight(value, keyPath, depth + 1)}</span>`;
        }

        if (index < keys.length - 1) result += ',';
        result += '\n';
    });
    result += `${indent}}`;
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

    // Find and highlight the target element
    setTimeout(() => {
        let targetElement = null;

        // Try different path combinations
        if (targetKey === 'evidence' && typeof index === 'string') {
            // For evidence items, look for the evidence ID
            targetElement = document.querySelector(`[data-path*="evidence.${index}"]`) ||
                document.querySelector(`[data-path*="${index}"]`);
        } else if (targetKey === 'explanation' && typeof index === 'number') {
            // For explanation items
            targetElement = document.querySelector(`[data-path*="explanation.${index}"]`);
        } else {
            // For simple keys
            const pathStr = [...basePath, targetKey].join('.');
            targetElement = document.querySelector(`[data-path="${pathStr}"]`) ||
                document.querySelector(`[data-path*="${targetKey}"]`);
        }

        if (targetElement) {
            targetElement.classList.add('json-highlight');
            targetElement.scrollIntoView({ behavior: 'smooth', block: 'center' });

            // Remove highlight after 3 seconds
            setTimeout(() => {
                targetElement.classList.remove('json-highlight');
            }, 3000);
        } else {
            // Fallback: search for any element containing the key
            const fallbackElement = document.querySelector(`[data-path*="${targetKey}"]`);
            if (fallbackElement) {
                fallbackElement.classList.add('json-highlight');
                fallbackElement.scrollIntoView({ behavior: 'smooth', block: 'center' });
                setTimeout(() => {
                    fallbackElement.classList.remove('json-highlight');
                }, 3000);
            } else {
                console.log('Target element not found for:', targetKey, index);
            }
        }
    }, 300); // Wait for panel to expand
}

// Close left panel
function closeLeftPanel() {
    const leftPanel = document.getElementById('leftPanel');
    if (leftPanel) {
        leftPanel.classList.remove('expanded');
        // Clear any inline width style that may have been set during resizing
        leftPanel.style.width = '';
        console.log('Left panel closed');
    }
}

// Close right panel
function closeRightPanel() {
    console.log('closeRightPanel called');
    const rightPanel = document.getElementById('rightPanel');
    if (rightPanel) {
        console.log('Current rightPanel classes:', rightPanel.className);
        rightPanel.classList.remove('expanded');
        // Clear any inline width style that may have been set during resizing
        rightPanel.style.width = '';
        console.log('Right panel closed, new classes:', rightPanel.className);

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
    } else {
        console.log('rightPanel element not found');
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

// Get OpenAI API key with user prompt
async function getOpenAIApiKey(forcePrompt = false) {
    // First check if we already have a key stored (unless forced to prompt)
    if (userApiKey && !forcePrompt) {
        return userApiKey;
    }

    // Check if key is configured in window object or config.js
    const configuredKey = (typeof window !== 'undefined' && window.OPENAI_TOKEN) ||
        (typeof OPENAI_TOKEN !== 'undefined' && OPENAI_TOKEN);

    if (configuredKey) {
        userApiKey = configuredKey;
        apiKeyRemembered = true;
        // Update button texts since API key is configured
        updateAIButtonTexts();
        return userApiKey;
    }

    // Prompt user for API key
    const apiKey = prompt(
        'To use AI features, please enter your OpenAI API key:\n\n' +
        'You can get your API key from: https://platform.openai.com/api-keys\n\n' +
        'Note: Your key will be stored temporarily in this session for convenience.'
    );

    if (!apiKey || !apiKey.trim()) {
        return null;
    }

    // Ask if user wants to remember the key for this session
    const remember = confirm(
        'Would you like to remember this API key for the rest of this session?\n\n' +
        'If you choose "OK", you won\'t need to enter it again until you refresh the page.\n' +
        'If you choose "Cancel", you\'ll be asked each time you use AI features.'
    );

    if (remember) {
        userApiKey = apiKey.trim();
        apiKeyRemembered = true;
        // Update button texts to reflect that API key is now remembered
        updateAIButtonTexts();
    }

    return apiKey.trim();
}

// Scrape web content using multiple proxy services and fallback strategies
async function scrapeWebContent(url) {
    try {
        console.log('Attempting to scrape content from:', url);

        // Try multiple CORS proxy services with different approaches
        const proxyServices = [
            // AllOrigins - good for many sites
            {
                url: `https://api.allorigins.win/get?url=${encodeURIComponent(url)}`,
                type: 'allorigins',
                timeout: 10000
            },
            // Proxy using a different service
            {
                url: `https://thingproxy.freeboard.io/fetch/${encodeURIComponent(url)}`,
                type: 'thingproxy',
                timeout: 10000
            },
            // Another CORS proxy
            {
                url: `https://api.codetabs.com/v1/proxy?quest=${encodeURIComponent(url)}`,
                type: 'codetabs',
                timeout: 10000
            }
        ];

        for (const proxy of proxyServices) {
            try {
                console.log(`Trying proxy service: ${proxy.type} - ${proxy.url}`);

                const controller = new AbortController();
                const timeoutId = setTimeout(() => controller.abort(), proxy.timeout);

                const response = await fetch(proxy.url, {
                    headers: {
                        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.0.0 Safari/537.36',
                        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
                        'Accept-Language': 'en-US,en;q=0.5',
                        'Accept-Encoding': 'gzip, deflate',
                        'DNT': '1',
                        'Connection': 'keep-alive',
                        'Upgrade-Insecure-Requests': '1'
                    },
                    signal: controller.signal
                });

                clearTimeout(timeoutId);

                if (response.ok) {
                    let content;
                    const contentType = response.headers.get('content-type') || '';

                    if (proxy.type === 'allorigins') {
                        const data = await response.json();
                        content = data.contents;
                    } else if (contentType.includes('application/json')) {
                        // Some proxies return JSON wrapped content
                        try {
                            const data = await response.json();
                            content = data.html || data.content || data.data || JSON.stringify(data);
                        } catch {
                            content = await response.text();
                        }
                    } else {
                        content = await response.text();
                    }

                    // Try to extract meaningful content
                    const extractedContent = extractTextContent(content, url);

                    if (extractedContent && extractedContent.length > 100) {
                        console.log(`✅ Successfully scraped content via ${proxy.type}, length:`, extractedContent.length);
                        return extractedContent;
                    } else {
                        console.log(`❌ ${proxy.type} returned insufficient content:`, extractedContent?.length || 0, 'chars');
                    }
                } else {
                    console.log(`❌ ${proxy.type} failed with status:`, response.status, response.statusText);
                }
            } catch (proxyError) {
                if (proxyError.name === 'AbortError') {
                    console.log(`⏰ ${proxy.type} timed out after ${proxy.timeout}ms`);
                } else {
                    console.log(`❌ ${proxy.type} error:`, proxyError.message);
                }
                continue;
            }
        }

        // If all proxies fail, try direct fetch (may fail due to CORS)
        console.log('🔄 All proxies failed, trying direct fetch...');
        try {
            const response = await fetch(url, {
                mode: 'cors',
                headers: {
                    'User-Agent': 'Mozilla/5.0 (compatible; ContentAnalyzer/1.0)'
                }
            });
            if (response.ok) {
                const content = await response.text();
                const extractedContent = extractTextContent(content, url);
                if (extractedContent && extractedContent.length > 100) {
                    console.log('✅ Direct fetch successful, length:', extractedContent.length);
                    return extractedContent;
                }
            }
        } catch (directError) {
            console.log('❌ Direct fetch failed:', directError.message);
        }

        // Final fallback - return a message explaining the issue
        console.log('❌ All content scraping methods failed for:', url);
        return generateFallbackContent(url);

    } catch (error) {
        console.error('💥 Critical error in scrapeWebContent:', error);
        return generateFallbackContent(url);
    }
}

// Extract text content from HTML with better parsing
function extractTextContent(htmlContent, sourceUrl) {
    try {
        if (!htmlContent || typeof htmlContent !== 'string') {
            return null;
        }

        // If it's already plain text, return as-is
        if (!htmlContent.includes('<') && !htmlContent.includes('>')) {
            return htmlContent.trim();
        }

        const parser = new DOMParser();
        const doc = parser.parseFromString(htmlContent, 'text/html');

        // Remove unwanted elements
        const unwantedSelectors = [
            'script', 'style', 'nav', 'header', 'footer', 'aside',
            '.advertisement', '.ads', '.social-share', '.comments',
            '.sidebar', '.menu', '.navigation', '.popup', '.modal',
            '[class*="cookie"]', '[class*="gdpr"]', '[id*="cookie"]'
        ];

        unwantedSelectors.forEach(selector => {
            doc.querySelectorAll(selector).forEach(el => el.remove());
        });

        // Try to find main content areas
        const contentSelectors = [
            'main', 'article', '[role="main"]',
            '.main-content', '.content', '.article-body', '.post-content',
            '#main', '#content', '#article', '.abstract', '.full-text',
            '.c-article-body', '.article-content', '.paper-content'
        ];

        let textContent = '';

        for (const selector of contentSelectors) {
            const mainContent = doc.querySelector(selector);
            if (mainContent) {
                textContent = mainContent.innerText || mainContent.textContent || '';
                if (textContent.length > 200) {
                    break;
                }
            }
        }

        // Fallback to body if no main content found
        if (!textContent || textContent.length < 200) {
            textContent = doc.body?.innerText || doc.body?.textContent || '';
        }

        // Clean up the text
        textContent = textContent
            .replace(/\s+/g, ' ')
            .replace(/\n\s*\n\s*/g, '\n')
            .replace(/^\s+|\s+$/g, '')
            .trim();

        // Filter out very short or repetitive content
        if (textContent.length < 100) {
            return null;
        }

        console.log(`📝 Extracted ${textContent.length} characters from ${sourceUrl}`);
        return textContent;

    } catch (error) {
        console.error('Error extracting text content:', error);
        return null;
    }
}

// Generate fallback content when scraping fails
function generateFallbackContent(url) {
    const domain = new URL(url).hostname;
    return `Unable to directly access content from ${domain}. This appears to be a restricted academic resource that requires institutional access or login credentials. 

The document URL is: ${url}

Common reasons for access restrictions:
• Paywall or subscription required
• Institutional login needed
• CORS/security policies blocking automated access
• Geographic restrictions

To analyze this content, you could:
1. Access the article manually and copy relevant excerpts
2. Use institutional access if available
3. Search for open access versions of the same research
4. Check if the authors have posted preprints elsewhere

The AI analysis will proceed with the available metadata and citation information.`;
}

// Analyze content with AI to find relevant passages
async function analyzeContentWithAI(content, explanationText, citation, apiKey, sourceUrl) {
    try {
        // Determine if we have actual content or fallback message
        const isAccessRestricted = content && content.includes('Unable to directly access content');

        // Truncate content if too long (GPT has token limits)
        const maxContentLength = 8000; // Approximately 2000 tokens
        const truncatedContent = content && content.length > maxContentLength
            ? content.substring(0, maxContentLength) + '...[content truncated]'
            : content || '';

        // Adjust system message based on whether we have content
        const systemMessage = isAccessRestricted
            ? `You are a scientific research assistant. When you cannot access the full content of a research paper, 
               provide helpful analysis based on the citation, URL, and explanation provided. 
               
               Return your response as a JSON object with this structure:
               {
                   "relevant_passages": [
                       {
                           "sentence": "Analysis based on citation and context",
                           "context_before": "Background information",
                           "context_after": "Implications and relevance",
                           "relevance_score": 0.8,
                           "note": "Analysis based on citation metadata due to access restrictions"
                       }
                   ],
                   "access_note": "Content could not be directly accessed due to restrictions"
               }`
            : `You are an expert at finding relevant passages in scientific documents. 
               Given a document content, an explanation, and a citation, find the most relevant sentences that support the explanation.
               
               Return your response as a JSON object with this structure:
               {
                   "relevant_passages": [
                       {
                           "sentence": "The exact sentence that supports the explanation",
                           "context_before": "1-2 sentences before for context", 
                           "context_after": "1-2 sentences after for context",
                           "relevance_score": 0.9
                       }
                   ]
               }
               
               Only include sentences with high relevance (score > 0.7). Maximum 5 passages.`;

        const userMessage = isAccessRestricted
            ? `I need help analyzing a research paper that I cannot directly access due to access restrictions.

**Explanation to support:** "${explanationText}"

**Citation:** ${citation || 'N/A'}

**Source URL:** ${sourceUrl}

**Access Status:** ${truncatedContent}

Based on the citation information and context, please provide an analysis of how this source likely supports the explanation. Include what kind of content or findings you would expect from this type of research.`
            : `Find passages in this document that support the given explanation:

**Explanation to support:** "${explanationText}"

**Expected citation:** ${citation || 'N/A'}

**Document content:**
${truncatedContent}

Please identify the most relevant sentences and provide context around them.`;

        const response = await fetch('https://api.openai.com/v1/chat/completions', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${apiKey}`
            },
            body: JSON.stringify({
                model: 'gpt-4o',
                messages: [
                    {
                        role: 'system',
                        content: systemMessage
                    },
                    {
                        role: 'user',
                        content: userMessage
                    }
                ],
                max_tokens: isAccessRestricted ? 800 : 1500,
                temperature: 0.1
            })
        });

        if (response.ok) {
            const result = await response.json();
            const analysisText = result.choices[0]?.message?.content;

            try {
                const analysisData = JSON.parse(analysisText);
                return analysisData;
            } catch (parseError) {
                console.error('Failed to parse AI response as JSON:', parseError);
                // Fallback: create a simple response
                return {
                    relevant_passages: [{
                        sentence: analysisText.substring(0, 200) + '...',
                        context_before: '',
                        context_after: '',
                        relevance_score: 0.8
                    }]
                };
            }
        } else {
            console.error('OpenAI API error:', response.status, response.statusText);
            return null;
        }
    } catch (error) {
        console.error('Error in AI analysis:', error);
        return null;
    }
}

// Display content in reading mode
function displayReadingMode(analysisResult, evidence, fromCache, extensionInfo = null) {
    const readingModeViewer = document.getElementById('readingModeViewer');

    let html = `
        <div class="source-info">
            <div class="source-title">Source Analysis</div>
            <a href="${evidence.source}" target="_blank" class="source-url">${evidence.source}</a>
            ${evidence.citation ? `<div style="margin-top: 0.5rem; font-style: italic; color: #666;">${evidence.citation}</div>` : ''}
        </div>
    `;

    if (fromCache) {
        html += '<div class="cache-indicator">Cached</div>';
    }

    if (extensionInfo && extensionInfo.useExtension) {
        html += `
            <div style="background: #d5f4e6; border: 1px solid #27ae60; border-radius: 6px; padding: 1rem; margin-bottom: 1.5rem;">
                <div style="font-weight: 600; color: #27ae60; margin-bottom: 0.5rem;">🎯 Chrome Extension Analysis</div>
                <div style="color: #155724; font-size: 0.9rem;">Text has been highlighted directly in the source document tab. The extension analyzed the full page content and highlighted relevant passages.</div>
            </div>
        `;
    }

    // Show access note if content was restricted
    if (analysisResult.access_note) {
        html += `
            <div style="background: #fff3cd; border: 1px solid #ffeaa7; border-radius: 6px; padding: 1rem; margin-bottom: 1.5rem;">
                <div style="font-weight: 600; color: #856404; margin-bottom: 0.5rem;">📋 Access Restricted</div>
                <div style="color: #856404; font-size: 0.9rem;">${analysisResult.access_note}</div>
                <div style="margin-top: 0.5rem;">
                    <button onclick="openInNewTab('${evidence.source}')" style="background: #ffc107; color: #000; border: none; padding: 0.5rem 1rem; border-radius: 4px; cursor: pointer; font-size: 0.85rem;">
                        📖 Open in New Tab
                    </button>
                </div>
            </div>
        `;
    }

    if (analysisResult.relevant_passages && analysisResult.relevant_passages.length > 0) {
        const title = analysisResult.access_note ? 'AI Analysis Based on Citation' : 'Relevant Passages';
        html += `<h3>${title}</h3>`;

        analysisResult.relevant_passages.forEach((passage, index) => {
            html += `
                <div class="relevant-passage">
                    ${passage.context_before ? `<div class="passage-context before">${passage.context_before}</div>` : ''}
                    <div class="highlighted-sentence">
                        ${passage.sentence}
                        ${passage.note ? `<div style="font-size: 0.85rem; color: #666; margin-top: 0.5rem; font-style: italic;">💡 ${passage.note}</div>` : ''}
                    </div>
                    ${passage.context_after ? `<div class="passage-context after">${passage.context_after}</div>` : ''}
                    ${passage.relevance_score ? `<div style="font-size: 0.8rem; color: #888; margin-top: 0.5rem;">Relevance: ${(passage.relevance_score * 100).toFixed(0)}%</div>` : ''}
                </div>
            `;
        });
    } else {
        html += '<div class="error-message">No relevant analysis could be generated for this source.</div>';
    }

    // Add helpful suggestions for restricted content
    if (analysisResult.access_note) {
        html += `
            <div style="background: #e8f4fd; border: 1px solid #3498db; border-radius: 6px; padding: 1rem; margin-top: 1.5rem;">
                <div style="font-weight: 600; color: #2c3e50; margin-bottom: 0.5rem;">💡 Suggestions</div>
                <ul style="color: #2c3e50; font-size: 0.9rem; margin: 0; padding-left: 1.2rem;">
                    <li>Try accessing through your institution's library</li>
                    <li>Search for open access versions on arXiv, ResearchGate, or author websites</li>
                    <li>Use Google Scholar to find related open access papers</li>
                    <li>Contact authors directly for preprints</li>
                </ul>
            </div>
        `;
    }

    readingModeViewer.innerHTML = html;
    readingModeViewer.style.display = 'block';
}

// Show error in reading mode
function showErrorInReadingMode(message) {
    const readingModeViewer = document.getElementById('readingModeViewer');
    readingModeViewer.innerHTML = `
        <div class="error-message">
            <strong>Unable to analyze source:</strong><br>
            ${message}
        </div>
    `;
    readingModeViewer.style.display = 'block';
}

// Open URL in new tab
function openInNewTab(url) {
    window.open(url, '_blank');
}

// View local evidence file
async function viewLocalEvidence(evidenceId) {
    try {
        // Load evidence source mapping
        const runType = getRunTypeFromCurrentFolder();
        const problemFolder = getProblemFolderFromCurrentFolder();

        if (!runType || !problemFolder) {
            showErrorInReadingMode('Cannot determine current folder context');
            return;
        }

        // Load evidence_source.json from the same folder as assessment.json
        const evidenceSourceData = await loadJSONFromAPI(problemFolder, 'evidence_source.json');

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

        const runType = getRunTypeFromCurrentFolder();
        const problemFolder = getProblemFolderFromCurrentFolder();
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

// Sample data functions for demo purposes
function getSampleAssessmentData() {
    return {
        "type": "assessment",
        "problem_id": "alloys_0003",
        "likert_score": -1,
        "continuous_score": 0.08,
        "confidence": 0.62,
        "explanation": [
            {
                "text": "The claim about WE43 magnesium alloy demonstrates several promising properties but requires verification of the specific oxygen concentration levels and fracture toughness values mentioned.",
                "evidence": ["ev_oxygen_study", "ev_fracture_test"]
            },
            {
                "text": "Laser powder bed fusion under high partial pressure oxygen is a novel approach that could theoretically achieve the claimed oxygen concentrations, though this requires careful validation.",
                "evidence": ["ev_lpbf_process", "ev_oxygen_uptake"]
            }
        ],
        "evidence": {
            "ev_oxygen_study": {
                "type": "literature review",
                "source": "https://example.com/oxygen-magnesium-study",
                "citation": "Smith et al. (2023). Oxygen uptake in magnesium alloys during additive manufacturing. Journal of Materials Science, 58(12), 1234-1245."
            },
            "ev_fracture_test": {
                "type": "experimental data",
                "source": "https://example.com/fracture-toughness-data",
                "citation": "Johnson, R. (2023). Fracture toughness measurements in WE43 alloy variants. Materials Testing, 45(8), 567-578."
            },
            "ev_lpbf_process": {
                "type": "simulation",
                "source": "https://example.com/lpbf-simulation",
                "citation": "Chen, L. et al. (2023). Computational modeling of oxygen transport in laser powder bed fusion. Additive Manufacturing, 34, 101234."
            },
            "ev_oxygen_uptake": {
                "type": "web search",
                "source": "https://example.com/oxygen-concentration-analysis",
                "citation": "Materials Research Database. High oxygen concentration effects in magnesium alloys. Retrieved 2023."
            }
        }
    };
}

function getSampleTrajectoryData() {
    return {
        "trajectory": {
            "thought_1": "I need to analyze the claim about WE43 magnesium alloy with specific oxygen concentrations and fracture toughness values. Let me start by understanding the baseline properties of WE43 alloy.",
            "tool_name_1": "search",
            "tool_args_1": {
                "query": "WE43 magnesium alloy baseline properties fracture toughness",
                "max_results": 5
            },
            "observation_1": "Found baseline WE43 properties: typical fracture toughness ranges from 14-18 MPa·m^0.5, density ~1.84 g/cm³, commonly used for biomedical applications.",
            "thought_2": "The claim mentions 10000 ppm oxygen concentration which is quite high for magnesium alloys. I should investigate if this level is achievable and what effects it might have.",
            "tool_name_2": "search",
            "tool_args_2": {
                "query": "magnesium alloy oxygen concentration 10000 ppm effects mechanical properties",
                "max_results": 3
            },
            "observation_2": "High oxygen concentrations in magnesium can increase hardness but may reduce ductility. 10000 ppm is at the upper end of typical ranges for processed magnesium alloys.",
            "thought_3": "Now I need to examine if laser powder bed fusion under high oxygen partial pressure can achieve these concentrations and maintain the claimed fracture toughness.",
            "tool_name_3": "search",
            "tool_args_3": {
                "query": "laser powder bed fusion magnesium oxygen partial pressure additive manufacturing",
                "max_results": 4
            },
            "observation_3": "LPBF of magnesium alloys typically requires controlled atmosphere. High oxygen partial pressure during processing could lead to increased oxygen pickup, but maintaining fracture toughness at 16 MPa·m^0.5 with high oxygen content needs verification."
        }
    };
}

function getSampleMappingData() {
    return {
        "explanation_to_trajectory": {
            "1": 1,
            "2": 3
        }
    };
}

