// Scientific Claim Assessment Interface - Main JavaScript

// Global variables
let currentData = null;
let evidenceData = null;
let trajectoryData = null;
let mappingData = null;
let userApiKey = null;
let apiKeyRemembered = false;
let evidenceCache = new Map(); // Cache for scraped content and AI analysis
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
            console.log('Trajectory data loaded successfully');
        }

        if (mappingDataResult) {
            mappingData = mappingDataResult;
            console.log('Mapping data loaded successfully');
        }

        // Load data into trajectory flow module
        if (window.trajectoryFlow && trajectoryDataResult && mappingDataResult) {
            window.trajectoryFlow.loadData(trajectoryDataResult, mappingDataResult);
            console.log('Trajectory flow data loaded');
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

// Load simulation files from sandbox folder for a given evidence item
async function loadSandboxFiles(evidenceId) {
    // The API server automatically handles sandbox folder mapping
    // So we can use the regular loadSimulationFiles function
    return loadSimulationFiles(evidenceId);
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
        // API server will automatically check sandbox folders
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
                // Use loadSandboxFiles which will try sandbox first, then fallback to regular simulation files
                const files = await loadSandboxFiles(evidenceId);
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
            <div class="explanation-header">
                <div class="explanation-text clickable-element" onclick="showInJsonPanel([${pathStr}], 'explanation', ${index})">${explanation.text}</div>
                <div class="explanation-actions">
                    ${hasEvidence ? `<button class="btn" onclick="toggleEvidenceSection(${index})">Evidence (${validEvidence.length})</button>` : ''}
                    ${mappingData ? `<button class="btn" onclick="showTrajectoryFlow(${index})">Show Trajectory</button>` : ''}
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
                            <div class="evidence-item" id="evidence-${index}-${evidenceId}" data-evidence-id="${evidenceId}" data-explanation-index="${index}" onclick="showInJsonPanel([${pathStr}], 'evidence', '${evidenceId}')">
                                <div class="evidence-header">
                                    <span class="evidence-type ${getEvidenceTypeClass(evidence.type)}">${evidence.type || 'Unknown'}</span>
                                </div>
                                <div class="evidence-citation-row">
                                    <div class="evidence-citation">${evidence.citation}</div>
                                    ${evidence.type && evidence.type.toLowerCase() === 'web search' ? `
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
                                    ` : ''}
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

// Show trajectory flow in left panel
function showTrajectoryFlow(explanationIndex = null) {
    // Check if trajectory is currently visible
    const isCurrentlyVisible = window.trajectoryFlow && window.trajectoryFlow.isVisible();

    if (window.trajectoryFlow && typeof window.trajectoryFlow.show === 'function') {
        window.trajectoryFlow.show(explanationIndex);

        // If trajectory was just shown (not hidden), expand evidence section
        if (!isCurrentlyVisible && explanationIndex !== null) {
            console.log(`Expanding evidence section for explanation ${explanationIndex}`);
            const evidenceSection = document.getElementById(`evidenceSection${explanationIndex}`);
            if (evidenceSection && !evidenceSection.classList.contains('expanded')) {
                toggleEvidenceSection(explanationIndex);
            }
        }
    } else {
        console.error('Trajectory flow module not loaded');
    }
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