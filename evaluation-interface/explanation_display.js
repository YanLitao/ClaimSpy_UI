// Explanation Display Module - Handles displaying explanations with evidence

// Display explanations
async function displayExplanations(explanations) {
    const container = document.getElementById('explanationsList');
    container.innerHTML = '';
    const assessmentPath = (typeof getAssessmentPath === 'function') ?
        getAssessmentPath() : [];

    // Handle case where explanations might be undefined or empty
    if (!explanations || !Array.isArray(explanations)) {
        container.innerHTML = '<p style="color: #666; padding: 1rem;">No explanations available.</p>';
        return;
    }

    // Pre-load all simulation files for all simulation evidence
    const simulationFilesCache = {};
    const evidenceData = window.evidenceManager ? window.evidenceManager.getEvidenceData() : null;

    if (evidenceData) {
        const simulationEvidence = Object.keys(evidenceData).filter(evidenceId => {
            const evidence = evidenceData[evidenceId];
            return evidence && evidence.type && evidence.type.toLowerCase().includes('simulation');
        });

        // Load all simulation files at once
        for (const evidenceId of simulationEvidence) {
            try {
                // Use loadSandboxFiles which will try sandbox first, then fallback to regular simulation files
                const files = window.simulationManager ?
                    await window.simulationManager.loadSandboxFiles(evidenceId) : [];
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
        const validEvidence = (explanation.evidence || []).filter(id => evidenceData && evidenceData[id] && evidenceData[id].citation);
        const hasEvidence = validEvidence.length > 0;

        // Collect all evidence types for this explanation
        const evidenceTypes = new Set();
        validEvidence.forEach(id => {
            const evidence = evidenceData[id];
            if (evidence && evidence.type) {
                evidenceTypes.add(evidence.type);
            }
        });

        const getEvidenceTypeClass = window.evidenceManager ?
            window.evidenceManager.getEvidenceTypeClass :
            (type) => type ? type.toLowerCase().replace(/\s+/g, '-') : 'unknown';

        const typesHTML = Array.from(evidenceTypes).map(type =>
            `<div><span class="evidence-type ${getEvidenceTypeClass(type)}">${type}</span></div>`
        ).join('');

        explanationEl.innerHTML = `
            <div class="explanation-header">
                <div class="explanation-text clickable-element" onclick="showInJsonPanel([${pathStr}], 'explanation', ${index})">${explanation.text}</div>
                <div class="explanation-actions">
                    ${hasEvidence ? `<button class="btn evidence-toggle-btn" onclick="toggleEvidenceSection(${index})" style="display: flex; flex-direction: column; align-items: flex-start; text-align: left; line-height: 1.3; min-width: 185px; width: auto;">
                        <div>${validEvidence.length} evidence:</div>
                        ${typesHTML ? `<div style="margin-top: 2px;">${typesHTML}</div>` : ''}
                    </button>` : ''}
                    ${(typeof mappingData !== 'undefined' && mappingData) ? `<button class="btn" onclick="showTrajectoryFlow(${index})">Show Trajectory</button>` : ''}
                </div>
            </div>
            <div class="explanation-content" id="explanationContent${index}">
                <textarea class="comment-box" placeholder="Add your comments about this explanation..."></textarea>
            </div>
            ${hasEvidence ? `
            <div class="evidence-section" id="evidenceSection${index}">
                <div class="evidence-header-toggle" onclick="toggleEvidenceSection(${index})">
                    <h4>${validEvidence.length} evidence:</h4>
                    <span class="collapse-icon">▼</span>
                </div>
                <div id="evidenceContent${index}">
                    ${(explanation.evidence || []).map(evidenceId => {
            const evidence = evidenceData ? evidenceData[evidenceId] : null;
            if (!evidence || !evidence.citation) return '';

            const isValidUrl = window.evidenceManager ?
                window.evidenceManager.isValidUrl :
                (url) => {
                    try { new URL(url); return true; } catch { return false; }
                };
            const hasValidUrl = isValidUrl(evidence.source);

            // Generate simulation files HTML if this is simulation evidence
            let simulationFilesHtml = '';
            if (evidence.type && evidence.type.toLowerCase().includes('simulation') && simulationFilesCache[evidenceId]) {
                const files = simulationFilesCache[evidenceId];
                if (files.length > 0) {
                    const formatFileSize = window.simulationManager ?
                        window.simulationManager.formatFileSize :
                        (bytes) => `${bytes} bytes`;

                    // Create file boxes (inline display)
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

                    // Create visualizations section (display all at bottom)
                    // Group files by visualization filename to avoid duplicates
                    const visualizationGroups = {};
                    files
                        .filter(file => file.hasVisualization && file.visualizationFile)
                        .forEach(file => {
                            if (!visualizationGroups[file.visualizationFile]) {
                                visualizationGroups[file.visualizationFile] = [];
                            }
                            visualizationGroups[file.visualizationFile].push(file.name);
                        });

                    const visualizationsHtml = Object.entries(visualizationGroups)
                        .map(([visualizationFile, associatedFiles]) => {
                            const runType = (typeof getRunTypeFromCurrentFolder === 'function') ?
                                getRunTypeFromCurrentFolder() : null;
                            const problemFolder = (typeof getProblemFolderFromCurrentFolder === 'function') ?
                                getProblemFolderFromCurrentFolder() : null;

                            // Create label showing all associated files
                            const filesLabel = associatedFiles.length > 1 ?
                                associatedFiles.join(', ') :
                                associatedFiles[0];

                            return `
                            <div class="visualization-item" style="margin-bottom: 1rem;">
                                <div style="font-size: 0.85rem; color: #666; margin-bottom: 0.5rem;">📊 Visualization for ${filesLabel}</div>
                                <div style="text-align: center;">
                                    <img src="/api/simulation-file/${runType}/${problemFolder}/${visualizationFile}" 
                                         alt="Visualization for ${filesLabel}"
                                         style="max-width: 100%; height: auto; border-radius: 4px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);"
                                         onerror="this.style.display='none'">
                                </div>
                            </div>
                        `;
                        }).join('');

                    simulationFilesHtml = `
                        <div class="simulation-files-section" id="simFiles${evidenceId}" style="margin-top: 0.5rem;">
                            <div class="simulation-files-header" style="font-size: 0.9rem; color: #666; margin-bottom: 0.5rem;">Simulation Files:</div>
                            <div class="simulation-files-content" id="simFilesContent${evidenceId}" style="margin-bottom: ${visualizationsHtml ? '1rem' : '0'};">${fileBoxesHtml}</div>
                            ${visualizationsHtml ? `<div class="visualizations-section" style="border-top: 1px solid #e9ecef; padding-top: 1rem;">${visualizationsHtml}</div>` : ''}
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
    if (window.evidenceManager && window.evidenceManager.updateLocalFileButtonTexts) {
        window.evidenceManager.updateLocalFileButtonTexts();
    }
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

// Export functions for global access
if (typeof window !== 'undefined') {
    window.explanationDisplay = {
        displayExplanations,
        toggleExplanationContent,
        markExplanation,
        formatTrajectoryText,
        formatSearchResults
    };
}