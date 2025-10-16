// Code Explanation Module - Handles displaying Python code with explanations

// Cache for explanation data
let explanationCache = new Map();

// Check if a file is an explanation file that should be hidden from evidence display
function isExplanationFile(filename) {
    return filename.endsWith('_explanation.json');
}

// Get the explanation filename for a Python file
function getExplanationFilename(pythonFilename) {
    if (!pythonFilename.endsWith('.py')) {
        return null;
    }
    const baseName = pythonFilename.replace(/\.py$/, '');
    return `${baseName}_explanation.json`;
}

// Load explanation data for a Python file
async function loadExplanationData(pythonFilename) {
    const explanationFilename = getExplanationFilename(pythonFilename);
    if (!explanationFilename) {
        return null;
    }

    // Check cache first
    if (explanationCache.has(explanationFilename)) {
        return explanationCache.get(explanationFilename);
    }

    try {
        const runType = getRunTypeFromCurrentFolder();
        const problemFolder = getProblemFolderFromCurrentFolder();

        if (!runType || !problemFolder) {
            console.warn('Cannot determine run type or problem folder for explanation loading');
            return null;
        }

        const response = await fetch(`/api/simulation-file/${runType}/${problemFolder}/${explanationFilename}`);
        if (!response.ok) {
            console.log(`No explanation file found: ${explanationFilename}`);
            return null;
        }

        const explanationData = await response.json();

        // Cache the data
        explanationCache.set(explanationFilename, explanationData);

        return explanationData;

    } catch (error) {
        console.warn(`Error loading explanation for ${pythonFilename}:`, error);
        return null;
    }
}

// Find matching code snippet for a line of code
function findMatchingSnippet(codeLine, snippets) {
    if (!snippets || !Array.isArray(snippets)) {
        return null;
    }

    for (const snippet of snippets) {
        if (snippet.code && snippet.code.includes(codeLine.trim())) {
            return snippet;
        }
    }
    return null;
}

// Create explanation overlay for a code snippet
function createExplanationOverlay(snippet, lineNumber, snippetNumber = null) {
    const overlay = document.createElement('div');
    overlay.className = 'code-explanation-overlay';
    overlay.style.cssText = `
        position: absolute;
        right: -320px;
        top: 0;
        width: 400px;
        background: rgba(255, 255, 255, 0.95);
        border: 1px solid #e1e8ed;
        border-radius: 8px;
        padding: 12px;
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        backdrop-filter: blur(5px);
        z-index: 1000;
        font-size: 0.85rem;
        line-height: 1.4;
        word-wrap: break-word;
        overflow-wrap: break-word;
        white-space: normal;
        visibility: hidden;
    `;

    // Create explanation content
    let explanationHtml = snippet.explanation;

    overlay.innerHTML = explanationHtml;
    return overlay;
}

// Apply explanations to Python code display
function applyCodeExplanations(codeElement, explanationData) {
    if (!explanationData || !explanationData.snippets) {
        return;
    }

    const snippets = explanationData.snippets;
    const fullCode = codeElement.textContent;
    const codeLines = fullCode.split('\n');

    // Create a wrapper for positioning overlays
    const wrapper = document.createElement('div');
    wrapper.style.cssText = 'position: relative; overflow: visible;';

    // Move the original code element inside the wrapper
    codeElement.parentNode.insertBefore(wrapper, codeElement);
    wrapper.appendChild(codeElement);

    // Calculate actual line height from the rendered code
    const computedStyle = window.getComputedStyle(codeElement);
    const fontSize = parseFloat(computedStyle.fontSize) || 14;
    const lineHeightStyle = computedStyle.lineHeight;

    let actualLineHeight;
    if (lineHeightStyle === 'normal') {
        actualLineHeight = fontSize * 1.4; // Match our CSS
    } else if (lineHeightStyle.endsWith('px')) {
        actualLineHeight = parseFloat(lineHeightStyle);
    } else if (!isNaN(parseFloat(lineHeightStyle))) {
        actualLineHeight = fontSize * parseFloat(lineHeightStyle);
    } else {
        actualLineHeight = fontSize * 1.4; // Fallback
    }

    // More precise line position calculation
    const measureLinePosition = (lineIndex) => {
        // Account for the padding-top of the pre element
        const preElement = codeElement.closest('pre');
        const preStyles = window.getComputedStyle(preElement);
        const paddingTop = parseFloat(preStyles.paddingTop) || 0;

        return paddingTop + (lineIndex * actualLineHeight);
    };

    // Process each snippet to find its position in the code
    snippets.forEach((snippet, snippetIndex) => {
        if (!snippet.code) {
            return;
        }

        // Normalize both the snippet and full code for better matching
        const normalizedSnippet = snippet.code.trim();
        const snippetLines = normalizedSnippet.split('\n');

        if (snippetLines.length === 0) return;

        // Try to find the snippet in the full code
        let matchStartIndex = -1;
        let matchEndIndex = -1;

        // Method 1: Look for exact multi-line match
        for (let i = 0; i <= codeLines.length - snippetLines.length; i++) {
            let isMatch = true;
            for (let j = 0; j < snippetLines.length; j++) {
                const codeLine = codeLines[i + j] ? codeLines[i + j].trim() : '';
                const snippetLine = snippetLines[j].trim();

                if (snippetLine && codeLine !== snippetLine) {
                    isMatch = false;
                    break;
                }
            }

            if (isMatch) {
                matchStartIndex = i;
                matchEndIndex = i + snippetLines.length - 1;
                break;
            }
        }

        // Method 2: If no exact match, try to find the first line and estimate range
        if (matchStartIndex === -1 && snippetLines.length > 0) {
            const firstLine = snippetLines[0].trim();
            if (firstLine) {
                for (let i = 0; i < codeLines.length; i++) {
                    if (codeLines[i].trim() === firstLine) {
                        matchStartIndex = i;
                        matchEndIndex = Math.min(i + snippetLines.length - 1, codeLines.length - 1);
                        break;
                    }
                }
            }
        }

        // Method 3: Use fuzzy matching for partial matches
        if (matchStartIndex === -1) {
            const firstLine = snippetLines[0].trim();
            if (firstLine.length > 10) { // Only for substantial lines
                for (let i = 0; i < codeLines.length; i++) {
                    const codeLine = codeLines[i].trim();
                    if (codeLine.includes(firstLine.substring(0, Math.min(20, firstLine.length)))) {
                        matchStartIndex = i;
                        matchEndIndex = Math.min(i + snippetLines.length - 1, codeLines.length - 1);
                        break;
                    }
                }
            }
        }

        if (matchStartIndex !== -1) {

            // Calculate precise positioning using actual line height
            const startLinePosition = measureLinePosition(matchStartIndex - 1);
            const endLinePosition = measureLinePosition(matchEndIndex);
            const snippetHeightPx = endLinePosition - startLinePosition;

            // Create and position the explanation overlay
            const overlay = createExplanationOverlay(snippet, matchStartIndex, snippetIndex + 1);

            // Position the overlay using pixel values for precision
            overlay.style.top = `${startLinePosition}px`;

            // Calculate optimal height: max of content height and snippet height
            // First, temporarily show overlay to measure content height
            overlay.style.visibility = 'visible';
            overlay.style.opacity = '0';
            wrapper.appendChild(overlay);

            const contentHeight = overlay.scrollHeight;
            const optimalHeight = Math.max(contentHeight, snippetHeightPx);

            overlay.style.height = `${optimalHeight}px`;
            overlay.style.display = 'flex';
            overlay.style.alignItems = 'center';

            // Hide overlay again
            overlay.style.visibility = 'hidden';
            wrapper.removeChild(overlay);

            // Add hover interaction
            let timeout;
            const showOverlay = () => {
                clearTimeout(timeout);
                overlay.style.opacity = '1';
                overlay.style.transform = 'translateX(0)';
                overlay.style.visibility = 'visible';
                overlay.style.right = '0px';
            };

            const hideOverlay = () => {
                timeout = setTimeout(() => {
                    overlay.style.opacity = '0';
                    overlay.style.transform = 'translateX(10px)';
                    setTimeout(() => {
                        overlay.style.visibility = 'hidden';
                    }, 300);
                }, 100);
            };

            // Initially hide the overlay
            overlay.style.opacity = '0';
            overlay.style.transform = 'translateX(10px)';
            overlay.style.visibility = 'hidden';
            overlay.style.transition = 'opacity 0.3s ease, transform 0.3s ease';

            // Add event listeners for hover
            overlay.addEventListener('mouseenter', showOverlay);
            overlay.addEventListener('mouseleave', hideOverlay);

            // Create a trigger area over the code lines using precise positioning
            const trigger = document.createElement('div');
            trigger.style.cssText = `
                position: absolute;
                left: 0;
                right: 0;
                height: ${Math.max(snippetHeightPx, actualLineHeight)}px;
                top: ${startLinePosition}px;
                z-index: 999;
                background: rgba(52, 152, 219, 0.05);
                border-left: 3px solid #3498db;
                cursor: help;
                opacity: 0;
                transition: opacity 0.2s ease;
                border-radius: 2px;
            `;

            // Add snippet indicator
            const indicator = document.createElement('div');
            indicator.style.cssText = `
                position: absolute;
                right: 8px;
                top: 2px;
                background: #3498db;
                color: white;
                border-radius: 10px;
                width: 20px;
                height: 20px;
                display: flex;
                align-items: center;
                justify-content: center;
                font-size: 11px;
                font-weight: bold;
                opacity: 0;
                transition: opacity 0.2s ease;
                pointer-events: none;
            `;
            indicator.textContent = snippetIndex + 1;

            trigger.addEventListener('mouseenter', () => {
                trigger.style.opacity = '0.8';
                indicator.style.opacity = '1';
                showOverlay();
            });

            trigger.addEventListener('mouseleave', () => {
                trigger.style.opacity = '0';
                indicator.style.opacity = '0';
                hideOverlay();
            });

            trigger.appendChild(indicator);
            wrapper.appendChild(trigger);
            wrapper.appendChild(overlay);
        }
    });
}

// Enhanced show simulation file function with explanation support
async function showSimulationFileWithExplanation(filename) {
    const runType = getRunTypeFromCurrentFolder();
    const problemFolder = getProblemFolderFromCurrentFolder();

    if (!runType || !problemFolder) {
        console.error('Cannot determine run type or problem folder');
        return;
    }

    // Check if this file has explanations available
    const hasExplanations = filename.endsWith('.py');
    let explanationData = null;

    if (hasExplanations) {
        explanationData = await loadExplanationData(filename);
    }

    // Show right panel with appropriate width
    const rightPanel = document.getElementById('rightPanel');
    if (!rightPanel.classList.contains('expanded')) {
        rightPanel.classList.add('expanded');
    }

    // Adjust panel width for files with explanations
    if (explanationData) {
        rightPanel.style.width = '800px';
    } else {
        rightPanel.style.width = ''; // Use default width
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

        // Check for visualization
        const baseName = filename.replace(/\.[^/.]+$/, '');
        const visualizationFilename = `${baseName}.png`;

        let visualizationHtml = '';
        try {
            const vizResponse = await fetch(`/api/simulation-file/${runType}/${problemFolder}/${visualizationFilename}`);
            if (vizResponse.ok) {
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
            }
        } catch (vizError) {
            // No visualization available
        }

        let contentHtml = '';

        if (filename.endsWith('.csv')) {
            contentHtml = formatCsvAsTable(content, filename);
        } else {
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
                <div class="content-display" style="position: relative; overflow: visible;">
                    <pre style="position: relative; overflow-x: auto; margin-right: 0;"><code class="language-${language}" id="codeContent">${displayContent}</code></pre>
                </div>
            `;
        }

        // Create header with explanation indicator
        let headerInfo = `Location: ${runType}/${problemFolder}/${filename}`;
        if (explanationData) {
            headerInfo += ` • 💡 ${explanationData.snippets ? explanationData.snippets.length : 0} explanations available`;
        }

        // Display in reading mode viewer
        readingModeViewer.innerHTML = `
            <div class="source-info">
                <div class="source-title">
                    ${explanationData ? '🐍📖' : '📄'} ${filename}
                    ${explanationData ? '<span style="color: #3498db; font-size: 0.8em; margin-left: 8px;">Interactive Explanations</span>' : ''}
                </div>
                <div style="margin-top: 0.5rem; font-style: italic; color: #666;">
                    ${headerInfo}
                </div>
                ${explanationData ? `
                    <div style="margin-top: 0.5rem; padding: 8px 12px; background: #e3f2fd; border-radius: 4px; font-size: 0.85rem; color: #1565c0;">
                        Hover over highlighted code sections to see explanations
                        ${explanationData.summary ? `<br><strong>Summary:</strong> ${explanationData.summary}` : ''}
                    </div>
                ` : ''}
            </div>
            <div style="margin-top: 1rem;">
                ${visualizationHtml}
                ${contentHtml}
            </div>
        `;

        readingModeViewer.style.display = 'block';

        // Apply explanations if available
        if (explanationData && filename.endsWith('.py')) {
            setTimeout(() => {
                const codeElement = document.getElementById('codeContent');
                if (codeElement) {
                    applyCodeExplanations(codeElement, explanationData);
                }
            }, 100);
        }

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

// Export functions for global use
window.isExplanationFile = isExplanationFile;
window.showSimulationFileWithExplanation = showSimulationFileWithExplanation;
window.loadExplanationData = loadExplanationData;