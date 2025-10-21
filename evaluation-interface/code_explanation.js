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
        z-index: 500;
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
function applyCodeExplanations(codeElement, explanationData, showByDefault = true) {
    if (!explanationData || !explanationData.snippets) {
        return;
    }

    const snippets = explanationData.snippets;
    const fullCode = codeElement.textContent;
    const codeLines = fullCode.split('\n');

    // Create a wrapper for positioning overlays
    const wrapper = document.createElement('div');
    wrapper.style.cssText = 'position: relative;';

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

            // Initially set overlay visibility based on mode
            if (showByDefault) {
                overlay.style.opacity = '1';
                overlay.style.transform = 'translateX(0)';
                overlay.style.visibility = 'visible';
                overlay.style.right = '0px';
                overlay.style.zIndex = '500';
            } else {
                overlay.style.opacity = '0';
                overlay.style.transform = 'translateX(10px)';
                overlay.style.visibility = 'hidden';
                overlay.style.zIndex = '500';
            }
            overlay.style.transition = 'opacity 0.3s ease, transform 0.3s ease';

            // Hover functions for z-index management
            const bringToFront = () => {
                overlay.style.zIndex = '600';
                trigger.style.zIndex = '599';
            };

            const sendToBack = () => {
                overlay.style.zIndex = '500';
                trigger.style.zIndex = '499';
            };

            // Add event listeners based on mode
            if (!showByDefault) {
                // Hover to show mode
                overlay.addEventListener('mouseenter', showOverlay);
                overlay.addEventListener('mouseleave', hideOverlay);
            } else {
                // Always show mode - just manage z-index
                overlay.addEventListener('mouseenter', bringToFront);
                overlay.addEventListener('mouseleave', sendToBack);
            }

            // Create a trigger area over the code lines using precise positioning
            const trigger = document.createElement('div');
            trigger.style.cssText = `
                position: absolute;
                left: -10px;
                right: 0;
                height: ${Math.max(snippetHeightPx, actualLineHeight)}px;
                top: ${startLinePosition}px;
                z-index: 499;
                background: rgba(52, 152, 219, 0.05);
                border-left: 3px solid #3498db;
                cursor: help;
                opacity: ${showByDefault ? '0.8' : '0'};
                transition: opacity 0.2s ease;
                border-radius: 2px;
            `;

            // Store reference for toggling
            trigger.explanationOverlay = overlay;
            trigger.showByDefault = showByDefault;

            if (!showByDefault) {
                // Hover to show mode
                trigger.addEventListener('mouseenter', () => {
                    trigger.style.opacity = '0.8';
                    showOverlay();
                });

                trigger.addEventListener('mouseleave', () => {
                    trigger.style.opacity = '0';
                    hideOverlay();
                });
            } else {
                // Always show mode - manage z-index on hover
                trigger.addEventListener('mouseenter', () => {
                    overlay.style.zIndex = '600';
                    trigger.style.zIndex = '599';
                });

                trigger.addEventListener('mouseleave', () => {
                    overlay.style.zIndex = '500';
                    trigger.style.zIndex = '499';
                });
            }

            wrapper.appendChild(trigger);
            wrapper.appendChild(overlay);
        }
    });
}

// Toggle explanation display mode
function toggleExplanationMode(showByDefault) {
    const codeElement = document.getElementById('codeContent');
    if (!codeElement) return;

    const wrapper = codeElement.closest('div[style*="position: relative"]');
    if (!wrapper) return;

    const triggers = wrapper.querySelectorAll('div[style*="z-index: 499"]');
    triggers.forEach(trigger => {
        const overlay = trigger.explanationOverlay;
        if (overlay) {
            // Remove all existing event listeners by cloning elements
            const newTrigger = trigger.cloneNode(true);
            const newOverlay = overlay.cloneNode(true);

            // Update references
            newTrigger.explanationOverlay = newOverlay;

            if (showByDefault) {
                // Show all overlays
                newTrigger.style.opacity = '0.8';
                newOverlay.style.opacity = '1';
                newOverlay.style.transform = 'translateX(0)';
                newOverlay.style.visibility = 'visible';
                newOverlay.style.right = '0px';
                newOverlay.style.zIndex = '500';
                newTrigger.style.zIndex = '499';

                // Add z-index management for always show mode
                newTrigger.addEventListener('mouseenter', () => {
                    newOverlay.style.zIndex = '600';
                    newTrigger.style.zIndex = '599';
                });

                newTrigger.addEventListener('mouseleave', () => {
                    newOverlay.style.zIndex = '500';
                    newTrigger.style.zIndex = '499';
                });

                newOverlay.addEventListener('mouseenter', () => {
                    newOverlay.style.zIndex = '600';
                    newTrigger.style.zIndex = '599';
                });

                newOverlay.addEventListener('mouseleave', () => {
                    newOverlay.style.zIndex = '500';
                    newTrigger.style.zIndex = '499';
                });
            } else {
                // Hide all overlays and enable hover to show
                newTrigger.style.opacity = '0';
                newOverlay.style.opacity = '0';
                newOverlay.style.transform = 'translateX(10px)';
                newOverlay.style.visibility = 'hidden';
                newOverlay.style.zIndex = '500';
                newTrigger.style.zIndex = '499';

                // Add hover to show functionality
                const showOverlay = () => {
                    newOverlay.style.opacity = '1';
                    newOverlay.style.transform = 'translateX(0)';
                    newOverlay.style.visibility = 'visible';
                    newOverlay.style.right = '0px';
                    newOverlay.style.zIndex = '600';
                    newTrigger.style.zIndex = '599';
                };

                const hideOverlay = () => {
                    setTimeout(() => {
                        newOverlay.style.opacity = '0';
                        newOverlay.style.transform = 'translateX(10px)';
                        newOverlay.style.zIndex = '500';
                        newTrigger.style.zIndex = '499';
                        setTimeout(() => {
                            newOverlay.style.visibility = 'hidden';
                        }, 300);
                    }, 100);
                };

                newTrigger.addEventListener('mouseenter', () => {
                    newTrigger.style.opacity = '0.8';
                    showOverlay();
                });

                newTrigger.addEventListener('mouseleave', () => {
                    newTrigger.style.opacity = '0';
                    hideOverlay();
                });

                newOverlay.addEventListener('mouseenter', showOverlay);
                newOverlay.addEventListener('mouseleave', hideOverlay);
            }

            // Replace old elements with new ones
            trigger.parentNode.replaceChild(newTrigger, trigger);
            overlay.parentNode.replaceChild(newOverlay, overlay);
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


        // Display in reading mode viewer with sticky header
        readingModeViewer.innerHTML = `
            <div class="source-info" style="position: sticky; top: 0; background: white; z-index: 1001; padding: 1rem; border-bottom: 1px solid #e1e8ed; margin: 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1);"
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <div class="source-title">
                        ${filename}
                    </div>
                    ${hasExplanations && explanationData ? `
                        <div style="display: flex; align-items: center; gap: 10px;">
                            <span style="font-size: 0.9rem; color: #666;">Code Explanations:</span>
                            <label style="display: flex; align-items: center; cursor: pointer;">
                                <input type="checkbox" id="explanationToggle" checked style="margin-right: 5px;">
                                <span style="font-size: 0.9rem;">Always Show</span>
                            </label>
                        </div>
                    ` : ''}
                </div>
            </div>
            <div style="padding: 1rem;">
                ${contentHtml}
            </div>
        `;

        readingModeViewer.style.display = 'block';

        // Apply explanations if available
        if (explanationData && filename.endsWith('.py')) {
            setTimeout(() => {
                const codeElement = document.getElementById('codeContent');
                if (codeElement) {
                    // Apply explanations with default show mode
                    applyCodeExplanations(codeElement, explanationData, true);

                    // Add toggle functionality
                    const toggleButton = document.getElementById('explanationToggle');
                    if (toggleButton) {
                        toggleButton.addEventListener('change', (e) => {
                            toggleExplanationMode(e.target.checked);
                        });
                    }
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
if (typeof window !== 'undefined') {
    window.codeExplanation = {
        isExplanationFile,
        showSimulationFileWithExplanation,
        loadExplanationData,
        getExplanationFilename,
        toggleExplanationMode
    };

    // Keep backward compatibility for existing functions
    window.isExplanationFile = isExplanationFile;
    window.showSimulationFileWithExplanation = showSimulationFileWithExplanation;
    window.loadExplanationData = loadExplanationData;
    window.toggleExplanationMode = toggleExplanationMode;
}