// Trajectory Flow Visualization
// This module handles the display of trajectory flow diagrams in the left panel

let flowTrajectoryData = null;
let flowMappingData = null;
let isTrajectoryFlowVisible = false;
let firstLoad = true;
let firstNodeX = -1;

// Common text formatting function
function formatText(text) {
    if (typeof text !== 'string') return '';
    return text.trim()
        // Convert \n to actual line breaks
        .replace(/\\n/g, '\n')
        // Convert line breaks to <br> tags
        .replace(/\n/g, '<br>')
        // Bold for **text**
        .replace(/\*\*(.*?)\*\*/g, '<b>$1</b>')
        // Bold for ## or ### headers (only the header line, not everything after)
        .replace(/^(###?)\s*(.+?)(<br>|$)/gm, '<b>$2</b>$3')
        // Convert "- text" at start of line to bullet point
        .replace(/^- (.*)$/gm, '• $1')
        // Convert XML-like tags to bold text format, but avoid already processed HTML tags
        .replace(/<(?!\/?(b|br|i|u|strong|em)\b)([A-Za-z][A-Za-z0-9]*?)>/g, '<b>$2</b>')
        .replace(/<\/(?!(b|br|i|u|strong|em)\b)([A-Za-z][A-Za-z0-9]*?)>/g, '<b>/$2</b>')
        // Clean up extra spaces but preserve <br> tags
        .replace(/\s+/g, ' ')
        .replace(/<br>\s+/g, '<br>')
        .replace(/\s+<br>/g, '<br>');
}

// Initialize trajectory flow system
function initTrajectoryFlow() {
    // This will be called when the page loads to set up event listeners
}

// Load trajectory and mapping data
function loadTrajectoryFlowData(trajectory, mapping) {
    flowTrajectoryData = trajectory;
    flowMappingData = mapping;
}

// Show trajectory flow in left panel
function showTrajectoryFlow(explanationIndex = null) {
    if (!flowTrajectoryData || !flowMappingData) {
        console.error('Trajectory or mapping data not loaded');
        return;
    }

    const leftPanel = document.getElementById('leftPanel');
    if (!leftPanel) {
        console.error('Left panel not found');
        return;
    }

    // Check if trajectory is currently visible and toggle if so
    if (isTrajectoryFlowVisible) {
        hideTrajectoryFlow();
        return;
    }

    // Expand left panel if not already expanded
    if (!leftPanel.classList.contains('expanded')) {
        leftPanel.classList.add('expanded');
    }

    // Hide JSON content if visible
    hideJsonContent();

    // Set flag to indicate trajectory flow is visible
    isTrajectoryFlowVisible = true;

    // Create trajectory flow content
    const flowContent = generateTrajectoryFlowHTML();

    // Display in left panel
    leftPanel.innerHTML = `
        <div class="trajectory-flow-container">
            <div class="trajectory-flow-header">
                <h3>Trajectory Flow</h3>
                <button class="close-btn" onclick="hideTrajectoryFlow()">✕</button>
            </div>
            <div class="trajectory-flow-content">
                ${flowContent}
            </div>
        </div>
    `;

    // Draw arrows after a short delay to ensure elements are rendered
    setTimeout(() => {
        // First recalculate positions based on actual content
        recalculateNodePositions();
        // Then draw arrows with correct positions
        drawFlowArrows();
        // Draw evidence connections if explanationIndex is provided
        if (explanationIndex !== null) {
            // Wait a bit more for evidence section to be expanded
            setTimeout(() => {
                drawEvidenceConnections(explanationIndex);
                // Auto-scroll to relevant positions after connections are drawn
                scrollToRelevantPositions(explanationIndex);
            }, 300);
        }
        // Set up scroll listeners for dynamic line updates
        setupScrollListeners();
    }, 100);
}

// Hide trajectory flow and restore previous state
function hideTrajectoryFlow() {
    isTrajectoryFlowVisible = false;

    // Clear evidence connections
    clearEvidenceConnections();

    // Remove scroll listeners
    if (window.evidenceScrollHandler) {
        window.removeEventListener('scroll', window.evidenceScrollHandler, true);
        const leftPanel = document.getElementById('leftPanel');
        if (leftPanel) {
            leftPanel.removeEventListener('scroll', window.evidenceScrollHandler);
        }
    }

    const leftPanel = document.getElementById('leftPanel');
    if (leftPanel) {
        // Collapse the left panel
        leftPanel.classList.remove('expanded');
        // Clear the content
        leftPanel.innerHTML = '<div class="welcome-message">Select an item to view details</div>';
    }
}

// Generate HTML for trajectory flow
function generateTrajectoryFlowHTML() {
    const trajectory = flowTrajectoryData.trajectory || flowTrajectoryData;
    const stepDependencies = flowMappingData.step_dependencies || {};

    // Extract all steps from trajectory
    const steps = extractStepsFromTrajectory(trajectory);

    // Create container with SVG for arrows
    let flowHTML = `
        <div class="flow-diagram-container">
            <svg class="flow-arrows" id="flowArrows"></svg>
            <div class="flow-diagram">
    `;

    // Generate step boxes with absolute positioning
    steps.forEach((step, index) => {
        flowHTML += generateStepNodeHTML(step, index, stepDependencies);
    });

    // Add reasoning box
    const reasoning = flowTrajectoryData.reasoning || (flowTrajectoryData.trajectory && flowTrajectoryData.trajectory.reasoning);
    const answer = flowTrajectoryData.answer || (flowTrajectoryData.trajectory && flowTrajectoryData.trajectory.answer);

    let nodeIndex = steps.length;
    if (reasoning) {
        let reasoningContent = generateReasoningNodeHTML(nodeIndex);
        flowHTML += reasoningContent;
        nodeIndex++;
    }

    // Add answer box
    if (answer) {
        let answerContent = generateAnswerNodeHTML(nodeIndex);
        flowHTML += answerContent;
    }

    flowHTML += `
            </div>
        </div>
    `;

    return flowHTML;
}

// Extract steps from trajectory data
function extractStepsFromTrajectory(trajectory) {
    const steps = [];
    const stepNumbers = new Set();

    // Find all step numbers
    Object.keys(trajectory).forEach(key => {
        const match = key.match(/^(thought|tool_name|tool_args|observation)_(\d+)$/);
        if (match) {
            stepNumbers.add(parseInt(match[2]));
        }
    });

    // Create step objects
    Array.from(stepNumbers).sort((a, b) => a - b).forEach(stepNum => {
        const step = {
            number: stepNum,
            thought: trajectory[`thought_${stepNum}`] || null,
            toolName: trajectory[`tool_name_${stepNum}`] || null,
            toolArgs: trajectory[`tool_args_${stepNum}`] || null,
            observation: trajectory[`observation_${stepNum}`] || null
        };
        steps.push(step);
    });

    return steps;
}

// Generate HTML for a single step node in flow chart
function generateStepNodeHTML(step, index, stepDependencies) {
    const stepNum = step.number;
    // Use temporary positioning, will be recalculated after rendering
    const yPosition = 50 + index * 200; // Initial rough positioning

    let html = `<div class="step-node" id="step-${stepNum}" data-step="${stepNum}" style="top: ${yPosition}px;">`;

    html += `<div class="step-header">Step ${stepNum}</div>`;

    // Collapsible thoughts display
    if (step.thought) {

        // Create 3-line summary (approximately 150 characters)
        const thoughtLines = step.thought.trim().split('\n').filter(line => line.trim());
        let thoughtSummary = thoughtLines.slice(0, 3).join('\n').trim();
        if (thoughtSummary.length > 150) {
            thoughtSummary = thoughtSummary.substring(0, 150) + '...';
        } else if (thoughtLines.length > 3) {
            thoughtSummary += '...';
        }

        const needsExpansion = step.thought.trim().length > 150 || thoughtLines.length > 3;

        html += `
            <div class="step-content thought-content">
                <span class="content-icon">Thought:</span>
                <div class="thought-container">
                    <div class="thought-summary" id="thought-summary-${stepNum}">
                        ${formatText(thoughtSummary)}
                    </div>
                    ${needsExpansion ? `
                        <div class="thought-full" id="thought-full-${stepNum}" style="display: none;">
                            ${formatText(step.thought)}
                        </div>
                        <button class="thought-toggle-btn" onclick="toggleThought(${stepNum})" id="thought-toggle-${stepNum}">
                            Show More
                        </button>
                    ` : ''}
                </div>
            </div>
        `;
    }

    if (step.toolName) {
        const toolDisplayText = step.toolName.toLowerCase() === 'finish' ?
            escapeHtml(step.toolName) :
            `${escapeHtml(step.toolName)}`;

        html += `
            <div class="step-content tool-content">
                <span class="content-icon">Tool: </span>
                <span class="content-text">${toolDisplayText}</span>
            </div>
        `;
    }

    if (step.observation !== null && step.observation) {
        if (Array.isArray(step.observation) && step.observation.length > 0) {
            // Show individual observation items as small boxes
            html += `
                <div class="step-content obs-content">
                    <span class="content-icon">Findings: </span>
                    <span class="content-text"></span>
                </div>
                <div class="observation-items">
            `;

            step.observation.forEach((item, itemIndex) => {
                const isSearchResult = item && typeof item === 'object' && item.source && item.citation;
                let itemPreview = '';

                if (isSearchResult) {
                    itemPreview = item.citation.length > 25 ?
                        item.citation.substring(0, 25) + '...' : item.citation;
                } else {
                    const itemStr = JSON.stringify(item);
                    itemPreview = itemStr.length > 25 ?
                        itemStr.substring(0, 25) + '...' : itemStr;
                }

                html += `
                    <div class="obs-item-box" id="obs-${stepNum}-${itemIndex}" data-step="${stepNum}" data-index="${itemIndex}" 
                         title="${escapeHtml(isSearchResult ? item.citation : JSON.stringify(item))}" onclick="toggleObsItem(${stepNum}, ${itemIndex})">
                        <span class="obs-item-index">${itemIndex}</span>
                        <div class="obs-item-content" style="cursor: pointer;">
                            <div class="obs-item-text obs-item-summary" id="obs-summary-${stepNum}-${itemIndex}" style="white-space: nowrap; overflow: hidden; text-overflow: ellipsis;">${escapeHtml(itemPreview)}</div>
                            <div class="obs-item-text obs-item-full" id="obs-full-${stepNum}-${itemIndex}" style="display: none; white-space: pre-wrap; word-wrap: break-word; line-height: 1.4;">${escapeHtml(isSearchResult ? item.citation : JSON.stringify(item))}</div>
                        </div>
                    </div>
                `;
            });

            html += '</div>';
        } else {
            // Single observation
            let obsPreview = '';
            if (typeof step.observation === 'string' && step.observation !== "Completed.") {
                obsPreview = step.observation.length > 40 ?
                    step.observation.substring(0, 40) + '...' : step.observation;
            } else {
                obsPreview = '';
            }

            html += `
                <div class="step-content obs-content">
                    <span class="content-icon"></span>
                    <span class="content-text">${escapeHtml(obsPreview)}</span>
                </div>
            `;
        }
    }

    html += '</div>';
    return html;
}

// Generate HTML for observation section (handles both list and single observations)
function generateObservationHTML(observation, stepNum) {
    let html = `
        <div class="step-section observation-section">
            <div class="section-title">Observation</div>
            <div class="section-content">
    `;

    if (Array.isArray(observation)) {
        // Handle list of observations (e.g., search results)
        observation.forEach((item, index) => {
            const isSearchResult = item && typeof item === 'object' && item.source && item.citation;
            html += `
                <div class="observation-item" id="obs-${stepNum}-${index}" data-step="${stepNum}" data-index="${index}">
                    <div class="observation-index">[${index}]</div>
                    <div class="observation-content">
            `;

            if (isSearchResult) {
                const citationPreview = item.citation.length > 80 ?
                    item.citation.substring(0, 80) + '...' : item.citation;
                html += `
                    <div class="search-result">
                        <div class="source-url">${escapeHtml(item.source)}</div>
                        <div class="citation" title="${escapeHtml(item.citation)}">
                            ${escapeHtml(citationPreview)}
                        </div>
                    </div>
                `;
            } else {
                const itemPreview = JSON.stringify(item).length > 80 ?
                    JSON.stringify(item).substring(0, 80) + '...' : JSON.stringify(item);
                html += `<div class="simple-observation">${escapeHtml(itemPreview)}</div>`;
            }

            html += '</div></div>';
        });
    } else {
        // Handle single observation
        if (typeof observation === 'string') {
            const obsPreview = observation.length > 100 ?
                observation.substring(0, 100) + '...' : observation;
            html += `<div class="single-observation">${escapeHtml(obsPreview)}</div>`;
        } else {
            const obsPreview = JSON.stringify(observation).length > 100 ?
                JSON.stringify(observation).substring(0, 100) + '...' : JSON.stringify(observation);
            html += `<div class="single-observation">${escapeHtml(obsPreview)}</div>`;
        }
    }

    html += '</div></div>';
    return html;
}



// Generate HTML for answer node
function generateAnswerNodeHTML(stepCount) {
    const answer = flowTrajectoryData.answer || (flowTrajectoryData.trajectory && flowTrajectoryData.trajectory.answer);
    return generateExpandableContentNodeHTML({
        nodeType: 'answer',
        nodeId: 'answer-node',
        title: 'Final Answer',
        content: answer || 'No answer available',
        yPosition: stepCount * 400,
        toggleFunction: 'toggleAnswer'
    });
}

// Generate HTML for reasoning node
function generateReasoningNodeHTML(stepCount) {
    const reasoning = flowTrajectoryData.reasoning || (flowTrajectoryData.trajectory && flowTrajectoryData.trajectory.reasoning);
    return generateExpandableContentNodeHTML({
        nodeType: 'reasoning',
        nodeId: 'reasoning-node',
        title: 'Final Reasoning',
        content: reasoning || 'No reasoning available',
        yPosition: stepCount * 400,
        toggleFunction: 'toggleReasoning'
    });
}

// Generic function to generate expandable content nodes (reasoning, answer, etc.)
function generateExpandableContentNodeHTML(config) {
    const { nodeType, nodeId, title, content, yPosition, toggleFunction } = config;

    if (!content || content === 'No reasoning available' || content === 'No answer available') {
        return `
            <div class="${nodeType}-node" id="${nodeId}" style="top: ${yPosition}px;">
                <div class="step-header">${title}</div>
                <div class="step-content ${nodeType}-content">
                    <span class="content-text">${content}</span>
                </div>
            </div>
        `;
    }

    // Create summary (approximately 150 characters)
    const contentLines = content.trim().split('\n').filter(line => line.trim());
    let contentSummary = contentLines.slice(0, 3).join('\n').trim();
    if (contentSummary.length > 150) {
        contentSummary = contentSummary.substring(0, 150) + '...';
    } else if (contentLines.length > 3) {
        contentSummary += '...';
    }

    const needsExpansion = content.trim().length > 150 || contentLines.length > 3;

    return `
        <div class="${nodeType}-node" id="${nodeId}" style="top: ${yPosition}px;">
            <div class="step-header">${title}</div>
            <div class="step-content ${nodeType}-content">
                <div class="${nodeType}-container">
                    <div class="${nodeType}-summary" id="${nodeType}-summary">
                        ${formatText(contentSummary)}
                    </div>
                    ${needsExpansion ? `
                        <div class="${nodeType}-full" id="${nodeType}-full" style="display: none;">
                            ${formatText(content)}
                        </div>
                        <button class="thought-toggle-btn" onclick="${toggleFunction}()" id="${nodeType}-toggle">
                            Show More
                        </button>
                    ` : ''}
                </div>
            </div>
        </div>
    `;
}

// Helper function to escape HTML
function escapeHtml(text) {
    if (typeof text !== 'string') return '';
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
}

// Hide JSON content if visible
function hideJsonContent() {
    const jsonPanel = document.querySelector('#leftPanel .json-viewer');
    if (jsonPanel) {
        jsonPanel.style.display = 'none';
    }
}

// Draw arrows between connected steps
function drawFlowArrows() {
    const svg = document.getElementById('flowArrows');
    if (!svg) return;

    // Clear existing arrows
    svg.innerHTML = '';

    const trajectory = flowTrajectoryData.trajectory || flowTrajectoryData;
    const stepDependencies = flowMappingData.step_dependencies || {};
    const steps = extractStepsFromTrajectory(trajectory);

    // Set SVG dimensions - calculate total height needed based on actual content
    const container = svg.parentElement;
    const flowDiagram = container.querySelector('.flow-diagram');

    // Calculate actual content height by finding the bottom-most element
    let actualHeight = 100; // Base padding
    const allNodes = container.querySelectorAll('.step-node, .reasoning-node, .answer-node');
    allNodes.forEach(node => {
        const nodeTop = parseInt(node.style.top) || node.offsetTop;
        const nodeHeight = node.offsetHeight;
        const nodeBottom = nodeTop + nodeHeight;
        actualHeight = Math.max(actualHeight, nodeBottom + 50); // Add some bottom padding
    });

    // Ensure minimum height for traditional calculation as fallback
    const totalSteps = steps.length;
    const reasoning = flowTrajectoryData.reasoning || (flowTrajectoryData.trajectory && flowTrajectoryData.trajectory.reasoning);
    const answer = flowTrajectoryData.answer || (flowTrajectoryData.trajectory && flowTrajectoryData.trajectory.answer);
    const extraNodes = (reasoning ? 1 : 0) + (answer ? 1 : 0);
    const minHeight = Math.max(actualHeight, (totalSteps + extraNodes) * 200 + 200);

    svg.setAttribute('width', container.offsetWidth);
    svg.setAttribute('height', minHeight);

    // Also set the flow diagram container height
    if (flowDiagram) {
        flowDiagram.style.minHeight = minHeight + 'px';
        flowDiagram.style.height = minHeight + 'px';
    }

    // Ensure SVG is properly positioned and sized
    const containerRect = container.getBoundingClientRect();
    svg.setAttribute('width', containerRect.width);

    // Force layout recalculation before drawing arrows
    container.offsetHeight; // Force reflow

    // Draw arrows for each dependency
    steps.forEach(step => {
        const stepNum = step.number;
        const dependency = stepDependencies[stepNum.toString()];
        if (dependency && dependency.depends_on_step !== -1) {
            drawArrow(svg, dependency.depends_on_step, stepNum, dependency.result_index);
        }
    });

    // Draw arrow to reasoning if it exists
    if (reasoning && steps.length > 0) {
        const lastStep = Math.max(...steps.map(s => s.number));
        drawArrowToReasoning(svg, lastStep);
    }

    // Draw arrow from reasoning to answer if both exist
    if (reasoning && answer) {
        drawArrowFromReasoningToAnswer(svg);
    }
}

// Draw a single arrow between two steps
function drawArrow(svg, fromStep, toStep, resultIndex) {

    let fromElement, toElement;

    // Handle connection from specific observation item if resultIndex is specified

    fromElement = document.getElementById(`step-${fromStep}`);

    toElement = document.getElementById(`step-${toStep}`);

    if (!fromElement || !toElement) return;

    // Force layout recalculation to ensure accurate positions
    fromElement.offsetHeight;
    toElement.offsetHeight;

    // Force another reflow to make sure positions are updated
    const container = svg.parentElement.querySelector('.flow-diagram');
    if (container) {
        container.offsetHeight; // Force container reflow
    }

    // Get the flow diagram container (parent of SVG)
    const flowContainer = svg.parentElement;
    const flowContainerRect = flowContainer.getBoundingClientRect();

    // Calculate positions relative to the flow container using offsetTop/offsetLeft
    let fromX, fromY;

    if (resultIndex !== null && resultIndex !== undefined && fromElement.id.includes('obs-')) {
        // Arrow from specific observation item - from right side
        fromX = fromElement.offsetLeft + fromElement.offsetWidth;
        if (firstLoad && firstNodeX === -1) {
            firstNodeX = fromX;
            firstLoad = false;
        }
        fromY = fromElement.offsetTop + fromElement.offsetHeight / 2;
    } else {
        // Arrow from step - from bottom center
        // Force reflow before reading positions

        fromX = fromElement.offsetLeft + fromElement.offsetWidth / 2;
        if (firstLoad && firstNodeX === -1) {
            firstNodeX = fromX - 25;
            firstLoad = false;
        }
        fromX = firstNodeX; // Keep all arrows aligned to first node's X position
        // Use parsed style.top instead of offsetTop for accurate positioning
        const fromTopValue = parseInt(fromElement.style.top) || fromElement.offsetTop;
        fromY = fromTopValue + fromElement.offsetHeight;
    }

    // Arrow to step - to top center
    // Force one more reflow before reading positions
    toElement.getBoundingClientRect(); // This forces a style recalculation

    const toX = firstNodeX;
    // Use parsed style.top instead of offsetTop for accurate positioning
    const toTopValue = parseInt(toElement.style.top) || toElement.offsetTop;
    const toY = toTopValue;

    // Create curved path for vertical flow
    let path;
    let midX; // Declare midX outside the condition

    // Straight vertical path from step bottom to next step top
    midX = (fromX + toX) / 2; // Set midX for vertical paths too
    const midY = (fromY + toY) / 2;

    // Ensure minimum path length for arrow visibility
    const pathLength = Math.abs(toY - fromY);
    if (pathLength < 20) {
        // For very short paths, make them straight lines
        path = `M ${fromX} ${fromY} L ${toX} ${toY}`;
    } else {
        path = `M ${fromX} ${fromY} Q ${fromX} ${midY} ${toX} ${midY} Q ${toX} ${midY} ${toX} ${toY}`;
    }

    // Create path element without marker
    const pathElement = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    pathElement.setAttribute('d', path);
    pathElement.setAttribute('stroke', '#007bff'); // All arrows blue
    pathElement.setAttribute('stroke-width', '2');
    pathElement.setAttribute('fill', 'none');
    pathElement.setAttribute('opacity', '1');

    // Create triangle arrowhead at the end point
    const triangle = document.createElementNS('http://www.w3.org/2000/svg', 'polygon');
    const arrowSize = 8;
    // Calculate triangle points for downward arrow (since most arrows go down)
    const trianglePoints = `${toX},${toY} ${toX - arrowSize},${toY - arrowSize} ${toX + arrowSize},${toY - arrowSize}`;
    triangle.setAttribute('points', trianglePoints);
    triangle.setAttribute('fill', '#007bff');
    triangle.setAttribute('stroke', '#007bff');
    triangle.setAttribute('stroke-width', '1');

    svg.appendChild(pathElement);
    svg.appendChild(triangle);
}

// Draw arrow to reasoning node
function drawArrowToReasoning(svg, fromStep) {
    const fromElement = document.getElementById(`step-${fromStep}`);
    const toElement = document.getElementById('reasoning-node');

    if (!fromElement || !toElement) return;

    // Use offset positions relative to the flow container
    const fromX = firstNodeX;
    const fromTopValue = parseInt(fromElement.style.top) || fromElement.offsetTop;
    const fromY = fromTopValue + fromElement.offsetHeight;
    const toX = fromX;
    const toTopValue = parseInt(toElement.style.top) || toElement.offsetTop;
    const toY = toTopValue;

    const midY = (fromY + toY) / 2;
    const path = `M ${fromX} ${fromY} Q ${fromX} ${midY} ${toX} ${midY} Q ${toX} ${midY} ${toX} ${toY}`;

    const pathElement = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    pathElement.setAttribute('d', path);
    pathElement.setAttribute('stroke', '#007bff'); // Blue arrow for reasoning too
    pathElement.setAttribute('stroke-width', '2');
    pathElement.setAttribute('fill', 'none');
    pathElement.setAttribute('opacity', '1');

    // Create triangle arrowhead at the end point
    const triangle = document.createElementNS('http://www.w3.org/2000/svg', 'polygon');
    const arrowSize = 8;
    const trianglePoints = `${toX},${toY} ${toX - arrowSize},${toY - arrowSize} ${toX + arrowSize},${toY - arrowSize}`;
    triangle.setAttribute('points', trianglePoints);
    triangle.setAttribute('fill', '#007bff');
    triangle.setAttribute('stroke', '#007bff');
    triangle.setAttribute('stroke-width', '1');

    svg.appendChild(pathElement);
    svg.appendChild(triangle);
}

// Draw arrow from reasoning to answer node
function drawArrowFromReasoningToAnswer(svg) {
    const fromElement = document.getElementById('reasoning-node');
    const toElement = document.getElementById('answer-node');

    if (!fromElement || !toElement) return;

    // Use offset positions relative to the flow container
    const fromX = firstNodeX;
    const fromTopValue = parseInt(fromElement.style.top) || fromElement.offsetTop;
    const fromY = fromTopValue + fromElement.offsetHeight;
    const toX = fromX
    const toTopValue = parseInt(toElement.style.top) || toElement.offsetTop;
    const toY = toTopValue;

    const midY = (fromY + toY) / 2;
    const path = `M ${fromX} ${fromY} Q ${fromX} ${midY} ${toX} ${midY} Q ${toX} ${midY} ${toX} ${toY}`;

    const pathElement = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    pathElement.setAttribute('d', path);
    pathElement.setAttribute('stroke', '#007bff'); // Blue arrow
    pathElement.setAttribute('stroke-width', '2');
    pathElement.setAttribute('fill', 'none');
    pathElement.setAttribute('opacity', '1');

    // Create triangle arrowhead at the end point
    const triangle = document.createElementNS('http://www.w3.org/2000/svg', 'polygon');
    const arrowSize = 8;
    const trianglePoints = `${toX},${toY} ${toX - arrowSize},${toY - arrowSize} ${toX + arrowSize},${toY - arrowSize}`;
    triangle.setAttribute('points', trianglePoints);
    triangle.setAttribute('fill', '#007bff');
    triangle.setAttribute('stroke', '#007bff');
    triangle.setAttribute('stroke-width', '1');

    svg.appendChild(pathElement);
    svg.appendChild(triangle);
}

// Check if trajectory flow is currently visible
function isTrajectoryFlowCurrentlyVisible() {
    return isTrajectoryFlowVisible;
}

// Draw connections between evidence items and trajectory steps
function drawEvidenceConnections(explanationIndex) {
    // Try to get mapping data from multiple sources
    const mappingData = flowMappingData || window.mappingData;
    if (!mappingData || !mappingData.evidence_to_trajectory) {
        return;
    }

    const evidenceToTrajectory = mappingData.evidence_to_trajectory;
    const evidenceSection = document.getElementById(`evidenceSection${explanationIndex}`);

    if (!evidenceSection) {
        console.error('Evidence section not found');
        return;
    }

    // Clear existing evidence connection lines
    const existingLines = document.querySelectorAll('.evidence-connection-line');
    existingLines.forEach(line => line.remove());

    // Find all evidence items in this explanation
    const evidenceItems = evidenceSection.querySelectorAll('.evidence-item[data-evidence-id]');

    evidenceItems.forEach(evidenceItem => {
        const evidenceId = evidenceItem.getAttribute('data-evidence-id');
        const explanationIdx = evidenceItem.getAttribute('data-explanation-index');
        const uniqueConnectionId = `${explanationIdx}-${evidenceId}`;

        // Get evidence type from data attribute
        const evidenceType = evidenceItem.getAttribute('data-evidence-type');

        // Check if evidence type is "reasoning"
        if (evidenceType === 'reasoning') {
            // Connect to reasoning node instead of trajectory steps
            const reasoningNode = document.getElementById('reasoning-node');
            if (reasoningNode) {
                drawEvidenceConnectionLine(evidenceItem, reasoningNode, uniqueConnectionId);
            }
        } else {
            // Regular evidence: connect to trajectory steps
            const trajectoryMapping = evidenceToTrajectory[evidenceId];

            if (trajectoryMapping && Array.isArray(trajectoryMapping) && trajectoryMapping.length > 0) {
                // Only connect to the first step, ignore additional steps
                const stepNumber = trajectoryMapping[0];

                // Connect to step node only, not to specific observation items
                const targetElement = document.querySelector(`#step-${stepNumber}`);

                if (targetElement) {
                    drawEvidenceConnectionLine(evidenceItem, targetElement, uniqueConnectionId);

                    // Add visual distinction to observation items if mapping has specific indices
                    if (trajectoryMapping.length > 1) {
                        for (let i = 1; i < trajectoryMapping.length; i++) {
                            const obsIndex = trajectoryMapping[i];
                            const obsElement = document.querySelector(`#obs-${stepNumber}-${obsIndex}`);
                            if (obsElement) {
                                obsElement.classList.add('highlighted-observation');
                            }
                        }
                    }
                }
            }
        }
    });
}

// Draw a connection line between evidence item and trajectory element
function drawEvidenceConnectionLine(evidenceElement, trajectoryElement, evidenceId) {
    // Create SVG for the connection line
    const connectionSvg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    connectionSvg.setAttribute('class', 'evidence-connection-line');
    connectionSvg.setAttribute('id', `connection-${evidenceId}`);
    connectionSvg.style.cssText = `
        position: fixed;
        top: 100px;
        left: 0;
        width: 100vw;
        height: calc(100vh - 100px);
        pointer-events: none;
        z-index: 900;
        overflow: visible;
    `;

    // Create defs section for the marker if not exists
    let defs = connectionSvg.querySelector('defs');
    if (!defs) {
        defs = document.createElementNS('http://www.w3.org/2000/svg', 'defs');

        // Blue thin arrowhead for evidence connections
        const arrowhead = document.createElementNS('http://www.w3.org/2000/svg', 'marker');
        arrowhead.setAttribute('id', `evidence-arrowhead-${evidenceId}`);
        arrowhead.setAttribute('markerWidth', '8');
        arrowhead.setAttribute('markerHeight', '6');
        arrowhead.setAttribute('refX', '7');
        arrowhead.setAttribute('refY', '3');
        arrowhead.setAttribute('orient', 'auto');

        const arrowPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        arrowPath.setAttribute('d', 'M0,0 L0,6 L8,3 z');
        arrowPath.setAttribute('fill', '#007bff');

        arrowhead.appendChild(arrowPath);
        defs.appendChild(arrowhead);
        connectionSvg.appendChild(defs);
    }

    // Add path element for right-angle turns instead of straight line
    const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    path.setAttribute('stroke', '#007bff');
    path.setAttribute('stroke-width', '2');
    path.setAttribute('fill', 'none');
    path.setAttribute('opacity', '1');

    connectionSvg.appendChild(path);

    // Store triangle reference for later creation
    connectionSvg.evidenceTriangle = document.createElementNS('http://www.w3.org/2000/svg', 'polygon');
    connectionSvg.evidenceTriangle.setAttribute('fill', '#007bff');
    connectionSvg.evidenceTriangle.setAttribute('stroke', '#007bff');
    connectionSvg.evidenceTriangle.setAttribute('stroke-width', '1');
    connectionSvg.appendChild(connectionSvg.evidenceTriangle);
    document.body.appendChild(connectionSvg);

    // Update path position with right-angle turns
    updateEvidenceConnectionPath(evidenceElement, trajectoryElement, path);

    // Store reference for updates
    if (!window.evidenceConnections) {
        window.evidenceConnections = new Map();
    }
    window.evidenceConnections.set(evidenceId, {
        svg: connectionSvg,
        line: path,
        evidenceElement: evidenceElement,
        trajectoryElement: trajectoryElement
    });
}

// Update the position of an evidence connection path with right-angle turns
function updateEvidenceConnectionPath(evidenceElement, trajectoryElement, pathElement) {
    try {
        const evidenceRect = evidenceElement.getBoundingClientRect();
        const trajectoryRect = trajectoryElement.getBoundingClientRect();

        // Calculate connection points - trajectory to evidence (from trajectory right to evidence left side)
        // Adjust coordinates relative to SVG position (which starts at top: 100px)
        const startX = trajectoryRect.right;
        const startY = trajectoryRect.top + trajectoryRect.height / 2 - 100; // Subtract header height
        const endX = evidenceRect.left; // Left edge of evidence item
        const endY = evidenceRect.top - 87.2; // Subtract header height

        // Create right-angle path: horizontal from trajectory, then vertical, then horizontal to evidence
        const midX = startX + (endX - startX) / 2;

        const pathData = `M ${startX} ${startY} L ${midX} ${startY} L ${midX} ${endY} L ${endX} ${endY}`;
        pathElement.setAttribute('d', pathData);

        // Make sure the path is visible with blue color
        pathElement.setAttribute('stroke', '#007bff');
        pathElement.setAttribute('stroke-width', '2');
        pathElement.setAttribute('opacity', '1');

        // Update triangle position at the end point (pointing right)
        const svg = pathElement.parentElement;
        if (svg && svg.evidenceTriangle) {
            const arrowSize = 6;
            const trianglePoints = `${endX},${endY} ${endX - arrowSize},${endY - arrowSize} ${endX - arrowSize},${endY + arrowSize}`;
            svg.evidenceTriangle.setAttribute('points', trianglePoints);
        }
    } catch (error) {
        console.error('Error updating evidence connection path:', error);
    }
}

// Legacy function for backward compatibility
function updateEvidenceConnectionLine(evidenceElement, trajectoryElement, lineElement) {
    // If it's actually a path element, use the new function
    if (lineElement.tagName === 'path') {
        updateEvidenceConnectionPath(evidenceElement, trajectoryElement, lineElement);
        return;
    }

    // Original line logic (kept for any remaining line elements)
    try {
        const evidenceRect = evidenceElement.getBoundingClientRect();
        const trajectoryRect = trajectoryElement.getBoundingClientRect();

        const startX = trajectoryRect.right;
        const startY = trajectoryRect.top + trajectoryRect.height / 2 - 100; // Subtract header height
        const endX = evidenceRect.left;
        const endY = evidenceRect.top + evidenceRect.height / 2 - 100; // Subtract header height

        lineElement.setAttribute('x1', startX);
        lineElement.setAttribute('y1', startY);
        lineElement.setAttribute('x2', endX);
        lineElement.setAttribute('y2', endY);
    } catch (error) {
        console.error('Error updating evidence connection line:', error);
    }
}

// Set up scroll listeners for dynamic line updates
function setupScrollListeners() {
    // Remove existing listeners to avoid duplicates
    if (window.evidenceScrollHandler) {
        window.removeEventListener('scroll', window.evidenceScrollHandler, true);
        const leftPanel = document.getElementById('leftPanel');
        if (leftPanel) {
            leftPanel.removeEventListener('scroll', window.evidenceScrollHandler);
        }
    }

    // Create new scroll handler
    window.evidenceScrollHandler = function () {
        if (window.evidenceConnections) {
            window.evidenceConnections.forEach((connection, evidenceId) => {
                updateEvidenceConnectionLine(
                    connection.evidenceElement,
                    connection.trajectoryElement,
                    connection.line
                );
            });
        }
    };

    // Add scroll listeners
    window.addEventListener('scroll', window.evidenceScrollHandler, true);
    const leftPanel = document.getElementById('leftPanel');
    if (leftPanel) {
        leftPanel.addEventListener('scroll', window.evidenceScrollHandler);
    }
}

// Clear all evidence connection lines
function clearEvidenceConnections() {
    if (window.evidenceConnections) {
        window.evidenceConnections.forEach((connection) => {
            connection.svg.remove();
        });
        window.evidenceConnections.clear();
    }
}

// Auto-scroll to relevant positions when trajectory flow is shown
function scrollToRelevantPositions(explanationIndex) {
    // Try to get mapping data from multiple sources
    const mappingData = flowMappingData || window.mappingData;
    if (!mappingData || !mappingData.evidence_to_trajectory) {
        return;
    }

    const evidenceToTrajectory = mappingData.evidence_to_trajectory;
    const evidenceSection = document.getElementById(`evidenceSection${explanationIndex}`);

    if (!evidenceSection) {
        return;
    }

    // Find the first evidence item in this explanation
    const firstEvidenceItem = evidenceSection.querySelector('.evidence-item[data-evidence-id]');
    if (!firstEvidenceItem) {
        return;
    }

    const evidenceId = firstEvidenceItem.getAttribute('data-evidence-id');
    const trajectoryMapping = evidenceToTrajectory[evidenceId];


    if (trajectoryMapping && Array.isArray(trajectoryMapping)) {
        const stepNumber = trajectoryMapping[0];
        const observationIndex = trajectoryMapping.length > 1 ? trajectoryMapping[1] : null;

        // Find the target step in trajectory flow
        let targetStepElement;
        if (observationIndex !== null) {
            targetStepElement = document.querySelector(`#obs-${stepNumber}-${observationIndex}`);
        }
        if (!targetStepElement) {
            targetStepElement = document.querySelector(`#step-${stepNumber}`);
        }

        if (targetStepElement) {
            // Scroll trajectory flow (left panel) to show the target step
            const leftPanel = document.getElementById('leftPanel');
            if (leftPanel) {
                // Use the correct scrollable container for trajectory flow
                const leftPanelContent = leftPanel.querySelector('.trajectory-flow-content') ||
                    leftPanel.querySelector('.left-panel-content') ||
                    leftPanel;
                const leftPanelRect = leftPanelContent.getBoundingClientRect();
                const targetRect = targetStepElement.getBoundingClientRect();

                // Calculate scroll position relative to the scrollable container
                const relativeTop = targetRect.top - leftPanelRect.top;
                const centerOffset = leftPanelRect.height / 2 - targetRect.height / 2;
                const scrollTop = leftPanelContent.scrollTop + relativeTop - centerOffset;

                // Force scroll even if smooth behavior fails
                leftPanelContent.scrollTop = Math.max(0, scrollTop);

                // Also try smooth scroll as fallback
                setTimeout(() => {
                    leftPanelContent.scrollTo({
                        top: Math.max(0, scrollTop),
                        behavior: 'smooth'
                    });
                }, 50);
            }
        }

        // Scroll main content area to show the evidence item
        setTimeout(() => {
            // Try multiple approaches for main content scrolling
            const mainContent = document.querySelector('.main-content') || document.querySelector('#contentArea');
            const evidenceRect = firstEvidenceItem.getBoundingClientRect();

            if (mainContent) {
                const mainContentRect = mainContent.getBoundingClientRect();
                const relativeTop = evidenceRect.top - mainContentRect.top;
                const centerOffset = window.innerHeight / 2 - evidenceRect.height / 2;
                const scrollTop = mainContent.scrollTop + relativeTop - centerOffset;

                // Force scroll
                mainContent.scrollTop = Math.max(0, scrollTop);

                // Try smooth scroll as well
                setTimeout(() => {
                    mainContent.scrollTo({
                        top: Math.max(0, scrollTop),
                        behavior: 'smooth'
                    });
                }, 50);
            } else {
                // Fallback to window scroll
                const viewportTop = evidenceRect.top;
                const centerOffset = window.innerHeight / 2 - evidenceRect.height / 2;
                const scrollTop = window.scrollY + viewportTop - centerOffset;

                window.scrollTo({
                    top: Math.max(0, scrollTop),
                    behavior: 'smooth'
                });
            }
        }, 300); // Increase delay to let left panel scroll first
    }
}



// Export functions for use in other scripts
window.trajectoryFlow = {
    init: initTrajectoryFlow,
    loadData: loadTrajectoryFlowData,
    show: showTrajectoryFlow,
    hide: hideTrajectoryFlow,
    isVisible: isTrajectoryFlowCurrentlyVisible,
    drawEvidenceConnections: drawEvidenceConnections,
    clearEvidenceConnections: clearEvidenceConnections,
    scrollToRelevantPositions: scrollToRelevantPositions,
    // Debug functions
    testDirectConnection: function () {
        // Create a simple test connection line
        const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        svg.style.cssText = `
            position: fixed;
            top: 100px;
            left: 0;
            width: 100vw;
            height: calc(100vh - 100px);
            pointer-events: none;
            z-index: 2000;
            background: rgba(255,0,0,0.1);
        `;

        const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        line.setAttribute('x1', '100');
        line.setAttribute('y1', '100');
        line.setAttribute('x2', '500');
        line.setAttribute('y2', '300');
        line.setAttribute('stroke', '#ff6b6b');
        line.setAttribute('stroke-width', '5');
        line.setAttribute('opacity', '1');

        svg.appendChild(line);
        document.body.appendChild(svg);

        // Remove after 5 seconds
        setTimeout(() => {
            svg.remove();
        }, 5000);
    }
};

// Toggle thought display and recalculate node positions
function toggleThought(stepNum) {
    const summaryEl = document.getElementById(`thought-summary-${stepNum}`);
    const fullEl = document.getElementById(`thought-full-${stepNum}`);
    const toggleBtn = document.getElementById(`thought-toggle-${stepNum}`);

    if (!summaryEl || !fullEl || !toggleBtn) return;

    const isExpanded = fullEl.style.display !== 'none';

    if (isExpanded) {
        // Collapse
        summaryEl.style.display = 'block';
        fullEl.style.display = 'none';
        toggleBtn.textContent = 'Show More';
    } else {
        // Expand
        summaryEl.style.display = 'none';
        fullEl.style.display = 'block';
        toggleBtn.textContent = 'Show Less';
    }

    // Recalculate positions after a short delay
    setTimeout(() => {
        recalculateNodePositions();
        // Clear and redraw arrows after repositioning
        const svg = document.getElementById('flowArrows');
        if (svg) {
            svg.innerHTML = '';
        }
        drawFlowArrows();
        // Redraw evidence connections if they exist
        if (document.querySelector('.evidence-connection-line')) {
            const activeExplanationIndex = getActiveExplanationIndex();
            if (activeExplanationIndex !== null) {
                clearEvidenceConnections();
                drawEvidenceConnections(activeExplanationIndex);
            }
        }
    }, 100);
}

// Recalculate node positions based on content height
function recalculateNodePositions() {
    const stepNodes = document.querySelectorAll('.step-node');
    let currentY = 50; // Starting position

    stepNodes.forEach((node, index) => {
        node.style.top = currentY + 'px';

        // Force reflow to get accurate height
        node.offsetHeight;

        // Get actual height of the node after reflow
        const nodeHeight = node.offsetHeight;
        currentY += nodeHeight + 30; // 30px spacing between nodes
    });

    // Update reasoning and answer node positions
    const reasoningNode = document.querySelector('.reasoning-node');
    if (reasoningNode) {
        reasoningNode.style.top = currentY + 'px';
        reasoningNode.offsetHeight; // Force reflow
        currentY += reasoningNode.offsetHeight + 30;
    }

    const answerNode = document.querySelector('.answer-node');
    if (answerNode) {
        answerNode.style.top = currentY + 'px';
        answerNode.offsetHeight; // Force reflow
        currentY += answerNode.offsetHeight + 30;
    }

    // Update the container height to accommodate all nodes
    const flowDiagram = document.querySelector('.flow-diagram');
    if (flowDiagram) {
        const totalHeight = currentY + 100; // Add some padding
        flowDiagram.style.minHeight = totalHeight + 'px';
        flowDiagram.style.height = totalHeight + 'px';

        // Update SVG height as well
        const svg = document.getElementById('flowArrows');
        if (svg) {
            svg.setAttribute('height', totalHeight);
        }

        // Force a complete reflow after all position updates
        flowDiagram.offsetHeight;

        // Force each step node to reflow again to ensure positions are committed
        stepNodes.forEach(node => {
            node.offsetTop; // Force individual reflow
        });
    }
}

// Helper function to get the currently active explanation index
function getActiveExplanationIndex() {
    const expandedEvidenceSection = document.querySelector('.evidence-section.expanded');
    if (expandedEvidenceSection) {
        const match = expandedEvidenceSection.id.match(/evidenceSection(\d+)/);
        return match ? parseInt(match[1]) : null;
    }
    return null;
}

// Toggle observation item display between summary and full text
function toggleObsItem(stepNum, itemIndex) {
    const summaryEl = document.getElementById(`obs-summary-${stepNum}-${itemIndex}`);
    const fullEl = document.getElementById(`obs-full-${stepNum}-${itemIndex}`);

    if (!summaryEl || !fullEl) return;

    const isExpanded = fullEl.style.display !== 'none';

    if (isExpanded) {
        // Collapse - show summary, hide full
        summaryEl.style.display = 'block';
        fullEl.style.display = 'none';
    } else {
        // Expand - hide summary, show full
        summaryEl.style.display = 'none';
        fullEl.style.display = 'block';
    }

    // Recalculate positions after content change
    setTimeout(() => {
        recalculateNodePositions();
        // Redraw arrows
        const svg = document.getElementById('flowArrows');
        if (svg) {
            svg.innerHTML = '';
        }
        drawFlowArrows();
        // Redraw evidence connections if they exist
        if (document.querySelector('.evidence-connection-line')) {
            const activeExplanationIndex = getActiveExplanationIndex();
            if (activeExplanationIndex !== null) {
                clearEvidenceConnections();
                drawEvidenceConnections(activeExplanationIndex);
            }
        }
    }, 50);
}

// Toggle reasoning display between summary and full text
function toggleReasoning() {
    const summaryEl = document.getElementById('reasoning-summary');
    const fullEl = document.getElementById('reasoning-full');
    const toggleBtn = document.getElementById('reasoning-toggle');

    if (!summaryEl || !fullEl || !toggleBtn) return;

    const isExpanded = fullEl.style.display !== 'none';

    if (isExpanded) {
        // Collapse
        summaryEl.style.display = 'block';
        fullEl.style.display = 'none';
        toggleBtn.textContent = 'Show More';
    } else {
        // Expand
        summaryEl.style.display = 'none';
        fullEl.style.display = 'block';
        toggleBtn.textContent = 'Show Less';
    }

    // Recalculate positions after a short delay
    setTimeout(() => {
        recalculateNodePositions();
        // Clear and redraw arrows after repositioning
        const svg = document.getElementById('flowArrows');
        if (svg) {
            svg.innerHTML = '';
        }
        drawFlowArrows();
        // Redraw evidence connections if they exist
        if (document.querySelector('.evidence-connection-line')) {
            const activeExplanationIndex = getActiveExplanationIndex();
            if (activeExplanationIndex !== null) {
                clearEvidenceConnections();
                drawEvidenceConnections(activeExplanationIndex);
            }
        }
    }, 100);
}

// Toggle answer display between summary and full text
function toggleAnswer() {
    const summaryEl = document.getElementById('answer-summary');
    const fullEl = document.getElementById('answer-full');
    const toggleBtn = document.getElementById('answer-toggle');

    if (!summaryEl || !fullEl || !toggleBtn) return;

    const isExpanded = fullEl.style.display !== 'none';

    if (isExpanded) {
        // Collapse
        summaryEl.style.display = 'block';
        fullEl.style.display = 'none';
        toggleBtn.textContent = 'Show More';
    } else {
        // Expand
        summaryEl.style.display = 'none';
        fullEl.style.display = 'block';
        toggleBtn.textContent = 'Show Less';
    }

    // Recalculate positions after a short delay
    setTimeout(() => {
        recalculateNodePositions();
        // Clear and redraw arrows after repositioning
        const svg = document.getElementById('flowArrows');
        if (svg) {
            svg.innerHTML = '';
        }
        drawFlowArrows();
        // Redraw evidence connections if they exist
        if (document.querySelector('.evidence-connection-line')) {
            const activeExplanationIndex = getActiveExplanationIndex();
            if (activeExplanationIndex !== null) {
                clearEvidenceConnections();
                drawEvidenceConnections(activeExplanationIndex);
            }
        }
    }, 100);
}