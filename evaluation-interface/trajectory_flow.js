// Trajectory Flow Visualization
// This module handles the display of trajectory flow diagrams in the left panel

let flowTrajectoryData = null;
let flowMappingData = null;
let isTrajectoryFlowVisible = false;

// Initialize trajectory flow system
function initTrajectoryFlow() {
    // This will be called when the page loads to set up event listeners
    console.log('Trajectory flow system initialized');
}

// Load trajectory and mapping data
function loadTrajectoryFlowData(trajectory, mapping) {
    flowTrajectoryData = trajectory;
    flowMappingData = mapping;
    console.log('Trajectory flow data loaded:', { trajectory: !!trajectory, mapping: !!mapping });
}

// Show trajectory flow in left panel
function showTrajectoryFlow(explanationIndex = null) {
    if (!flowTrajectoryData || !flowMappingData) {
        console.error('Trajectory or mapping data not available');
        return;
    }

    const leftPanel = document.getElementById('leftPanel');
    if (!leftPanel) {
        console.error('Left panel not found');
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
                <h3>🔍 Trajectory Flow</h3>
                <button class="close-btn" onclick="hideTrajectoryFlow()">✕</button>
            </div>
            <div class="trajectory-flow-content">
                ${flowContent}
            </div>
        </div>
    `;

    // Add CSS for trajectory flow if not already added
    addTrajectoryFlowCSS();

    // Draw arrows after a short delay to ensure elements are rendered
    setTimeout(() => {
        drawFlowArrows();
        // Draw evidence connections if explanationIndex is provided
        if (explanationIndex !== null) {
            // Wait a bit more for evidence section to be expanded
            setTimeout(() => {
                drawEvidenceConnections(explanationIndex);
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
        // You may want to restore previous content here
        // For now, just clear the content
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
    const yPosition = index * 400; // Further increased vertical spacing for observation items

    let html = `<div class="step-node" id="step-${stepNum}" data-step="${stepNum}" style="top: ${yPosition}px;">`;

    html += `<div class="step-header">Step ${stepNum}</div>`;

    // Compact content display
    if (step.thought) {
        const thoughtPreview = step.thought.length > 60 ?
            step.thought.substring(0, 60) + '...' : step.thought;
        html += `
            <div class="step-content thought-content">
                <span class="content-icon">💭</span>
                <span class="content-text" title="${escapeHtml(step.thought)}">
                    ${escapeHtml(thoughtPreview)}
                </span>
            </div>
        `;
    }

    if (step.toolName) {
        html += `
            <div class="step-content tool-content">
                <span class="content-icon">🔧</span>
                <span class="content-text">${escapeHtml(step.toolName)}</span>
            </div>
        `;
    }

    if (step.observation !== null) {
        if (Array.isArray(step.observation)) {
            // Show individual observation items as small boxes
            html += `
                <div class="step-content obs-content">
                    <span class="content-icon">👁️</span>
                    <span class="content-text">Observations:</span>
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
                         title="${escapeHtml(isSearchResult ? item.citation : JSON.stringify(item))}">
                        <span class="obs-item-index">[${itemIndex}]</span>
                        <span class="obs-item-text">${escapeHtml(itemPreview)}</span>
                    </div>
                `;
            });

            html += '</div>';
        } else {
            // Single observation
            let obsPreview = '';
            if (typeof step.observation === 'string') {
                obsPreview = step.observation.length > 40 ?
                    step.observation.substring(0, 40) + '...' : step.observation;
            } else {
                obsPreview = 'Data result';
            }

            html += `
                <div class="step-content obs-content">
                    <span class="content-icon">👁️</span>
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
            <div class="section-title">👁️ Observation</div>
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

// Generate HTML for reasoning node
function generateReasoningNodeHTML(stepCount) {
    const reasoning = flowTrajectoryData.reasoning || (flowTrajectoryData.trajectory && flowTrajectoryData.trajectory.reasoning);
    const reasoningPreview = reasoning && reasoning.length > 80 ?
        reasoning.substring(0, 80) + '...' : reasoning;
    const yPosition = stepCount * 400; // Updated to match step spacing

    return `
        <div class="reasoning-node" id="reasoning-node" style="top: ${yPosition}px;">
            <div class="step-header">📋 Final Reasoning</div>
            <div class="step-content reasoning-content">
                <span class="content-text" title="${escapeHtml(reasoning || '')}">
                    ${escapeHtml(reasoningPreview || 'No reasoning available')}
                </span>
            </div>
        </div>
    `;
}

// Generate HTML for answer node
function generateAnswerNodeHTML(stepCount) {
    const answer = flowTrajectoryData.answer || (flowTrajectoryData.trajectory && flowTrajectoryData.trajectory.answer);
    const answerPreview = answer && answer.length > 80 ?
        answer.substring(0, 80) + '...' : answer;
    const yPosition = stepCount * 400; // Updated to match step spacing

    return `
        <div class="answer-node" id="answer-node" style="top: ${yPosition}px;">
            <div class="step-header">📝 Final Answer</div>
            <div class="step-content answer-content">
                <span class="content-text" title="${escapeHtml(answer || '')}">
                    ${escapeHtml(answerPreview || 'No answer available')}
                </span>
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

// Add CSS for trajectory flow
function addTrajectoryFlowCSS() {
    if (document.getElementById('trajectory-flow-css')) {
        return; // CSS already added
    }

    const style = document.createElement('style');
    style.id = 'trajectory-flow-css';
    style.textContent = `
        .trajectory-flow-container {
            height: 100%;
            display: flex;
            flex-direction: column;
            background: #f8f9fa;
        }
        
        .trajectory-flow-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 1rem;
            background: #fff;
            border-bottom: 1px solid #e9ecef;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }
        
        .trajectory-flow-header h3 {
            margin: 0;
            color: #2c3e50;
            font-size: 1.1rem;
        }
        
        .close-btn {
            background: none;
            border: none;
            font-size: 1.2rem;
            cursor: pointer;
            padding: 0.25rem 0.5rem;
            color: #6c757d;
            border-radius: 4px;
        }
        
        .close-btn:hover {
            background: #f8f9fa;
            color: #495057;
        }
        
        .trajectory-flow-content {
            flex: 1;
            overflow: auto;
            padding: 1rem;
            position: relative;
        }
        
        .flow-diagram-container {
            position: relative;
            min-height: 100%;
            width: 100%;
        }
        
        .flow-arrows {
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            pointer-events: none;
            z-index: 5;
        }
        
        .flow-diagram {
            position: relative;
            z-index: 2;
            padding: 0 50px;
            min-height: 100%;
            display: flex;
            flex-direction: column;
            align-items: center;
        }
        
        .step-node, .reasoning-node, .answer-node {
            position: absolute;
            left: 50%;
            transform: translateX(-50%);
            width: 300px;
            background: #fff;
            border: 2px solid #dee2e6;
            border-radius: 8px;
            padding: 0.75rem;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            transition: all 0.2s ease;
            z-index: 3;
        }
        
        .step-node:hover, .reasoning-node:hover, .answer-node:hover {
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
            transform: translateX(-50%) translateY(-2px);
        }
        
        .step-node {
            border-left: 4px solid #007bff;
        }
        
        .reasoning-node {
            border-left: 4px solid #6f42c1;
        }
        
        .answer-node {
            border-left: 4px solid #28a745;
        }
        
        .step-header {
            font-weight: bold;
            color: #2c3e50;
            margin-bottom: 0.5rem;
            font-size: 0.9rem;
            text-align: center;
            padding-bottom: 0.25rem;
            border-bottom: 1px solid #e9ecef;
        }
        
        .step-content {
            display: flex;
            align-items: flex-start;
            gap: 0.5rem;
            margin-bottom: 0.5rem;
            font-size: 0.75rem;
            line-height: 1.3;
        }
        
        .step-content:last-child {
            margin-bottom: 0;
        }
        
        .content-icon {
            flex-shrink: 0;
            font-size: 0.8rem;
        }
        
        .content-text {
            flex: 1;
            color: #6c757d;
            word-wrap: break-word;
        }
        
        .thought-content .content-text {
            color: #007bff;
        }
        
        .tool-content .content-text {
            color: #28a745;
            font-family: 'Courier New', monospace;
            font-weight: bold;
        }
        
        .obs-content .content-text {
            color: #ffc107;
            font-style: italic;
        }
        
        .reasoning-content .content-text {
            color: #6f42c1;
        }
        
        .answer-content .content-text {
            color: #28a745;
        }
        
        .observation-items {
            margin-top: 0.5rem;
            display: flex;
            flex-wrap: wrap;
            gap: 0.25rem;
        }
        
        .obs-item-box {
            display: flex;
            align-items: center;
            gap: 0.25rem;
            background: #fff3cd;
            border: 1px solid #ffc107;
            border-radius: 4px;
            padding: 0.25rem 0.5rem;
            font-size: 0.7rem;
            cursor: pointer;
            transition: all 0.2s ease;
            max-width: 120px;
        }
        
        .obs-item-box:hover {
            background: #fff3cd;
            border-color: #ffcd39;
            transform: translateY(-1px);
            box-shadow: 0 2px 4px rgba(255, 193, 7, 0.3);
        }
        
        .obs-item-index {
            background: #ffc107;
            color: #212529;
            padding: 0.1rem 0.3rem;
            border-radius: 2px;
            font-weight: bold;
            font-size: 0.65rem;
            flex-shrink: 0;
        }
        
        .obs-item-text {
            color: #856404;
            font-weight: 500;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            flex: 1;
        }
        
        /* Arrow markers */
        .flow-arrows defs {
            display: none;
        }
        
        /* Evidence connection lines */
        .evidence-connection-line {
            z-index: 1000 !important;
            position: fixed !important;
            overflow: visible !important;
        }
        
        .evidence-connection-line line {
            stroke: #ff6b6b !important;
            stroke-width: 3 !important;
            stroke-dasharray: 8,4 !important;
            opacity: 0.9 !important;
            animation: dash 2s linear infinite !important;
        }
        
        @keyframes dash {
            to {
                stroke-dashoffset: -12;
            }
        }
        
        /* Responsive design */
        @media (max-width: 768px) {
            .trajectory-flow-content {
                padding: 0.5rem;
            }
            
            .step-node, .reasoning-node {
                width: 180px;
            }
            
            .flow-diagram {
                padding-left: 30px;
            }
        }
    `;

    // Add SVG arrow markers
    const defs = document.createElementNS('http://www.w3.org/2000/svg', 'defs');

    // Blue arrowhead for regular connections
    const arrowhead = document.createElementNS('http://www.w3.org/2000/svg', 'marker');
    arrowhead.setAttribute('id', 'arrowhead');
    arrowhead.setAttribute('markerWidth', '10');
    arrowhead.setAttribute('markerHeight', '7');
    arrowhead.setAttribute('refX', '9');
    arrowhead.setAttribute('refY', '3.5');
    arrowhead.setAttribute('orient', 'auto');
    const arrowPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    arrowPath.setAttribute('d', 'M0,0 L0,7 L10,3.5 z');
    arrowPath.setAttribute('fill', '#007bff');
    arrowhead.appendChild(arrowPath);
    defs.appendChild(arrowhead);

    // Purple arrowhead for reasoning connections
    const arrowheadPurple = document.createElementNS('http://www.w3.org/2000/svg', 'marker');
    arrowheadPurple.setAttribute('id', 'arrowhead-purple');
    arrowheadPurple.setAttribute('markerWidth', '10');
    arrowheadPurple.setAttribute('markerHeight', '7');
    arrowheadPurple.setAttribute('refX', '9');
    arrowheadPurple.setAttribute('refY', '3.5');
    arrowheadPurple.setAttribute('orient', 'auto');
    const arrowPathPurple = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    arrowPathPurple.setAttribute('d', 'M0,0 L0,7 L10,3.5 z');
    arrowPathPurple.setAttribute('fill', '#6f42c1');
    arrowheadPurple.appendChild(arrowPathPurple);
    defs.appendChild(arrowheadPurple);

    // Yellow arrowhead for specific result connections
    const arrowheadYellow = document.createElementNS('http://www.w3.org/2000/svg', 'marker');
    arrowheadYellow.setAttribute('id', 'arrowhead-yellow');
    arrowheadYellow.setAttribute('markerWidth', '10');
    arrowheadYellow.setAttribute('markerHeight', '7');
    arrowheadYellow.setAttribute('refX', '9');
    arrowheadYellow.setAttribute('refY', '3.5');
    arrowheadYellow.setAttribute('orient', 'auto');
    const arrowPathYellow = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    arrowPathYellow.setAttribute('d', 'M0,0 L0,7 L10,3.5 z');
    arrowPathYellow.setAttribute('fill', '#ffc107');
    arrowheadYellow.appendChild(arrowPathYellow);
    defs.appendChild(arrowheadYellow);

    // Red arrowhead for evidence connections
    const arrowheadRed = document.createElementNS('http://www.w3.org/2000/svg', 'marker');
    arrowheadRed.setAttribute('id', 'evidence-arrowhead');
    arrowheadRed.setAttribute('markerWidth', '10');
    arrowheadRed.setAttribute('markerHeight', '7');
    arrowheadRed.setAttribute('refX', '9');
    arrowheadRed.setAttribute('refY', '3.5');
    arrowheadRed.setAttribute('orient', 'auto');
    const arrowPathRed = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    arrowPathRed.setAttribute('d', 'M0,0 L0,7 L10,3.5 z');
    arrowPathRed.setAttribute('fill', '#ff6b6b');
    arrowheadRed.appendChild(arrowPathRed);
    defs.appendChild(arrowheadRed);

    // Add defs to document
    document.body.appendChild(defs);
    setTimeout(() => {
        const svg = document.getElementById('flowArrows');
        if (svg && !svg.querySelector('defs')) {
            svg.appendChild(defs.cloneNode(true));
        }
    }, 50);

    document.head.appendChild(style);
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

    // Set SVG dimensions - calculate total height needed
    const container = svg.parentElement;
    const totalSteps = steps.length;
    const reasoning = flowTrajectoryData.reasoning || (flowTrajectoryData.trajectory && flowTrajectoryData.trajectory.reasoning);
    const answer = flowTrajectoryData.answer || (flowTrajectoryData.trajectory && flowTrajectoryData.trajectory.answer);
    const extraNodes = (reasoning ? 1 : 0) + (answer ? 1 : 0);
    const minHeight = (totalSteps + extraNodes) * 400 + 100; // Add padding

    svg.setAttribute('width', container.offsetWidth);
    svg.setAttribute('height', Math.max(container.offsetHeight, minHeight));

    // Also set the flow diagram container height
    const flowDiagram = container.querySelector('.flow-diagram');
    if (flowDiagram) {
        flowDiagram.style.minHeight = minHeight + 'px';
    }

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
    if (resultIndex !== null && resultIndex !== undefined) {
        fromElement = document.getElementById(`obs-${fromStep}-${resultIndex}`);
        console.log(`Looking for obs element: obs-${fromStep}-${resultIndex}`, fromElement);
        // If specific observation item not found, fall back to step
        if (!fromElement) {
            console.log(`Obs element not found, falling back to step-${fromStep}`);
            fromElement = document.getElementById(`step-${fromStep}`);
        }
    } else {
        fromElement = document.getElementById(`step-${fromStep}`);
    }

    toElement = document.getElementById(`step-${toStep}`);

    if (!fromElement || !toElement) return;

    const fromRect = fromElement.getBoundingClientRect();
    const toRect = toElement.getBoundingClientRect();
    const containerRect = svg.parentElement.getBoundingClientRect();

    // Calculate positions relative to container
    let fromX, fromY;
    if (resultIndex !== null && resultIndex !== undefined && fromElement.id.includes('obs-')) {
        // Arrow from specific observation item - from right side
        fromX = fromRect.right - containerRect.left;
        fromY = fromRect.top + fromRect.height / 2 - containerRect.top;
    } else {
        // Arrow from step - from bottom center
        fromX = fromRect.left + fromRect.width / 2 - containerRect.left;
        fromY = fromRect.bottom - containerRect.top;
    }

    // Arrow to step - to top center
    const toX = toRect.left + toRect.width / 2 - containerRect.left;
    const toY = toRect.top - containerRect.top;

    // Create curved path for vertical flow
    let path;
    let midX; // Declare midX outside the condition

    if (resultIndex !== null && resultIndex !== undefined && fromElement.id.includes('obs-')) {
        // Curved path from observation item to step top
        midX = (fromX + toX) / 2 + 30;
        path = `M ${fromX} ${fromY} Q ${midX} ${fromY} ${midX} ${(fromY + toY) / 2} Q ${midX} ${toY} ${toX} ${toY}`;
    } else {
        // Straight vertical path from step bottom to next step top
        midX = (fromX + toX) / 2; // Set midX for vertical paths too
        const midY = (fromY + toY) / 2;
        path = `M ${fromX} ${fromY} Q ${fromX} ${midY} ${toX} ${midY} Q ${toX} ${midY} ${toX} ${toY}`;
    }

    // Create path element
    const pathElement = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    pathElement.setAttribute('d', path);
    pathElement.setAttribute('stroke', '#007bff'); // All arrows blue
    pathElement.setAttribute('stroke-width', '2');
    pathElement.setAttribute('fill', 'none');
    pathElement.setAttribute('marker-end', 'url(#arrowhead)'); // All arrows use blue arrowhead

    // Add dash pattern for specific result connections
    if (resultIndex !== null && resultIndex !== undefined) {
        pathElement.setAttribute('stroke-dasharray', '5,5');
    }

    svg.appendChild(pathElement);

    // Add label text if result index specified
    if (resultIndex !== null && resultIndex !== undefined) {
        const textElement = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        textElement.setAttribute('x', midX);
        textElement.setAttribute('y', (fromY + toY) / 2 - 5);
        textElement.setAttribute('text-anchor', 'middle');
        textElement.setAttribute('fill', '#007bff'); // Blue label text
        textElement.setAttribute('font-size', '10px');
        textElement.setAttribute('font-weight', 'bold');
        textElement.textContent = `[${resultIndex}]`;
        svg.appendChild(textElement);
    }
}

// Draw arrow to reasoning node
function drawArrowToReasoning(svg, fromStep) {
    const fromElement = document.getElementById(`step-${fromStep}`);
    const toElement = document.getElementById('reasoning-node');

    if (!fromElement || !toElement) return;

    const fromRect = fromElement.getBoundingClientRect();
    const toRect = toElement.getBoundingClientRect();
    const containerRect = svg.parentElement.getBoundingClientRect();

    const fromX = fromRect.left + fromRect.width / 2 - containerRect.left;
    const fromY = fromRect.bottom - containerRect.top;
    const toX = toRect.left + toRect.width / 2 - containerRect.left;
    const toY = toRect.top - containerRect.top;

    const midY = (fromY + toY) / 2;
    const path = `M ${fromX} ${fromY} Q ${fromX} ${midY} ${toX} ${midY} Q ${toX} ${midY} ${toX} ${toY}`;

    const pathElement = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    pathElement.setAttribute('d', path);
    pathElement.setAttribute('stroke', '#007bff'); // Blue arrow for reasoning too
    pathElement.setAttribute('stroke-width', '2');
    pathElement.setAttribute('fill', 'none');
    pathElement.setAttribute('marker-end', 'url(#arrowhead)'); // Use blue arrowhead

    svg.appendChild(pathElement);
}

// Draw arrow from reasoning to answer node
function drawArrowFromReasoningToAnswer(svg) {
    const fromElement = document.getElementById('reasoning-node');
    const toElement = document.getElementById('answer-node');

    if (!fromElement || !toElement) return;

    const fromRect = fromElement.getBoundingClientRect();
    const toRect = toElement.getBoundingClientRect();
    const containerRect = svg.parentElement.getBoundingClientRect();

    const fromX = fromRect.left + fromRect.width / 2 - containerRect.left;
    const fromY = fromRect.bottom - containerRect.top;
    const toX = toRect.left + toRect.width / 2 - containerRect.left;
    const toY = toRect.top - containerRect.top;

    const midY = (fromY + toY) / 2;
    const path = `M ${fromX} ${fromY} Q ${fromX} ${midY} ${toX} ${midY} Q ${toX} ${midY} ${toX} ${toY}`;

    const pathElement = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    pathElement.setAttribute('d', path);
    pathElement.setAttribute('stroke', '#007bff'); // Blue arrow
    pathElement.setAttribute('stroke-width', '2');
    pathElement.setAttribute('fill', 'none');
    pathElement.setAttribute('marker-end', 'url(#arrowhead)'); // Use blue arrowhead

    svg.appendChild(pathElement);
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
        console.log('No evidence to trajectory mapping available');
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
        const trajectoryMapping = evidenceToTrajectory[evidenceId];

        if (trajectoryMapping && Array.isArray(trajectoryMapping)) {
            const stepNumber = trajectoryMapping[0];
            const observationIndex = trajectoryMapping.length > 1 ? trajectoryMapping[1] : null;

            // Find the target element in trajectory flow
            let targetElement;
            if (observationIndex !== null) {
                // Connect to specific observation item
                targetElement = document.querySelector(`#obs-${stepNumber}-${observationIndex}`);
            } else {
                // Connect to step node
                targetElement = document.querySelector(`#step-${stepNumber}`);
            }

            if (targetElement) {
                drawEvidenceConnectionLine(evidenceItem, targetElement, uniqueConnectionId);
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
        top: 0;
        left: 0;
        width: 100vw;
        height: 100vh;
        pointer-events: none;
        z-index: 1000;
        overflow: visible;
    `;

    // Create defs section for the marker if not exists
    let defs = connectionSvg.querySelector('defs');
    if (!defs) {
        defs = document.createElementNS('http://www.w3.org/2000/svg', 'defs');

        // Red arrowhead for evidence connections
        const arrowhead = document.createElementNS('http://www.w3.org/2000/svg', 'marker');
        arrowhead.setAttribute('id', `evidence-arrowhead-${evidenceId}`);
        arrowhead.setAttribute('markerWidth', '10');
        arrowhead.setAttribute('markerHeight', '7');
        arrowhead.setAttribute('refX', '9');
        arrowhead.setAttribute('refY', '3.5');
        arrowhead.setAttribute('orient', 'auto');

        const arrowPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        arrowPath.setAttribute('d', 'M0,0 L0,7 L10,3.5 z');
        arrowPath.setAttribute('fill', '#ff6b6b');

        arrowhead.appendChild(arrowPath);
        defs.appendChild(arrowhead);
        connectionSvg.appendChild(defs);
    }

    // Add line element
    const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    line.setAttribute('stroke', '#ff6b6b');
    line.setAttribute('stroke-width', '3');
    line.setAttribute('stroke-dasharray', '8,4');
    line.setAttribute('marker-end', `url(#evidence-arrowhead-${evidenceId})`);
    line.setAttribute('opacity', '0.9');

    connectionSvg.appendChild(line);
    document.body.appendChild(connectionSvg);

    // Update line position
    updateEvidenceConnectionLine(evidenceElement, trajectoryElement, line);

    console.log(`Drawing evidence connection from ${evidenceId}:`, {
        evidenceElement,
        trajectoryElement,
        connectionSvg,
        line,
        evidenceRect: evidenceElement.getBoundingClientRect(),
        trajectoryRect: trajectoryElement.getBoundingClientRect()
    });

    // Store reference for updates
    if (!window.evidenceConnections) {
        window.evidenceConnections = new Map();
    }
    window.evidenceConnections.set(evidenceId, {
        svg: connectionSvg,
        line: line,
        evidenceElement: evidenceElement,
        trajectoryElement: trajectoryElement
    });
}

// Update the position of an evidence connection line
function updateEvidenceConnectionLine(evidenceElement, trajectoryElement, lineElement) {
    try {
        const evidenceRect = evidenceElement.getBoundingClientRect();
        const trajectoryRect = trajectoryElement.getBoundingClientRect();

        // Calculate connection points - evidence left to trajectory right
        const startX = evidenceRect.left;
        const startY = evidenceRect.top + evidenceRect.height / 2;
        const endX = trajectoryRect.right;
        const endY = trajectoryRect.top + trajectoryRect.height / 2;

        console.log(`Updating connection line:`, {
            evidenceRect,
            trajectoryRect,
            startX, startY, endX, endY
        });

        lineElement.setAttribute('x1', startX);
        lineElement.setAttribute('y1', startY);
        lineElement.setAttribute('x2', endX);
        lineElement.setAttribute('y2', endY);

        // Make sure the line is visible
        lineElement.setAttribute('stroke', '#ff6b6b');
        lineElement.setAttribute('stroke-width', '3');
        lineElement.setAttribute('opacity', '0.9');
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



// Export functions for use in other scripts
window.trajectoryFlow = {
    init: initTrajectoryFlow,
    loadData: loadTrajectoryFlowData,
    show: showTrajectoryFlow,
    hide: hideTrajectoryFlow,
    isVisible: isTrajectoryFlowCurrentlyVisible,
    drawEvidenceConnections: drawEvidenceConnections,
    clearEvidenceConnections: clearEvidenceConnections,
    // Debug functions
    testDirectConnection: function () {
        // Create a simple test connection line
        const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        svg.style.cssText = `
            position: fixed;
            top: 0;
            left: 0;
            width: 100vw;
            height: 100vh;
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

        console.log('Test connection created:', svg, line);

        // Remove after 5 seconds
        setTimeout(() => {
            svg.remove();
            console.log('Test connection removed');
        }, 5000);
    }
};