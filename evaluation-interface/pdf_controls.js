// PDF Control Functions for Enhanced PDF Viewer
let pdfZoomLevels = {}; // Store zoom levels for each PDF

// Zoom PDF function using PDF.js
function zoomPDF(evidenceId, direction) {
    const currentZoom = pdfZoomLevels[evidenceId] || 100;
    let newZoom = currentZoom;

    if (direction === 'in') {
        newZoom = Math.min(currentZoom + 25, 300); // Max 300%
    } else if (direction === 'out') {
        newZoom = Math.max(currentZoom - 25, 50); // Min 50%
    }

    pdfZoomLevels[evidenceId] = newZoom;

    // Update zoom level display
    const zoomLevelSpan = document.getElementById(`zoomLevel_${evidenceId}`);
    if (zoomLevelSpan) {
        zoomLevelSpan.textContent = `${newZoom}%`;
    }

    // Send zoom command to PDF.js viewer
    const pdfFrame = document.getElementById(`pdfFrame_${evidenceId}`);
    if (pdfFrame && pdfFrame.contentWindow) {
        try {
            // Use PDF.js viewer commands
            const zoomValue = newZoom / 100; // PDF.js expects decimal values
            pdfFrame.contentWindow.postMessage({
                type: 'setZoom',
                zoom: zoomValue
            }, '*');
        } catch (error) {
            console.warn('Could not communicate with PDF.js viewer:', error);
            // Fallback: update iframe src with new zoom
            const currentSrc = pdfFrame.src;
            const basePath = currentSrc.split('#')[0];
            const pageMatch = currentSrc.match(/page=(\d+)/);
            const currentPage = pageMatch ? pageMatch[1] : '1';
            pdfFrame.src = `${basePath}#page=${currentPage}&zoom=${newZoom}`;
        }
    }
}

// Navigate PDF pages
function navigatePage(evidenceId, direction) {
    const pageInput = document.getElementById(`pageInput_${evidenceId}`);
    if (!pageInput) return;

    let currentPage = parseInt(pageInput.value) || 1;

    if (direction === 'next') {
        currentPage += 1;
    } else if (direction === 'prev') {
        currentPage = Math.max(1, currentPage - 1);
    }

    pageInput.value = currentPage;
    jumpToPage(evidenceId, currentPage);
}

// Jump to specific page
function jumpToPage(evidenceId, pageNumber) {
    const pdfFrame = document.getElementById(`pdfFrame_${evidenceId}`);
    const pageInput = document.getElementById(`pageInput_${evidenceId}`);

    if (!pdfFrame || !pageInput) return;

    const page = Math.max(1, parseInt(pageNumber) || 1);
    pageInput.value = page;

    // Send page navigation command to PDF.js viewer
    if (pdfFrame.contentWindow) {
        try {
            pdfFrame.contentWindow.postMessage({
                type: 'setPage',
                page: page
            }, '*');
        } catch (error) {
            console.warn('Could not communicate with PDF.js viewer:', error);
            // Fallback: update iframe src
            const currentZoom = pdfZoomLevels[evidenceId] || 100;
            const currentSrc = pdfFrame.src;
            const basePath = currentSrc.split('#')[0];
            pdfFrame.src = `${basePath}#page=${page}&zoom=${currentZoom}`;
        }
    }
}

// Jump to page with specific highlight using exact text matching
function jumpToPageWithHighlight(evidenceId, pageNumber, highlightText) {
    console.log(`jumpToPageWithHighlight called: evidenceId=${evidenceId}, page=${pageNumber}, text=${highlightText}`);

    // First jump to the specified page
    jumpToPage(evidenceId, pageNumber);

    // Then perform exact text search on that page
    const pdfFrame = document.getElementById(`pdfFrame_${evidenceId}`);
    console.log('PDF frame found:', !!pdfFrame);

    if (pdfFrame && pdfFrame.contentWindow) {
        try {
            // Perform exact text search on the specified page
            const performExactSearch = () => {
                try {
                    const pdfViewer = pdfFrame.contentWindow.PDFViewerApplication;
                    if (pdfViewer && pdfViewer.findController) {
                        console.log('PDF.js findController found, executing exact text search...');
                        console.log('Available methods:', Object.keys(pdfViewer.findController));

                        // Clear previous search
                        if (typeof pdfViewer.findController.reset === 'function') {
                            pdfViewer.findController.reset();
                        }

                        // Method 1: Try the most direct approach - setting state and executing
                        try {
                            // Set the search query directly
                            pdfViewer.findController.state = {
                                query: highlightText,
                                caseSensitive: false,
                                entireWord: false,
                                highlightAll: false,
                                findPrevious: false,
                                matchDiacritics: false,
                                phraseSearch: true
                            };

                            // Execute the search
                            if (typeof pdfViewer.findController.executeCommand === 'function') {
                                pdfViewer.findController.executeCommand('find');
                                console.log('Search executed via state + executeCommand');
                                return;
                            } else if (typeof pdfViewer.findController._search === 'function') {
                                pdfViewer.findController._search();
                                console.log('Search executed via _search()');
                                return;
                            }
                        } catch (e) {
                            console.warn('Direct state method failed:', e);
                        }

                        // Method 2: Try using the findbar programmatically
                        try {
                            const findInput = pdfViewer.findBar?.findField;
                            if (findInput) {
                                findInput.value = highlightText;
                                // Trigger the search event
                                const event = new Event('input', { bubbles: true });
                                findInput.dispatchEvent(event);
                                console.log('Search executed via findBar input');
                                return;
                            }
                        } catch (e) {
                            console.warn('FindBar method failed:', e);
                        }

                        // Method 3: Try keyboard shortcut simulation
                        try {
                            // Open find bar with Ctrl+F
                            const ctrlF = new KeyboardEvent('keydown', {
                                key: 'f',
                                ctrlKey: true,
                                bubbles: true
                            });
                            pdfFrame.contentDocument.dispatchEvent(ctrlF);

                            setTimeout(() => {
                                const findInput = pdfViewer.findBar?.findField;
                                if (findInput) {
                                    findInput.value = highlightText;
                                    const enterEvent = new KeyboardEvent('keydown', {
                                        key: 'Enter',
                                        bubbles: true
                                    });
                                    findInput.dispatchEvent(enterEvent);
                                    console.log('Search executed via keyboard simulation');
                                }
                            }, 100);
                            return;
                        } catch (e) {
                            console.warn('Keyboard simulation failed:', e);
                        }

                        console.warn('All direct methods failed, using URL fallback');
                        performUrlSearch();
                    } else {
                        console.warn('PDF.js findController not available, using URL fallback');
                        performUrlSearch();
                    }
                } catch (error) {
                    console.warn('Direct API search failed:', error);
                    performUrlSearch();
                }
            };

            // URL-based search fallback with exact phrase matching
            const performUrlSearch = () => {
                try {
                    console.log('Using URL-based exact phrase search...');
                    const currentSrc = pdfFrame.src;
                    const basePath = currentSrc.split('#')[0];

                    // Simple approach: add search term directly to URL
                    const searchQuery = encodeURIComponent(highlightText);
                    const newSrc = `${basePath}#page=${pageNumber}&search=${searchQuery}&zoom=125`;

                    console.log('New URL with search:', newSrc);
                    console.log('Original text to match:', highlightText);

                    // Force reload with new search
                    pdfFrame.src = newSrc;
                    console.log('URL updated with search parameters');
                } catch (error) {
                    console.warn('URL search approach failed:', error);
                }
            };

            // Execute search when PDF is fully ready
            const executeWhenReady = () => {
                const viewer = pdfFrame.contentWindow?.PDFViewerApplication;
                if (viewer && viewer.pdfDocument && viewer.pdfViewer && viewer.findController) {
                    console.log('PDF.js fully loaded, executing search...');
                    performExactSearch();
                } else {
                    console.log('PDF.js not fully ready, waiting...');
                    setTimeout(executeWhenReady, 300);
                }
            };

            // Start checking when iframe content is ready
            if (pdfFrame.contentDocument && pdfFrame.contentDocument.readyState === 'complete') {
                executeWhenReady();
            } else {
                pdfFrame.onload = executeWhenReady;
            }

        } catch (error) {
            console.warn('Error in exact text search:', error);
        }
    } else {
        console.error('PDF frame not found for evidenceId:', evidenceId);
    }
}

// Override the original showLocalPDF function with enhanced version
async function showLocalPDF(evidenceMapping, evidenceId) {
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
    loadingOverlay.textContent = 'Loading file...';

    try {
        // Construct the file path
        const runType = getRunTypeFromCurrentFolder();
        const problemFolder = getProblemFolderFromCurrentFolder();

        // The file should be accessible via the API server
        const filePath = `/api/data/${runType}/${problemFolder}/evidences/${evidenceMapping.filename}`;

        // Determine if it's a PDF file for enhanced features
        const isPDF = evidenceMapping.filename.toLowerCase().endsWith('.pdf');

        // Extract highlight information from new structure
        let highlightPage = 1;
        let highlightInfo = [];

        if (evidenceMapping.highlights && Array.isArray(evidenceMapping.highlights)) {
            highlightInfo = evidenceMapping.highlights;
            // Find the first highlight page or use page 1
            const firstHighlight = highlightInfo.find(h => h.page);
            highlightPage = firstHighlight ? firstHighlight.page : 1;
        } else if (evidenceMapping.highlight) {
            // Fallback to old highlight format
            highlightPage = evidenceMapping.highlight;
        }

        // Create enhanced file viewer HTML
        const fileViewerHtml = `
            <div class="pdf-viewer-container" id="pdfContainer_${evidenceId}">
                <div class="pdf-header">
                    <h4>📄 ${evidenceMapping.filename}</h4>
                    <div class="pdf-info">
                        <small>Evidence ID: ${evidenceId}</small>
                        ${evidenceMapping.original_url ? `<br><small>Original: <a href="${evidenceMapping.original_url}" target="_blank">${evidenceMapping.original_url}</a></small>` : ''}
                        ${isPDF && highlightPage > 1 ? `<br><small>Auto-jump to page: ${highlightPage}</small>` : ''}
                        ${highlightInfo.length > 0 ? `<br><small>Highlights: ${highlightInfo.length} items</small>` : ''}
                    </div>
                    <div class="pdf-actions">
                        <button onclick="openPDFInNewTab('${filePath}')" class="btn btn-small">
                            🔗 Open in New Tab
                        </button>
                    </div>
                </div>
                <div class="pdf-content" id="pdfContent_${evidenceId}">
                    ${isPDF ?
                `<iframe src="/pdfjs/web/viewer.html?file=${encodeURIComponent(filePath)}#page=${highlightPage}&zoom=page-fit" 
                        id="pdfFrame_${evidenceId}" 
                        width="100%" 
                        height="600px" 
                        frameborder="0"
                        onload="handlePDFFrameLoad('${evidenceId}')"
                        onerror="handlePDFFrameError('${evidenceId}')"></iframe>`
                : evidenceMapping.filename.toLowerCase().endsWith('.txt') ?
                    `<iframe src="${filePath}" width="100%" height="600px" style="border: 1px solid #ddd; background: white;"></iframe>` :
                    ['png', 'jpg', 'jpeg', 'gif'].includes(evidenceMapping.filename.split('.').pop().toLowerCase()) ?
                        `<img src="${filePath}" style="max-width: 100%; height: auto;" alt="${evidenceMapping.filename}" />` :
                        `<embed src="${filePath}" type="application/pdf" width="100%" height="600px" />`
            }
                </div>
                <div class="citation-info">
                    <h5>Citation Information:</h5>
                    <p>${evidenceMapping.citation || 'No citation available'}</p>
                    ${highlightInfo.length > 0 ? `
                        <div class="highlight-info">
                            <h6>Highlights:</h6>
                            <div class="highlight-list">
                                ${highlightInfo.map((highlight, index) => {
                const safeText = highlight.text.replace(/'/g, "\\'").replace(/"/g, '\\"').replace(/\n/g, '\\n');
                return `
                                    <div class="highlight-item" onclick="jumpToPageWithHighlight('${evidenceId}', ${highlight.page}, '${safeText}')" style="cursor: pointer; padding: 10px; margin: 6px 0; border: 2px solid #ffc107; border-radius: 6px; background: #fff9c4; transition: all 0.2s ease;">
                                        <span class="highlight-page" style="font-weight: bold; color: #856404; font-size: 0.9rem;">📄 Page ${highlight.page}:</span>
                                        <span class="highlight-text" style="background-color: #ffeb3b; padding: 4px 8px; border-radius: 4px; margin-left: 8px; border: 1px solid #ffc107; display: inline-block; margin-top: 4px;">${highlight.text}</span>
                                        <div style="font-size: 0.8rem; color: #6c757d; margin-top: 4px;">🎯 Click to jump and highlight in PDF</div>
                                    </div>
                                `;
            }).join('')}
                            </div>
                        </div>
                    ` : ''}
                </div>
            </div>
        `;

        // Display the file viewer
        readingModeViewer.innerHTML = fileViewerHtml;
        readingModeViewer.style.display = 'block';

        // Initialize zoom level and setup PDF.js communication
        if (isPDF) {
            pdfZoomLevels[evidenceId] = 100;

            // Setup PDF.js viewer communication after iframe loads
            const pdfFrame = document.getElementById(`pdfFrame_${evidenceId}`);
            if (pdfFrame) {
                pdfFrame.onload = function () {
                    // Apply initial highlights if available
                    if (highlightInfo.length > 0) {
                        const firstHighlight = highlightInfo[0];
                        if (firstHighlight.text && firstHighlight.page) {
                            console.log('Applying initial highlight on load...');

                            // Check if PDF.js is ready, if not, wait for it
                            const checkAndApplyHighlight = () => {
                                if (pdfFrame.contentWindow && pdfFrame.contentWindow.PDFViewerApplication) {
                                    // PDF.js is ready, apply highlight immediately
                                    jumpToPageWithHighlight(evidenceId, firstHighlight.page, firstHighlight.text);
                                } else {
                                    // PDF.js not ready yet, wait a bit and try again
                                    setTimeout(checkAndApplyHighlight, 200);
                                }
                            };

                            checkAndApplyHighlight();
                        }
                    }
                };
            }
        }

    } catch (error) {
        console.error('Error displaying file:', error);
        showErrorInReadingMode('Failed to load evidence file');
    } finally {
        loadingOverlay.style.display = 'none';
    }
}

// Listen for messages from PDF.js viewer
window.addEventListener('message', function (event) {
    if (event.data && event.data.type === 'pdfViewerUpdate') {
        // Handle updates from PDF.js viewer (page changes, zoom changes, etc.)
        const { evidenceId, page, zoom } = event.data;

        if (evidenceId && page) {
            const pageInput = document.getElementById(`pageInput_${evidenceId}`);
            if (pageInput) {
                pageInput.value = page;
            }
        }

        if (evidenceId && zoom) {
            pdfZoomLevels[evidenceId] = Math.round(zoom * 100);
            const zoomLevelSpan = document.getElementById(`zoomLevel_${evidenceId}`);
            if (zoomLevelSpan) {
                zoomLevelSpan.textContent = `${Math.round(zoom * 100)}%`;
            }
        }
    }
}, false);

// Error handling for PDF.js iframe
function handlePDFFrameLoad(evidenceId) {
    console.log(`PDF frame loaded for evidence ${evidenceId}`);

    const pdfFrame = document.getElementById(`pdfFrame_${evidenceId}`);
    if (pdfFrame && pdfFrame.contentWindow) {
        try {
            // Enhanced error suppression for PDF.js iframe
            pdfFrame.contentWindow.addEventListener('error', function (e) {
                if (e.message && (e.message.includes('injectLinkAnnotations') ||
                    e.message.includes('annotation') ||
                    e.message.includes('render') && e.message.includes('method must be called'))) {
                    console.warn('Suppressed PDF.js iframe annotation error:', e.message);
                    e.preventDefault();
                    return false;
                }
            });

            // Also suppress console.error in the iframe if possible
            if (pdfFrame.contentWindow.console) {
                const originalIframeError = pdfFrame.contentWindow.console.error;
                pdfFrame.contentWindow.console.error = function (...args) {
                    const message = args.join(' ');
                    if (message.includes('injectLinkAnnotations') ||
                        message.includes('AnnotationLayerBuilder') ||
                        message.includes('pdf_page_view.js')) {
                        console.warn('Suppressed PDF.js iframe console error:', message);
                        return;
                    }
                    originalIframeError.apply(pdfFrame.contentWindow.console, args);
                };
            }
        } catch (error) {
            // Ignore cross-origin errors - this is expected
            console.log('Could not access PDF.js iframe content (cross-origin) - this is normal');
        }
    }
}

function handlePDFFrameError(evidenceId) {
    console.error(`PDF frame failed to load for evidence ${evidenceId}`);
    const pdfContent = document.getElementById(`pdfContent_${evidenceId}`);
    if (pdfContent) {
        pdfContent.innerHTML = `
            <div style="padding: 2rem; text-align: center; color: #666;">
                <p>⚠️ Failed to load PDF viewer</p>
                <p>Try opening the file in a new tab instead.</p>
            </div>
        `;
    }
}

// Open PDF in new tab
function openPDFInNewTab(filePath) {
    console.log('Opening PDF in new tab:', filePath);
    const newTabUrl = `/pdfjs/web/viewer.html?file=${encodeURIComponent(filePath)}`;
    window.open(newTabUrl, '_blank');
}

// Make functions globally accessible
window.handlePDFFrameLoad = handlePDFFrameLoad;
window.handlePDFFrameError = handlePDFFrameError;
window.showLocalPDF = showLocalPDF;
window.jumpToPageWithHighlight = jumpToPageWithHighlight;
window.jumpToPage = jumpToPage;
window.zoomPDF = zoomPDF;
window.navigatePage = navigatePage;
window.openPDFInNewTab = openPDFInNewTab;

// Global error handler to suppress PDF.js annotation errors
window.addEventListener('error', function (e) {
    if (e.message && (e.message.includes('injectLinkAnnotations') ||
        e.message.includes('annotation') ||
        e.message.includes('render') && e.message.includes('method must be called'))) {
        console.warn('Suppressed PDF.js annotation error:', e.message);
        e.preventDefault();
        return false;
    }
});

// Enhanced console error suppression for PDF.js
const originalConsoleError = console.error;
console.error = function (...args) {
    const message = args.join(' ');
    if (message.includes('injectLinkAnnotations') ||
        message.includes('render') && message.includes('annotation') ||
        message.includes('AnnotationLayerBuilder') ||
        message.includes('pdf_page_view.js') ||
        message.includes('method must be called before')) {
        console.warn('Suppressed PDF.js error:', message);
        return;
    }
    originalConsoleError.apply(console, args);
};

// Additional suppression for unhandled promise rejections from PDF.js
window.addEventListener('unhandledrejection', function (e) {
    if (e.reason && e.reason.message &&
        (e.reason.message.includes('injectLinkAnnotations') ||
            e.reason.message.includes('annotation'))) {
        console.warn('Suppressed PDF.js promise rejection:', e.reason.message);
        e.preventDefault();
        return false;
    }
});