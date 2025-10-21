// PDF Control Functions for Enhanced PDF Viewer
let pdfZoomLevels = {}; // Store zoom levels for each PDF
let activeHighlights = {}; // Store active highlights for repositioning on zoom

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

            // Reposition highlights after zoom change
            setTimeout(() => {
                repositionHighlights(evidenceId);
            }, 500); // Wait for zoom to complete
        } catch (error) {
            console.warn('Could not communicate with PDF.js viewer:', error);
            // Fallback: update iframe src with new zoom
            const currentSrc = pdfFrame.src;
            const basePath = currentSrc.split('#')[0];
            const pageMatch = currentSrc.match(/page=(\d+)/);
            const currentPage = pageMatch ? pageMatch[1] : '1';
            pdfFrame.src = `${basePath}#page=${currentPage}&zoom=${newZoom}`;

            // Reposition highlights after iframe reload
            setTimeout(() => {
                repositionHighlights(evidenceId);
            }, 1000);
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

// Apply all highlights to the PDF at once
async function applyAllHighlights(evidenceId, highlightInfo) {
    const pdfFrame = document.getElementById(`pdfFrame_${evidenceId}`);
    if (!pdfFrame) {
        console.error('PDF frame not found for evidenceId:', evidenceId);
        return;
    }

    // Find the first highlight page and jump to it
    const firstHighlight = highlightInfo.find(h => h.page);
    const startPage = firstHighlight ? firstHighlight.page : 1;

    jumpToPage(evidenceId, startPage);

    // Apply text layer-based precise highlighting and wait for completion
    const highlightingSuccess = await applyHighlightsViaPDFJS(evidenceId, highlightInfo);

    // Only scroll to highlights if highlighting was successful
    if (highlightingSuccess) {
        scrollToFirstHighlight(evidenceId);
    }
} async function applyHighlightsViaPDFJS(evidenceId, highlightInfo) {
    const pdfFrame = document.getElementById(`pdfFrame_${evidenceId}`);
    if (!pdfFrame || !pdfFrame.contentWindow) return;

    // Store highlight info for repositioning on zoom
    activeHighlights[evidenceId] = highlightInfo;

    // Wait for PDF.js to be fully loaded
    const waitForPDFReady = () => {
        return new Promise((resolve, reject) => {
            const maxAttempts = 30;
            let attempts = 0;

            const checkReady = () => {
                try {
                    const pdfViewer = pdfFrame.contentWindow.PDFViewerApplication;
                    if (pdfViewer && pdfViewer.pdfDocument) {
                        resolve(pdfViewer);
                    } else if (attempts < maxAttempts) {
                        attempts++;
                        setTimeout(checkReady, 200);
                    } else {
                        reject(new Error('PDF.js not ready after maximum attempts - possible server error'));
                    }
                } catch (error) {
                    // Handle iframe access errors (cross-origin, etc.)
                    if (attempts < maxAttempts) {
                        attempts++;
                        setTimeout(checkReady, 200);
                    } else {
                        reject(new Error(`PDF.js access error: ${error.message}`));
                    }
                }
            };

            checkReady();
        });
    };

    try {
        const pdfViewer = await waitForPDFReady();

        // Process each highlight and wait for completion
        const highlightPromises = [];
        for (const highlight of highlightInfo) {
            if (!highlight.text || !highlight.page) continue;

            try {
                const promise = highlightTextOnPage(pdfViewer, highlight.page, highlight.text, evidenceId);
                highlightPromises.push(promise);
            } catch (error) {
                console.warn(`Failed to highlight "${highlight.text}" on page ${highlight.page}:`, error);
            }
        }

        // Wait for all highlights to be processed
        await Promise.all(highlightPromises);

        return true; // Indicate successful completion

    } catch (error) {
        console.warn('Error in text layer highlighting:', error);
        return false;
    }
}// 在指定页面上精确高亮文本
async function highlightTextOnPage(pdfViewer, pageNumber, searchText, evidenceId) {
    try {
        // 获取页面对象
        const page = await pdfViewer.pdfDocument.getPage(pageNumber);
        const textContent = await page.getTextContent();

        // 提取所有文本内容
        const textItems = textContent.items;
        let fullText = '';
        let itemPositions = [];

        // 构建完整文本并记录每个字符的位置信息
        textItems.forEach((item, itemIndex) => {
            const startPos = fullText.length;
            fullText += item.str;
            const endPos = fullText.length;

            itemPositions.push({
                itemIndex,
                startPos,
                endPos,
                transform: item.transform,
                str: item.str,
                fontName: item.fontName,
                height: item.height,
                width: item.width
            });

            // 添加空格分隔（如果需要）
            if (itemIndex < textItems.length - 1) {
                fullText += ' ';
                itemPositions.push({
                    itemIndex: -1, // 空格标记
                    startPos: fullText.length - 1,
                    endPos: fullText.length,
                    isSpace: true
                });
            }
        });

        // 查找匹配的文本
        const searchIndex = fullText.toLowerCase().indexOf(searchText.toLowerCase());
        if (searchIndex === -1) {
            console.warn(`Text "${searchText}" not found on page ${pageNumber}`);
            return;
        }

        // 找到匹配文本的位置信息
        const matchEndIndex = searchIndex + searchText.length;
        const matchingItems = itemPositions.filter(pos =>
            !pos.isSpace &&
            pos.startPos < matchEndIndex &&
            pos.endPos > searchIndex
        );

        if (matchingItems.length > 0) {
            // 创建高亮覆盖层
            await createHighlightOverlay(pdfViewer, pageNumber, matchingItems, evidenceId);
        }

    } catch (error) {
        console.error(`Error highlighting text on page ${pageNumber}:`, error);
    }
}

// Reposition highlights after zoom changes
async function repositionHighlights(evidenceId) {
    const highlightInfo = activeHighlights[evidenceId];
    if (!highlightInfo || highlightInfo.length === 0) {
        return;
    }

    // Clear existing highlights
    clearHighlights(evidenceId);

    // Reapply highlights with new zoom level
    const pdfFrame = document.getElementById(`pdfFrame_${evidenceId}`);
    if (!pdfFrame || !pdfFrame.contentWindow) return;

    try {
        const pdfViewer = pdfFrame.contentWindow.PDFViewerApplication;
        if (!pdfViewer || !pdfViewer.pdfDocument) {
            console.warn('PDF.js document not available for repositioning');
            return;
        }

        // Reprocess each highlight
        for (const highlight of highlightInfo) {
            if (!highlight.text || !highlight.page) continue;

            try {
                await highlightTextOnPage(pdfViewer, highlight.page, highlight.text, evidenceId);
            } catch (error) {
                console.warn(`Failed to reposition highlight "${highlight.text}" on page ${highlight.page}:`, error);
            }
        }

    } catch (error) {
        console.warn('Error repositioning highlights:', error);
    }
}

// Clear existing highlights
function clearHighlights(evidenceId) {
    const highlightElements = document.querySelectorAll(`.custom-highlight-${evidenceId}`);
    highlightElements.forEach(element => {
        element.remove();
    });
}

// 创建高亮覆盖层
async function createHighlightOverlay(pdfViewer, pageNumber, textItems, evidenceId) {
    try {
        // 获取页面容器
        const pageView = pdfViewer.pdfViewer.getPageView(pageNumber - 1);
        if (!pageView || !pageView.div) {
            console.warn(`Page view not found for page ${pageNumber}`);
            return;
        }

        const pageContainer = pageView.div;
        const viewport = pageView.viewport;

        // 为每个文本项创建高亮矩形
        textItems.forEach((item, index) => {
            const transform = item.transform;
            if (!transform || transform.length < 6) return;

            // PDF坐标转换为屏幕坐标
            const x = transform[4];
            const y = transform[5];
            const width = item.width || 50; // 默认宽度
            const height = item.height || 12; // 默认高度

            // 转换坐标系（PDF坐标系Y轴向上，DOM坐标系Y轴向下）
            const screenCoords = viewport.convertToViewportPoint(x, y);
            const screenX = screenCoords[0];
            const screenY = screenCoords[1] - height; // 调整Y坐标

            // 创建高亮元素
            const highlightDiv = document.createElement('div');
            highlightDiv.className = `custom-highlight-${evidenceId}`;
            highlightDiv.style.cssText = `
                position: absolute;
                left: ${screenX}px;
                top: ${screenY - 10}px;
                width: ${width * viewport.scale}px;
                height: ${height * viewport.scale}px;
                background-color: rgba(255, 255, 0, 0.3);
                border: 1px solid rgba(255, 193, 7, 0.8);
                pointer-events: none;
                z-index: 100;
                border-radius: 2px;
            `;

            pageContainer.appendChild(highlightDiv);
        });

    } catch (error) {
        console.error('Error creating highlight overlay:', error);
    }
}

function jumpToPageWithHighlight(evidenceId, pageNumber, highlightText) {
    // First jump to the specified page
    jumpToPage(evidenceId, pageNumber);

    // Then perform exact text search on that page
    const pdfFrame = document.getElementById(`pdfFrame_${evidenceId}`);

    if (pdfFrame && pdfFrame.contentWindow) {
        try {
            // Perform exact text search on the specified page
            const performExactSearch = () => {
                try {
                    const pdfViewer = pdfFrame.contentWindow.PDFViewerApplication;
                    if (pdfViewer && pdfViewer.findController) {
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
                                return;
                            } else if (typeof pdfViewer.findController._search === 'function') {
                                pdfViewer.findController._search();
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
                    const currentSrc = pdfFrame.src;
                    const basePath = currentSrc.split('#')[0];

                    // Simple approach: add search term directly to URL
                    const searchQuery = encodeURIComponent(highlightText);
                    const newSrc = `${basePath}#page=${pageNumber}&search=${searchQuery}&zoom=125`;

                    // Force reload with new search
                    pdfFrame.src = newSrc;
                } catch (error) {
                    console.warn('URL search approach failed:', error);
                }
            };

            // Execute search when PDF is fully ready
            const executeWhenReady = () => {
                const viewer = pdfFrame.contentWindow?.PDFViewerApplication;
                if (viewer && viewer.pdfDocument && viewer.pdfViewer && viewer.findController) {
                    performExactSearch();
                } else {
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
        // Server expects format: /api/data/{folder}/evidences/{filename}
        const filePath = `/api/data/${problemFolder}/evidences/${evidenceMapping.filename}`;
        console.log(`Loading evidence file: ${filePath}`);

        // Log the file path for debugging - PDF.js will handle the actual loading
        console.log(`📄 Preparing to load PDF: ${filePath}`);

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
                    <div class="pdf-actions">
                        <button onclick="openPDFInNewTab('${filePath}')" class="btn btn-small">
                            Open in New Tab
                        </button>
                    </div>
                </div>
                <div class="pdf-content" id="pdfContent_${evidenceId}">
                    ${isPDF ?
                `<iframe src="/pdfjs/web/viewer.html?file=${encodeURIComponent(filePath)}#page=${highlightPage}&zoom=200" 
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
            </div>
        `;

        // Display the file viewer
        readingModeViewer.innerHTML = fileViewerHtml;
        readingModeViewer.style.display = 'block';

        // Initialize zoom level and setup PDF.js communication
        if (isPDF) {
            pdfZoomLevels[evidenceId] = 200; // Set to 200% zoom

            // Setup PDF.js viewer communication after iframe loads
            const pdfFrame = document.getElementById(`pdfFrame_${evidenceId}`);
            if (pdfFrame) {
                pdfFrame.onload = function () {
                    // Apply all highlights if available
                    if (highlightInfo.length > 0) {
                        // Check if PDF.js is ready, if not, wait for it
                        const checkAndApplyAllHighlights = async () => {
                            if (pdfFrame.contentWindow && pdfFrame.contentWindow.PDFViewerApplication) {
                                // PDF.js is ready, apply all highlights
                                await applyAllHighlights(evidenceId, highlightInfo);
                            } else {
                                // PDF.js not ready yet, wait a bit and try again
                                setTimeout(checkAndApplyAllHighlights, 200);
                            }
                        };

                        checkAndApplyAllHighlights();
                    }
                };
            }
        }

    } catch (error) {
        console.error('Error displaying file:', error);
        console.error('File path attempted:', filePath);
        showErrorInReadingMode(`Failed to load evidence file: ${error.message || error}`);
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

            // Reposition highlights when zoom changes from PDF.js controls
            setTimeout(() => {
                repositionHighlights(evidenceId);
            }, 300);
        }
    }
}, false);

// Error handling for PDF.js iframe
function handlePDFFrameLoad(evidenceId) {
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
                } else if (e.message && (e.message.includes('Unexpected server response') ||
                    e.message.includes('400') || e.message.includes('404') ||
                    e.message.includes('Invalid data path format') || e.message.includes('not allowed'))) {
                    console.error('PDF server error:', e.message);
                    const pdfContent = document.getElementById(`pdfContent_${evidenceId}`);
                    if (pdfContent) {
                        pdfContent.innerHTML = `
                            <div style="padding: 2rem; text-align: center; color: #d32f2f;">
                                <p>⚠️ Failed to load PDF file</p>
                                <p>Server error: ${e.message.includes('400') ? 'File not found or not allowed' : 'Server error'}</p>
                                <p>Please check if the file exists on the server.</p>
                                <button onclick="location.reload()" style="margin-top: 1rem; padding: 0.5rem 1rem; background: #1976d2; color: white; border: none; border-radius: 4px; cursor: pointer;">
                                    Retry
                                </button>
                            </div>
                        `;
                    }
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
        }
    }

    // No longer need setTimeout since we wait for highlights to complete
    // The scroll will be triggered automatically after highlights are applied
}

// Enhanced scroll to first highlight with better timing and error handling
function scrollToFirstHighlight(evidenceId) {
    const highlightInfo = activeHighlights[evidenceId];
    if (!highlightInfo || highlightInfo.length === 0) {
        return;
    }

    const pdfFrame = document.getElementById(`pdfFrame_${evidenceId}`);
    if (!pdfFrame || !pdfFrame.contentWindow) {
        return;
    }

    // Enhanced waiting mechanism with better error handling
    const waitForHighlightsAndScroll = (attempts = 0) => {
        const maxAttempts = 20; // Increased max attempts

        try {
            const pdfViewer = pdfFrame.contentWindow.PDFViewerApplication;
            if (!pdfViewer || !pdfViewer.pdfDocument || !pdfViewer.pdfViewer) {
                if (attempts < maxAttempts) {
                    setTimeout(() => waitForHighlightsAndScroll(attempts + 1), 400);
                    return;
                } else {
                    console.warn('PDF.js not ready after maximum attempts');
                    return;
                }
            }

            // Find the first highlight
            const firstHighlight = highlightInfo.find(h => h.text && h.page);
            if (!firstHighlight) {
                return;
            }

            // Get the page view for the first highlight
            const pageView = pdfViewer.pdfViewer.getPageView(firstHighlight.page - 1);
            if (!pageView || !pageView.div) {
                if (attempts < maxAttempts) {
                    setTimeout(() => waitForHighlightsAndScroll(attempts + 1), 400);
                    return;
                } else {
                    console.warn(`Page view not found after maximum attempts`);
                    return;
                }
            }

            // Look for highlight elements in the page
            const highlightElements = pageView.div.querySelectorAll(`.custom-highlight-${evidenceId}`);

            if (highlightElements.length > 0) {
                // Found highlights - scroll to the first one
                const firstHighlightElement = highlightElements[0];

                // Try multiple scrolling approaches for better compatibility
                const performScroll = () => {
                    // Method 1: Use PDF viewer container scroll (most reliable)
                    const viewerContainer = pdfViewer.pdfViewer.container;
                    if (viewerContainer) {
                        try {
                            const containerRect = viewerContainer.getBoundingClientRect();
                            const highlightRect = firstHighlightElement.getBoundingClientRect();

                            // Calculate scroll position to center the highlight
                            const scrollTop = viewerContainer.scrollTop +
                                (highlightRect.top - containerRect.top) -
                                (containerRect.height / 3); // Position at top third instead of center

                            viewerContainer.scrollTo({
                                top: Math.max(0, scrollTop),
                                behavior: 'smooth'
                            });

                            return true;
                        } catch (scrollError) {
                            console.warn('Container scroll failed:', scrollError);
                        }
                    }

                    // Method 2: Fallback to element scrollIntoView
                    try {
                        firstHighlightElement.scrollIntoView({
                            behavior: 'smooth',
                            block: 'start',
                            inline: 'nearest'
                        });
                        return true;
                    } catch (scrollError) {
                        console.warn('scrollIntoView failed:', scrollError);
                        return false;
                    }
                };

                // Perform the scroll
                if (!performScroll()) {
                    console.warn('All scroll methods failed');
                }

            } else if (attempts < maxAttempts) {
                // No highlights found yet - wait longer
                setTimeout(() => waitForHighlightsAndScroll(attempts + 1), 400);
            } else {
                // No highlights found after waiting - scroll to page
                try {
                    const viewerContainer = pdfViewer.pdfViewer.container;
                    if (viewerContainer && pageView.div) {
                        const containerRect = viewerContainer.getBoundingClientRect();
                        const pageRect = pageView.div.getBoundingClientRect();
                        const scrollTop = viewerContainer.scrollTop + (pageRect.top - containerRect.top);

                        viewerContainer.scrollTo({
                            top: Math.max(0, scrollTop),
                            behavior: 'smooth'
                        });
                    }
                } catch (error) {
                    console.warn('Page scroll failed:', error);
                }
            }

        } catch (error) {
            console.warn('Error in enhanced scroll function:', error);
            if (attempts < maxAttempts) {
                setTimeout(() => waitForHighlightsAndScroll(attempts + 1), 400);
            }
        }
    };

    // Start the enhanced waiting process
    waitForHighlightsAndScroll();
} function handlePDFFrameError(evidenceId) {
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
    const newTabUrl = `/pdfjs/web/viewer.html?file=${encodeURIComponent(filePath)}`;
    window.open(newTabUrl, '_blank');
}

// Make functions globally accessible
window.handlePDFFrameLoad = handlePDFFrameLoad;
window.handlePDFFrameError = handlePDFFrameError;
window.scrollToFirstHighlight = scrollToFirstHighlight;
window.showLocalPDF = showLocalPDF;
window.applyAllHighlights = applyAllHighlights;
window.applyHighlightsViaPDFJS = applyHighlightsViaPDFJS;
window.highlightTextOnPage = highlightTextOnPage;
window.createHighlightOverlay = createHighlightOverlay;
window.repositionHighlights = repositionHighlights;
window.clearHighlights = clearHighlights;
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
    // Don't suppress server response errors - let them show for debugging
    if (message.includes('Unexpected server response') ||
        message.includes('Invalid data path format') ||
        message.includes('ResponseException')) {
        console.error('PDF Server Error:', message);
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
    // Handle server response errors
    if (e.reason && e.reason.message &&
        (e.reason.message.includes('Unexpected server response') ||
            e.reason.message.includes('ResponseException'))) {
        console.error('PDF server error (unhandled promise):', e.reason.message);
        // Don't suppress these - they are important for debugging
        return;
    }
});