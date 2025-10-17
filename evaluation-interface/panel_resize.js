// Panel resizing functionality
document.addEventListener('DOMContentLoaded', function () {
    let isResizing = false;
    let currentResizer = null;
    let startX = 0;
    let startWidth = 0;

    // Initialize resizers
    const leftResizer = document.getElementById('leftResizer');
    const rightResizer = document.getElementById('rightResizer');
    const leftPanel = document.getElementById('leftPanel');
    const rightPanel = document.getElementById('rightPanel');
    const container = document.querySelector('.container');

    // Update resize handle visibility
    function updateResizeHandleVisibility() {
        if (leftPanel.classList.contains('expanded')) {
            leftResizer.style.display = 'block';
            // Position the left resizer at the right edge of left panel
            const leftWidth = leftPanel.offsetWidth || 500;
            leftResizer.style.left = leftWidth + 'px';
        } else {
            leftResizer.style.display = 'none';
        }

        if (rightPanel.classList.contains('expanded')) {
            container.classList.add('right-panel-expanded');
            rightResizer.style.display = 'block';
            // Position the right resizer at the left edge of right panel
            // Use a timeout to ensure the panel has finished expanding
            setTimeout(() => {
                const rightWidth = rightPanel.offsetWidth;
                if (rightWidth > 0) {
                    rightResizer.style.right = rightWidth + 'px';
                } else {
                    // Fallback to CSS default width
                    rightResizer.style.right = '400px';
                }
            }, 300);  // Increased timeout for better reliability
        } else {
            container.classList.remove('right-panel-expanded');
            rightResizer.style.display = 'none';
        }
    }

    // Watch for panel changes
    const observer = new MutationObserver((mutations) => {
        updateResizeHandleVisibility();
    });
    if (leftPanel) observer.observe(leftPanel, { attributes: true, attributeFilter: ['class'] });
    if (rightPanel) observer.observe(rightPanel, { attributes: true, attributeFilter: ['class'] });

    // Initial visibility check
    updateResizeHandleVisibility();

    if (leftResizer) {
        leftResizer.addEventListener('mousedown', initResize);
    }
    if (rightResizer) {
        rightResizer.addEventListener('mousedown', initResize);
    }

    function initResize(e) {
        isResizing = true;
        currentResizer = e.target;
        startX = e.clientX;

        if (currentResizer.id === 'leftResizer') {
            // Left resizer controls left panel
            startWidth = leftPanel.offsetWidth;
        } else if (currentResizer.id === 'rightResizer') {
            // Right resizer controls right panel
            startWidth = rightPanel.offsetWidth;
        }

        currentResizer.classList.add('dragging');
        document.addEventListener('mousemove', doResize);
        document.addEventListener('mouseup', stopResize);
        // Add window listeners as backup to catch mouse events outside the document
        window.addEventListener('mousemove', doResize);
        window.addEventListener('mouseup', stopResize);
        e.preventDefault();
    }

    function doResize(e) {
        if (!isResizing || !currentResizer) return;

        const dx = e.clientX - startX;
        let newWidth;

        if (currentResizer.id === 'leftResizer') {
            // Left resizer controls left panel width
            newWidth = Math.max(200, Math.min(800, startWidth + dx));
            leftPanel.style.width = newWidth + 'px';
            // Update resizer position to follow the new panel edge
            currentResizer.style.left = newWidth + 'px';
        } else if (currentResizer.id === 'rightResizer') {
            // Right resizer controls right panel width
            newWidth = Math.max(200, Math.min(800, startWidth - dx));
            rightPanel.style.width = newWidth + 'px';
            // Update resizer position to follow the new panel edge
            currentResizer.style.right = newWidth + 'px';
        }

        // Prevent text selection during drag
        e.preventDefault();
    }

    function stopResize(e) {
        if (!isResizing) return;

        isResizing = false;
        if (currentResizer) {
            currentResizer.classList.remove('dragging');
        }
        currentResizer = null;

        // Remove event listeners from both document and window for better reliability
        document.removeEventListener('mousemove', doResize);
        document.removeEventListener('mouseup', stopResize);
        window.removeEventListener('mousemove', doResize);
        window.removeEventListener('mouseup', stopResize);

        // Prevent any remaining event bubbling
        if (e) {
            e.preventDefault();
            e.stopPropagation();
        }
    }
});