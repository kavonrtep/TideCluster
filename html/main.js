// Main JavaScript file for TRC Comparative Analysis Report

// Tab switching functionality
function showTab(tabName) {
    // Hide all tab contents
    const tabContents = document.querySelectorAll(".tab-content");
    tabContents.forEach(content => {
        content.classList.remove("active");
    });
    
    // Remove active class from all tab buttons
    const tabButtons = document.querySelectorAll(".tab-button");
    tabButtons.forEach(button => {
        button.classList.remove("active");
    });
    
    // Show selected tab content
    document.getElementById(tabName).classList.add("active");
    
    // Add active class to clicked button
    event.target.classList.add("active");
    
    // Initialize the appropriate visualization
    switch(tabName) {
        case "overview":
            if (typeof initOverviewTable === 'function') {
                initOverviewTable();
            } else {
                console.error('initOverviewTable function not found');
            }
            break;
        case "shared-families":
            if (typeof initSharedMatrix === 'function') {
                initSharedMatrix();
            } else {
                console.error('initSharedMatrix function not found');
            }
            break;
        case "detailed-families":
            if (typeof initDetailedMatrix === 'function') {
                initDetailedMatrix();
            } else {
                console.error('initDetailedMatrix function not found');
            }
            break;
        case "plot":
            if (typeof initPlotTab === 'function') {
                initPlotTab();
            } else {
                console.error('initPlotTab function not found - plot.js may not be loaded');
                console.log('plotJsLoaded flag:', window.plotJsLoaded);
                console.log('Available functions:', Object.keys(window).filter(k => k.includes('init')));

                const plotContainer = document.getElementById('plot-container');
                if (plotContainer) {
                    plotContainer.innerHTML =
                        '<div style="color: red; padding: 20px; border: 1px solid red; margin: 20px;">' +
                        '<h3>Plot functionality not available</h3>' +
                        '<p>The plot.js file may not be loaded correctly.</p>' +
                        '<p>Please check the browser console for more details.</p>' +
                        '</div>';
                }
            }
            break;
    }
}

// Initialize the report
document.addEventListener("DOMContentLoaded", function() {
    console.log("DOM loaded, checking available functions:");
    console.log("initOverviewTable:", typeof initOverviewTable);
    console.log("initSharedMatrix:", typeof initSharedMatrix);
    console.log("initDetailedMatrix:", typeof initDetailedMatrix);
    console.log("initPlotTab:", typeof initPlotTab);
    console.log("plotJsLoaded flag:", window.plotJsLoaded);
    console.log("data object:", typeof data);

    if (typeof initOverviewTable === 'function') {
        initOverviewTable();
    } else {
        console.error('initOverviewTable not available');
    }

    // Add event listeners for view mode radio buttons
    const viewModeInputs = document.querySelectorAll('input[name="view-mode"]');
    viewModeInputs.forEach(input => {
        input.addEventListener("change", function() {
            if (document.getElementById("detailed-families").classList.contains("active")) {
                if (typeof updateDetailedMatrix === 'function') {
                    updateDetailedMatrix();
                } else {
                    console.error('updateDetailedMatrix not available');
                }
            }
        });
    });
});