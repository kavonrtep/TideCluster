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
            initOverviewTable();
            break;
        case "shared-families":
            initSharedMatrix();
            break;
        case "detailed-families":
            initDetailedMatrix();
            break;
        case "plot":
            initPlotTab();
            break;
    }
}

// Initialize the report
document.addEventListener("DOMContentLoaded", function() {
    initOverviewTable();
    
    // Add event listeners for view mode radio buttons
    const viewModeInputs = document.querySelectorAll('input[name="view-mode"]');
    viewModeInputs.forEach(input => {
        input.addEventListener("change", function() {
            if (document.getElementById("detailed-families").classList.contains("active")) {
                initDetailedMatrix();
            }
        });
    });
});