// Utility functions for TRC Comparative Analysis Report

// Format numbers with appropriate units
function formatLength(length) {
    if (length >= 1000000) {
        return (length / 1000000).toFixed(2) + " Mbp";
    } else if (length >= 1000) {
        return (length / 1000).toFixed(2) + " kbp";
    } else {
        return length + " bp";
    }
}

// Create color scale for heatmaps
function getColorScale(values, colorScheme = "Blues") {
    const minVal = Math.min(...values);
    const maxVal = Math.max(...values);
    
    if (minVal === maxVal) {
        return (val) => "#f0f0f0";
    }
    
    return function(val) {
        const intensity = (val - minVal) / (maxVal - minVal);
        
        if (colorScheme === "Blues") {
            const blue = Math.floor(255 - intensity * 200);
            return `rgb(${blue}, ${blue}, 255)`;
        } else if (colorScheme === "Reds") {
            const red = Math.floor(255 - intensity * 200);
            return `rgb(255, ${red}, ${red})`;
        } else {
            // Default grayscale
            const gray = Math.floor(255 - intensity * 200);
            return `rgb(${gray}, ${gray}, ${gray})`;
        }
    };
}

// Create tooltip
function createTooltip() {
    const tooltip = document.createElement("div");
    tooltip.id = "tooltip";
    tooltip.style.cssText = `
        position: absolute;
        background: rgba(0, 0, 0, 0.8);
        color: white;
        padding: 8px;
        border-radius: 4px;
        font-size: 12px;
        pointer-events: none;
        z-index: 1000;
        display: none;
    `;
    document.body.appendChild(tooltip);
    return tooltip;
}

// Show tooltip
function showTooltip(event, text) {
    const tooltip = document.getElementById("tooltip") || createTooltip();
    tooltip.innerHTML = text;
    tooltip.style.display = "block";
    tooltip.style.left = (event.pageX + 10) + "px";
    tooltip.style.top = (event.pageY - 10) + "px";
}

// Hide tooltip
function hideTooltip() {
    const tooltip = document.getElementById("tooltip");
    if (tooltip) {
        tooltip.style.display = "none";
    }
}