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

// Sizing for the sample matrices (Shared / Detailed families), scaled to the
// number of samples: readable text with few samples, the original compact
// 7-8px layout once there are many. Linear between FEW and MANY samples.
function matrixDensity(nSamples) {
    const FEW = 8, MANY = 48;
    const t = Math.min(1, Math.max(0, (nSamples - FEW) / (MANY - FEW)));
    const lerp = (a, b) => Math.round(a + (b - a) * t);
    // Beyond MANY keep the original 1000/n rule; lerp(64, 21) meets it there.
    const col = nSamples > MANY ? Math.max(20, 1000 / nSamples) : lerp(64, 21);
    return {
        font: lerp(14, 8),          // matrix cells, row headers
        headerFont: lerp(14, 7),    // vertical column headers, family IDs
        col: Math.round(col),       // sample column width (px)
        rowHeader: lerp(140, 80),   // row-header column width (px)
        familyCol: lerp(90, 60),    // Family ID column width (px)
    };
}

// Inline style that hands matrixDensity() to the CSS variables in styles.css.
function matrixDensityStyle(d) {
    return `--mx-font:${d.font}px;--mx-header-font:${d.headerFont}px;` +
           `--mx-col:${d.col}px;--mx-row-header:${d.rowHeader}px;` +
           `--mx-family-col:${d.familyCol}px;`;
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
        font-size: 14px;
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