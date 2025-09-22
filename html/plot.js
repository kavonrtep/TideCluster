// Karyotype Plot functionality for single chromosome visualization
console.log("plot.js loaded successfully");

let selectedFamilies = [];
let maxSelectedFamilies = 10;

// Qualitative color palette for up to 10 families
const plotColors = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", 
    "#ffff33", "#a65628", "#f781bf", "#999999", "#1f78b4"
];

function initPlotTab() {
    console.log("Initializing plot tab");
    console.log("Data structure:", data);
    console.log("Karyotype data available:", data && data.karyotype_data ? "Yes" : "No");
    if (data && data.karyotype_data) {
        console.log("Karyotype data samples:", Object.keys(data.karyotype_data));
    }
    console.log("Detailed families available:", data && data.detailed_families ? data.detailed_families.length : "No");

    populateChromosomeSelect();
    populateFamilyList();
    updateSelectedFamiliesDisplay();
    updatePlotButton();
    addChromosomeSelectListener();
}

function populateChromosomeSelect() {
    const select = document.getElementById("chromosome-select");
    if (!select) {
        console.error("Chromosome select element not found");
        return;
    }

    // Clear existing options except the first one
    while (select.children.length > 1) {
        select.removeChild(select.lastChild);
    }

    const chromosomesWithSize = [];

    // Collect all unique chromosomes with their maximum size across samples
    if (data.karyotype_data) {
        console.log("Found karyotype data, processing chromosomes...");
        const chromosomeMaxSizes = {};

        Object.values(data.karyotype_data).forEach(sampleData => {
            if (sampleData.contig_info) {
                sampleData.contig_info.forEach(contig => {
                    const chrName = contig.contig;
                    const chrLength = contig.length;

                    if (!chromosomeMaxSizes[chrName] || chromosomeMaxSizes[chrName] < chrLength) {
                        chromosomeMaxSizes[chrName] = chrLength;
                    }
                });
            }
        });

        // Convert to array for sorting
        Object.keys(chromosomeMaxSizes).forEach(chrName => {
            chromosomesWithSize.push({
                name: chrName,
                size: chromosomeMaxSizes[chrName]
            });
        });

        console.log("Found", chromosomesWithSize.length, "chromosomes");
    } else {
        console.log("No karyotype data found");
    }

    // Sort chromosomes by size (largest first)
    chromosomesWithSize.sort((a, b) => b.size - a.size);

    chromosomesWithSize.forEach(chr => {
        const option = document.createElement("option");
        option.value = chr.name;
        option.textContent = chr.name;
        select.appendChild(option);
    });

    console.log("Populated chromosome select with", chromosomesWithSize.length, "options");
}

function populateFamilyList() {
    const container = document.getElementById("family-list");
    
    if (!data || !data.detailed_families) {
        container.innerHTML = "<p>No family data available</p>";
        console.error("No detailed_families data found");
        return;
    }
    
    const families = data.detailed_families;
    container.innerHTML = "";
    
    console.log("Populating family list with", families.length, "families");
    
    families.forEach((family, index) => {
        const familyId = `SF_${String(family.family_id).padStart(4, "0")}`;
        const annotation = family.prevalent_annot || "No annotation";
        
        const familyItem = document.createElement("div");
        familyItem.className = "family-item";
        familyItem.dataset.familyId = familyId.toLowerCase();
        familyItem.dataset.annotation = String(annotation || "").toLowerCase();
        
        familyItem.innerHTML = `
            <div class="family-checkbox">
                <input type="checkbox" id="family-${familyId}" 
                       onchange="toggleFamily('${familyId}')" 
                       ${selectedFamilies.includes(familyId) ? 'checked' : ''}>
                <label for="family-${familyId}">
                    <span class="family-id">${familyId}</span>
                    <span class="family-annotation"> - ${annotation}</span>
                </label>
            </div>
        `;
        
        container.appendChild(familyItem);
    });
    
    console.log("Family list populated with", container.children.length, "items");
}

function filterFamilies() {
    const searchTerm = document.getElementById("family-search").value.toLowerCase();
    const familyItems = document.querySelectorAll(".family-item");
    
    console.log("Filtering families with term:", searchTerm);
    console.log("Found", familyItems.length, "family items");
    
    let visibleCount = 0;
    familyItems.forEach(item => {
        const familyId = item.dataset.familyId;
        const annotation = item.dataset.annotation;
        
        const matches = familyId.includes(searchTerm) || annotation.includes(searchTerm);
        item.style.display = matches ? "block" : "none";
        
        if (matches) visibleCount++;
    });
    
    console.log("Showing", visibleCount, "families after filter");
}

function toggleFamily(familyId) {
    const checkbox = document.getElementById(`family-${familyId}`);
    
    if (checkbox.checked) {
        if (selectedFamilies.length >= maxSelectedFamilies) {
            checkbox.checked = false;
            alert(`Maximum ${maxSelectedFamilies} families can be selected.`);
            return;
        }
        selectedFamilies.push(familyId);
    } else {
        const index = selectedFamilies.indexOf(familyId);
        if (index > -1) {
            selectedFamilies.splice(index, 1);
        }
    }
    
    updateSelectedFamiliesDisplay();
    updatePlotButton();
}

function updateSelectedFamiliesDisplay() {
    const container = document.getElementById("selected-families");
    
    if (selectedFamilies.length === 0) {
        container.innerHTML = "<p>No families selected</p>";
        return;
    }
    
    container.innerHTML = selectedFamilies.map((familyId, index) => {
        const color = plotColors[index % plotColors.length];
        return `
            <div class="selected-family" style="border-left: 4px solid ${color};">
                <span>${familyId}</span>
                <button onclick="removeFamily('${familyId}')" class="remove-btn">×</button>
            </div>
        `;
    }).join("");
}

function removeFamily(familyId) {
    const index = selectedFamilies.indexOf(familyId);
    if (index > -1) {
        selectedFamilies.splice(index, 1);
        
        // Uncheck the checkbox
        const checkbox = document.getElementById(`family-${familyId}`);
        if (checkbox) {
            checkbox.checked = false;
        }
        
        updateSelectedFamiliesDisplay();
        updatePlotButton();
    }
}

function clearSelection() {
    selectedFamilies = [];
    
    // Uncheck all checkboxes
    document.querySelectorAll(".family-item input[type='checkbox']").forEach(cb => {
        cb.checked = false;
    });
    
    updateSelectedFamiliesDisplay();
    updatePlotButton();
}

function updatePlotButton() {
    const button = document.getElementById("plot-button");
    const chromosome = document.getElementById("chromosome-select").value;
    
    button.disabled = !chromosome || selectedFamilies.length === 0;
}

// Add event listener for chromosome selection change
function addChromosomeSelectListener() {
    const chromosomeSelect = document.getElementById("chromosome-select");
    if (chromosomeSelect) {
        chromosomeSelect.addEventListener("change", updatePlotButton);
        console.log("Chromosome select listener added");
    } else {
        console.log("Chromosome select element not found");
    }
}

function generatePlot() {
    const chromosome = document.getElementById("chromosome-select").value;
    const container = document.getElementById("plot-container");
    
    if (!chromosome || selectedFamilies.length === 0) {
        container.innerHTML = "<p>Please select a chromosome and at least one family.</p>";
        return;
    }
    
    if (!data.karyotype_data || Object.keys(data.karyotype_data).length === 0) {
        container.innerHTML = "<p>No karyotype data available for plotting.</p>";
        return;
    }
    
    // Generate SVG plot
    const plotHtml = generateKaryotypePlotSVG(chromosome, selectedFamilies);
    container.innerHTML = plotHtml;
    
    // Initialize zoom functionality after DOM is updated
    setTimeout(() => {
        initZoomFunctionality(container.querySelector('svg'));
    }, 100);
}

function generateKaryotypePlotSVG(chromosome, families) {
    const samples = data.samples;
    const svgWidth = 1040; // 30% wider (800 * 1.3)
    const svgHeight = samples.length * 40 + 150; // Space for each sample + legend
    
    // Find maximum chromosome length for scaling
    let maxLength = 0;
    samples.forEach(sample => {
        if (data.karyotype_data[sample] && data.karyotype_data[sample].contig_info) {
            const contigInfo = data.karyotype_data[sample].contig_info.find(c => c.contig === chromosome);
            if (contigInfo) {
                maxLength = Math.max(maxLength, contigInfo.length);
            }
        }
    });
    
    if (maxLength === 0) {
        return '<p>No data found for chromosome ' + chromosome + '</p>';
    }
    
    const chrWidth = 780; // 30% wider (600 * 1.3)
    const chrHeight = 20;
    const sampleSpacing = 35;
    const leftMargin = 150;
    const topMargin = 50;
    
    const plotId = `plot-svg-${Date.now()}`;
    
    let svg = `
        <div class="plot-wrapper">
            <div class="zoom-controls">
                <button onclick="resetZoom('${plotId}')">Reset Zoom</button>
                <span class="zoom-info">Click and drag to select area for zooming</span>
            </div>
            <svg id="${plotId}" width="${svgWidth}" height="${svgHeight}" style="border: 1px solid #ddd; background: white; cursor: grab; user-select: none;">
                <defs>
                    <style>
                        .chromosome-rect { fill: lightgrey; stroke: darkgrey; stroke-width: 1; }
                        .sample-label { font-family: Arial, sans-serif; font-size: 12px; text-anchor: end; }
                        .chr-label { font-family: Arial, sans-serif; font-size: 14px; font-weight: bold; text-anchor: middle; }
                        .length-label { font-family: Arial, sans-serif; font-size: 10px; fill: #666; }
                        .legend-text { font-family: Arial, sans-serif; font-size: 12px; }
                        .legend-title { font-family: Arial, sans-serif; font-size: 14px; font-weight: bold; }
                    </style>
                </defs>
                
                <g id="main-content">
                    <!-- Chromosome header -->
                    <text x="${leftMargin + chrWidth/2}" y="30" class="chr-label">${chromosome}</text>
    `;
    
    // Plot each sample
    samples.forEach((sample, sampleIdx) => {
        const y = topMargin + sampleIdx * sampleSpacing;
        
        // Sample label
        svg += '<text x="' + (leftMargin - 10) + '" y="' + (y + chrHeight/2 + 4) + '" class="sample-label">' + sample + '</text>';
        
        // Check if sample has this chromosome
        if (data.karyotype_data[sample] && data.karyotype_data[sample].contig_info) {
            const contigInfo = data.karyotype_data[sample].contig_info.find(c => c.contig === chromosome);
            
            if (contigInfo) {
                const chrLength = contigInfo.length;
                const proportionalWidth = (chrLength / maxLength) * chrWidth;
                
                // Draw chromosome rectangle
                svg += '<rect x="' + leftMargin + '" y="' + y + '" width="' + proportionalWidth + '" height="' + chrHeight + '" class="chromosome-rect"/>';
                
                // Add length label
                const lengthLabel = formatLength(chrLength);
                svg += '<text x="' + (leftMargin + proportionalWidth + 5) + '" y="' + (y + chrHeight/2 + 4) + '" class="length-label">' + lengthLabel + '</text>';
                
                // Plot satellite families
                if (data.karyotype_data[sample].contigs && data.karyotype_data[sample].contigs[chromosome]) {
                    const satellites = data.karyotype_data[sample].contigs[chromosome].satellites;
                    
                    families.forEach((familyId, familyIdx) => {
                        const color = plotColors[familyIdx % plotColors.length];
                        
                        satellites.forEach(satellite => {
                            if (satellite.family === familyId) {
                                const startPos = (satellite.start / chrLength) * proportionalWidth;
                                const width = Math.max(((satellite.end - satellite.start) / chrLength) * proportionalWidth, 2);
                                
                                svg += '<rect x="' + (leftMargin + startPos) + '" y="' + y + '" width="' + width + '" height="' + chrHeight + '" ' +
                                        'fill="' + color + '" stroke="' + color + '" stroke-width="1" opacity="0.8">' +
                                        '<title>' + familyId + ' in ' + sample + ': ' + satellite.start.toLocaleString() + '-' + satellite.end.toLocaleString() + ' bp</title>' +
                                        '</rect>';
                            }
                        });
                    });
                }
            }
        }
    });
    
    // Add legend
    const legendY = topMargin + samples.length * sampleSpacing + 30;
    svg += '<text x="' + leftMargin + '" y="' + legendY + '" class="legend-title">Satellite Families:</text>';
    
    families.forEach((familyId, index) => {
        const color = plotColors[index % plotColors.length];
        const x = leftMargin + (index % 5) * 120;
        const y = legendY + 20 + Math.floor(index / 5) * 20;

        svg += '<rect x="' + x + '" y="' + (y - 10) + '" width="15" height="15" fill="' + color + '" stroke="' + color + '"/>';
        svg += '<text x="' + (x + 20) + '" y="' + y + '" class="legend-text">' + familyId + '</text>';
    });

    svg += '</g></svg></div>';
    
    return svg;
}

function initZoomFunctionality(svg) {
    if (!svg) {
        console.error('SVG element not found for zoom initialization');
        return;
    }
    
    console.log('Initializing zoom functionality for SVG:', svg.id);
    
    const svgWidth = parseInt(svg.getAttribute('width'));
    const svgHeight = parseInt(svg.getAttribute('height'));
    let viewBox = { x: 0, y: 0, width: svgWidth, height: svgHeight };
    let isSelecting = false;
    let selectionBox = null;
    let startX = 0;
    let startY = 0;
    
    function updateViewBox() {
        svg.setAttribute('viewBox', `${viewBox.x} ${viewBox.y} ${viewBox.width} ${viewBox.height}`);
        console.log('ViewBox updated:', viewBox);
    }
    
    function resetZoom() {
        viewBox = { x: 0, y: 0, width: svgWidth, height: svgHeight };
        updateViewBox();
    }
    
    // Make reset function globally accessible
    window.resetZoom = function(plotId) {
        if (plotId === svg.id) {
            resetZoom();
        }
    };
    
    function createSelectionBox() {
        if (selectionBox) {
            selectionBox.remove();
        }
        selectionBox = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        selectionBox.setAttribute('fill', 'rgba(0, 123, 255, 0.2)');
        selectionBox.setAttribute('stroke', '#007bff');
        selectionBox.setAttribute('stroke-width', '2');
        selectionBox.setAttribute('stroke-dasharray', '5,5');
        selectionBox.setAttribute('pointer-events', 'none');
        svg.appendChild(selectionBox);
        console.log('Selection box created');
    }
    
    function getSVGCoords(e) {
        const rect = svg.getBoundingClientRect();
        const x = ((e.clientX - rect.left) / rect.width) * viewBox.width + viewBox.x;
        const y = ((e.clientY - rect.top) / rect.height) * viewBox.height + viewBox.y;
        return { x, y };
    }
    
    function updateSelectionBox(x1, y1, x2, y2) {
        if (!selectionBox) return;
        const x = Math.min(x1, x2);
        const y = Math.min(y1, y2);
        const width = Math.abs(x2 - x1);
        const height = Math.abs(y2 - y1);
        
        selectionBox.setAttribute('x', x);
        selectionBox.setAttribute('y', y);
        selectionBox.setAttribute('width', width);
        selectionBox.setAttribute('height', height);
    }
    
    function zoomToSelection(x1, y1, x2, y2) {
        const selectionX = Math.min(x1, x2);
        const selectionY = Math.min(y1, y2);
        const selectionWidth = Math.abs(x2 - x1);
        const selectionHeight = Math.abs(y2 - y1);
        
        console.log('Zoom selection:', { selectionX, selectionY, selectionWidth, selectionHeight });
        
        if (selectionWidth < 10 || selectionHeight < 10) {
            console.log('Selection too small, ignoring');
            return;
        }
        
        viewBox.x = selectionX;
        viewBox.y = selectionY;
        viewBox.width = selectionWidth;
        viewBox.height = selectionHeight;
        
        updateViewBox();
    }
    
    svg.addEventListener('mousedown', function(e) {
        console.log('mousedown event triggered on SVG');
        e.preventDefault();
        e.stopPropagation();
        const coords = getSVGCoords(e);
        console.log('start coords:', coords);
        startX = coords.x;
        startY = coords.y;
        isSelecting = true;
        createSelectionBox();
        svg.style.cursor = 'crosshair';
        document.body.style.userSelect = 'none';
    });
    
    svg.addEventListener('mousemove', function(e) {
        if (isSelecting) {
            console.log('mousemove during selection');
            e.preventDefault();
            const coords = getSVGCoords(e);
            updateSelectionBox(startX, startY, coords.x, coords.y);
        }
    });
    
    svg.addEventListener('mouseup', function(e) {
        if (isSelecting) {
            console.log('mouseup event triggered on SVG');
            e.preventDefault();
            const coords = getSVGCoords(e);
            console.log('end coords:', coords);
            
            zoomToSelection(startX, startY, coords.x, coords.y);
            
            if (selectionBox) {
                selectionBox.remove();
                selectionBox = null;
            }
            
            isSelecting = false;
            svg.style.cursor = 'grab';
            document.body.style.userSelect = '';
        }
    });
    
    document.addEventListener('mouseup', function() {
        if (isSelecting) {
            console.log('document mouseup - cleaning up selection');
            if (selectionBox) {
                selectionBox.remove();
                selectionBox = null;
            }
            isSelecting = false;
            svg.style.cursor = 'grab';
            document.body.style.userSelect = '';
        }
    });
    
    console.log('Zoom event listeners attached to SVG:', svg.id);
}

// Mark plot.js as loaded
window.plotJsLoaded = true;