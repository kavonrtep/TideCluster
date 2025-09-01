// Karyotype visualization functionality

function showKaryotype(familyId) {
    // Create new window for karyotype visualization
    const karyotypeWindow = window.open('', '_blank', 'width=1200,height=800,scrollbars=yes,resizable=yes');
    
    if (!data.karyotype_data || Object.keys(data.karyotype_data).length === 0) {
        karyotypeWindow.document.write(`
            <html>
                <head><title>Karyotype View - ${familyId}</title></head>
                <body><p>No karyotype data available.</p></body>
            </html>
        `);
        return;
    }
    
    // Generate HTML content for the new window
    const karyotypeHtml = generateKaryotypeHTML(familyId, '', '', false);
    
    karyotypeWindow.document.write(karyotypeHtml);
    karyotypeWindow.document.close();
}

function generateKaryotypeHTML(familyId, chromosomeFilter = '', sampleFilter = '', useRegex = false) {
    const filteredSamples = data.samples.filter(sample => {
        if (sampleFilter === '') return true;
        if (useRegex) {
            try {
                const regex = new RegExp(sampleFilter, 'i');
                return regex.test(sample);
            } catch (e) {
                console.warn('Invalid regex pattern for sample filter:', sampleFilter, e);
                return false;
            }
        } else {
            return sample.toLowerCase().includes(sampleFilter.toLowerCase());
        }
    });
    
    let contentHtml = "";
    
    // Generate karyotype visualization for each filtered sample
    filteredSamples.forEach(sample => {
        const sampleData = data.karyotype_data[sample];
        if (!sampleData) return;
        
        contentHtml += `<div class="sample-karyotype">`;
        contentHtml += `<h4>${sample}</h4>`;
        contentHtml += `<div class="chromosomes-container">`;
        
        // Sort contigs by order from GFF3 file and apply chromosome filter
        const sortedContigs = sampleData.contig_info
            .filter(contigInfo => {
                if (chromosomeFilter === '') return true;
                if (useRegex) {
                    try {
                        const regex = new RegExp(chromosomeFilter, 'i');
                        return regex.test(contigInfo.contig);
                    } catch (e) {
                        console.warn('Invalid regex pattern for chromosome filter:', chromosomeFilter, e);
                        return false;
                    }
                } else {
                    return contigInfo.contig.toLowerCase().includes(chromosomeFilter.toLowerCase());
                }
            })
            .sort((a, b) => a.order - b.order);
        
        if (sortedContigs.length === 0) {
            contentHtml += `<p>No chromosomes match the filter "${chromosomeFilter}"</p>`;
        } else {
            // Get all contig lengths to create scale - use fixed width for better comparison
            const maxLength = Math.max(...sortedContigs.map(c => c.length));
            const baseScale = 1000; // Base scale for visualization
            
            sortedContigs.forEach(contigInfo => {
                const contigName = contigInfo.contig;
                const contigLength = contigInfo.length;
                const contigData = sampleData.contigs[contigName];
                
                if (!contigData) return;
                
                // Calculate proportional width based on actual length
                const chromosomeWidth = Math.max((contigLength / maxLength) * baseScale, 100); // Minimum 100px
                
                contentHtml += `<div class="chromosome-linear">`;
                contentHtml += `<div class="chromosome-header">`;
                contentHtml += `<div class="chromosome-name-linear">${contigName}</div>`;
                contentHtml += `<div class="chromosome-length-linear">${formatLength(contigLength)}</div>`;
                contentHtml += `</div>`;
                contentHtml += `<div class="chromosome-bar-linear" style="width: ${chromosomeWidth}px;">`;
                
                // Add scale markings
                const numTicks = Math.min(Math.floor(chromosomeWidth / 100), 10); // Max 10 ticks
                if (numTicks > 1) {
                    for (let i = 0; i <= numTicks; i++) {
                        const tickPos = (i / numTicks) * chromosomeWidth;
                        const genomicPos = Math.floor((i / numTicks) * contigLength);
                        contentHtml += `<div class="scale-tick" style="left: ${tickPos}px;" title="${formatLength(genomicPos)}"></div>`;
                    }
                }
                
                // Add satellite family positions
                contigData.satellites.forEach(satellite => {
                    if (satellite.family === familyId) {
                        const startPos = (satellite.start / contigLength) * chromosomeWidth;
                        const width = Math.max(((satellite.end - satellite.start) / contigLength) * chromosomeWidth, 3);
                        
                        contentHtml += `<div class="satellite-region-linear" style="left: ${startPos}px; width: ${width}px;" 
                                title="${contigName}: ${satellite.start.toLocaleString()}-${satellite.end.toLocaleString()} bp (${formatLength(satellite.end - satellite.start)})"></div>`;
                    }
                });
                
                contentHtml += `</div>`;
                
                // Add position scale below the bar
                contentHtml += `<div class="position-scale" style="width: ${chromosomeWidth}px;">`;
                contentHtml += `<span class="scale-start">0</span>`;
                contentHtml += `<span class="scale-end">${formatLength(contigLength)}</span>`;
                contentHtml += `</div>`;
                
                contentHtml += `</div>`;
            });
        }
        
        contentHtml += `</div></div>`;
    });
    
    if (contentHtml === "") {
        contentHtml = `<p>No ${familyId} satellite families found in the filtered karyotype data.</p>`;
    }
    
    // Return complete HTML document
    return `
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Karyotype View - ${familyId}</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        .karyotype-header {
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .filters {
            display: flex;
            gap: 20px;
            align-items: center;
            margin-top: 15px;
        }
        .filter-group {
            display: flex;
            flex-direction: column;
            gap: 5px;
        }
        .filter-group label {
            font-weight: bold;
            font-size: 14px;
        }
        .filter-group input {
            padding: 8px;
            border: 1px solid #ddd;
            border-radius: 4px;
            width: 200px;
        }
        .apply-button {
            background: #007bff;
            color: white;
            border: none;
            padding: 10px 20px;
            border-radius: 4px;
            cursor: pointer;
            height: fit-content;
            align-self: end;
        }
        .apply-button:hover {
            background: #0056b3;
        }
        .sample-karyotype {
            background: white;
            margin: 20px 0;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .sample-karyotype h4 {
            margin: 0 0 15px 0;
            color: #333;
            border-bottom: 2px solid #007bff;
            padding-bottom: 5px;
        }
        .chromosome-linear {
            margin: 15px 0;
            padding: 10px;
            border: 1px solid #eee;
            border-radius: 5px;
            background: #fafafa;
        }
        .chromosome-header {
            display: flex;
            justify-content: space-between;
            margin-bottom: 10px;
        }
        .chromosome-name-linear {
            font-weight: bold;
            color: #333;
        }
        .chromosome-length-linear {
            color: #666;
            font-size: 0.9em;
        }
        .chromosome-bar-linear {
            position: relative;
            height: 20px;
            background: linear-gradient(to right, #e8f4fd, #b8e0ff);
            border: 1px solid #007bff;
            border-radius: 10px;
            margin: 10px 0;
        }
        .scale-tick {
            position: absolute;
            top: -5px;
            height: 30px;
            width: 1px;
            background: #666;
        }
        .satellite-region-linear {
            position: absolute;
            top: 0;
            height: 100%;
            background: #ff6b6b;
            border-radius: 3px;
            opacity: 0.8;
        }
        .satellite-region-linear:hover {
            opacity: 1;
            box-shadow: 0 2px 4px rgba(255,107,107,0.3);
        }
        .position-scale {
            display: flex;
            justify-content: space-between;
            font-size: 0.8em;
            color: #666;
        }
    </style>
    <script>
        const data = ${JSON.stringify(data, null, 2)};
        
        function formatLength(length) {
            if (length >= 1000000) {
                return (length / 1000000).toFixed(2) + ' Mbp';
            } else if (length >= 1000) {
                return (length / 1000).toFixed(2) + ' kbp';
            } else {
                return length + ' bp';
            }
        }
        
        
        function applyFilters() {
            const chromosomeFilter = document.getElementById('chromosome-filter').value;
            const sampleFilter = document.getElementById('sample-filter').value;
            const useRegex = document.getElementById('use-regex').checked;
            const container = document.getElementById('karyotype-content');
            
            const filteredHtml = generateFilteredContent('${familyId}', chromosomeFilter, sampleFilter, useRegex);
            container.innerHTML = filteredHtml;
        }
        
        function generateFilteredContent(familyId, chromosomeFilter, sampleFilter, useRegex) {
            const filteredSamples = data.samples.filter(sample => {
                if (sampleFilter === '') return true;
                if (useRegex) {
                    try {
                        const regex = new RegExp(sampleFilter, 'i');
                        return regex.test(sample);
                    } catch (e) {
                        console.warn('Invalid regex pattern for sample filter:', sampleFilter, e);
                        return false;
                    }
                } else {
                    return sample.toLowerCase().includes(sampleFilter.toLowerCase());
                }
            });
            
            let contentHtml = "";
            
            filteredSamples.forEach(sample => {
                const sampleData = data.karyotype_data[sample];
                if (!sampleData) return;
                
                contentHtml += '<div class="sample-karyotype">';
                contentHtml += '<h4>' + sample + '</h4>';
                contentHtml += '<div class="chromosomes-container">';
                
                const sortedContigs = sampleData.contig_info
                    .filter(contigInfo => {
                        if (chromosomeFilter === '') return true;
                        if (useRegex) {
                            try {
                                const regex = new RegExp(chromosomeFilter, 'i');
                                return regex.test(contigInfo.contig);
                            } catch (e) {
                                console.warn('Invalid regex pattern for chromosome filter:', chromosomeFilter, e);
                                return false;
                            }
                        } else {
                            return contigInfo.contig.toLowerCase().includes(chromosomeFilter.toLowerCase());
                        }
                    })
                    .sort((a, b) => a.order - b.order);
                
                if (sortedContigs.length === 0) {
                    contentHtml += '<p>No chromosomes match the filter "' + chromosomeFilter + '"</p>';
                } else {
                    const maxLength = Math.max(...sortedContigs.map(c => c.length));
                    const baseScale = 1000;
                    
                    sortedContigs.forEach(contigInfo => {
                        const contigName = contigInfo.contig;
                        const contigLength = contigInfo.length;
                        const contigData = sampleData.contigs[contigName];
                        
                        if (!contigData) return;
                        
                        const chromosomeWidth = Math.max((contigLength / maxLength) * baseScale, 100);
                        
                        contentHtml += '<div class="chromosome-linear">';
                        contentHtml += '<div class="chromosome-header">';
                        contentHtml += '<div class="chromosome-name-linear">' + contigName + '</div>';
                        contentHtml += '<div class="chromosome-length-linear">' + formatLength(contigLength) + '</div>';
                        contentHtml += '</div>';
                        contentHtml += '<div class="chromosome-bar-linear" style="width: ' + chromosomeWidth + 'px;">';
                        
                        const numTicks = Math.min(Math.floor(chromosomeWidth / 100), 10);
                        if (numTicks > 1) {
                            for (let i = 0; i <= numTicks; i++) {
                                const tickPos = (i / numTicks) * chromosomeWidth;
                                const genomicPos = Math.floor((i / numTicks) * contigLength);
                                contentHtml += '<div class="scale-tick" style="left: ' + tickPos + 'px;" title="' + formatLength(genomicPos) + '"></div>';
                            }
                        }
                        
                        contigData.satellites.forEach(satellite => {
                            if (satellite.family === familyId) {
                                const startPos = (satellite.start / contigLength) * chromosomeWidth;
                                const width = Math.max(((satellite.end - satellite.start) / contigLength) * chromosomeWidth, 3);
                                
                                contentHtml += '<div class="satellite-region-linear" style="left: ' + startPos + 'px; width: ' + width + 'px;" title="' + 
                                    contigName + ': ' + satellite.start.toLocaleString() + '-' + satellite.end.toLocaleString() + 
                                    ' bp (' + formatLength(satellite.end - satellite.start) + ')"></div>';
                            }
                        });
                        
                        contentHtml += '</div>';
                        contentHtml += '<div class="position-scale" style="width: ' + chromosomeWidth + 'px;">';
                        contentHtml += '<span class="scale-start">0</span>';
                        contentHtml += '<span class="scale-end">' + formatLength(contigLength) + '</span>';
                        contentHtml += '</div>';
                        contentHtml += '</div>';
                    });
                }
                
                contentHtml += '</div></div>';
            });
            
            if (contentHtml === "") {
                contentHtml = '<p>No ' + familyId + ' satellite families found in the filtered karyotype data.</p>';
            }
            
            return contentHtml;
        }
    </script>
</head>
<body>
    <div class="karyotype-header">
        <h2>Karyotype View - ${familyId}</h2>
        <div class="filters">
            <div class="filter-group">
                <label for="chromosome-filter">Chromosome/Contig Filter:</label>
                <input type="text" id="chromosome-filter" placeholder="e.g., chr1, chr1$, chr[14]">
            </div>
            <div class="filter-group">
                <label for="sample-filter">Sample Filter:</label>
                <input type="text" id="sample-filter" placeholder="e.g., sample1, wild$, (mut|ctrl)">
            </div>
            <div class="filter-group">
                <label>
                    <input type="checkbox" id="use-regex"> Use Regular Expressions
                </label>
                <small style="color: #666;">Enable for patterns like chr1$, chr[14], (mut|ctrl)</small>
            </div>
            <button class="apply-button" onclick="applyFilters()">Apply Filters</button>
        </div>
    </div>
    <div id="karyotype-content">
        ${contentHtml}
    </div>
    <script>
        // Initialize with all data when page loads
        window.onload = function() {
            const container = document.getElementById('karyotype-content');
            if (container.innerHTML.trim() === '') {
                const filteredHtml = generateFilteredContent('${familyId}', '', '', false);
                container.innerHTML = filteredHtml;
            }
        };
    </script>
</body>
</html>`;
}

