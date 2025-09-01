// Detailed families matrix functionality

function initDetailedMatrix() {
    const container = document.getElementById("detailed-matrix-container");
    const families = data.detailed_families;
    const samples = data.samples;
    
    // Get selected view mode
    const viewMode = document.querySelector('input[name="view-mode"]:checked').value;
    const isLengthView = viewMode === "length";
    const isArrayView = viewMode === "array-count";
    
    // Get all values for color scaling
    const allValues = [];
    families.forEach(family => {
        samples.forEach(sample => {
            let key;
            if (isLengthView) {
                key = `${sample}_length`;
            } else if (isArrayView) {
                key = `${sample}_array_count`;
            } else {
                key = `${sample}_trc_count`;
            }
            allValues.push(family[key] || 0);
        });
    });
    
    const colorScale = getColorScale(allValues.filter(v => v > 0), "Reds");
    
    let html = `
        <div class="matrix-container">
            <table class="detailed-matrix">
                <thead>
                    <tr>
                        <th class="corner-cell">Family ID</th>
    `;
    
    // Sample headers
    samples.forEach(sample => {
        html += `<th class="sample-header">${sample}</th>`;
    });
    
    html += `<th class="annotation-header">Prevalent Annotation</th>`;
    html += `</tr></thead><tbody>`;
    
    // Family rows
    families.forEach(family => {
        html += `<tr>`;
        const familyIdStr = `SF_${String(family.family_id).padStart(4, "0")}`;
        const hasKaryotypeData = data.karyotype_data && Object.keys(data.karyotype_data).length > 0;
        
        if (hasKaryotypeData) {
            html += `<td class="family-id clickable" onclick="showKaryotype('${familyIdStr}')">${familyIdStr}</td>`;
        } else {
            html += `<td class="family-id">${familyIdStr}</td>`;
        }
        
        samples.forEach(sample => {
            const trcKey = `${sample}_trc_count`;
            const lengthKey = `${sample}_length`;
            const arrayKey = `${sample}_array_count`;
            const trcCount = family[trcKey] || 0;
            const length = family[lengthKey] || 0;
            const arrayCount = family[arrayKey] || 0;
            
            let displayValue, colorValue;
            if (isLengthView) {
                displayValue = formatLength(length);
                colorValue = length;
            } else if (isArrayView) {
                displayValue = arrayCount;
                colorValue = arrayCount;
            } else {
                displayValue = trcCount;
                colorValue = trcCount;
            }
            const bgColor = colorValue > 0 ? colorScale(colorValue) : "#f9f9f9";
            
            // Get contigs/chromosomes for this family in this sample
            const familyIdStr = `SF_${String(family.family_id).padStart(4, "0")}`;
            let contigs = [];
            if (data.karyotype_data && data.karyotype_data[sample]) {
                const sampleData = data.karyotype_data[sample];
                if (sampleData.contigs) {
                    Object.keys(sampleData.contigs).forEach(contigName => {
                        const contigData = sampleData.contigs[contigName];
                        if (contigData.satellites && contigData.satellites.some(sat => sat.family === familyIdStr)) {
                            contigs.push(contigName);
                        }
                    });
                }
            }
            
            const contigsList = contigs.length > 0 ? contigs.join(', ') : 'None';
            const tooltipText = `Family ${family.family_id} in ${sample}<br>` +
                              `TRCs: ${trcCount}<br>` +
                              `Arrays: ${arrayCount}<br>` +
                              `Length: ${formatLength(length)}<br>` +
                              `Contigs: ${contigsList}`;
            
            html += `
                <td class="matrix-cell ${colorValue > 0 ? 'has-value' : 'empty'}" 
                    style="background-color: ${bgColor}"
                    onmouseover="showTooltip(event, '${tooltipText}')"
                    onmouseout="hideTooltip()">
                    ${colorValue > 0 ? displayValue : ""}
                </td>
            `;
        });
        
        // Prevalent annotation column
        const annotation = family.prevalent_annot || "";
        html += `<td class="annotation-cell">${annotation}</td>`;
        html += `</tr>`;
    });
    
    html += `
                </tbody>
            </table>
        </div>
        <div class="matrix-legend">
            <p><strong>Interpretation:</strong></p>
            <ul>
                <li>Rows are satellite families sorted by family index</li>
                <li>Columns are samples ordered by hierarchical clustering</li>
                <li>Cells show ${isLengthView ? "total genomic length" : isArrayView ? "number of annotated arrays" : "number of TRCs"} for each family-sample combination</li>
                <li>Empty cells indicate the family is absent in that sample</li>
                <li>Darker colors indicate ${isLengthView ? "longer total length" : isArrayView ? "more arrays" : "more TRCs"}</li>
            </ul>
        </div>
    `;
    
    container.innerHTML = html;
}