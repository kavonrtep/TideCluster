// Shared families matrix functionality

function initSharedMatrix() {
    const container = document.getElementById("shared-matrix-container");
    const matrix = data.shared_matrix;
    const samples = data.samples;
    
    // Get all values for color scaling
    const allValues = [];
    for (let i = 0; i < samples.length; i++) {
        for (let j = 0; j < samples.length; j++) {
            allValues.push(matrix[i][j]);
        }
    }
    
    const colorScale = getColorScale(allValues, "Blues");
    
    let html = `
        <div class="matrix-container">
            <table class="shared-matrix">
                <thead>
                    <tr>
                        <th class="corner-cell"></th>
    `;
    
    // Column headers
    samples.forEach(sample => {
        html += `<th class="sample-header">${sample}</th>`;
    });
    
    html += `
                    </tr>
                </thead>
                <tbody>
    `;
    
    // Matrix rows
    for (let i = 0; i < samples.length; i++) {
        html += `<tr>`;
        html += `<th class="sample-header row-header">${samples[i]}</th>`;
        
        for (let j = 0; j < samples.length; j++) {
            const value = matrix[i][j];
            const bgColor = colorScale(value);
            const isDiagonal = i === j;
            
            html += `
                <td class="matrix-cell ${isDiagonal ? 'diagonal' : ''}" 
                    style="background-color: ${bgColor}"
                    onmouseover="showTooltip(event, '${samples[i]} vs ${samples[j]}: ${value} shared families')"
                    onmouseout="hideTooltip()">
                    ${value}
                </td>
            `;
        }
        
        html += `</tr>`;
    }
    
    html += `
                </tbody>
            </table>
        </div>
        <div class="matrix-legend">
            <p><strong>Interpretation:</strong></p>
            <ul>
                <li>Rows and columns are samples ordered by hierarchical clustering</li>
                <li>Diagonal cells show total families present in each sample</li>
                <li>Off-diagonal cells show shared families between sample pairs</li>
                <li>Darker colors indicate higher numbers of shared families</li>
            </ul>
        </div>
    `;
    
    container.innerHTML = html;
}