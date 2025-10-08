// Detailed families matrix functionality

let detailedMatrixData = {
    families: [],
    samples: [],
    colorScale: null,
    viewMode: 'trc-count',
    currentPage: 0,
    rowsPerPage: 100,
    filteredFamilies: [],
    searchTerm: ''
};

function initDetailedMatrix() {
    const container = document.getElementById("detailed-matrix-container");
    detailedMatrixData.families = data.detailed_families;
    detailedMatrixData.samples = data.samples;
    detailedMatrixData.filteredFamilies = [...detailedMatrixData.families];

    updateDetailedMatrix();
}

function updateDetailedMatrix() {
    const container = document.getElementById("detailed-matrix-container");
    const families = detailedMatrixData.filteredFamilies;
    const samples = detailedMatrixData.samples;

    // Get selected view mode
    const viewModeElement = document.querySelector('input[name="view-mode"]:checked');
    const viewMode = viewModeElement ? viewModeElement.value : 'trc-count';
    detailedMatrixData.viewMode = viewMode;

    const isLengthView = viewMode === "length";
    const isArrayView = viewMode === "array-count";

    // Calculate color scale for current page only to improve performance
    const startIdx = detailedMatrixData.currentPage * detailedMatrixData.rowsPerPage;
    const endIdx = Math.min(startIdx + detailedMatrixData.rowsPerPage, families.length);
    const currentPageFamilies = families.slice(startIdx, endIdx);

    const allValues = [];
    currentPageFamilies.forEach(family => {
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
    detailedMatrixData.colorScale = colorScale;

    // Calculate pagination info
    const totalPages = Math.ceil(families.length / detailedMatrixData.rowsPerPage);
    const startRow = startIdx + 1;
    const endRow = Math.min(endIdx, families.length);

    // Dynamic column width based on number of samples
    let columnWidth = 35;
    if (samples.length > 20) {
        columnWidth = Math.max(20, Math.min(35, 1000 / samples.length));
    }

    let html = `
        <div class="detailed-matrix-controls">
            <div class="pagination-info">
                Showing families ${startRow}-${endRow} of ${families.length}
                ${detailedMatrixData.searchTerm ? `(filtered from ${detailedMatrixData.families.length} total)` : ''}
                (Page ${detailedMatrixData.currentPage + 1} of ${totalPages})
            </div>
            <div class="pagination-controls">
                <input type="number" id="rows-per-page" value="${detailedMatrixData.rowsPerPage}" min="10" max="500" step="10">
                <label for="rows-per-page">rows per page</label>
                <button onclick="changeDetailedPage(-1)" ${detailedMatrixData.currentPage === 0 ? 'disabled' : ''}>Previous</button>
                <button onclick="changeDetailedPage(1)" ${detailedMatrixData.currentPage >= totalPages - 1 ? 'disabled' : ''}>Next</button>
                <select id="page-select" onchange="goToDetailedPage(this.value)">
    `;

    for (let i = 0; i < totalPages; i++) {
        html += `<option value="${i}" ${i === detailedMatrixData.currentPage ? 'selected' : ''}>Page ${i + 1}</option>`;
    }

    html += `
                </select>
            </div>
        </div>
        <div class="matrix-container">
            <table class="detailed-matrix" style="table-layout: fixed;">
                <thead>
                    <tr>
                        <th class="corner-cell" style="width: 60px;">Family ID</th>
    `;

    // Sample headers with dynamic width
    samples.forEach(sample => {
        html += `<th class="sample-header" style="width: ${columnWidth}px; max-width: ${columnWidth}px;">${sample}</th>`;
    });

    html += `<th class="annotation-header" style="width: 200px;">Prevalent Annotation</th>`;
    html += `</tr></thead><tbody>`;

    // Family rows - only current page
    currentPageFamilies.forEach(family => {
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

            // Simplified tooltip to improve performance
            const tooltipText = `Family ${family.family_id} in ${sample}: ${displayValue}`;

            html += `
                <td class="matrix-cell ${colorValue > 0 ? 'has-value' : 'empty'}"
                    style="background-color: ${bgColor}; width: ${columnWidth}px; max-width: ${columnWidth}px;"
                    onmouseover="showTooltip(event, '${tooltipText}')"
                    onmouseout="hideTooltip()">
                    ${colorValue > 0 ? displayValue : ""}
                </td>
            `;
        });

        // Prevalent annotation column
        const annotation = family.prevalent_annot || "";
        html += `<td class="annotation-cell" style="width: 200px;">${annotation}</td>`;
        html += `</tr>`;
    });

    html += `
                </tbody>
            </table>
        </div>
        <div class="matrix-legend">
            <p><strong>Interpretation:</strong></p>
            <ul>
                <li>Rows are satellite families sorted by family index (showing ${startRow}-${endRow} of ${families.length}${detailedMatrixData.searchTerm ? ` filtered` : ''})</li>
                <li>Columns are samples ordered by hierarchical clustering</li>
                <li>Cells show ${isLengthView ? "total genomic length" : isArrayView ? "number of annotated arrays" : "number of TRCs"} for each family-sample combination</li>
                <li>Empty cells indicate the family is absent in that sample</li>
                <li>Darker colors indicate ${isLengthView ? "longer total length" : isArrayView ? "more arrays" : "more TRCs"}</li>
                ${detailedMatrixData.searchTerm ? `<li><strong>Filter active:</strong> Showing families matching "${detailedMatrixData.searchTerm}"</li>` : ''}
            </ul>
        </div>
    `;

    container.innerHTML = html;

    // Update rows per page listener
    const rowsPerPageInput = document.getElementById('rows-per-page');
    if (rowsPerPageInput) {
        rowsPerPageInput.addEventListener('change', function() {
            const newRowsPerPage = parseInt(this.value);
            if (newRowsPerPage >= 10 && newRowsPerPage <= 500) {
                detailedMatrixData.rowsPerPage = newRowsPerPage;
                detailedMatrixData.currentPage = 0;
                updateDetailedMatrix();
            }
        });
    }

    // Restore search input value (this is in the matrix container, not affected by updateDetailedMatrix)
    setTimeout(() => {
        const searchInput = document.getElementById('detailed-family-search');
        if (searchInput && detailedMatrixData.searchTerm) {
            searchInput.value = detailedMatrixData.searchTerm;
        }
    }, 10);
}

function changeDetailedPage(direction) {
    const totalPages = Math.ceil(detailedMatrixData.filteredFamilies.length / detailedMatrixData.rowsPerPage);
    const newPage = detailedMatrixData.currentPage + direction;

    if (newPage >= 0 && newPage < totalPages) {
        detailedMatrixData.currentPage = newPage;
        updateDetailedMatrix();
    }
}

function goToDetailedPage(pageIndex) {
    const pageNum = parseInt(pageIndex);
    const totalPages = Math.ceil(detailedMatrixData.filteredFamilies.length / detailedMatrixData.rowsPerPage);

    if (pageNum >= 0 && pageNum < totalPages) {
        detailedMatrixData.currentPage = pageNum;
        updateDetailedMatrix();
    }
}

// Filter families based on search term
function filterDetailedFamilies() {
    const searchInput = document.getElementById('detailed-family-search');
    if (!searchInput) {
        console.error('Search input not found');
        return;
    }

    const searchTerm = searchInput.value.toLowerCase();
    detailedMatrixData.searchTerm = searchTerm;

    console.log('Filtering detailed families with term:', searchTerm);
    console.log('Total families to filter:', detailedMatrixData.families.length);

    if (!detailedMatrixData.families || detailedMatrixData.families.length === 0) {
        console.error('No families data available for filtering');
        return;
    }

    if (searchTerm === '') {
        // No filter, show all families
        detailedMatrixData.filteredFamilies = [...detailedMatrixData.families];
    } else {
        // Filter families based on family ID and annotation
        detailedMatrixData.filteredFamilies = detailedMatrixData.families.filter((family, index) => {
            try {
                const familyIdStr = `SF_${String(family.family_id).padStart(4, "0")}`;
                const annotation = String(family.prevalent_annot || '');

                const matches = familyIdStr.toLowerCase().includes(searchTerm) ||
                               annotation.toLowerCase().includes(searchTerm);

                if (index < 3) { // Debug first few families
                    console.log(`Family ${index}:`, {
                        family_id: family.family_id,
                        familyIdStr,
                        prevalent_annot: family.prevalent_annot,
                        annotation,
                        matches
                    });
                }

                return matches;
            } catch (error) {
                console.error('Error filtering family at index', index, ':', error, family);
                return false;
            }
        });
    }

    // Reset to first page when filtering
    detailedMatrixData.currentPage = 0;

    console.log('Filtered families:', detailedMatrixData.filteredFamilies.length, 'of', detailedMatrixData.families.length);

    // Update the matrix display
    updateDetailedMatrix();
}

// Clear the filter
function clearDetailedFilter() {
    const searchInput = document.getElementById('detailed-family-search');
    if (searchInput) {
        searchInput.value = '';
    }

    detailedMatrixData.searchTerm = '';
    detailedMatrixData.filteredFamilies = [...detailedMatrixData.families];
    detailedMatrixData.currentPage = 0;

    console.log('Cleared detailed families filter');

    updateDetailedMatrix();
}