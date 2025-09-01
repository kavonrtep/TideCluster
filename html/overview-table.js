// Overview table functionality

function initOverviewTable() {
    const container = document.getElementById("overview-table-container");
    const stats = data.overview_stats;
    
    let html = `
        <table class="overview-table">
            <thead>
                <tr>
                    <th>Sample Name</th>
                    <th>Number of TRCs</th>
                    <th>Number of Families</th>
                    <th>Number of Unique Families</th>
                    <th>Total Length</th>
                </tr>
            </thead>
            <tbody>
    `;
    
    stats.forEach(row => {
        html += `
            <tr>
                <td class="sample-name">${row.Sample}</td>
                <td class="numeric">${row.Number_of_TRCs.toLocaleString()}</td>
                <td class="numeric">${row.Number_of_Families.toLocaleString()}</td>
                <td class="numeric">${row.Number_of_Unique_Families.toLocaleString()}</td>
                <td class="numeric">${row.Total_Length}</td>
            </tr>
        `;
    });
    
    html += `
            </tbody>
        </table>
    `;
    
    container.innerHTML = html;
}