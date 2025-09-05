// TideCluster Visualizer - Main Application
class TideClusterViz {
    constructor() {
        this.data = null;
        this.scales = { x: null };
        this.selectedAnnotations = [];
        this.colorMap = new Map();
        this.availableColors = [];
        this.filteredAnnotations = [];
        this.allAnnotations = [];
        
        // Selection state
        this.selectedRegion = null; // { startBp, endBp }
        this.isSelecting = false;
        this.selectionStart = 0;
        this.overviewWidth = 0;
        this.overviewHeight = 180;
        this.detailWidth = 0;
        this.detailHeight = 300;
        
        // DOM elements
        this.elements = {};
        
        this.init();
    }
    
    async init() {
        try {
            await this.loadData();
            this.setupDOM();
            this.setupEventListeners();
            this.processAnnotations();
            this.setupScales();
            this.render();
            this.applyInitialSelection();
        } catch (error) {
            console.error('Failed to initialize visualizer:', error);
            document.body.innerHTML = '<div style="padding: 20px; color: red;">Error loading visualization: ' + error.message + '</div>';
        }
    }
    
    async loadData() {
        if (window.TIDECLUSTER_DATA) {
            this.data = window.TIDECLUSTER_DATA;
            console.log('Loaded data:', this.data.features.length, 'features');
        } else {
            throw new Error('No data found. Data should be embedded in the HTML.');
        }
    }
    
    setupDOM() {
        // Cache DOM elements
        this.elements = {
            title: document.getElementById('title'),
            search: document.getElementById('annotation-search'),
            allList: document.getElementById('all-annotations'),
            selectedList: document.getElementById('selected-annotations'),
            allCount: document.getElementById('all-count'),
            selectedCount: document.getElementById('selected-count'),
            addBtn: document.getElementById('add-btn'),
            removeBtn: document.getElementById('remove-btn'),
            resetBtn: document.getElementById('reset-btn'),
            selectTopBtn: document.getElementById('select-top-btn'),
            legendItems: document.getElementById('legend-items'),
            
            // Overview elements
            overviewSvg: document.getElementById('overview-svg'),
            overviewGenomeTrack: document.getElementById('overview-genome-track'),
            overviewFeaturesLayer: document.getElementById('overview-features-layer'),
            overviewLabelsLayer: document.getElementById('overview-labels-layer'),
            zoomRegionHighlight: document.getElementById('zoom-region-highlight'),
            
            // Detail elements
            detailSvg: document.getElementById('detail-svg'),
            detailGenomeTrack: document.getElementById('detail-genome-track'),
            detailFeaturesLayer: document.getElementById('detail-features-layer'),
            detailLabelsLayer: document.getElementById('detail-labels-layer'),
            detailRegionInfo: document.getElementById('detail-region-info'),
            
            tooltip: document.getElementById('tooltip'),
            zoomReset: document.getElementById('zoom-reset')
        };
        
        // Set title
        this.elements.title.textContent = this.data.title;
        
        // Setup color palette
        this.availableColors = [...this.data.palette];
    }
    
    setupEventListeners() {
        // Search functionality
        this.elements.search.addEventListener('input', (e) => {
            this.filterAnnotations(e.target.value);
        });
        
        // List controls
        this.elements.addBtn.addEventListener('click', () => {
            this.addSelectedAnnotation();
        });
        
        this.elements.removeBtn.addEventListener('click', () => {
            this.removeSelectedAnnotation();
        });
        
        // Action buttons
        this.elements.resetBtn.addEventListener('click', () => {
            this.resetSelection();
        });
        
        this.elements.selectTopBtn.addEventListener('click', () => {
            this.selectTopN(5);
        });
        
        // Reset selection
        this.elements.zoomReset.addEventListener('click', () => {
            this.resetSelection();
        });
        
        // Overview SVG interactions for selection
        this.setupOverviewInteractions();
    }
    
    setupOverviewInteractions() {
        const overviewSvg = this.elements.overviewSvg;
        
        // Mouse selection on overview
        overviewSvg.addEventListener('mousedown', (e) => {
            this.startSelection(e);
        });
        
        overviewSvg.addEventListener('mousemove', (e) => {
            this.updateSelection(e);
            this.updateTooltip(e);
        });
        
        overviewSvg.addEventListener('mouseup', (e) => {
            this.endSelection(e);
        });
        
        overviewSvg.addEventListener('mouseleave', () => {
            this.cancelSelection();
            this.hideTooltip();
        });
        
        // Tooltip on detail view
        const detailSvg = this.elements.detailSvg;
        detailSvg.addEventListener('mousemove', (e) => {
            this.updateTooltip(e);
        });
        
        detailSvg.addEventListener('mouseleave', () => {
            this.hideTooltip();
        });
        
        // Prevent context menu
        overviewSvg.addEventListener('contextmenu', (e) => e.preventDefault());
        detailSvg.addEventListener('contextmenu', (e) => e.preventDefault());
    }
    
    processAnnotations() {
        // Count annotations
        const annotationCounts = new Map();
        
        this.data.features.forEach(feature => {
            const annotation = feature.annotation_main === null ? 'NA' : feature.annotation_main;
            annotationCounts.set(annotation, (annotationCounts.get(annotation) || 0) + 1);
        });
        
        // Sort by count (descending) then alphabetically
        this.allAnnotations = Array.from(annotationCounts.entries())
            .map(([annotation, count]) => ({ annotation, count }))
            .sort((a, b) => {
                if (b.count !== a.count) return b.count - a.count;
                return a.annotation.localeCompare(b.annotation);
            });
        
        this.filteredAnnotations = [...this.allAnnotations];
        this.updateAnnotationLists();
    }
    
    filterAnnotations(searchTerm) {
        const term = searchTerm.toLowerCase();
        this.filteredAnnotations = this.allAnnotations.filter(item => 
            item.annotation.toLowerCase().includes(term)
        );
        this.updateAllAnnotationsList();
    }
    
    updateAnnotationLists() {
        this.updateAllAnnotationsList();
        this.updateSelectedAnnotationsList();
        this.updateCounts();
        this.updateButtons();
    }
    
    updateAllAnnotationsList() {
        const list = this.elements.allList;
        list.innerHTML = '';
        
        this.filteredAnnotations.forEach(item => {
            const li = document.createElement('li');
            li.dataset.annotation = item.annotation;
            
            const info = document.createElement('div');
            info.className = 'annotation-info';
            
            const name = document.createElement('span');
            name.className = 'annotation-name';
            name.textContent = item.annotation;
            
            const count = document.createElement('span');
            count.className = 'annotation-count';
            count.textContent = `(${item.count})`;
            
            info.appendChild(name);
            info.appendChild(count);
            li.appendChild(info);
            
            li.addEventListener('click', () => {
                this.selectAnnotationInList(li, 'all');
            });
            
            li.addEventListener('dblclick', () => {
                this.addAnnotation(item.annotation);
            });
            
            list.appendChild(li);
        });
    }
    
    updateSelectedAnnotationsList() {
        const list = this.elements.selectedList;
        list.innerHTML = '';
        
        this.selectedAnnotations.forEach(annotation => {
            const li = document.createElement('li');
            li.dataset.annotation = annotation;
            
            const info = document.createElement('div');
            info.className = 'annotation-info';
            
            const chip = document.createElement('div');
            chip.className = 'color-chip';
            chip.style.backgroundColor = this.colorMap.get(annotation);
            
            const name = document.createElement('span');
            name.className = 'annotation-name';
            name.textContent = annotation;
            
            const removeBtn = document.createElement('button');
            removeBtn.className = 'remove-annotation';
            removeBtn.innerHTML = '&times;';
            removeBtn.addEventListener('click', (e) => {
                e.stopPropagation();
                this.removeAnnotation(annotation);
            });
            
            info.appendChild(chip);
            info.appendChild(name);
            li.appendChild(info);
            li.appendChild(removeBtn);
            
            li.addEventListener('click', () => {
                this.selectAnnotationInList(li, 'selected');
            });
            
            list.appendChild(li);
        });
        
        this.updateLegend();
    }
    
    updateLegend() {
        const container = this.elements.legendItems;
        container.innerHTML = '';
        
        this.selectedAnnotations.forEach(annotation => {
            const item = document.createElement('div');
            item.className = 'legend-item';
            
            const chip = document.createElement('div');
            chip.className = 'color-chip';
            chip.style.backgroundColor = this.colorMap.get(annotation);
            
            const name = document.createElement('span');
            name.textContent = annotation;
            
            item.appendChild(chip);
            item.appendChild(name);
            container.appendChild(item);
        });
    }
    
    selectAnnotationInList(li, listType) {
        // Clear previous selection in this list
        const list = listType === 'all' ? this.elements.allList : this.elements.selectedList;
        list.querySelectorAll('li.selected').forEach(item => {
            item.classList.remove('selected');
        });
        
        li.classList.add('selected');
        this.updateButtons();
    }
    
    updateButtons() {
        const allSelected = this.elements.allList.querySelector('li.selected');
        const selectedSelected = this.elements.selectedList.querySelector('li.selected');
        
        this.elements.addBtn.disabled = !allSelected || this.selectedAnnotations.length >= this.data.maxSelected;
        this.elements.removeBtn.disabled = !selectedSelected;
    }
    
    addSelectedAnnotation() {
        const selected = this.elements.allList.querySelector('li.selected');
        if (selected) {
            const annotation = selected.dataset.annotation;
            this.addAnnotation(annotation);
        }
    }
    
    removeSelectedAnnotation() {
        const selected = this.elements.selectedList.querySelector('li.selected');
        if (selected) {
            const annotation = selected.dataset.annotation;
            this.removeAnnotation(annotation);
        }
    }
    
    addAnnotation(annotation) {
        if (this.selectedAnnotations.includes(annotation)) return;
        if (this.selectedAnnotations.length >= this.data.maxSelected) {
            alert(`Maximum ${this.data.maxSelected} annotations can be selected.`);
            return;
        }
        
        this.selectedAnnotations.push(annotation);
        this.assignColor(annotation);
        this.updateAnnotationLists();
        this.render(); // Re-render both overview and detail views
    }
    
    removeAnnotation(annotation) {
        const index = this.selectedAnnotations.indexOf(annotation);
        if (index > -1) {
            this.selectedAnnotations.splice(index, 1);
            this.releaseColor(annotation);
            this.updateAnnotationLists();
            this.render(); // Re-render both overview and detail views
        }
    }
    
    assignColor(annotation) {
        if (!this.colorMap.has(annotation) && this.availableColors.length > 0) {
            const color = this.availableColors.shift();
            this.colorMap.set(annotation, color);
        }
    }
    
    releaseColor(annotation) {
        if (this.colorMap.has(annotation)) {
            const color = this.colorMap.get(annotation);
            this.colorMap.delete(annotation);
            this.availableColors.unshift(color);
        }
    }
    
    resetSelection() {
        this.selectedAnnotations = [];
        this.colorMap.clear();
        this.availableColors = [...this.data.palette];
        this.updateAnnotationLists();
        this.render(); // Re-render both overview and detail views
    }
    
    selectTopN(n) {
        this.resetSelection();
        const toSelect = this.allAnnotations.slice(0, Math.min(n, this.data.maxSelected));
        toSelect.forEach(item => {
            this.selectedAnnotations.push(item.annotation);
            this.assignColor(item.annotation);
        });
        this.updateAnnotationLists();
        this.render(); // Re-render both overview and detail views
    }
    
    updateCounts() {
        this.elements.allCount.textContent = `(${this.filteredAnnotations.length})`;
        this.elements.selectedCount.textContent = `(${this.selectedAnnotations.length}/${this.data.maxSelected})`;
    }
    
    applyInitialSelection() {
        if (this.data.initialAnnotations && this.data.initialAnnotations.length > 0) {
            this.data.initialAnnotations.forEach(annotation => {
                if (this.allAnnotations.some(item => item.annotation === annotation)) {
                    this.addAnnotation(annotation);
                }
            });
        }
    }
    
    setupScales() {
        this.overviewWidth = 2400; // Fixed width
        this.detailWidth = 2400; // Fixed width
        
        // Add margin to overview for easier selection at genome ends
        this.overviewMargin = 30; // 30px margin on each side
        this.overviewPlotWidth = this.overviewWidth - (2 * this.overviewMargin);
        
        // Fixed overview scale with margin
        this.overviewScale = (bp) => this.overviewMargin + (bp / this.data.genome_length) * this.overviewPlotWidth;
    }
    
    getDetailScale(startBp, endBp) {
        // Dynamic detail scale based on selected region
        const regionSize = endBp - startBp;
        return (bp) => ((bp - startBp) / regionSize) * this.detailWidth;
    }
    
    render() {
        this.renderOverview();
        this.renderDetail();
        this.updateZoomRegionHighlight();
    }
    
    // Selection interaction methods
    startSelection(e) {
        const rect = this.elements.overviewSvg.getBoundingClientRect();
        const x = e.clientX - rect.left;
        // Convert x position to genome coordinates considering margin
        this.selectionStart = Math.max(0, Math.min(this.data.genome_length, 
            ((x - this.overviewMargin) / this.overviewPlotWidth) * this.data.genome_length));
        this.isSelecting = true;
        
        // Clear existing selection box
        this.clearSelectionBox();
    }
    
    updateSelection(e) {
        if (!this.isSelecting) return;
        
        const rect = this.elements.overviewSvg.getBoundingClientRect();
        const x = e.clientX - rect.left;
        // Convert x position to genome coordinates considering margin
        const currentBp = Math.max(0, Math.min(this.data.genome_length,
            ((x - this.overviewMargin) / this.overviewPlotWidth) * this.data.genome_length));
        
        const startBp = Math.min(this.selectionStart, currentBp);
        const endBp = Math.max(this.selectionStart, currentBp);
        
        this.showSelectionBox(startBp, endBp);
    }
    
    endSelection(e) {
        if (!this.isSelecting) return;
        
        const rect = this.elements.overviewSvg.getBoundingClientRect();
        const x = e.clientX - rect.left;
        // Convert x position to genome coordinates considering margin
        const endBp = Math.max(0, Math.min(this.data.genome_length,
            ((x - this.overviewMargin) / this.overviewPlotWidth) * this.data.genome_length));
        
        const startBp = Math.min(this.selectionStart, endBp);
        const finalEndBp = Math.max(this.selectionStart, endBp);
        
        // Minimum selection size
        const minSelection = this.data.genome_length * 0.001; // 0.1% of genome
        if (finalEndBp - startBp > minSelection) {
            this.selectedRegion = { 
                startBp: Math.max(0, startBp), 
                endBp: Math.min(this.data.genome_length, finalEndBp) 
            };
            this.render();
        }
        
        this.isSelecting = false;
        this.clearSelectionBox();
    }
    
    cancelSelection() {
        this.isSelecting = false;
        this.clearSelectionBox();
    }
    
    resetSelection() {
        this.selectedRegion = null;
        this.render();
    }
    
    showSelectionBox(startBp, endBp) {
        this.clearSelectionBox();
        
        const x1 = this.overviewScale(startBp);
        const x2 = this.overviewScale(endBp);
        const width = x2 - x1;
        
        const selectionBox = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        selectionBox.setAttribute('x', x1);
        selectionBox.setAttribute('y', 10);
        selectionBox.setAttribute('width', width);
        selectionBox.setAttribute('height', this.overviewHeight - 20);
        selectionBox.classList.add('selection-box');
        selectionBox.id = 'temp-selection-box';
        
        this.elements.overviewSvg.appendChild(selectionBox);
    }
    
    clearSelectionBox() {
        const existing = document.getElementById('temp-selection-box');
        if (existing) {
            existing.remove();
        }
    }
    
    updateZoomRegionHighlight() {
        const highlight = this.elements.zoomRegionHighlight;
        highlight.innerHTML = '';
        
        if (!this.selectedRegion) return;
        
        const x1 = this.overviewScale(this.selectedRegion.startBp);
        const x2 = this.overviewScale(this.selectedRegion.endBp);
        const width = x2 - x1;
        
        const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        rect.setAttribute('x', x1);
        rect.setAttribute('y', 10);
        rect.setAttribute('width', width);
        rect.setAttribute('height', this.overviewHeight - 20);
        rect.classList.add('zoom-region');
        
        highlight.appendChild(rect);
        
        // Update region info
        const regionSizeMb = ((this.selectedRegion.endBp - this.selectedRegion.startBp) / 1000000).toFixed(2);
        this.elements.detailRegionInfo.textContent = 
            `(${this.selectedRegion.startBp.toLocaleString()} - ${this.selectedRegion.endBp.toLocaleString()} bp, ${regionSizeMb} Mb)`;
    }
    
    renderOverview() {
        // Render genome track
        this.renderOverviewGenomeTrack();
        
        // Render features
        this.renderOverviewFeatures();
        
        // Render contig labels (minimal)
        this.renderOverviewLabels();
    }
    
    renderOverviewGenomeTrack() {
        const track = this.elements.overviewGenomeTrack;
        track.innerHTML = '';
        
        const y = this.overviewHeight / 2 + 15; // Move genome line down to give more space above
        
        // Main genome line (only in the plot area, not in margins)
        const genomeLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        genomeLine.setAttribute('x1', this.overviewMargin);
        genomeLine.setAttribute('y1', y);
        genomeLine.setAttribute('x2', this.overviewMargin + this.overviewPlotWidth);
        genomeLine.setAttribute('y2', y);
        genomeLine.classList.add('genome-line');
        track.appendChild(genomeLine);
        
        // Contig separators
        this.data.contigs.forEach((contig, i) => {
            if (i > 0) {
                const x = this.overviewScale(contig.offset_bp);
                
                const separator = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                separator.setAttribute('x1', x);
                separator.setAttribute('y1', y - 20);
                separator.setAttribute('x2', x);
                separator.setAttribute('y2', y);
                separator.classList.add('contig-separator');
                track.appendChild(separator);
            }
        });
    }
    
    renderOverviewFeatures() {
        const layer = this.elements.overviewFeaturesLayer;
        layer.innerHTML = '';
        
        if (this.selectedAnnotations.length === 0) return;
        
        const y = this.overviewHeight / 2 + 15; // Match genome line position
        const rectHeight = 25;
        const selectedSet = new Set(this.selectedAnnotations);
        
        // Render features at overview scale
        this.data.features.forEach(feature => {
            const annotation = feature.annotation_main === null ? 'NA' : feature.annotation_main;
            if (!selectedSet.has(annotation)) return;
            
            const color = this.colorMap.get(annotation);
            if (!color) return;
            
            const contig = this.data.contigs.find(c => c.name === feature.contig);
            if (!contig) return;
            
            const startBp = contig.offset_bp + feature.start - 1;
            const endBp = contig.offset_bp + feature.end - 1;
            
            const x1 = this.overviewScale(startBp);
            const x2 = this.overviewScale(endBp);
            const width = Math.max(x2 - x1, 1);
            
            const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
            rect.setAttribute('x', x1);
            rect.setAttribute('y', y + 2); // Position below genome line
            rect.setAttribute('width', width);
            rect.setAttribute('height', rectHeight);
            rect.setAttribute('fill', color);
            rect.setAttribute('stroke', color); // Same color as fill
            rect.classList.add('feature-rect');
            
            layer.appendChild(rect);
        });
    }
    
    renderOverviewLabels() {
        const labelsLayer = this.elements.overviewLabelsLayer;
        labelsLayer.innerHTML = '';
        
        // Only show major contigs to avoid clutter
        const majorContigs = this.data.contigs.filter(contig => 
            contig.length > this.data.genome_length * 0.05 // > 5% of genome
        );
        
        const y = this.overviewHeight / 2 + 15 - 50; // Position relative to moved genome line
        
        majorContigs.forEach(contig => {
            const startX = Math.round(this.overviewScale(contig.offset_bp));
            
            const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            label.setAttribute('x', startX + 10);
            label.setAttribute('y', Math.round(y));
            label.setAttribute('transform', `rotate(-90 ${startX + 10} ${Math.round(y)})`);
            label.setAttribute('text-anchor', 'start');
            label.setAttribute('font-size', '10px'); // Overview font size
            label.classList.add('contig-label');
            label.textContent = contig.name;
            labelsLayer.appendChild(label);
        });
    }
    
    renderDetail() {
        if (!this.selectedRegion) {
            // Clear detail view
            this.elements.detailGenomeTrack.innerHTML = '';
            this.elements.detailFeaturesLayer.innerHTML = '';
            this.elements.detailLabelsLayer.innerHTML = '';
            this.elements.detailRegionInfo.textContent = '(No region selected)';
            return;
        }
        
        this.renderDetailGenomeTrack();
        this.renderDetailFeatures();
        this.renderDetailLabels();
    }
    
    renderDetailGenomeTrack() {
        const track = this.elements.detailGenomeTrack;
        track.innerHTML = '';
        
        const y = this.detailHeight / 2;
        const detailScale = this.getDetailScale(this.selectedRegion.startBp, this.selectedRegion.endBp);
        
        // Main genome line
        const genomeLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        genomeLine.setAttribute('x1', 0);
        genomeLine.setAttribute('y1', y);
        genomeLine.setAttribute('x2', this.detailWidth);
        genomeLine.setAttribute('y2', y);
        genomeLine.classList.add('genome-line');
        track.appendChild(genomeLine);
        
        // Contig separators within selected region
        this.data.contigs.forEach((contig, i) => {
            const contigStart = contig.offset_bp;
            const contigEnd = contig.offset_bp + contig.length;
            
            // Check if contig boundary is in selected region
            if (contigStart > this.selectedRegion.startBp && contigStart < this.selectedRegion.endBp) {
                const x = detailScale(contigStart);
                
                const separator = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                separator.setAttribute('x1', x);
                separator.setAttribute('y1', y - 25);
                separator.setAttribute('x2', x);
                separator.setAttribute('y2', y);
                separator.classList.add('contig-separator');
                track.appendChild(separator);
            }
        });
    }
    
    renderDetailFeatures() {
        const layer = this.elements.detailFeaturesLayer;
        layer.innerHTML = '';
        
        if (this.selectedAnnotations.length === 0 || !this.selectedRegion) return;
        
        const y = this.detailHeight / 2;
        const rectHeight = 35;
        const selectedSet = new Set(this.selectedAnnotations);
        const detailScale = this.getDetailScale(this.selectedRegion.startBp, this.selectedRegion.endBp);
        
        // Filter features in selected region
        const visibleFeatures = this.data.features.filter(feature => {
            const annotation = feature.annotation_main === null ? 'NA' : feature.annotation_main;
            if (!selectedSet.has(annotation)) return false;
            
            const contig = this.data.contigs.find(c => c.name === feature.contig);
            if (!contig) return false;
            
            const featureStartBp = contig.offset_bp + feature.start - 1;
            const featureEndBp = contig.offset_bp + feature.end - 1;
            
            // Check overlap with selected region
            return featureEndBp >= this.selectedRegion.startBp && 
                   featureStartBp <= this.selectedRegion.endBp;
        });
        
        // Sort by length (shorter on top)
        visibleFeatures.sort((a, b) => b.length - a.length);
        
        visibleFeatures.forEach(feature => {
            const annotation = feature.annotation_main === null ? 'NA' : feature.annotation_main;
            const color = this.colorMap.get(annotation);
            if (!color) return;
            
            const contig = this.data.contigs.find(c => c.name === feature.contig);
            const featureStartBp = contig.offset_bp + feature.start - 1;
            const featureEndBp = contig.offset_bp + feature.end - 1;
            
            // Clip to selected region
            const clippedStart = Math.max(featureStartBp, this.selectedRegion.startBp);
            const clippedEnd = Math.min(featureEndBp, this.selectedRegion.endBp);
            
            const x1 = detailScale(clippedStart);
            const x2 = detailScale(clippedEnd);
            const width = Math.max(x2 - x1, 1);
            
            const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
            rect.setAttribute('x', x1);
            rect.setAttribute('y', y + 2); // Position below genome line
            rect.setAttribute('width', width);
            rect.setAttribute('height', rectHeight);
            rect.setAttribute('fill', color);
            rect.setAttribute('stroke', color); // Same color as fill
            rect.classList.add('feature-rect');
            
            // Store feature data for tooltip
            rect._featureData = feature;
            
            layer.appendChild(rect);
        });
    }
    
    renderDetailLabels() {
        const labelsLayer = this.elements.detailLabelsLayer;
        labelsLayer.innerHTML = '';
        
        if (!this.selectedRegion) return;
        
        const y = this.detailHeight / 2 - 20;
        const detailScale = this.getDetailScale(this.selectedRegion.startBp, this.selectedRegion.endBp);
        
        // Find contigs that overlap with selected region
        this.data.contigs.forEach(contig => {
            const contigStart = contig.offset_bp;
            const contigEnd = contig.offset_bp + contig.length;
            
            // Check if contig overlaps with selected region
            if (contigEnd < this.selectedRegion.startBp || contigStart > this.selectedRegion.endBp) {
                return;
            }
            
            // Position label at start of contig within visible region
            const labelBp = Math.max(contigStart, this.selectedRegion.startBp);
            const x = Math.round(detailScale(labelBp) + 10); // Round to whole pixels
            
            const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            label.setAttribute('x', x);
            label.setAttribute('y', Math.round(y)); // Round Y position too
            label.setAttribute('transform', `rotate(-90 ${x} ${Math.round(y)})`);
            label.setAttribute('text-anchor', 'start');
            label.setAttribute('font-size', '14px'); // Bigger font for detail view
            label.classList.add('contig-label');
            label.textContent = contig.name;
            labelsLayer.appendChild(label);
        });
    }
    
    // Tooltip functionality
    updateTooltip(e) {
        const target = e.target;
        if (target._featureData) {
            this.showTooltip(e, target._featureData);
        } else {
            this.hideTooltip();
        }
    }
    
    showTooltip(e, feature) {
        const tooltip = this.elements.tooltip;
        
        const annotation = feature.annotation_main === null ? 'NA' : feature.annotation_main;
        
        let content = `
            <div class="tooltip-field"><span class="tooltip-label">Name:</span> ${feature.name || 'N/A'}</div>
            <div class="tooltip-field"><span class="tooltip-label">Annotation:</span> ${annotation}</div>
            <div class="tooltip-field"><span class="tooltip-label">Raw:</span> ${feature.annotation_raw}</div>
            <div class="tooltip-field"><span class="tooltip-label">Type:</span> ${feature.repeat_type || 'N/A'}</div>
            <div class="tooltip-field"><span class="tooltip-label">Contig:</span> ${feature.contig}</div>
            <div class="tooltip-field"><span class="tooltip-label">Position:</span> ${feature.start.toLocaleString()}-${feature.end.toLocaleString()}</div>
            <div class="tooltip-field"><span class="tooltip-label">Length:</span> ${feature.length.toLocaleString()} bp</div>
        `;
        
        if (feature.ssr) {
            content += `<div class="tooltip-field"><span class="tooltip-label">SSR:</span> ${feature.ssr}</div>`;
        }
        
        tooltip.innerHTML = content;
        tooltip.style.left = (e.pageX + 10) + 'px';
        tooltip.style.top = (e.pageY - 10) + 'px';
        tooltip.classList.add('visible');
    }
    
    hideTooltip() {
        this.elements.tooltip.classList.remove('visible');
    }
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    new TideClusterViz();
});