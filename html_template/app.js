// TideCluster Visualizer - Main Application
class TideClusterViz {
    constructor() {
        this.data = null;
        this.scales = { x: null };
        this.selectedAnnotations = [];
        this.annotationVisibility = new Map(); // Track visibility state for each annotation
        this.colorMap = new Map();
        this.availableColors = [];
        this.filteredAnnotations = [];
        this.allAnnotations = [];
        
        // Color picker state
        this.currentColorPickerAnnotation = null;
        this.allAvailableColors = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
            '#ff6b6b', '#4ecdc4', '#45b7d1', '#f39c12', '#e74c3c',
            '#9b59b6', '#3498db', '#2ecc71', '#f1c40f', '#e67e22',
            '#34495e', '#95a5a6', '#16a085', '#27ae60', '#8e44ad',
            '#2980b9', '#c0392b', '#d35400', '#7f8c8d', '#2c3e50'
        ];
        
        // Selection state - three levels
        this.selectedRegion = null; // Overview to detail { startBp, endBp }
        this.detailSelectedRegion = null; // Detail to zoom { startBp, endBp }
        this.isSelectingOverview = false;
        this.isSelectingDetail = false;
        this.selectionStartOverview = 0;
        this.selectionStartDetail = 0;
        this.overviewWidth = 0;
        this.overviewHeight = 180;
        this.detailWidth = 0;
        this.detailHeight = 180;
        this.zoomWidth = 0;
        this.zoomHeight = 180;
        
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
            detailZoomRegionHighlight: document.getElementById('detail-zoom-region-highlight'),
            
            // Zoom elements
            zoomSvg: document.getElementById('zoom-svg'),
            zoomGenomeTrack: document.getElementById('zoom-genome-track'),
            zoomFeaturesLayer: document.getElementById('zoom-features-layer'),
            zoomLabelsLayer: document.getElementById('zoom-labels-layer'),
            zoomRegionInfo: document.getElementById('zoom-region-info'),
            
            tooltip: document.getElementById('tooltip'),
            zoomReset: document.getElementById('zoom-reset'),
            detailZoomReset: document.getElementById('detail-zoom-reset'),
            
            // Color picker elements
            colorPicker: document.getElementById('color-picker'),
            colorPickerGrid: document.getElementById('color-picker-grid'),
            colorPickerClose: document.getElementById('color-picker-close')
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
        
        // Reset selections
        this.elements.zoomReset.addEventListener('click', () => {
            this.resetAllSelections();
        });
        
        this.elements.detailZoomReset.addEventListener('click', () => {
            this.resetDetailSelection();
        });
        
        // SVG interactions for selection
        this.setupOverviewInteractions();
        this.setupDetailInteractions();
        
        // Color picker interactions
        this.setupColorPicker();
    }
    
    setupOverviewInteractions() {
        const overviewSvg = this.elements.overviewSvg;
        
        // Mouse selection on overview
        overviewSvg.addEventListener('mousedown', (e) => {
            this.startOverviewSelection(e);
        });
        
        overviewSvg.addEventListener('mousemove', (e) => {
            this.updateOverviewSelection(e);
            this.updateTooltip(e);
        });
        
        overviewSvg.addEventListener('mouseup', (e) => {
            this.endOverviewSelection(e);
        });
        
        overviewSvg.addEventListener('mouseleave', () => {
            this.cancelOverviewSelection();
            this.hideTooltip();
        });
        
        // Prevent context menu
        overviewSvg.addEventListener('contextmenu', (e) => e.preventDefault());
    }
    
    setupDetailInteractions() {
        const detailSvg = this.elements.detailSvg;
        
        // Mouse selection on detail view
        detailSvg.addEventListener('mousedown', (e) => {
            this.startDetailSelection(e);
        });
        
        detailSvg.addEventListener('mousemove', (e) => {
            this.updateDetailSelection(e);
            this.updateTooltip(e);
        });
        
        detailSvg.addEventListener('mouseup', (e) => {
            this.endDetailSelection(e);
        });
        
        detailSvg.addEventListener('mouseleave', () => {
            this.cancelDetailSelection();
            this.hideTooltip();
        });
        
        // Tooltip on zoom view
        const zoomSvg = this.elements.zoomSvg;
        zoomSvg.addEventListener('mousemove', (e) => {
            this.updateTooltip(e);
        });
        
        zoomSvg.addEventListener('mouseleave', () => {
            this.hideTooltip();
        });
        
        // Prevent context menu
        detailSvg.addEventListener('contextmenu', (e) => e.preventDefault());
        zoomSvg.addEventListener('contextmenu', (e) => e.preventDefault());
    }
    
    setupColorPicker() {
        // Close button
        this.elements.colorPickerClose.addEventListener('click', () => {
            this.hideColorPicker();
        });
        
        // Close when clicking outside
        document.addEventListener('click', (e) => {
            if (this.elements.colorPicker.classList.contains('visible') && 
                !this.elements.colorPicker.contains(e.target) &&
                !e.target.classList.contains('color-chip')) {
                this.hideColorPicker();
            }
        });
        
        // Initialize color picker grid
        this.initializeColorPickerGrid();
    }
    
    initializeColorPickerGrid() {
        const grid = this.elements.colorPickerGrid;
        grid.innerHTML = '';
        
        this.allAvailableColors.forEach(color => {
            const colorOption = document.createElement('div');
            colorOption.className = 'color-option';
            colorOption.style.backgroundColor = color;
            colorOption.dataset.color = color;
            colorOption.title = color;
            
            colorOption.addEventListener('click', () => {
                if (this.currentColorPickerAnnotation) {
                    this.changeAnnotationColor(this.currentColorPickerAnnotation, color);
                }
                this.hideColorPicker();
            });
            
            grid.appendChild(colorOption);
        });
    }
    
    showColorPicker(annotation) {
        this.currentColorPickerAnnotation = annotation;
        
        // Update current color indication
        const currentColor = this.colorMap.get(annotation);
        this.elements.colorPickerGrid.querySelectorAll('.color-option').forEach(option => {
            if (option.dataset.color === currentColor) {
                option.classList.add('current');
            } else {
                option.classList.remove('current');
            }
        });
        
        this.elements.colorPicker.classList.add('visible');
    }
    
    hideColorPicker() {
        this.elements.colorPicker.classList.remove('visible');
        this.currentColorPickerAnnotation = null;
    }
    
    changeAnnotationColor(annotation, newColor) {
        const oldColor = this.colorMap.get(annotation);
        
        // Check if new color is already in use by another annotation
        const colorInUse = Array.from(this.colorMap.entries()).find(([ann, color]) => 
            ann !== annotation && color === newColor
        );
        
        if (colorInUse) {
            // Swap colors
            const otherAnnotation = colorInUse[0];
            this.colorMap.set(otherAnnotation, oldColor);
            this.colorMap.set(annotation, newColor);
        } else {
            // Simple color change
            this.colorMap.set(annotation, newColor);
            
            // Return old color to available colors if it's not the new one
            if (oldColor && oldColor !== newColor && !this.availableColors.includes(oldColor)) {
                this.availableColors.push(oldColor);
            }
            
            // Remove new color from available colors
            const availableIndex = this.availableColors.indexOf(newColor);
            if (availableIndex > -1) {
                this.availableColors.splice(availableIndex, 1);
            }
        }
        
        // Update displays
        this.updateAnnotationLists();
        this.render();
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
            chip.title = 'Click to change color';
            chip.addEventListener('click', (e) => {
                e.stopPropagation();
                this.showColorPicker(annotation);
            });
            
            const name = document.createElement('span');
            name.className = 'annotation-name';
            name.textContent = annotation;
            
            info.appendChild(chip);
            info.appendChild(name);
            
            // Create controls container (visibility + remove button)
            const controlsContainer = document.createElement('div');
            controlsContainer.className = 'annotation-controls';
            
            // Visibility checkbox
            const visibilityContainer = document.createElement('div');
            visibilityContainer.className = 'visibility-control';
            
            const visibilityCheckbox = document.createElement('input');
            visibilityCheckbox.type = 'checkbox';
            visibilityCheckbox.checked = this.annotationVisibility.get(annotation) !== false; // Default to visible
            visibilityCheckbox.addEventListener('change', (e) => {
                e.stopPropagation();
                this.toggleAnnotationVisibility(annotation, e.target.checked);
            });
            
            const visibilityLabel = document.createElement('label');
            visibilityLabel.appendChild(visibilityCheckbox);
            visibilityLabel.appendChild(document.createTextNode('Show'));
            visibilityLabel.title = 'Toggle visibility in plots';
            
            visibilityContainer.appendChild(visibilityLabel);
            
            const removeBtn = document.createElement('button');
            removeBtn.className = 'remove-annotation';
            removeBtn.innerHTML = '&times;';
            removeBtn.addEventListener('click', (e) => {
                e.stopPropagation();
                this.removeAnnotation(annotation);
            });
            
            controlsContainer.appendChild(visibilityContainer);
            controlsContainer.appendChild(removeBtn);
            
            li.appendChild(info);
            li.appendChild(controlsContainer);
            
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
            
            // Add visual indication for hidden annotations
            const isVisible = this.annotationVisibility.get(annotation) !== false;
            if (!isVisible) {
                item.style.opacity = '0.5';
                item.title = 'Hidden annotation (click Show checkbox to display)';
            }
            
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
        
        // Set default visibility to true
        this.annotationVisibility.set(annotation, true);
        
        this.assignColor(annotation);
        this.updateAnnotationLists();
        this.render(); // Re-render both overview and detail views
    }
    
    toggleAnnotationVisibility(annotation, visible) {
        this.annotationVisibility.set(annotation, visible);
        this.render(); // Re-render all views
        this.updateLegend();
    }
    
    removeAnnotation(annotation) {
        const index = this.selectedAnnotations.indexOf(annotation);
        if (index > -1) {
            this.selectedAnnotations.splice(index, 1);
            this.releaseColor(annotation);
            // Clean up visibility state
            this.annotationVisibility.delete(annotation);
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
        this.zoomWidth = 2400; // Fixed width
        
        // Add margin to overview for easier selection at genome ends
        this.overviewMargin = 30; // 30px margin on each side
        this.overviewPlotWidth = this.overviewWidth - (2 * this.overviewMargin);
        
        // Detail margins
        this.detailMargin = 30;
        this.detailPlotWidth = this.detailWidth - (2 * this.detailMargin);
        
        // Zoom margins
        this.zoomMargin = 30;
        this.zoomPlotWidth = this.zoomWidth - (2 * this.zoomMargin);
        
        // Fixed overview scale with margin
        this.overviewScale = (bp) => this.overviewMargin + (bp / this.data.genome_length) * this.overviewPlotWidth;
    }
    
    getDetailScale() {
        if (!this.selectedRegion) return null;
        // Dynamic detail scale based on selected region
        const startBp = this.selectedRegion.startBp;
        const endBp = this.selectedRegion.endBp;
        const regionSize = endBp - startBp;
        return (bp) => ((bp - startBp) / regionSize) * this.detailPlotWidth + this.detailMargin;
    }
    
    render() {
        this.renderOverview();
        this.renderDetail();
        this.renderZoom();
        this.updateZoomRegionHighlight();
        this.updateDetailZoomRegionHighlight();
    }
    
    // Overview selection interaction methods
    startOverviewSelection(e) {
        const rect = this.elements.overviewSvg.getBoundingClientRect();
        const x = e.clientX - rect.left;
        // Convert x position to genome coordinates considering margin
        this.selectionStartOverview = Math.max(0, Math.min(this.data.genome_length, 
            ((x - this.overviewMargin) / this.overviewPlotWidth) * this.data.genome_length));
        this.isSelectingOverview = true;
        
        // Clear existing selection box
        this.clearOverviewSelectionBox();
    }
    
    updateOverviewSelection(e) {
        if (!this.isSelectingOverview) return;
        
        const rect = this.elements.overviewSvg.getBoundingClientRect();
        const x = e.clientX - rect.left;
        // Convert x position to genome coordinates considering margin
        const currentBp = Math.max(0, Math.min(this.data.genome_length,
            ((x - this.overviewMargin) / this.overviewPlotWidth) * this.data.genome_length));
        
        const startBp = Math.min(this.selectionStartOverview, currentBp);
        const endBp = Math.max(this.selectionStartOverview, currentBp);
        
        this.showOverviewSelectionBox(startBp, endBp);
    }
    
    endOverviewSelection(e) {
        if (!this.isSelectingOverview) return;
        
        const rect = this.elements.overviewSvg.getBoundingClientRect();
        const x = e.clientX - rect.left;
        // Convert x position to genome coordinates considering margin
        const endBp = Math.max(0, Math.min(this.data.genome_length,
            ((x - this.overviewMargin) / this.overviewPlotWidth) * this.data.genome_length));
        
        const startBp = Math.min(this.selectionStartOverview, endBp);
        const finalEndBp = Math.max(this.selectionStartOverview, endBp);
        
        // Minimum selection size
        const minSelection = this.data.genome_length * 0.001; // 0.1% of genome
        if (finalEndBp - startBp > minSelection) {
            this.selectedRegion = { 
                startBp: Math.max(0, startBp), 
                endBp: Math.min(this.data.genome_length, finalEndBp) 
            };
            // Reset detail selection when overview changes
            this.detailSelectedRegion = null;
            this.render();
        }
        
        this.isSelectingOverview = false;
        this.clearOverviewSelectionBox();
    }
    
    cancelOverviewSelection() {
        this.isSelectingOverview = false;
        this.clearOverviewSelectionBox();
    }
    
    resetAllSelections() {
        this.selectedRegion = null;
        this.detailSelectedRegion = null;
        this.render();
    }
    
    resetDetailSelection() {
        this.detailSelectedRegion = null;
        this.render();
    }
    
    // Detail selection interaction methods
    startDetailSelection(e) {
        if (!this.selectedRegion) return; // Can only select on detail if overview region is selected
        
        const rect = this.elements.detailSvg.getBoundingClientRect();
        const x = e.clientX - rect.left;
        const detailScale = this.getDetailScale();
        
        // Convert detail x position back to genome coordinates
        const regionSize = this.selectedRegion.endBp - this.selectedRegion.startBp;
        const relativeX = (x - this.detailMargin) / this.detailPlotWidth;
        this.selectionStartDetail = this.selectedRegion.startBp + (relativeX * regionSize);
        this.selectionStartDetail = Math.max(this.selectedRegion.startBp, 
            Math.min(this.selectedRegion.endBp, this.selectionStartDetail));
        this.isSelectingDetail = true;
        
        this.clearDetailSelectionBox();
    }
    
    updateDetailSelection(e) {
        if (!this.isSelectingDetail || !this.selectedRegion) return;
        
        const rect = this.elements.detailSvg.getBoundingClientRect();
        const x = e.clientX - rect.left;
        
        // Convert detail x position back to genome coordinates
        const regionSize = this.selectedRegion.endBp - this.selectedRegion.startBp;
        const relativeX = (x - this.detailMargin) / this.detailPlotWidth;
        const currentBp = this.selectedRegion.startBp + (relativeX * regionSize);
        const clampedBp = Math.max(this.selectedRegion.startBp, 
            Math.min(this.selectedRegion.endBp, currentBp));
        
        const startBp = Math.min(this.selectionStartDetail, clampedBp);
        const endBp = Math.max(this.selectionStartDetail, clampedBp);
        
        this.showDetailSelectionBox(startBp, endBp);
    }
    
    endDetailSelection(e) {
        if (!this.isSelectingDetail || !this.selectedRegion) return;
        
        const rect = this.elements.detailSvg.getBoundingClientRect();
        const x = e.clientX - rect.left;
        
        // Convert detail x position back to genome coordinates
        const regionSize = this.selectedRegion.endBp - this.selectedRegion.startBp;
        const relativeX = (x - this.detailMargin) / this.detailPlotWidth;
        const endBp = this.selectedRegion.startBp + (relativeX * regionSize);
        const clampedEndBp = Math.max(this.selectedRegion.startBp, 
            Math.min(this.selectedRegion.endBp, endBp));
        
        const startBp = Math.min(this.selectionStartDetail, clampedEndBp);
        const finalEndBp = Math.max(this.selectionStartDetail, clampedEndBp);
        
        // Minimum selection size relative to detail view
        const minSelection = (this.selectedRegion.endBp - this.selectedRegion.startBp) * 0.01; // 1% of detail region
        if (finalEndBp - startBp > minSelection) {
            this.detailSelectedRegion = {
                startBp: Math.max(this.selectedRegion.startBp, startBp),
                endBp: Math.min(this.selectedRegion.endBp, finalEndBp)
            };
            this.render();
        }
        
        this.isSelectingDetail = false;
        this.clearDetailSelectionBox();
    }
    
    cancelDetailSelection() {
        this.isSelectingDetail = false;
        this.clearDetailSelectionBox();
    }
    
    showOverviewSelectionBox(startBp, endBp) {
        this.clearOverviewSelectionBox();
        
        const x1 = this.overviewScale(startBp);
        const x2 = this.overviewScale(endBp);
        const width = x2 - x1;
        
        const selectionBox = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        selectionBox.setAttribute('x', x1);
        selectionBox.setAttribute('y', 10);
        selectionBox.setAttribute('width', width);
        selectionBox.setAttribute('height', this.overviewHeight - 20);
        selectionBox.classList.add('selection-box');
        selectionBox.id = 'temp-overview-selection-box';
        
        this.elements.overviewSvg.appendChild(selectionBox);
    }
    
    clearOverviewSelectionBox() {
        const existing = document.getElementById('temp-overview-selection-box');
        if (existing) {
            existing.remove();
        }
    }
    
    showDetailSelectionBox(startBp, endBp) {
        this.clearDetailSelectionBox();
        
        if (!this.selectedRegion) return;
        
        const detailScale = this.getDetailScale();
        const x1 = detailScale(startBp);
        const x2 = detailScale(endBp);
        const width = x2 - x1;
        
        const selectionBox = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        selectionBox.setAttribute('x', x1);
        selectionBox.setAttribute('y', 10);
        selectionBox.setAttribute('width', width);
        selectionBox.setAttribute('height', this.detailHeight - 20);
        selectionBox.classList.add('selection-box');
        selectionBox.id = 'temp-detail-selection-box';
        
        this.elements.detailSvg.appendChild(selectionBox);
    }
    
    clearDetailSelectionBox() {
        const existing = document.getElementById('temp-detail-selection-box');
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
    
    updateDetailZoomRegionHighlight() {
        const highlight = this.elements.detailZoomRegionHighlight;
        highlight.innerHTML = '';
        
        if (!this.detailSelectedRegion || !this.selectedRegion) return;
        
        const detailScale = this.getDetailScale();
        const x1 = detailScale(this.detailSelectedRegion.startBp);
        const x2 = detailScale(this.detailSelectedRegion.endBp);
        const width = x2 - x1;
        
        const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        rect.setAttribute('x', x1);
        rect.setAttribute('y', 10);
        rect.setAttribute('width', width);
        rect.setAttribute('height', this.detailHeight - 20);
        rect.classList.add('zoom-region');
        
        highlight.appendChild(rect);
        
        // Update zoom region info
        const zoomRegionSizeKb = ((this.detailSelectedRegion.endBp - this.detailSelectedRegion.startBp) / 1000).toFixed(2);
        this.elements.zoomRegionInfo.textContent = 
            `(${this.detailSelectedRegion.startBp.toLocaleString()} - ${this.detailSelectedRegion.endBp.toLocaleString()} bp, ${zoomRegionSizeKb} kb)`;
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
            if (this.annotationVisibility.get(annotation) === false) return; // Skip hidden annotations
            
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
        const detailScale = this.getDetailScale();
        
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
        const detailScale = this.getDetailScale();
        
        // Filter features in selected region
        const visibleFeatures = this.data.features.filter(feature => {
            const annotation = feature.annotation_main === null ? 'NA' : feature.annotation_main;
            if (!selectedSet.has(annotation)) return false;
            if (this.annotationVisibility.get(annotation) === false) return false; // Skip hidden annotations
            
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
        const detailScale = this.getDetailScale();
        
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
    
    renderZoom() {
        if (!this.detailSelectedRegion) {
            // Clear zoom view if no detail region selected
            this.elements.zoomGenomeTrack.innerHTML = '';
            this.elements.zoomFeaturesLayer.innerHTML = '';
            this.elements.zoomLabelsLayer.innerHTML = '';
            this.elements.zoomRegionInfo.textContent = '';
            return;
        }
        
        // Render genome track
        this.renderZoomGenomeTrack();
        
        // Render features
        this.renderZoomFeatures();
        
        // Render contig labels
        this.renderZoomLabels();
    }
    
    renderZoomGenomeTrack() {
        const track = this.elements.zoomGenomeTrack;
        track.innerHTML = '';
        
        if (!this.detailSelectedRegion) return;
        
        const y = this.zoomHeight / 2 + 15;
        
        // Main genome line
        const genomeLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        genomeLine.setAttribute('x1', this.zoomMargin);
        genomeLine.setAttribute('y1', y);
        genomeLine.setAttribute('x2', this.zoomMargin + this.zoomPlotWidth);
        genomeLine.setAttribute('y2', y);
        genomeLine.classList.add('genome-line');
        track.appendChild(genomeLine);
        
        // Contig separators within zoom region
        this.data.contigs.forEach(contig => {
            const contigStart = contig.offset_bp;
            const contigEnd = contig.offset_bp + contig.length;
            
            // Check if contig boundary is within zoom region
            if (contigStart >= this.detailSelectedRegion.startBp && contigStart <= this.detailSelectedRegion.endBp) {
                const zoomScale = this.getZoomScale();
                const x = zoomScale(contigStart);
                
                const separator = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                separator.setAttribute('x1', x);
                separator.setAttribute('y1', y - 10);
                separator.setAttribute('x2', x);
                separator.setAttribute('y2', y + 10);
                separator.classList.add('contig-separator');
                track.appendChild(separator);
            }
        });
    }
    
    renderZoomFeatures() {
        const layer = this.elements.zoomFeaturesLayer;
        layer.innerHTML = '';
        
        if (!this.detailSelectedRegion || this.selectedAnnotations.length === 0) return;
        
        const zoomScale = this.getZoomScale();
        const y = this.zoomHeight / 2 + 15;
        const rectHeight = 45; // Tallest for maximum zoom level (overview: 25, detail: 35, zoom: 45)
        const selectedSet = new Set(this.selectedAnnotations);
        
        
        // Filter features that overlap with zoom region
        const visibleFeatures = this.data.features.filter(feature => {
            const annotation = feature.annotation_main === null ? 'NA' : feature.annotation_main;
            if (!selectedSet.has(annotation)) return false;
            if (this.annotationVisibility.get(annotation) === false) return false; // Skip hidden annotations
            
            const contig = this.data.contigs.find(c => c.name === feature.contig);
            if (!contig) return false;
            
            const featureStartBp = contig.offset_bp + feature.start - 1;
            const featureEndBp = contig.offset_bp + feature.end - 1;
            
            const overlaps = featureEndBp >= this.detailSelectedRegion.startBp && 
                   featureStartBp <= this.detailSelectedRegion.endBp;
            
            return overlaps;
        });
        
        
        // Sort by length (shorter on top) to match detail view
        visibleFeatures.sort((a, b) => b.length - a.length);
        
        visibleFeatures.forEach(feature => {
            const annotation = feature.annotation_main === null ? 'NA' : feature.annotation_main;
            const color = this.colorMap.get(annotation);
            if (!color) return;
            
            const contig = this.data.contigs.find(c => c.name === feature.contig);
            if (!contig) return;
            
            const featureStartBp = contig.offset_bp + feature.start - 1;
            const featureEndBp = contig.offset_bp + feature.end - 1;
            
            const x1 = Math.max(this.zoomMargin, zoomScale(Math.max(featureStartBp, this.detailSelectedRegion.startBp)));
            const x2 = Math.min(this.zoomMargin + this.zoomPlotWidth, zoomScale(Math.min(featureEndBp, this.detailSelectedRegion.endBp)));
            const width = Math.max(1, x2 - x1);
            
            const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
            rect.setAttribute('x', x1);
            rect.setAttribute('y', y + 2); // Position below genome line, consistent with other panels
            rect.setAttribute('width', width);
            rect.setAttribute('height', rectHeight);
            rect.setAttribute('fill', color);
            rect.setAttribute('stroke', color);
            rect.classList.add('feature-rect');
            
            // Store feature data for tooltip
            rect._featureData = feature;
            
            layer.appendChild(rect);
        });
    }
    
    renderZoomLabels() {
        const labelsLayer = this.elements.zoomLabelsLayer;
        labelsLayer.innerHTML = '';
        
        if (!this.detailSelectedRegion) return;
        
        const y = this.zoomHeight / 2 - 20;
        const zoomScale = this.getZoomScale();
        
        // Find contigs that overlap with zoom region
        this.data.contigs.forEach(contig => {
            const contigStart = contig.offset_bp;
            const contigEnd = contig.offset_bp + contig.length;
            
            // Check if contig overlaps with zoom region
            if (contigEnd < this.detailSelectedRegion.startBp || contigStart > this.detailSelectedRegion.endBp) {
                return;
            }
            
            // Position label at start of contig within visible region
            const labelBp = Math.max(contigStart, this.detailSelectedRegion.startBp);
            const x = Math.round(zoomScale(labelBp) + 10);
            
            const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            label.setAttribute('x', x);
            label.setAttribute('y', Math.round(y));
            label.setAttribute('transform', `rotate(-90 ${x} ${Math.round(y)})`);
            label.setAttribute('text-anchor', 'start');
            label.setAttribute('font-size', '16px'); // Largest font for zoom view
            label.classList.add('contig-label');
            label.textContent = contig.name;
            labelsLayer.appendChild(label);
        });
    }
    
    getZoomScale() {
        if (!this.detailSelectedRegion) return null;
        const startBp = this.detailSelectedRegion.startBp;
        const endBp = this.detailSelectedRegion.endBp;
        const regionSize = endBp - startBp;
        return (bp) => ((bp - startBp) / regionSize) * this.zoomPlotWidth + this.zoomMargin;
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