#!/usr/bin/env python3

"""
TideCluster Results Visualizer
Converts TideCluster GFF3 + genome .fai into interactive HTML visualization
"""

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path
from urllib.parse import unquote


def parse_gff3_attributes(attr_string):
    """Parse GFF3 attributes string into dictionary"""
    attributes = {}
    if not attr_string or attr_string == '.':
        return attributes
    
    # Split by semicolon and parse key=value pairs
    for pair in attr_string.split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            # URL decode the value
            attributes[key] = unquote(value)
    
    return attributes


def url_decode(text):
    """URL decode text with common encodings"""
    if not text or text == "NA":
        return text
    return unquote(text)


def derive_annotation_main(annotation_raw):
    """
    Derive main annotation from raw annotation string.
    Returns tuple: (annotation_main, is_ambiguous)
    """
    if not annotation_raw or annotation_raw == "NA":
        return None, False
    
    # Split by commas and parse percentages
    labels = [label.strip() for label in annotation_raw.split(',')]
    
    max_pct = 0
    best_label = None
    
    for label in labels:
        # Extract percentage from parentheses
        pct_match = re.search(r'\((\d+(?:\.\d+)?)%\)', label)
        if pct_match:
            pct_val = float(pct_match.group(1))
            if pct_val > max_pct:
                max_pct = pct_val
                # Remove percentage part for annotation_main
                best_label = re.sub(r'\s*\([^)]*\)\s*$', '', label).strip()
    
    if max_pct > 50:
        return best_label, False
    else:
        return "ambiguous", True


def load_fai(fai_path):
    """Load genome FAI file and return contig information"""
    contigs = []
    offset_bp = 0
    
    with open(fai_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                name = parts[0]
                length = int(parts[1])
                
                contigs.append({
                    'name': name,
                    'length': length,
                    'offset_bp': offset_bp
                })
                
                offset_bp += length
    
    return contigs, offset_bp


def load_gff3(gff3_path):
    """Load and process GFF3 file"""
    features = []
    
    with open(gff3_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) != 9:
                continue
            
            seqname, source, feature_type, start, end, score, strand, phase, attributes = parts
            
            # Only process tandem_repeat features
            if feature_type != 'tandem_repeat':
                continue
            
            # Parse attributes
            attr_dict = parse_gff3_attributes(attributes)
            
            # Extract and decode relevant attributes
            name_val = attr_dict.get('Name')
            repeat_type = attr_dict.get('repeat_type')
            annotation_raw = attr_dict.get('annotation', 'NA')
            ssr_val = attr_dict.get('ssr')
            
            # Derive annotation_main
            annotation_main, is_ambiguous = derive_annotation_main(annotation_raw)
            
            # Create feature record
            feature_record = {
                'contig': seqname,
                'start': int(start),
                'end': int(end),
                'length': int(end) - int(start) + 1,
                'name': name_val,
                'repeat_type': repeat_type,
                'annotation_main': annotation_main,
                'annotation_raw': annotation_raw,
                'is_ambiguous': is_ambiguous
            }
            
            # Add SSR if present
            if ssr_val:
                feature_record['ssr'] = ssr_val
            
            features.append(feature_record)
    
    return features


def main():
    parser = argparse.ArgumentParser(
        description='TideCluster Results Visualizer - Convert GFF3 + FAI to interactive HTML'
    )
    parser.add_argument('--gff3', required=True, help='Path to TideCluster GFF3 file')
    parser.add_argument('--fai', required=True, help='Path to genome FASTA index (.fai) file')
    parser.add_argument('--out', default='output.html', help='Output HTML file path (default: output.html)')
    parser.add_argument('--title', default='TideCluster Results', help='Visualization title')
    parser.add_argument('--maxSelected', type=int, default=10, help='Maximum selected annotations (default: 10)')
    parser.add_argument('--palette', 
                       default='#1f77b4,#ff7f0e,#2ca02c,#d62728,#9467bd,#8c564b,#e377c2,#7f7f7f,#bcbd22,#17becf',
                       help='Color palette (comma-separated hex colors)')
    parser.add_argument('--initialAnnotations', default='', 
                       help='Preselected annotations (comma-separated)')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.gff3):
        print(f"Error: GFF3 file not found: {args.gff3}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.fai):
        print(f"Error: FAI file not found: {args.fai}", file=sys.stderr)
        sys.exit(1)
    
    # Get script directory for template files
    script_dir = Path(__file__).parent
    template_dir = script_dir / 'html_template'
    
    if not template_dir.exists():
        print(f"Error: Template directory not found: {template_dir}", file=sys.stderr)
        sys.exit(1)
    
    print("Loading data...")
    
    # Load genome FAI
    contigs, genome_length = load_fai(args.fai)
    print(f"Genome length: {genome_length:,} bp across {len(contigs)} contigs")
    
    # Load GFF3
    features = load_gff3(args.gff3)
    print(f"Loaded {len(features)} tandem repeat features")
    
    if len(features) == 0:
        print("Warning: No tandem repeat features found in GFF3")
    
    # Create output directory
    output_path = Path(args.out)
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_files_dir = output_dir / 'output_files'
    output_files_dir.mkdir(exist_ok=True)
    
    print(f"Creating output in: {output_files_dir}")
    
    # Copy template files
    shutil.copy2(template_dir / 'app.js', output_files_dir)
    shutil.copy2(template_dir / 'styles.css', output_files_dir)
    
    # Create data.json
    data_json = {
        'title': args.title,
        'genome_length': genome_length,
        'maxSelected': args.maxSelected,
        'palette': args.palette.split(','),
        'initialAnnotations': args.initialAnnotations.split(',') if args.initialAnnotations else [],
        'contigs': contigs,
        'features': features
    }
    
    # Write JSON
    json_file = output_files_dir / 'data.json'
    with open(json_file, 'w') as f:
        json.dump(data_json, f, indent=2)
    
    # Create HTML file with embedded data
    html_template = template_dir / 'index.html'
    with open(html_template, 'r') as f:
        html_content = f.read()
    
    # Embed data in HTML
    data_json_str = json.dumps(data_json)
    html_content = html_content.replace('{{DATA_JSON}}', data_json_str)
    
    # Write output HTML
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    print("Visualization created successfully!")
    print(f"Open {args.out} in your browser")


if __name__ == '__main__':
    main()