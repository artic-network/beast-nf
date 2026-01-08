#!/usr/bin/env python3
"""
Generate HTML report for BEAST analysis
"""

import argparse
import sys
import os
import re
from datetime import datetime
from pathlib import Path
try:
    from Bio import SeqIO
except ImportError:
    print("Warning: Biopython not found. Install with: pip install biopython", file=sys.stderr)
    SeqIO = None


def parse_fasta_info(fasta_file):
    """Extract information from FASTA file"""
    if SeqIO is None:
        return None
    
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    info = {
        'num_taxa': len(sequences),
        'seq_length': len(sequences[0].seq) if sequences else 0,
        'taxa': []
    }
    
    # Extract dates from sequence names
    date_pattern = r'(\d{4}-\d{2}-\d{2}|\d{4})'
    for seq in sequences:
        date_match = re.search(date_pattern, seq.id)
        date = date_match.group(1) if date_match else 'Unknown'
        info['taxa'].append({
            'name': seq.id,
            'date': date,
            'length': len(seq.seq)
        })
    
    return info


def parse_loganalyser_output(loganalyser_file):
    """Parse the loganalyser output file"""
    results = []
    
    with open(loganalyser_file, 'r') as f:
        lines = f.readlines()
    
    # Find the table section
    in_table = False
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        # Split by whitespace
        parts = line.split()
        if len(parts) >= 5:
            # Assume format: statistic mean stderr median hpdLower hpdUpper ESS
            results.append({
                'statistic': parts[0],
                'mean': parts[1],
                'stderr': parts[2] if len(parts) > 2 else '-',
                'median': parts[3] if len(parts) > 3 else '-',
                'hpd_lower': parts[4] if len(parts) > 4 else '-',
                'hpd_upper': parts[5] if len(parts) > 5 else '-',
                'ess': parts[6] if len(parts) > 6 else '-'
            })
    
    return results


def read_file_content(filepath):
    """Safely read file content"""
    try:
        with open(filepath, 'r') as f:
            return f.read()
    except:
        return None


def generate_html_report(fasta_file, template_file, log_file, loganalyser_file, 
                        svg_file, output_file, run_info):
    """Generate comprehensive HTML report"""
    
    # Parse input information
    fasta_info = parse_fasta_info(fasta_file) if SeqIO else None
    
    # Parse loganalyser results
    log_results = parse_loganalyser_output(loganalyser_file)
    
    # Read SVG content
    svg_content = read_file_content(svg_file)
    
    # Start HTML
    html = []
    html.append('<!DOCTYPE html>')
    html.append('<html lang="en">')
    html.append('<head>')
    html.append('    <meta charset="UTF-8">')
    html.append('    <meta name="viewport" content="width=device-width, initial-scale=1.0">')
    html.append('    <title>BEAST Analysis Report</title>')
    html.append('    <style>')
    html.append('''
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .header h1 {
            margin: 0 0 10px 0;
            font-size: 2.5em;
        }
        .header .subtitle {
            opacity: 0.9;
            font-size: 1.1em;
        }
        .section {
            background: white;
            padding: 25px;
            margin-bottom: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .section h2 {
            color: #667eea;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
            margin-top: 0;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }
        th {
            background-color: #667eea;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }
        td {
            padding: 10px 12px;
            border-bottom: 1px solid #e0e0e0;
        }
        tr:hover {
            background-color: #f8f9ff;
        }
        .info-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 15px;
            margin: 15px 0;
        }
        .info-item {
            background: #f8f9ff;
            padding: 15px;
            border-radius: 5px;
            border-left: 4px solid #667eea;
        }
        .info-item .label {
            font-weight: 600;
            color: #666;
            font-size: 0.9em;
            text-transform: uppercase;
            margin-bottom: 5px;
        }
        .info-item .value {
            font-size: 1.3em;
            color: #333;
        }
        .tree-container {
            text-align: center;
            margin: 20px 0;
            padding: 20px;
            background: #fafafa;
            border-radius: 8px;
        }
        .tree-container svg {
            max-width: 100%;
            height: auto;
        }
        .footer {
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            color: #666;
            font-size: 0.9em;
        }
        .badge {
            display: inline-block;
            padding: 4px 8px;
            background: #4caf50;
            color: white;
            border-radius: 3px;
            font-size: 0.85em;
            font-weight: 600;
        }
        .badge.warning {
            background: #ff9800;
        }
        .badge.error {
            background: #f44336;
        }
    ''')
    html.append('    </style>')
    html.append('</head>')
    html.append('<body>')
    
    # Header
    html.append('    <div class="header">')
    html.append('        <h1>🧬 BEAST Analysis Report</h1>')
    html.append(f'        <div class="subtitle">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</div>')
    html.append('    </div>')
    
    # Input Data Section
    html.append('    <div class="section">')
    html.append('        <h2>📊 Input Data</h2>')
    
    if fasta_info:
        html.append('        <div class="info-grid">')
        html.append('            <div class="info-item">')
        html.append('                <div class="label">Number of Taxa</div>')
        html.append(f'                <div class="value">{fasta_info["num_taxa"]}</div>')
        html.append('            </div>')
        html.append('            <div class="info-item">')
        html.append('                <div class="label">Sequence Length</div>')
        html.append(f'                <div class="value">{fasta_info["seq_length"]} bp</div>')
        html.append('            </div>')
        html.append('            <div class="info-item">')
        html.append('                <div class="label">Template</div>')
        html.append(f'                <div class="value">{os.path.basename(template_file)}</div>')
        html.append('            </div>')
        html.append('        </div>')
        
        # Taxa table (if not too many)
        if fasta_info['num_taxa'] <= 50:
            html.append('        <h3>Taxa and Sampling Dates</h3>')
            html.append('        <table>')
            html.append('            <tr><th>Taxon</th><th>Date</th><th>Sequence Length</th></tr>')
            for taxon in fasta_info['taxa']:
                html.append(f'            <tr><td>{taxon["name"]}</td><td>{taxon["date"]}</td><td>{taxon["length"]} bp</td></tr>')
            html.append('        </table>')
    
    html.append('    </div>')
    
    # Run Information Section
    html.append('    <div class="section">')
    html.append('        <h2>⚙️ Analysis Details</h2>')
    html.append('        <div class="info-grid">')
    html.append('            <div class="info-item">')
    html.append('                <div class="label">Chain Length</div>')
    html.append(f'                <div class="value">{run_info.get("chain_length", "N/A"):,}</div>')
    html.append('            </div>')
    html.append('            <div class="info-item">')
    html.append('                <div class="label">Log Every</div>')
    html.append(f'                <div class="value">{run_info.get("log_every", "N/A"):,}</div>')
    html.append('            </div>')
    html.append('            <div class="info-item">')
    html.append('                <div class="label">Burn-in</div>')
    html.append(f'                <div class="value">{run_info.get("burnin", "N/A")}%</div>')
    html.append('            </div>')
    html.append('            <div class="info-item">')
    html.append('                <div class="label">Runtime</div>')
    html.append(f'                <div class="value">{run_info.get("runtime", "N/A")}</div>')
    html.append('            </div>')
    html.append('        </div>')
    html.append('    </div>')
    
    # Parameter Estimates Section
    html.append('    <div class="section">')
    html.append('        <h2>📈 Parameter Estimates</h2>')
    
    if log_results:
        html.append('        <table>')
        html.append('            <tr>')
        html.append('                <th>Parameter</th>')
        html.append('                <th>Mean</th>')
        html.append('                <th>Std Error</th>')
        html.append('                <th>Median</th>')
        html.append('                <th>95% HPD Lower</th>')
        html.append('                <th>95% HPD Upper</th>')
        html.append('                <th>ESS</th>')
        html.append('            </tr>')
        
        for result in log_results:
            # Check ESS and add warning badge if needed
            ess_cell = result['ess']
            try:
                ess_val = float(result['ess'])
                if ess_val < 200:
                    ess_cell = f'{result["ess"]} <span class="badge error">Low</span>'
                elif ess_val < 500:
                    ess_cell = f'{result["ess"]} <span class="badge warning">Fair</span>'
                else:
                    ess_cell = f'{result["ess"]} <span class="badge">Good</span>'
            except:
                pass
            
            html.append('            <tr>')
            html.append(f'                <td><strong>{result["statistic"]}</strong></td>')
            html.append(f'                <td>{result["mean"]}</td>')
            html.append(f'                <td>{result["stderr"]}</td>')
            html.append(f'                <td>{result["median"]}</td>')
            html.append(f'                <td>{result["hpd_lower"]}</td>')
            html.append(f'                <td>{result["hpd_upper"]}</td>')
            html.append(f'                <td>{ess_cell}</td>')
            html.append('            </tr>')
        
        html.append('        </table>')
    else:
        html.append('        <p>No parameter estimates available.</p>')
    
    html.append('    </div>')
    
    # Tree Visualization Section
    html.append('    <div class="section">')
    html.append('        <h2>🌳 Maximum Clade Credibility Tree</h2>')
    
    if svg_content:
        html.append('        <div class="tree-container">')
        html.append(svg_content)
        html.append('        </div>')
    else:
        html.append('        <p>Tree visualization not available.</p>')
    
    html.append('    </div>')
    
    # Footer
    html.append('    <div class="footer">')
    html.append('        <p>Generated by BEAST-NF Pipeline</p>')
    html.append('        <p>For more information, visit <a href="https://www.beast2.org/">beast2.org</a></p>')
    html.append('    </div>')
    
    html.append('</body>')
    html.append('</html>')
    
    # Write HTML file
    with open(output_file, 'w') as f:
        f.write('\n'.join(html))
    
    print(f"Report generated: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Generate HTML report for BEAST analysis')
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--template', required=True, help='BEAST template file')
    parser.add_argument('--log', required=True, help='BEAST log file')
    parser.add_argument('--loganalyser', required=True, help='Loganalyser output file')
    parser.add_argument('--svg', required=True, help='Tree SVG file')
    parser.add_argument('--output', required=True, help='Output HTML file')
    parser.add_argument('--chain-length', type=int, default=0, help='MCMC chain length')
    parser.add_argument('--log-every', type=int, default=0, help='Logging frequency')
    parser.add_argument('--burnin', type=int, default=0, help='Burn-in percentage')
    parser.add_argument('--runtime', default='N/A', help='Runtime')
    
    args = parser.parse_args()
    
    run_info = {
        'chain_length': args.chain_length,
        'log_every': args.log_every,
        'burnin': args.burnin,
        'runtime': args.runtime
    }
    
    generate_html_report(
        args.fasta,
        args.template,
        args.log,
        args.loganalyser,
        args.svg,
        args.output,
        run_info
    )


if __name__ == '__main__':
    main()
