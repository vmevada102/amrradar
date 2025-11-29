# Generate interactive HTML dashboard for results/Identified_organisms.tsv

# ==============================
#   AUTO-INSTALL DEPENDENCIES
# ==============================
required_packages = ["pandas", "numpy", "scipy"]
import subprocess
import sys

for pkg in required_packages:
    try:
        __import__(pkg)
    except ImportError:
        print("Installing missing package: {}".format(pkg))
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])

# ==============================
#   IMPORTS
# ==============================
import os
import json
import re
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram as scipy_dendrogram

# ==============================
#   FILE PATHS / CONFIG
# ==============================
INPUT_TSV = os.path.join("results", "Identified_organisms.tsv")
OUTPUT_DIR = os.path.join("results", "Charts")
OUTPUT_HTML = os.path.join(OUTPUT_DIR, "Identified_organisms.html")

# Number of organisms to include in heatmap and clustering (top by mean percentage)
HEATMAP_TOP_N = 50


# ==============================
#   PARSING FUNCTION
# ==============================
def parse_organisms(detail):
    """
    Parse 'identified_organism_detail' strings like:
    'Homo sapiens (1613, 0.23%) Pseudomonas aeruginosa (xxx, yyy%) ...'
    into a dict: { 'Homo sapiens': 0.23, 'Pseudomonas aeruginosa': yyy, ... }
    """
    organisms = {}
    if pd.isna(detail):
        return organisms

    pattern = r'([A-Za-z\s\.\-]+)\s*\(([\d,\.]+\s*,\s*([\d\.]+)%)\)'
    matches = re.findall(pattern, str(detail))

    for name, _, pct in matches:
        name = name.strip()
        try:
            organisms[name] = float(pct)
        except ValueError:
            # Skip values that cannot be parsed as float
            pass

    return organisms


# ==============================
#   MAIN PROCESS
# ==============================
def main():
    # Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load TSV using UTF-8 (we already tested that file is UTF-8)
    print("Reading TSV file as UTF-8:", INPUT_TSV)
    df = pd.read_csv(INPUT_TSV, sep="\t", encoding="utf-8", engine="python")

    # Parse organism details
    df["parsed"] = df["identified_organism_detail"].apply(parse_organisms)

    # Preserve sample order as in file
    samples = list(dict.fromkeys(df["sample"].tolist()))

    # Sample -> group mapping
    sample_groups = {}
    for _, row in df.iterrows():
        s = row["sample"]
        g = row["group"]
        if s not in sample_groups:
            sample_groups[s] = g

    # Organism -> sample -> percentage mapping
    organism_sample_pct = {}
    for _, row in df.iterrows():
        sample = row["sample"]
        for org, pct in row["parsed"].items():
            organism_sample_pct.setdefault(org, {})[sample] = pct

    organisms_all = sorted(organism_sample_pct.keys())

    # Compute mean percentage per organism for ranking
    num_samples = len(samples)
    org_mean = {}
    for org, sdict in organism_sample_pct.items():
        total_pct = sum(sdict.values())
        mean_pct = total_pct / float(num_samples) if num_samples > 0 else 0.0
        org_mean[org] = mean_pct

    # Organisms sorted by mean percentage descending
    organisms_sorted = [
        org for org, _ in sorted(org_mean.items(), key=lambda x: x[1], reverse=True)
    ]

    # Group values for dropdown
    groups = ["ALL"] + sorted(set(sample_groups.values()))

    # ==============================
    #   HEATMAP + CLUSTERING
    # ==============================
    heatmap_orgs_base = organisms_sorted[:HEATMAP_TOP_N]
    heatmap_matrix_base = [
        [organism_sample_pct.get(org, {}).get(s, 0.0) for s in samples]
        for org in heatmap_orgs_base
    ]

    heatmap_organisms = heatmap_orgs_base
    heatmap_matrix = heatmap_matrix_base
    dendro_js = None

    if len(heatmap_orgs_base) >= 2:
        # Perform hierarchical clustering
        X = np.array(heatmap_matrix_base)
        Z = linkage(X, method="average", metric="euclidean")
        ddata = scipy_dendrogram(
            Z, labels=heatmap_orgs_base, orientation="right", no_plot=True
        )

        icoord = ddata["icoord"]
        dcoord = ddata["dcoord"]
        leaves = ddata["leaves"]

        # Map leaf index to y-position (same logic as scipy dendrogram)
        leaf_y = {}
        for i, leaf_idx in enumerate(leaves):
            y = 5 + 10 * i
            label = heatmap_orgs_base[leaf_idx]
            leaf_y[label] = y

        # Sort by y to get tick order
        tick_items = sorted(leaf_y.items(), key=lambda kv: kv[1])
        tickvals = [v for _, v in tick_items]
        ticktext = [k for k, _ in tick_items]

        # Reorder heatmap rows to match dendrogram leaf order
        org_index = {org: i for i, org in enumerate(heatmap_orgs_base)}
        heatmap_organisms = ticktext
        heatmap_matrix = [heatmap_matrix_base[org_index[o]] for o in heatmap_organisms]

        dendro_js = {
            "icoord": icoord,
            "dcoord": dcoord,
            "tickvals": tickvals,
            "ticktext": ticktext,
        }

    # Pack data for JavaScript
    data_js = {
        "samples": samples,
        "data": organism_sample_pct,
        "organisms_all": organisms_all,
        "organisms_sorted": organisms_sorted,
        "sample_groups": sample_groups,
        "groups": groups,
        "heatmap": {
            "organisms": heatmap_organisms,
            "samples": samples,
            "matrix": heatmap_matrix,
        },
        "dendrogram": dendro_js,
        "heatmap_top_n": HEATMAP_TOP_N,
    }

    data_json = json.dumps(data_js, ensure_ascii=True)

    # ==============================
    #   BUILD HTML CONTENT
    # ==============================
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>Interactive Organism Percentage Dashboard</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
        }}
        .controls {{
            margin-bottom: 20px;
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
            align-items: center;
        }}
        .control-block {{
            display: flex;
            align-items: center;
            gap: 6px;
        }}
        select, input[type="range"] {{
            padding: 4px 8px;
            font-size: 14px;
        }}
        #organismSelect {{
            min-width: 260px;
        }}
        #barContainer, #heatmapContainer {{
            width: 100%;
        }}
        #chart {{
            width: 100%;
            height: 650px;
        }}
        #heatmapDiv {{
            width: 100%;
            height: 700px;
        }}
        #dendrogramDiv {{
            width: 100%;
            height: 400px;
        }}
        small {{
            color: #555;
        }}
    </style>
</head>
<body>
    <h2>Interactive Organism Percentage per Sample</h2>

    <div class="controls">
        <div class="control-block">
            <label for="groupSelect">Group:</label>
            <select id="groupSelect"></select>
        </div>

        <div class="control-block">
            <label for="viewSelect">View:</label>
            <select id="viewSelect">
                <option value="bar" selected>Bar chart</option>
                <option value="heatmap">Heatmap + dendrogram (top {HEATMAP_TOP_N})</option>
            </select>
        </div>

        <div class="control-block">
            <label for="topNSelect">Top N organisms (bar view):</label>
            <select id="topNSelect">
                <option value="5">5</option>
                <option value="10" selected>10</option>
                <option value="20">20</option>
                <option value="30">30</option>
                <option value="50">50</option>
                <option value="ALL">ALL</option>
            </select>
        </div>

        <div class="control-block">
            <label for="organismSelect">Organisms:</label>
            <select id="organismSelect" multiple size="5"></select>
            <small>(Use Ctrl/Shift to select multiple)</small>
        </div>

        <div class="control-block">
            <label for="thresholdSlider">Min %:</label>
            <input id="thresholdSlider" type="range" min="0" max="5" step="0.01" value="0.00">
            <span id="thresholdValue">0.00%</span>
        </div>

        <div class="control-block">
            <label><input type="checkbox" id="logToggle"> Log10 Y-axis (bar view)</label>
        </div>
    </div>

    <div id="barContainer">
        <div id="chart"></div>
    </div>

    <div id="heatmapContainer" style="display:none;">
        <h3>Hierarchical clustering (dendrogram)</h3>
        <div id="dendrogramDiv"></div>
        <h3>Heatmap (organism x sample)</h3>
        <div id="heatmapDiv"></div>
    </div>

    <script>
    const chartData = {data_json};

    const samples = chartData.samples;
    const organismsAll = chartData.organisms_all;
    const organismsSorted = chartData.organisms_sorted;
    const data = chartData.data;
    const sampleGroups = chartData.sample_groups;
    const groups = chartData.groups;
    const heatmapData = chartData.heatmap;
    const dendrogramData = chartData.dendrogram;
    const heatmapTopN = chartData.heatmap_top_n;

    const groupSelect = document.getElementById('groupSelect');
    const viewSelect = document.getElementById('viewSelect');
    const topNSelect = document.getElementById('topNSelect');
    const organismSelect = document.getElementById('organismSelect');
    const thresholdSlider = document.getElementById('thresholdSlider');
    const thresholdValueSpan = document.getElementById('thresholdValue');
    const logToggle = document.getElementById('logToggle');

    const barContainer = document.getElementById('barContainer');
    const heatmapContainer = document.getElementById('heatmapContainer');

    // Populate group dropdown
    groups.forEach(function(g) {{
        const opt = document.createElement('option');
        opt.value = g;
        opt.textContent = (g === 'ALL') ? 'ALL groups' : g;
        groupSelect.appendChild(opt);
    }});

    function populateOrganismList() {{
        organismSelect.innerHTML = "";
        const N = topNSelect.value;
        let list;
        if (N === "ALL") {{
            list = organismsSorted.slice();
        }} else {{
            const nInt = parseInt(N);
            list = organismsSorted.slice(0, nInt);
        }}
        list.forEach(function(org) {{
            const opt = document.createElement('option');
            opt.value = org;
            opt.textContent = org;
            organismSelect.appendChild(opt);
        }});
    }}

    function getSelectedOrganisms() {{
        return Array.from(organismSelect.selectedOptions).map(function(o) {{ return o.value; }});
    }}

    function getFilteredSamples() {{
        const group = groupSelect.value;
        if (group === "ALL") {{
            return samples.slice();
        }}
        return samples.filter(function(s) {{ return sampleGroups[s] === group; }});
    }}

    function updateThresholdLabel() {{
        const v = parseFloat(thresholdSlider.value);
        thresholdValueSpan.textContent = v.toFixed(2) + "%";
    }}

    function updateBarChart() {{
        const selectedOrgs = getSelectedOrganisms();
        const filteredSamples = getFilteredSamples();
        const threshold = parseFloat(thresholdSlider.value);
        const isLog = logToggle.checked;

        const traces = [];

        selectedOrgs.forEach(function(org) {{
            const orgData = data[org] || {{}};
            const x = filteredSamples;
            const y = filteredSamples.map(function(s) {{
                let v = orgData[s] || 0;
                if (v < threshold) {{
                    v = 0;
                }}
                if (isLog && v === 0) {{
                    return null;
                }}
                return v;
            }});

            traces.push({{
                x: x,
                y: y,
                type: "bar",
                name: org
            }});
        }});

        const layout = {{
            title: "Organism Percentage per Sample (stacked)",
            barmode: "stack",
            xaxis: {{
                title: "Samples",
                tickangle: -45,
                automargin: true
            }},
            yaxis: {{
                title: "Percentage",
                type: isLog ? "log" : "linear",
                rangemode: "tozero"
            }},
            margin: {{
                b: 160
            }},
            legend: {{
                orientation: "h",
                y: -0.3
            }}
        }};

        Plotly.react("chart", traces, layout, {{responsive: true}});
    }}

    function renderDendrogram() {{
        if (!dendrogramData) {{
            Plotly.react("dendrogramDiv", [], {{
                title: "Dendrogram not available (insufficient organisms)"
            }});
            return;
        }}

        const icoord = dendrogramData.icoord;
        const dcoord = dendrogramData.dcoord;
        const tickvals = dendrogramData.tickvals;
        const ticktext = dendrogramData.ticktext;

        const traces = [];
        for (let i = 0; i < icoord.length; i++) {{
            traces.push({{
                x: dcoord[i],
                y: icoord[i],
                mode: "lines",
                type: "scatter",
                line: {{width: 1}}
            }});
        }}

        const layout = {{
            xaxis: {{
                showticklabels: false,
                zeroline: false,
                showgrid: false
            }},
            yaxis: {{
                tickvals: tickvals,
                ticktext: ticktext,
                showgrid: false,
                zeroline: false,
                automargin: true
            }},
            margin: {{
                l: 150,
                r: 20,
                t: 20,
                b: 20
            }}
        }};

        Plotly.react("dendrogramDiv", traces, layout, {{responsive: true}});
    }}

    function renderHeatmap() {{
        const threshold = parseFloat(thresholdSlider.value);
        const hmOrgs = heatmapData.organisms;
        const hmSamples = heatmapData.samples;
        const baseMatrix = heatmapData.matrix;

        const z = baseMatrix.map(function(row) {{
            return row.map(function(v) {{
                return (v >= threshold) ? v : 0;
            }});
        }});

        const trace = {{
            x: hmSamples,
            y: hmOrgs,
            z: z,
            type: "heatmap",
            colorbar: {{ title: "%" }}
        }};

        const layout = {{
            xaxis: {{
                tickangle: -45,
                automargin: true
            }},
            yaxis: {{
                automargin: true
            }},
            margin: {{
                l: 150,
                r: 20,
                t: 20,
                b: 160
            }}
        }};

        Plotly.react("heatmapDiv", [trace], layout, {{responsive: true}});
    }}

    function updateView() {{
        const view = viewSelect.value;
        if (view === "bar") {{
            barContainer.style.display = "";
            heatmapContainer.style.display = "none";
            updateBarChart();
        }} else {{
            barContainer.style.display = "none";
            heatmapContainer.style.display = "";
            renderDendrogram();
            renderHeatmap();
        }}
    }}

    // Event listeners
    topNSelect.addEventListener("change", function() {{
        populateOrganismList();
        // Auto-select first up to 3 organisms when list changes
        for (let i = 0; i < organismSelect.options.length; i++) {{
            organismSelect.options[i].selected = i < 3;
        }}
        updateView();
    }});

    organismSelect.addEventListener("change", updateView);
    groupSelect.addEventListener("change", updateView);

    thresholdSlider.addEventListener("input", function() {{
        updateThresholdLabel();
        updateView();
    }});

    logToggle.addEventListener("change", function() {{
        if (viewSelect.value === "bar") {{
            updateBarChart();
        }}
    }});

    viewSelect.addEventListener("change", updateView);

    // Initial setup
    function initialize() {{
        groupSelect.value = "ALL";
        populateOrganismList();
        // Select first up to 3 organisms
        for (let i = 0; i < Math.min(3, organismSelect.options.length); i++) {{
            organismSelect.options[i].selected = true;
        }}
        updateThresholdLabel();
        updateView();
    }}

    initialize();
    </script>
</body>
</html>
"""

    # Write HTML file
    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html_content)

    print("HTML file created at: {}".format(OUTPUT_HTML))


if __name__ == "__main__":
    main()
