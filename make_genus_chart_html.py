# Generate interactive HTML dashboard for results/Identified_organisms_genus.tsv

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
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram as scipy_dendrogram

# ==============================
#   FILE PATHS / CONFIG
# ==============================
INPUT_TSV = os.path.join("results", "Identified_organisms_genus.tsv")
OUTPUT_DIR = os.path.join("results", "Charts")
OUTPUT_HTML = os.path.join(OUTPUT_DIR, "Identified_organisms_genus.html")

# Max number of genera to include in heatmap and clustering (top by mean percentage)
HEATMAP_TOP_N = 50


# ==============================
#   MAIN PROCESS
# ==============================
def main():
    # Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load TSV (UTF-8; we already know your files are UTF-8)
    print("Reading genus TSV file as UTF-8:", INPUT_TSV)
    df = pd.read_csv(INPUT_TSV, sep="\t", encoding="utf-8", engine="python")

    # Optional: if percent_basis has multiple types, keep only "total"
    if "percent_basis" in df.columns:
        if df["percent_basis"].nunique() > 1:
            print("Filtering rows where percent_basis == 'total'")
            df = df[df["percent_basis"] == "total"].copy()

    # Preserve sample order as in file
    samples = list(dict.fromkeys(df["sample"].tolist()))

    # Sample -> group mapping
    sample_groups = {}
    for _, row in df.iterrows():
        s = row["sample"]
        g = row["group"]
        if s not in sample_groups:
            sample_groups[s] = g

    # Genus -> sample -> percentage mapping
    genus_sample_pct = {}
    for _, row in df.iterrows():
        sample = row["sample"]
        genus = row["genus"]
        pct = float(row["percent"])
        if genus not in genus_sample_pct:
            genus_sample_pct[genus] = {}
        # Sum if multiple rows per genus/sample exist
        genus_sample_pct[genus][sample] = genus_sample_pct[genus].get(sample, 0.0) + pct

    genera_all = sorted(genus_sample_pct.keys())

    # Compute mean percentage per genus for ranking
    num_samples = len(samples)
    genus_mean = {}
    for genus, sdict in genus_sample_pct.items():
        total_pct = sum(sdict.values())
        mean_pct = total_pct / float(num_samples) if num_samples > 0 else 0.0
        genus_mean[genus] = mean_pct

    # Genera sorted by mean percentage descending
    genera_sorted = [
        genus for genus, _ in sorted(genus_mean.items(), key=lambda x: x[1], reverse=True)
    ]

    # Group values for dropdown
    groups = ["ALL"] + sorted(set(sample_groups.values()))

    # ==============================
    #   HEATMAP + CLUSTERING
    # ==============================
    heatmap_genera_base = genera_sorted[:HEATMAP_TOP_N]
    heatmap_matrix_base = [
        [genus_sample_pct.get(genus, {}).get(s, 0.0) for s in samples]
        for genus in heatmap_genera_base
    ]

    heatmap_genera = heatmap_genera_base
    heatmap_matrix = heatmap_matrix_base
    dendro_js = None

    if len(heatmap_genera_base) >= 2:
        # Perform hierarchical clustering
        X = np.array(heatmap_matrix_base)
        Z = linkage(X, method="average", metric="euclidean")
        ddata = scipy_dendrogram(
            Z, labels=heatmap_genera_base, orientation="right", no_plot=True
        )

        icoord = ddata["icoord"]
        dcoord = ddata["dcoord"]
        leaves = ddata["leaves"]

        # Map leaf index to y-position (same logic as scipy dendrogram)
        leaf_y = {}
        for i, leaf_idx in enumerate(leaves):
            y = 5 + 10 * i
            label = heatmap_genera_base[leaf_idx]
            leaf_y[label] = y

        # Sort by y to get tick order
        tick_items = sorted(leaf_y.items(), key=lambda kv: kv[1])
        tickvals = [v for _, v in tick_items]
        ticktext = [k for k, _ in tick_items]

        # Reorder heatmap rows to match dendrogram leaf order
        genus_index = {genus: i for i, genus in enumerate(heatmap_genera_base)}
        heatmap_genera = ticktext
        heatmap_matrix = [heatmap_matrix_base[genus_index[g]] for g in heatmap_genera]

        dendro_js = {
            "icoord": icoord,
            "dcoord": dcoord,
            "tickvals": tickvals,
            "ticktext": ticktext,
        }

    # Pack data for JavaScript
    data_js = {
        "samples": samples,
        "data": genus_sample_pct,
        "genera_all": genera_all,
        "genera_sorted": genera_sorted,
        "sample_groups": sample_groups,
        "groups": groups,
        "heatmap": {
            "genera": heatmap_genera,
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
    <title>Interactive Genus Percentage Dashboard</title>
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
        #genusSelect {{
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
    <h2>Interactive Genus Percentage per Sample</h2>

    <div class="controls">
        <div class="control-block">
            <label for="groupSelect">Group:</label>
            <select id="groupSelect"></select>
        </div>

        <div class="control-block">
            <label for="viewSelect">View:</label>
            <select id="viewSelect">
                <option value="bar" selected>Bar chart</option>
                <option value="heatmap">Heatmap + dendrogram (up to top {HEATMAP_TOP_N})</option>
            </select>
        </div>

        <div class="control-block">
            <label for="topNSelect">Top N genera (bar view):</label>
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
            <label for="genusSelect">Genera:</label>
            <select id="genusSelect" multiple size="5"></select>
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

        <div class="control-block">
            <label for="heatmapNSelect">Heatmap top N genera:</label>
            <select id="heatmapNSelect">
                <option value="10">10</option>
                <option value="20" selected>20</option>
                <option value="30">30</option>
                <option value="50">50</option>
                <option value="ALL">ALL (max {HEATMAP_TOP_N})</option>
            </select>
        </div>
    </div>

    <div id="barContainer">
        <div id="chart"></div>
    </div>

    <div id="heatmapContainer" style="display:none;">
        <h3>Hierarchical clustering (dendrogram)</h3>
        <div id="dendrogramDiv"></div>
        <h3>Heatmap (genus x sample)</h3>
        <div id="heatmapDiv"></div>
    </div>

    <script>
    const chartData = {data_json};

    const samples = chartData.samples;
    const generaAll = chartData.genera_all;
    const generaSorted = chartData.genera_sorted;
    const data = chartData.data;
    const sampleGroups = chartData.sample_groups;
    const groups = chartData.groups;
    const heatmapData = chartData.heatmap;
    const dendrogramData = chartData.dendrogram;
    const heatmapTopNMax = chartData.heatmap_top_n;

    const groupSelect = document.getElementById('groupSelect');
    const viewSelect = document.getElementById('viewSelect');
    const topNSelect = document.getElementById('topNSelect');
    const genusSelect = document.getElementById('genusSelect');
    const thresholdSlider = document.getElementById('thresholdSlider');
    const thresholdValueSpan = document.getElementById('thresholdValue');
    const logToggle = document.getElementById('logToggle');
    const heatmapNSelect = document.getElementById('heatmapNSelect');

    const barContainer = document.getElementById('barContainer');
    const heatmapContainer = document.getElementById('heatmapContainer');

    // Populate group dropdown
    groups.forEach(function(g) {{
        const opt = document.createElement('option');
        opt.value = g;
        opt.textContent = (g === 'ALL') ? 'ALL groups' : g;
        groupSelect.appendChild(opt);
    }});

    function populateGenusList() {{
        genusSelect.innerHTML = "";
        const N = topNSelect.value;
        let list;
        if (N === "ALL") {{
            list = generaSorted.slice();
        }} else {{
            const nInt = parseInt(N);
            list = generaSorted.slice(0, nInt);
        }}
        list.forEach(function(genus) {{
            const opt = document.createElement('option');
            opt.value = genus;
            opt.textContent = genus;
            genusSelect.appendChild(opt);
        }});
    }}

    function getSelectedGenera() {{
        return Array.from(genusSelect.selectedOptions).map(function(o) {{ return o.value; }});
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
        const selectedGenera = getSelectedGenera();
        const filteredSamples = getFilteredSamples();
        const threshold = parseFloat(thresholdSlider.value);
        const isLog = logToggle.checked;

        const traces = [];

        selectedGenera.forEach(function(genus) {{
            const genusData = data[genus] || {{}};
            const x = filteredSamples;
            const y = filteredSamples.map(function(s) {{
                let v = genusData[s] || 0;
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
                name: genus
            }});
        }});

        const layout = {{
            title: "Genus Percentage per Sample (stacked)",
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
                title: "Dendrogram not available (insufficient genera)"
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
        const Nopt = heatmapNSelect.value;

        const hmGeneraAll = heatmapData.genera;
        const hmSamples = heatmapData.samples;
        const baseMatrixAll = heatmapData.matrix;

        let N;
        if (Nopt === "ALL") {{
            N = hmGeneraAll.length;
        }} else {{
            N = Math.min(parseInt(Nopt), hmGeneraAll.length);
        }}

        const hmGenera = hmGeneraAll.slice(0, N);
        const baseMatrix = baseMatrixAll.slice(0, N);

        const z = baseMatrix.map(function(row) {{
            return row.map(function(v) {{
                return (v >= threshold) ? v : 0;
            }});
        }});

        const trace = {{
            x: hmSamples,
            y: hmGenera,
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
        populateGenusList();
        // Auto-select first up to 3 genera when list changes
        for (let i = 0; i < genusSelect.options.length; i++) {{
            genusSelect.options[i].selected = i < 3;
        }}
        updateView();
    }});

    genusSelect.addEventListener("change", updateView);
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
    heatmapNSelect.addEventListener("change", function() {{
        if (viewSelect.value === "heatmap") {{
            renderHeatmap();
        }}
    }});

    // Initial setup
    function initialize() {{
        groupSelect.value = "ALL";
        populateGenusList();
        // Select first up to 3 genera
        for (let i = 0; i < Math.min(3, genusSelect.options.length); i++) {{
            genusSelect.options[i].selected = true;
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

    print("Genus HTML file created at: {}".format(OUTPUT_HTML))


if __name__ == "__main__":
    main()
