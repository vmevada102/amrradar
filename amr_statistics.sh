#!/usr/bin/env bash
#
# amr_statistics.sh  -- simple, robust runner (no conda activation)
#
# Usage:
#   ./amr_statistics.sh --run-python [path/to/One_Health_AMR_analysis.py]
#   ./amr_statistics.sh --run-python          # auto-detect script in same dir
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUMMARY_XLSX="${SCRIPT_DIR}/results/OneHealthAMR_AMRFinder_summary.xlsx"
PY_SCRIPT=""
RUN_PY=false

# parse args (allow --run-python with or without an argument)
if [[ $# -gt 0 ]]; then
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --run-python) RUN_PY=true; shift; if [[ $# -gt 0 && ! "$1" =~ ^- ]]; then PY_SCRIPT="$1"; shift; fi ;;
      --help|-h) cat <<'USAGE'
Usage:
  ./amr_statistics.sh --run-python [path/to/One_Health_AMR_analysis.py]
USAGE
        exit 0 ;;
      *) echo "Unknown argument: $1" >&2; shift ;;
    esac
  done
fi

# input file check
if [[ ! -f "$SUMMARY_XLSX" ]]; then
  echo "ERROR: Summary Excel not found at: $SUMMARY_XLSX" >&2
  exit 2
fi
echo "Using summary file: $SUMMARY_XLSX"

# determine next numbered output folder
BASE_DIR="${SCRIPT_DIR}/results/tools"
mkdir -p "$BASE_DIR"
max_num=0
while IFS= read -r -d '' d; do
  name="$(basename "$d")"
  if [[ $name =~ ^([0-9]+)_ ]]; then
    num=${BASH_REMATCH[1]}
    if (( num > max_num )); then max_num=$num; fi
  fi
done < <(find "$BASE_DIR" -maxdepth 1 -mindepth 1 -type d -print0 2>/dev/null || true)
N=$((max_num + 1))
OUTDIR="${BASE_DIR}/${N}_AMRStatistics"
mkdir -p "$OUTDIR"
echo "Created output directory: $OUTDIR"

# pick python
PY_BIN="$(command -v python3 || command -v python || true)"
if [[ -z "$PY_BIN" ]]; then
  echo "ERROR: No python or python3 found in PATH." >&2
  exit 4
fi
echo "Using python: $PY_BIN"

# create log
LOG="$OUTDIR/analysis.log"
: > "$LOG"
echo "Wrapper started at: $(date -u)" >> "$LOG"
echo "Python: $PY_BIN" >> "$LOG"
echo "OUTDIR: $OUTDIR" >> "$LOG"
echo "" >> "$LOG"

if [[ "$RUN_PY" != true ]]; then
  echo "Prepared output folder: $OUTDIR"
  echo "To run: ./amr_statistics.sh --run-python [path/to/One_Health_AMR_analysis.py]"
  exit 0
fi

# resolve python script
if [[ -z "$PY_SCRIPT" ]]; then
  if [[ -f "${SCRIPT_DIR}/One_Health_AMR_analysis.py" ]]; then
    PY_SCRIPT="${SCRIPT_DIR}/One_Health_AMR_analysis.py"
  else
    echo "ERROR: --run-python requested but Python script not found." | tee -a "$LOG"
    exit 6
  fi
fi
if [[ ! -f "$PY_SCRIPT" ]]; then
  echo "ERROR: Python analysis script not found at: $PY_SCRIPT" | tee -a "$LOG"
  exit 6
fi

# run unbuffered
echo "Running command:" | tee -a "$LOG"
echo "\"$PY_BIN\" -u \"$PY_SCRIPT\" \"$SUMMARY_XLSX\" \"$OUTDIR\"" | tee -a "$LOG"
if "$PY_BIN" -u "$PY_SCRIPT" "$SUMMARY_XLSX" "$OUTDIR" >> "$LOG" 2>&1; then
  echo "Python analysis finished successfully at: $(date -u)" | tee -a "$LOG"
  echo "" | tee -a "$LOG"
  echo "Files created:" | tee -a "$LOG"
  ls -1 "$OUTDIR" | tee -a "$LOG"
  echo "" | tee -a "$LOG"
  if [[ -f "$OUTDIR/run_summary.txt" ]]; then
    echo "---- run_summary.txt ----" | tee -a "$LOG"
    sed -n '1,200p' "$OUTDIR/run_summary.txt" | tee -a "$LOG"
    echo "-------------------------" | tee -a "$LOG"
    # print it to terminal as well (automatic summary)
    echo ""
    echo "SUMMARY (from run_summary.txt):"
    cat "$OUTDIR/run_summary.txt"
    echo ""
  fi
  exit 0
else
  rc=$?
  echo "Python analysis failed with exit code $rc at: $(date -u)" | tee -a "$LOG"
  echo "---- last 200 lines of $LOG ----"
  tail -n 200 "$LOG" || true
  exit $rc
fi
