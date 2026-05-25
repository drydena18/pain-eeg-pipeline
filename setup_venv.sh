#!/usr/bin/env bash
# =============================================================================
# setup_venv.sh — pain-alpha-dynamics Python environment setup
#
# Usage:
#   bash setup_venv.sh                  # standard GPU/desktop install
#   bash setup_venv.sh --headless       # HPC without GPU (Mesa offscreen)
#   bash setup_venv.sh --python /path/to/python3.11
#
# What this does:
#   1. Creates a venv at ~/envs/thesis  (override with THESIS_VENV)
#   2. Upgrades pip + setuptools
#   3. Installs requirements.txt
#   4. Installs Mesa pyvista variant on --headless
#   5. Smoke-tests the key imports
#   6. Prints the python path to paste into expXX.json
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
VENV_DIR="${THESIS_VENV:-$HOME/envs/thesis}"
PYTHON_BIN="python3"
HEADLESS=0
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REQ_FILE="$SCRIPT_DIR/requirements.txt"

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --headless)   HEADLESS=1; shift ;;
        --python)     PYTHON_BIN="$2"; shift 2 ;;
        --venv)       VENV_DIR="$2"; shift 2 ;;
        *)            echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# ── Preflight checks ──────────────────────────────────────────────────────────
echo ""
echo "=========================================="
echo "  pain-alpha-dynamics venv setup"
echo "=========================================="

if ! command -v "$PYTHON_BIN" &>/dev/null; then
    echo "ERROR: Python not found at '$PYTHON_BIN'"
    echo "  Try: bash setup_venv.sh --python /usr/bin/python3.11"
    exit 1
fi

PY_VERSION=$("$PYTHON_BIN" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PY_MAJOR=$(echo "$PY_VERSION" | cut -d. -f1)
PY_MINOR=$(echo "$PY_VERSION" | cut -d. -f2)

echo "Python binary  : $("$PYTHON_BIN" --version)"
echo "Venv target    : $VENV_DIR"
echo "Requirements   : $REQ_FILE"
echo "Headless/Mesa  : $([[ $HEADLESS -eq 1 ]] && echo yes || echo no)"
echo ""

if [[ $PY_MAJOR -lt 3 || ($PY_MAJOR -eq 3 && $PY_MINOR -lt 10) ]]; then
    echo "ERROR: Python >= 3.10 required (got $PY_VERSION)"
    exit 1
fi

if [[ ! -f "$REQ_FILE" ]]; then
    echo "ERROR: requirements.txt not found at $REQ_FILE"
    exit 1
fi

# ── Create venv ───────────────────────────────────────────────────────────────
if [[ -d "$VENV_DIR" ]]; then
    echo "Venv already exists at $VENV_DIR"
    read -rp "  Re-use existing venv? [Y/n] " REPLY
    REPLY="${REPLY:-Y}"
    if [[ "$REPLY" =~ ^[Nn] ]]; then
        echo "  Removing existing venv..."
        rm -rf "$VENV_DIR"
        "$PYTHON_BIN" -m venv "$VENV_DIR"
    fi
else
    echo "Creating venv at $VENV_DIR ..."
    "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

VENV_PYTHON="$VENV_DIR/bin/python"
VENV_PIP="$VENV_DIR/bin/pip"

# ── Upgrade pip ───────────────────────────────────────────────────────────────
echo ""
echo "Upgrading pip + setuptools..."
"$VENV_PIP" install --quiet --upgrade pip setuptools wheel

# ── Install requirements ──────────────────────────────────────────────────────
echo "Installing requirements.txt..."
if [[ $HEADLESS -eq 1 ]]; then
    # On headless HPC: replace plain pyvista with the OSMesa-enabled build,
    # and skip pyvistaqt which requires a Qt display.
    echo "  [headless] Filtering out pyvistaqt; installing pyvista[osmesa] instead..."
    grep -v "^pyvistaqt" "$REQ_FILE" > /tmp/req_headless.txt
    "$VENV_PIP" install --quiet -r /tmp/req_headless.txt
    "$VENV_PIP" install --quiet "pyvista[osmesa]>=0.44"
    rm /tmp/req_headless.txt
else
    "$VENV_PIP" install --quiet -r "$REQ_FILE"
fi

# ── Smoke tests ───────────────────────────────────────────────────────────────
echo ""
echo "Running smoke tests..."

SMOKE_RESULT=$("$VENV_PYTHON" - <<'PYEOF'
import sys
failures = []

tests = [
    ("numpy",      "import numpy as np; print(f'  numpy      {np.__version__}')"),
    ("scipy",      "import scipy; print(f'  scipy      {scipy.__version__}')"),
    ("pandas",     "import pandas as pd; print(f'  pandas     {pd.__version__}')"),
    ("matplotlib", "import matplotlib; print(f'  matplotlib {matplotlib.__version__}')"),
    ("mne",        "import mne; print(f'  mne        {mne.__version__}')"),
    ("specparam",  "import specparam; print(f'  specparam  {specparam.__version__}')"),
    ("imageio",    "import imageio; print(f'  imageio    {imageio.__version__}')"),
    ("pyvista",    "import pyvista; print(f'  pyvista    {pyvista.__version__}')"),
]

for name, code in tests:
    try:
        exec(code)
    except Exception as e:
        failures.append(f"  FAIL {name}: {e}")

if failures:
    print("FAILURES:")
    for f in failures:
        print(f)
    sys.exit(1)
else:
    print("  All imports OK")
PYEOF
)
echo "$SMOKE_RESULT"
SMOKE_EXIT=$?

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=========================================="
if [[ $SMOKE_EXIT -eq 0 ]]; then
    echo "  Setup complete."
else
    echo "  Setup finished with import failures above."
    echo "  Check the error messages and re-run if needed."
fi
echo ""
echo "  Python path (paste into exp01.json):"
echo "  \"python_exe\": \"$VENV_PYTHON\""
echo ""
echo "  To activate:"
echo "  source $VENV_DIR/bin/activate"
echo ""
if [[ $HEADLESS -eq 1 ]]; then
    echo "  Headless HPC reminder:"
    echo "  Set cfg.source.render.use_mesa = true in your JSON config."
    echo "  If OSMesa is not installed system-wide, ask your sysadmin:"
    echo "    apt install libosmesa6-dev   (Debian/Ubuntu)"
    echo "    module load mesa             (if using environment modules)"
    echo ""
fi
echo "=========================================="