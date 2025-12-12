#!/bin/bash
# WASP2 Paper - Generate All Figures
#
# Usage: ./run_all.sh [figure_number]
#   ./run_all.sh      # Generate all figures
#   ./run_all.sh 1    # Generate only Figure 1
#   ./run_all.sh 2    # Generate only Figure 2
#   ./run_all.sh 3    # Generate only Figure 3

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "WASP2 Paper Figure Generation"
echo "=========================================="

generate_figure1() {
    echo ""
    echo "Generating Figure 1: Read Mapping..."
    python figure1/scripts/generate_figure1.py
}

generate_figure2() {
    echo ""
    echo "Generating Figure 2: Counting Performance..."
    python figure2/scripts/generate_figure2.py
}

generate_figure3() {
    echo ""
    echo "Generating Figure 3: Statistical Analysis..."
    python figure3/scripts/generate_figure3.py
}

generate_figure4() {
    echo ""
    echo "Figure 4: Single Cell (not yet implemented)"
    # python figure4/scripts/generate_figure4.py
}

# Main
if [ $# -eq 0 ]; then
    # Generate all figures
    generate_figure1
    generate_figure2
    generate_figure3
    generate_figure4
    echo ""
    echo "=========================================="
    echo "All figures generated!"
    echo "=========================================="
else
    # Generate specific figure
    case $1 in
        1) generate_figure1 ;;
        2) generate_figure2 ;;
        3) generate_figure3 ;;
        4) generate_figure4 ;;
        *) echo "Unknown figure: $1. Use 1-4." ;;
    esac
fi

echo ""
echo "Output locations:"
echo "  figure1/plots/figure1.{png,pdf}"
echo "  figure2/plots/figure2.{png,pdf}"
echo "  figure3/plots/figure3.{png,pdf}"
