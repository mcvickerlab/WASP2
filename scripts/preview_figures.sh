#!/bin/bash
# Quick preview script for WASP2 publication figures
# Usage: ./preview_figures.sh [figure_number]

FIGURES_DIR="$(dirname "$0")/../figures"

show_usage() {
    echo "WASP2 Figure Preview"
    echo "===================="
    echo ""
    echo "Usage: $0 [figure_number|all|list]"
    echo ""
    echo "Options:"
    echo "  1-5     View specific figure (opens PNG)"
    echo "  all     View all figures"
    echo "  list    List all available figures"
    echo "  stats   Show supplementary table"
    echo ""
    echo "Examples:"
    echo "  $0 1        # View Figure 1 (runtime comparison)"
    echo "  $0 all      # View all figures"
    echo "  $0 list     # List all figures with descriptions"
}

list_figures() {
    echo "WASP2 Publication Figures"
    echo "========================="
    echo ""
    echo "Figure 1: Runtime Performance Comparison"
    echo "  Files: figure1_runtime_comparison.{pdf,png}"
    echo "  Description: Bar chart comparing WASP1, STAR+WASP, WASP2 runtime"
    echo ""
    echo "Figure 2: Pipeline Step Breakdown"
    echo "  Files: figure2_pipeline_breakdown.{pdf,png}"
    echo "  Description: Horizontal bar chart of pipeline step timings"
    echo ""
    echo "Figure 3: Read Processing Throughput"
    echo "  Files: figure3_throughput_comparison.{pdf,png}"
    echo "  Description: Throughput comparison (M reads/second)"
    echo ""
    echo "Figure 4: Variant and Read Statistics"
    echo "  Files: figure4_variant_coverage.{pdf,png}"
    echo "  Description: 2-panel figure (read distribution + variant stats)"
    echo ""
    echo "Figure 5: GATK Comparison"
    echo "  Files: figure5_gatk_comparison.{pdf,png}"
    echo "  Description: 3-panel scatter plot (SNP, INS, DEL accuracy)"
    echo ""
    echo "Supplementary Table S1:"
    echo "  File: supplementary_table_S1.csv"
    echo "  Description: Detailed pipeline statistics"
    echo ""
}

show_stats() {
    echo "Supplementary Table S1 - Pipeline Statistics"
    echo "============================================="
    echo ""
    if [ -f "$FIGURES_DIR/supplementary_table_S1.csv" ]; then
        column -t -s, "$FIGURES_DIR/supplementary_table_S1.csv"
    else
        echo "Error: supplementary_table_S1.csv not found"
        exit 1
    fi
}

open_figure() {
    local fig_num=$1
    local fig_file="$FIGURES_DIR/figure${fig_num}_*.png"

    # Find the actual file
    local actual_file=$(ls $fig_file 2>/dev/null | head -1)

    if [ -z "$actual_file" ]; then
        echo "Error: Figure $fig_num not found"
        exit 1
    fi

    echo "Opening: $(basename $actual_file)"

    # Try different viewers
    if command -v xdg-open &> /dev/null; then
        xdg-open "$actual_file"
    elif command -v eog &> /dev/null; then
        eog "$actual_file"
    elif command -v display &> /dev/null; then
        display "$actual_file"
    else
        echo "No image viewer found. Please open manually:"
        echo "$actual_file"
    fi
}

# Main logic
case "${1:-list}" in
    1|2|3|4|5)
        open_figure "$1"
        ;;
    all)
        echo "Opening all figures..."
        for i in 1 2 3 4 5; do
            open_figure "$i"
        done
        ;;
    list)
        list_figures
        ;;
    stats)
        show_stats
        ;;
    -h|--help)
        show_usage
        ;;
    *)
        echo "Error: Invalid option '$1'"
        echo ""
        show_usage
        exit 1
        ;;
esac
