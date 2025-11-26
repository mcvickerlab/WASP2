#!/bin/bash
# Run WASP2 indel simulation validation

set -e

# Activate wasp2-dev conda environment (requires Python 3.10+ for dataclass slots)
eval "$(conda shell.bash hook)"
conda activate wasp2-dev

echo "=========================================="
echo "WASP2 Indel Validation Simulation"
echo "=========================================="
echo ""

# Check if tier specified
TIER=${1:-minimum}

if [[ "$TIER" != "minimum" && "$TIER" != "moderate" && "$TIER" != "comprehensive" ]]; then
    echo "Usage: $0 [minimum|moderate|comprehensive]"
    echo ""
    echo "Tiers:"
    echo "  minimum       -  90 tests (~10 min)  - Proves algorithm works"
    echo "  moderate      - 270 tests (~30 min)  - Shows robustness across coverage"
    echo "  comprehensive - 810 tests (~2 hrs)   - Edge cases and large indels"
    echo ""
    exit 1
fi

echo "Running $TIER tier simulation..."
echo ""

# Create output directory with timestamp
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="simulation_results_${TIER}_${TIMESTAMP}"
mkdir -p "$OUTDIR"

echo "Output directory: $OUTDIR"
echo ""

# Run simulation
python simulate_indel_ase_v2.py \
    --tier "$TIER" \
    --workdir "$OUTDIR" \
    --keep

# Check if successful
if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Simulation complete!"
    echo ""
    echo "Results saved to:"
    echo "  - Summary: $OUTDIR/simulation_results.csv"
    echo "  - Reference: $OUTDIR/reference.fa"
    echo "  - VCF: $OUTDIR/variants.vcf.gz"
    echo "  - Final BAM: $OUTDIR/wasp2_output/keep.merged.bam"
    echo ""
    echo "To view results:"
    echo "  cat $OUTDIR/simulation_results.csv | column -t -s,"
else
    echo ""
    echo "❌ Simulation failed"
    echo "Check logs in: $OUTDIR"
    exit 1
fi
