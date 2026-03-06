#!/usr/bin/env bash
# download_chr21.sh -- Download chr21 real data for WASP2 pipeline validation
#
# This script downloads publicly available genomic data from the 1000 Genomes
# Project and UCSC to create a real-data test set for validating all 4 WASP2
# Nextflow pipelines beyond synthetic test data.
#
# Data sources:
#   - Reference FASTA: UCSC hg38 chr21 (GRCh38)
#   - VCF: 1000 Genomes NYGC 30x, 2022 SHAPEIT5-phased panel (3,202 samples)
#   - CRAM/BAM: 1000 Genomes NYGC 30x high-coverage alignments
#
# Samples available (choose via SAMPLE variable):
#   - NA12878 (CEU, commonly used benchmark, GIAB reference)
#   - HG00731 (PUR, Puerto Rican trio child, used in WASP2 benchmarks)
#
# Estimated disk space:
#   - chr21 reference FASTA:       ~40 MB (uncompressed)
#   - chr21 BAM (one sample):      ~1-2 GB
#   - chr21 VCF (het SNPs only):   ~5 MB
#   - BWA index for chr21:         ~50 MB
#   - STAR index for chr21:        ~1.5 GB (RNA-seq pipeline only)
#   - Total minimum:               ~2 GB (ATAC/scATAC/OUTRIDER)
#   - Total with STAR index:       ~3.5 GB (RNA-seq pipeline)
#
# Prerequisites:
#   samtools >= 1.10, bcftools, tabix, bgzip, bwa, curl
#   Optional: STAR (for RNA-seq STAR index generation)
#
# Usage:
#   bash download_chr21.sh                    # Default: NA12878
#   SAMPLE=HG00731 bash download_chr21.sh     # Use HG00731 instead
#   SKIP_STAR_INDEX=1 bash download_chr21.sh  # Skip STAR index build
#   DRY_RUN=1 bash download_chr21.sh          # Show what would be downloaded

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WASP2_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"

# Sample selection (override with SAMPLE env var)
SAMPLE="${SAMPLE:-NA12878}"
DRY_RUN="${DRY_RUN:-0}"
SKIP_STAR_INDEX="${SKIP_STAR_INDEX:-0}"
SKIP_FASTQ="${SKIP_FASTQ:-0}"
THREADS="${THREADS:-4}"

# ── URLs ─────────────────────────────────────────────────────────────────────

# Reference genome: UCSC hg38 chr21 (~12 MB compressed, ~40 MB uncompressed)
CHR21_FASTA_URL="https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr21.fa.gz"

# Full GRCh38 reference (needed for CRAM-to-BAM conversion, ~3 GB)
# This is the same reference used by 1000 Genomes NYGC for alignment
GRCH38_FULL_REF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# VCF: 1000 Genomes 2022 SHAPEIT5-phased panel (3,202 samples, chr21)
# Two VCF sources available:
#   - 2022 panel (SHAPEIT5, preferred): contains SNV + INDEL + SV, better phasing
#   - 2020 panel (SHAPEIT2): SNV + INDEL only
ONEKG_VCF_2022_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
ONEKG_VCF_2020_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.filtered.shapeit2-duohmm-phased.vcf.gz"

# Use 2022 panel by default (better phasing quality)
ONEKG_VCF_URL="${ONEKG_VCF_2022_URL}"

# CRAM alignment files (NYGC 30x high-coverage, GRCh38)
# Pattern: ftp.sra.ebi.ac.uk/vol1/run/ERR.../ERR.../SAMPLE.final.cram
# These are whole-genome CRAMs; we subset to chr21 during download.
#
# NA12878 (CEU):
#   CRAM: ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram
#   FASTQ R1: ftp.sra.ebi.ac.uk/vol1/fastq/ERR323/004/ERR3239334/ERR3239334_1.fastq.gz
#   FASTQ R2: ftp.sra.ebi.ac.uk/vol1/fastq/ERR323/004/ERR3239334/ERR3239334_2.fastq.gz
#   ENA Run: ERR3239334
#
# HG00731 (PUR):
#   CRAM: ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/PUR/HG00731/high_cov_alignment/HG00731.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.cram
#   Note: HG00731 CRAM is in the HGSV SV discovery collection, not the main 2504 collection

declare -A CRAM_URLS=(
    ["NA12878"]="https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram"
    ["HG00731"]="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/PUR/HG00731/high_cov_alignment/HG00731.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.cram"
)

declare -A FASTQ_R1_URLS=(
    ["NA12878"]="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR323/004/ERR3239334/ERR3239334_1.fastq.gz"
)
declare -A FASTQ_R2_URLS=(
    ["NA12878"]="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR323/004/ERR3239334/ERR3239334_2.fastq.gz"
)

# Note: HG00731 FASTQ files are from the HGSV collection and are very large
# (whole-genome paired-end). For pipeline testing, it's more practical to
# extract FASTQ from the chr21 BAM using samtools fastq.

# ── Logging ──────────────────────────────────────────────────────────────────

log_info()  { echo "[INFO]  $(date +%H:%M:%S) $*"; }
log_warn()  { echo "[WARN]  $(date +%H:%M:%S) $*" >&2; }
log_error() { echo "[ERROR] $(date +%H:%M:%S) $*" >&2; }

# ── Dependency check ─────────────────────────────────────────────────────────

check_deps() {
    local missing=()
    for cmd in samtools bcftools tabix bgzip bwa curl; do
        if ! command -v "$cmd" &>/dev/null; then
            missing+=("$cmd")
        fi
    done

    if [[ ${#missing[@]} -gt 0 ]]; then
        log_error "Missing required tools: ${missing[*]}"
        log_error "Install via: conda activate WASP2_dev2  (or install individually)"
        exit 1
    fi

    # Check samtools version >= 1.10 (needed for CRAM support with remote refs)
    local st_ver
    st_ver=$(samtools --version | head -1 | grep -oE '[0-9]+\.[0-9]+')
    local st_major st_minor
    st_major=$(echo "$st_ver" | cut -d. -f1)
    st_minor=$(echo "$st_ver" | cut -d. -f2)
    if [[ "$st_major" -lt 1 ]] || { [[ "$st_major" -eq 1 ]] && [[ "$st_minor" -lt 10 ]]; }; then
        log_error "samtools >= 1.10 required (found ${st_ver}). Upgrade for CRAM support."
        exit 1
    fi
    log_info "Dependencies OK (samtools ${st_ver})"
}

# ── Disk space estimate ─────────────────────────────────────────────────────

estimate_disk_space() {
    log_info ""
    log_info "Estimated disk space requirements:"
    log_info "  chr21 reference FASTA (uncompressed):  ~40 MB"
    log_info "  GRCh38 full reference (for CRAM conv): ~3 GB (temp, deleted after)"
    log_info "  chr21 BAM (30x coverage):              ~1-2 GB"
    log_info "  chr21 het VCF (one sample):             ~5 MB"
    log_info "  BWA index for chr21:                   ~50 MB"
    if [[ "${SKIP_STAR_INDEX}" != "1" ]]; then
        log_info "  STAR index for chr21:                  ~1.5 GB"
        log_info "  ─────────────────────────────────────"
        log_info "  TOTAL (peak, during GRCh38 download):  ~6 GB"
        log_info "  TOTAL (final, after cleanup):          ~3.5 GB"
    else
        log_info "  ─────────────────────────────────────"
        log_info "  TOTAL (peak, during GRCh38 download):  ~5 GB"
        log_info "  TOTAL (final, after cleanup):          ~2 GB"
    fi
    log_info ""
}

# ── Step 1: Download chr21 reference FASTA ───────────────────────────────────

download_reference() {
    local chr21_fa="${DATA_DIR}/chr21.fa"
    local chr21_fa_gz="${DATA_DIR}/chr21.fa.gz"

    if [[ -f "${chr21_fa}" ]] && [[ -f "${chr21_fa}.fai" ]]; then
        log_info "chr21 reference already present -- skipping"
        return 0
    fi

    log_info "Downloading chr21 reference FASTA from UCSC (~12 MB compressed)..."
    log_info "  URL: ${CHR21_FASTA_URL}"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log_info "  [DRY RUN] Would download chr21.fa.gz"
        return 0
    fi

    curl -L --progress-bar -o "${chr21_fa_gz}" "${CHR21_FASTA_URL}"

    # Decompress
    log_info "Decompressing..."
    gunzip -f "${chr21_fa_gz}"

    # Create FASTA index
    log_info "Indexing reference..."
    samtools faidx "${chr21_fa}"

    # Create sequence dictionary (needed by some tools)
    samtools dict "${chr21_fa}" > "${DATA_DIR}/chr21.dict"

    log_info "Reference ready: ${chr21_fa} ($(du -h "${chr21_fa}" | cut -f1))"
}

# ── Step 2: Build BWA index ─────────────────────────────────────────────────

build_bwa_index() {
    local chr21_fa="${DATA_DIR}/chr21.fa"
    local bwa_dir="${DATA_DIR}/bwa_index"

    if [[ -f "${bwa_dir}/chr21.fa.bwt" ]]; then
        log_info "BWA index already present -- skipping"
        return 0
    fi

    log_info "Building BWA index for chr21..."

    if [[ "${DRY_RUN}" == "1" ]]; then
        log_info "  [DRY RUN] Would build BWA index"
        return 0
    fi

    mkdir -p "${bwa_dir}"

    # Symlink reference into bwa_index dir (nf-core convention)
    ln -sf "${chr21_fa}" "${bwa_dir}/chr21.fa"
    ln -sf "${chr21_fa}.fai" "${bwa_dir}/chr21.fa.fai"

    bwa index "${bwa_dir}/chr21.fa"
    log_info "BWA index ready at ${bwa_dir}/"
}

# ── Step 3: Download and subset VCF ─────────────────────────────────────────

download_vcf() {
    local vcf_out="${DATA_DIR}/${SAMPLE}_chr21_het.vcf.gz"

    if [[ -f "${vcf_out}" ]] && [[ -f "${vcf_out}.tbi" ]]; then
        log_info "VCF already present -- skipping"
        return 0
    fi

    local full_vcf="${DATA_DIR}/1kGP_chr21_full.vcf.gz"

    log_info "Downloading 1000 Genomes chr21 phased VCF (~407 MB)..."
    log_info "  URL: ${ONEKG_VCF_URL}"
    log_info "  Source: 1kGP 2022 SHAPEIT5-phased panel (3,202 samples, GRCh38)"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log_info "  [DRY RUN] Would download and subset VCF"
        return 0
    fi

    # Download full chr21 VCF
    if [[ ! -f "${full_vcf}" ]]; then
        curl -L --progress-bar -o "${full_vcf}" "${ONEKG_VCF_URL}"
        curl -sL -o "${full_vcf}.tbi" "${ONEKG_VCF_URL}.tbi"
    fi

    # Verify download
    if [[ ! -f "${full_vcf}" ]] || [[ $(wc -c < "${full_vcf}") -lt 1000000 ]]; then
        log_error "VCF download failed or file too small"
        rm -f "${full_vcf}" "${full_vcf}.tbi"
        exit 1
    fi

    # Verify sample exists in VCF
    log_info "Verifying ${SAMPLE} is in the VCF..."
    local sample_list
    sample_list="$(bcftools query -l "${full_vcf}" 2>/dev/null)"
    if ! echo "${sample_list}" | grep -q "^${SAMPLE}$"; then
        log_error "${SAMPLE} not found in VCF!"
        log_error "Available samples (first 10):"
        echo "${sample_list}" | head -10 >&2
        exit 1
    fi
    log_info "  Confirmed: ${SAMPLE} present"

    # Subset to sample, SNPs only, het sites only
    log_info "Subsetting to ${SAMPLE} het SNPs..."
    bcftools view -s "${SAMPLE}" -v snps "${full_vcf}" 2>/dev/null \
        | bcftools view -i 'GT="het"' -Oz -o "${vcf_out}" 2>/dev/null
    tabix -p vcf "${vcf_out}"

    # Also create a VCF with ALL variants for the sample (not just het SNPs)
    # Some pipelines may want indels and hom-alt sites too
    local vcf_all="${DATA_DIR}/${SAMPLE}_chr21_all.vcf.gz"
    log_info "Creating full variant VCF (all genotypes, SNPs+indels)..."
    bcftools view -s "${SAMPLE}" -i 'GT="alt"' "${full_vcf}" 2>/dev/null \
        | bcftools view -Oz -o "${vcf_all}" 2>/dev/null
    tabix -p vcf "${vcf_all}"

    # Count variants
    local n_het n_all
    n_het=$(bcftools view -H "${vcf_out}" 2>/dev/null | wc -l | tr -d ' ')
    n_all=$(bcftools view -H "${vcf_all}" 2>/dev/null | wc -l | tr -d ' ')
    log_info "  ${n_het} het SNPs, ${n_all} total variants for ${SAMPLE} on chr21"

    # Clean up full VCF to save disk space
    log_info "Removing full VCF to save space..."
    rm -f "${full_vcf}" "${full_vcf}.tbi"

    log_info "VCF ready: ${vcf_out}"
}

# ── Step 4: Download and subset CRAM → chr21 BAM ────────────────────────────

download_bam() {
    local bam_out="${DATA_DIR}/${SAMPLE}_chr21.bam"

    if [[ -f "${bam_out}" ]] && [[ -f "${bam_out}.bai" ]]; then
        log_info "BAM already present -- skipping"
        return 0
    fi

    local cram_url="${CRAM_URLS[${SAMPLE}]:-}"
    if [[ -z "${cram_url}" ]]; then
        log_error "No CRAM URL configured for sample ${SAMPLE}"
        log_error "Supported samples: ${!CRAM_URLS[*]}"
        exit 1
    fi

    log_info "Downloading chr21 reads from CRAM (streaming, no full download needed)..."
    log_info "  CRAM: ${cram_url}"
    log_info "  This streams the CRAM, extracts chr21 reads, and writes BAM."
    log_info "  Requires the GRCh38 reference for CRAM decoding."

    if [[ "${DRY_RUN}" == "1" ]]; then
        log_info "  [DRY RUN] Would stream CRAM and extract chr21 BAM"
        return 0
    fi

    # We need the full GRCh38 reference for CRAM decoding
    local grch38_ref="${DATA_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    if [[ ! -f "${grch38_ref}" ]]; then
        log_info "Downloading GRCh38 full reference for CRAM decoding (~3 GB)..."
        log_info "  URL: ${GRCH38_FULL_REF_URL}"
        log_info "  (This will be deleted after BAM extraction to save space)"
        curl -L --progress-bar -o "${grch38_ref}" "${GRCH38_FULL_REF_URL}"
        samtools faidx "${grch38_ref}"
    fi

    # Stream CRAM, subset to chr21, output as BAM
    # The CRAM index (.crai) should be alongside the CRAM on the server
    log_info "Extracting chr21 reads from CRAM → BAM..."
    log_info "  (This may take 10-30 minutes depending on bandwidth)"

    samtools view -b -@ "${THREADS}" \
        -T "${grch38_ref}" \
        "${cram_url}" \
        chr21 \
        -o "${bam_out}.unsorted"

    # Sort and index
    log_info "Sorting BAM..."
    samtools sort -@ "${THREADS}" -o "${bam_out}" "${bam_out}.unsorted"
    rm -f "${bam_out}.unsorted"

    log_info "Indexing BAM..."
    samtools index "${bam_out}"

    local n_reads
    n_reads=$(samtools view -c "${bam_out}")
    local bam_size
    bam_size=$(du -h "${bam_out}" | cut -f1)
    log_info "BAM ready: ${bam_out} (${bam_size}, ${n_reads} reads)"

    # Clean up full reference (keep only chr21)
    log_info "Removing full GRCh38 reference to save space..."
    rm -f "${grch38_ref}" "${grch38_ref}.fai"
}

# ── Step 5: Extract FASTQ from chr21 BAM ────────────────────────────────────

extract_fastq() {
    local bam_in="${DATA_DIR}/${SAMPLE}_chr21.bam"
    local fq_r1="${DATA_DIR}/${SAMPLE}_chr21_R1.fastq.gz"
    local fq_r2="${DATA_DIR}/${SAMPLE}_chr21_R2.fastq.gz"

    if [[ "${SKIP_FASTQ}" == "1" ]]; then
        log_info "Skipping FASTQ extraction (SKIP_FASTQ=1)"
        return 0
    fi

    if [[ -f "${fq_r1}" ]] && [[ -f "${fq_r2}" ]]; then
        log_info "FASTQ files already present -- skipping"
        return 0
    fi

    if [[ ! -f "${bam_in}" ]]; then
        log_error "BAM not found -- run download_bam first"
        exit 1
    fi

    log_info "Extracting paired-end FASTQ from chr21 BAM..."

    if [[ "${DRY_RUN}" == "1" ]]; then
        log_info "  [DRY RUN] Would extract FASTQ from BAM"
        return 0
    fi

    # Sort by name first (required for proper paired-end extraction)
    local namesorted="${DATA_DIR}/${SAMPLE}_chr21.namesorted.bam"
    samtools sort -n -@ "${THREADS}" -o "${namesorted}" "${bam_in}"

    # Extract paired FASTQ
    samtools fastq -@ "${THREADS}" \
        -1 "${fq_r1}" \
        -2 "${fq_r2}" \
        -0 /dev/null \
        -s /dev/null \
        -n \
        "${namesorted}"

    rm -f "${namesorted}"

    local r1_size r2_size
    r1_size=$(du -h "${fq_r1}" | cut -f1)
    r2_size=$(du -h "${fq_r2}" | cut -f1)
    log_info "FASTQ ready: R1=${r1_size}, R2=${r2_size}"
}

# ── Step 6: Build STAR index (optional, for nf-rnaseq) ──────────────────────

build_star_index() {
    local star_dir="${DATA_DIR}/star_index"

    if [[ "${SKIP_STAR_INDEX}" == "1" ]]; then
        log_info "Skipping STAR index build (SKIP_STAR_INDEX=1)"
        return 0
    fi

    if [[ -d "${star_dir}" ]] && [[ -f "${star_dir}/SA" ]]; then
        log_info "STAR index already present -- skipping"
        return 0
    fi

    if ! command -v STAR &>/dev/null; then
        log_warn "STAR not found -- skipping STAR index build"
        log_warn "Install STAR or set SKIP_STAR_INDEX=1 to suppress this warning"
        return 0
    fi

    log_info "Building STAR index for chr21 (this takes ~10 min and ~8 GB RAM)..."

    if [[ "${DRY_RUN}" == "1" ]]; then
        log_info "  [DRY RUN] Would build STAR index"
        return 0
    fi

    # We need a GTF for STAR. Use Ensembl chr21 annotation or a subset.
    # For simplicity, build without annotation (still usable, just no splice junctions)
    local chr21_fa="${DATA_DIR}/chr21.fa"
    mkdir -p "${star_dir}"

    STAR --runMode genomeGenerate \
        --genomeDir "${star_dir}" \
        --genomeFastaFiles "${chr21_fa}" \
        --genomeSAindexNbases 11 \
        --runThreadN "${THREADS}"

    log_info "STAR index ready at ${star_dir}/"
}

# ── Step 7: Generate pipeline samplesheets ───────────────────────────────────

generate_samplesheets() {
    log_info "Generating pipeline samplesheets..."

    if [[ "${DRY_RUN}" == "1" ]]; then
        log_info "  [DRY RUN] Would generate samplesheets"
        return 0
    fi

    local ss_dir="${SCRIPT_DIR}/samplesheets"
    mkdir -p "${ss_dir}"

    # Paths (use absolute paths for Nextflow compatibility)
    local bam="${DATA_DIR}/${SAMPLE}_chr21.bam"
    local bai="${DATA_DIR}/${SAMPLE}_chr21.bam.bai"
    local fq_r1="${DATA_DIR}/${SAMPLE}_chr21_R1.fastq.gz"
    local fq_r2="${DATA_DIR}/${SAMPLE}_chr21_R2.fastq.gz"

    # ── nf-atacseq samplesheet (FASTQ input) ──
    cat > "${ss_dir}/atacseq_samplesheet.csv" <<EOF
sample,fastq_1,fastq_2,sample_name
${SAMPLE}_chr21,${fq_r1},${fq_r2},${SAMPLE}
EOF
    log_info "  Created: atacseq_samplesheet.csv"

    # ── nf-rnaseq samplesheet (FASTQ input) ──
    cat > "${ss_dir}/rnaseq_samplesheet.csv" <<EOF
sample,fastq_1,fastq_2
${SAMPLE}_chr21,${fq_r1},${fq_r2}
EOF
    log_info "  Created: rnaseq_samplesheet.csv"

    # ── nf-outrider samplesheet (BAM input) ──
    # OUTRIDER needs multiple samples for outlier detection.
    # With a single sample, you can test the pipeline but not the
    # statistical analysis. For proper OUTRIDER testing, download
    # additional samples by running this script multiple times with
    # different SAMPLE values.
    cat > "${ss_dir}/outrider_samplesheet.csv" <<EOF
sample,bam,bai
${SAMPLE}_chr21,${bam},${bai}
EOF
    log_info "  Created: outrider_samplesheet.csv"
    log_info "  NOTE: OUTRIDER requires multiple samples. Re-run with different"
    log_info "        SAMPLE values and append rows for proper testing."

    # ── nf-scatac samplesheet (fragments input) ──
    # scATAC needs fragments.tsv.gz from CellRanger-ATAC, which cannot be
    # derived from a bulk BAM. This samplesheet is a placeholder template.
    cat > "${ss_dir}/scatac_samplesheet.csv" <<EOF
sample,fragments,cellranger_dir,barcode_tag,chemistry
${SAMPLE}_chr21,PLACEHOLDER_fragments.tsv.gz,,CB,10x-atac-v2
EOF
    log_info "  Created: scatac_samplesheet.csv"
    log_info "  NOTE: scATAC requires CellRanger-ATAC fragments.tsv.gz output."
    log_info "        This samplesheet is a template; update paths to real"
    log_info "        10x scATAC-seq data for the sample."

    log_info "Samplesheets saved to ${ss_dir}/"
}

# ── Step 8: Generate Nextflow config for real data ───────────────────────────

generate_configs() {
    log_info "Generating Nextflow test configs..."

    if [[ "${DRY_RUN}" == "1" ]]; then
        log_info "  [DRY RUN] Would generate Nextflow configs"
        return 0
    fi

    local cfg_dir="${SCRIPT_DIR}/configs"
    mkdir -p "${cfg_dir}"

    local vcf_het="${DATA_DIR}/${SAMPLE}_chr21_het.vcf.gz"
    local chr21_fa="${DATA_DIR}/chr21.fa"
    local bwa_dir="${DATA_DIR}/bwa_index"
    local star_dir="${DATA_DIR}/star_index"
    local ss_dir="${SCRIPT_DIR}/samplesheets"

    # ── nf-atacseq config ──
    cat > "${cfg_dir}/test_real_atacseq.config" <<EOF
/*
 * Real data test config for nf-atacseq
 * Sample: ${SAMPLE} chr21 (1000 Genomes 30x)
 * Usage: nextflow run pipelines/nf-atacseq -c tests/real_data/configs/test_real_atacseq.config
 */
params {
    config_profile_name        = 'Real data test (${SAMPLE} chr21)'
    config_profile_description = 'Real chr21 data from 1000 Genomes for pipeline validation'

    max_cpus   = 8
    max_memory = '32.GB'
    max_time   = '12.h'

    input     = '${ss_dir}/atacseq_samplesheet.csv'
    fasta     = '${chr21_fa}'
    bwa_index = '${bwa_dir}'
    vcf       = '${vcf_het}'

    macs_gsize = '4.6e7'   // chr21 effective genome size (~46 MB)
}
EOF

    # ── nf-rnaseq config ──
    cat > "${cfg_dir}/test_real_rnaseq.config" <<EOF
/*
 * Real data test config for nf-rnaseq
 * Sample: ${SAMPLE} chr21 (1000 Genomes 30x)
 * Usage: nextflow run pipelines/nf-rnaseq -c tests/real_data/configs/test_real_rnaseq.config
 *
 * NOTE: This uses WGS data aligned to chr21, not actual RNA-seq data.
 *       Useful for testing pipeline mechanics but not for biological interpretation.
 *       For real RNA-seq validation, use GTEx or ENCODE chr21-subsetted BAMs.
 */
params {
    config_profile_name        = 'Real data test (${SAMPLE} chr21)'
    config_profile_description = 'Real chr21 data from 1000 Genomes for pipeline validation'

    max_cpus   = 8
    max_memory = '32.GB'
    max_time   = '12.h'

    input      = '${ss_dir}/rnaseq_samplesheet.csv'
    vcf        = '${vcf_het}'
    star_index = '${star_dir}'
    gtf        = null  // TODO: Download Ensembl chr21 GTF for splice-aware alignment
}
EOF

    # ── nf-outrider config ──
    cat > "${cfg_dir}/test_real_outrider.config" <<EOF
/*
 * Real data test config for nf-outrider
 * Sample: ${SAMPLE} chr21 (1000 Genomes 30x)
 * Usage: nextflow run pipelines/nf-outrider -c tests/real_data/configs/test_real_outrider.config
 */
params {
    config_profile_name        = 'Real data test (${SAMPLE} chr21)'
    config_profile_description = 'Real chr21 data from 1000 Genomes for pipeline validation'

    max_cpus   = 8
    max_memory = '64.GB'
    max_time   = '24.h'

    input = '${ss_dir}/outrider_samplesheet.csv'
    vcf   = '${vcf_het}'
    gtf   = null  // TODO: Download Ensembl chr21 GTF
}
EOF

    # ── nf-scatac config (template) ──
    cat > "${cfg_dir}/test_real_scatac.config" <<EOF
/*
 * Real data test config for nf-scatac (TEMPLATE)
 * Usage: nextflow run pipelines/nf-scatac -c tests/real_data/configs/test_real_scatac.config
 *
 * NOTE: scATAC-seq requires CellRanger-ATAC output (fragments.tsv.gz).
 *       The bulk WGS data from 1000 Genomes cannot be used directly.
 *       Update the samplesheet with real 10x scATAC-seq data paths.
 */
params {
    config_profile_name        = 'Real data test (${SAMPLE} chr21)'
    config_profile_description = 'Real chr21 data for scATAC pipeline validation'

    max_cpus   = 8
    max_memory = '64.GB'
    max_time   = '24.h'

    input = '${ss_dir}/scatac_samplesheet.csv'
    vcf   = '${vcf_het}'

    min_fragments_per_cell = 1000
    min_cells_per_snp      = 3
    min_count              = 10
}
EOF

    log_info "Configs saved to ${cfg_dir}/"
}

# ── Summary ──────────────────────────────────────────────────────────────────

print_summary() {
    echo ""
    log_info "============================================================"
    log_info "  WASP2 Real Data Download Complete"
    log_info "============================================================"
    echo ""
    echo "  Sample:     ${SAMPLE}"
    echo "  Region:     chr21 (GRCh38/hg38)"
    echo "  Source:     1000 Genomes Project NYGC 30x high-coverage"
    echo ""

    if [[ "${DRY_RUN}" == "1" ]]; then
        echo "  [DRY RUN -- no files were downloaded]"
        echo ""
        return
    fi

    echo "  Files:"
    local chr21_fa="${DATA_DIR}/chr21.fa"
    local bam="${DATA_DIR}/${SAMPLE}_chr21.bam"
    local vcf="${DATA_DIR}/${SAMPLE}_chr21_het.vcf.gz"
    local vcf_all="${DATA_DIR}/${SAMPLE}_chr21_all.vcf.gz"
    local fq_r1="${DATA_DIR}/${SAMPLE}_chr21_R1.fastq.gz"
    local fq_r2="${DATA_DIR}/${SAMPLE}_chr21_R2.fastq.gz"

    [[ -f "${chr21_fa}" ]]  && echo "    Reference:  ${chr21_fa}  ($(du -h "${chr21_fa}" | cut -f1))"
    [[ -f "${bam}" ]]       && echo "    BAM:        ${bam}  ($(du -h "${bam}" | cut -f1))"
    [[ -f "${vcf}" ]]       && echo "    VCF (het):  ${vcf}  ($(du -h "${vcf}" | cut -f1))"
    [[ -f "${vcf_all}" ]]   && echo "    VCF (all):  ${vcf_all}  ($(du -h "${vcf_all}" | cut -f1))"
    [[ -f "${fq_r1}" ]]     && echo "    FASTQ R1:   ${fq_r1}  ($(du -h "${fq_r1}" | cut -f1))"
    [[ -f "${fq_r2}" ]]     && echo "    FASTQ R2:   ${fq_r2}  ($(du -h "${fq_r2}" | cut -f1))"

    if [[ -f "${bam}" ]]; then
        local n_reads
        n_reads=$(samtools view -c "${bam}" 2>/dev/null || echo "?")
        echo "    Reads:      ${n_reads}"
    fi
    if [[ -f "${vcf}" ]]; then
        local n_vars
        n_vars=$(bcftools view -H "${vcf}" 2>/dev/null | wc -l | tr -d ' ')
        echo "    Het SNPs:   ${n_vars}"
    fi

    echo ""
    echo "  Samplesheets:  ${SCRIPT_DIR}/samplesheets/"
    echo "  Configs:       ${SCRIPT_DIR}/configs/"
    echo ""
    echo "  Total disk:    $(du -sh "${DATA_DIR}" | cut -f1)"
    echo ""
    echo "  Run pipelines with:"
    echo "    nextflow run pipelines/nf-atacseq  -c tests/real_data/configs/test_real_atacseq.config"
    echo "    nextflow run pipelines/nf-rnaseq   -c tests/real_data/configs/test_real_rnaseq.config"
    echo "    nextflow run pipelines/nf-outrider -c tests/real_data/configs/test_real_outrider.config"
    echo ""
}

# ── Main ─────────────────────────────────────────────────────────────────────

main() {
    log_info "WASP2 Real Data Downloader"
    log_info "Sample: ${SAMPLE} | Region: chr21 | Assembly: GRCh38/hg38"
    echo ""

    check_deps
    estimate_disk_space

    mkdir -p "${DATA_DIR}"

    # Confirm before proceeding (unless DRY_RUN)
    if [[ "${DRY_RUN}" != "1" ]]; then
        read -rp "Proceed with download? [y/N] " confirm
        if [[ "${confirm}" != "y" && "${confirm}" != "Y" ]]; then
            log_info "Aborted."
            exit 0
        fi
    fi

    download_reference
    build_bwa_index
    download_vcf
    download_bam
    extract_fastq
    build_star_index
    generate_samplesheets
    generate_configs
    print_summary
}

main "$@"
