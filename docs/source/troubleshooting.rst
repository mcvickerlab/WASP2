Troubleshooting Guide
=====================

Solutions to common issues organized by module and error type.

.. contents:: Categories
   :local:
   :depth: 2

Installation Issues
-------------------

Rust extension fails to build
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Error**:

.. code-block:: text

   error: failed to run custom build command for `wasp2-rust`

**Causes and solutions**:

1. **Missing Rust compiler**:

   .. code-block:: bash

      curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
      source $HOME/.cargo/env

2. **Missing libclang**:

   .. code-block:: bash

      # Ubuntu/Debian
      sudo apt-get install libclang-dev

      # macOS
      brew install llvm
      export LIBCLANG_PATH=$(brew --prefix llvm)/lib

      # CentOS/RHEL
      sudo yum install clang-devel

3. **Incompatible maturin version**:

   .. code-block:: bash

      pip install --upgrade maturin
      pip install wasp2 --no-cache-dir

ImportError: No module named 'wasp2_rust'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Rust extension didn't build. Options:

1. **Rebuild with verbose output**:

   .. code-block:: bash

      pip uninstall wasp2
      pip install wasp2 -v

2. **Use Python fallback** (slower):

   .. code-block:: bash

      export WASP2_DISABLE_RUST=1

Counting Module Issues
----------------------

No output SNPs
~~~~~~~~~~~~~~

**Symptoms**: Output file empty or has only header

**Diagnostic steps**:

.. code-block:: bash

   # 1. Check sample exists in VCF
   bcftools query -l variants.vcf.gz | grep -w "sample1"

   # 2. Check for heterozygous SNPs
   bcftools view -s sample1 -g het variants.vcf.gz | head -20

   # 3. Check BAM has reads
   samtools view -c sample.bam

   # 4. Check coordinate overlap
   # Get first SNP position
   bcftools view -H variants.vcf.gz | head -1 | cut -f1,2
   # Check BAM coverage there
   samtools view sample.bam chr1:1000000-1000100 | head

**Solutions**:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Problem
     - Solution
   * - Sample not in VCF
     - Use exact sample name from ``bcftools query -l``
   * - No het SNPs for sample
     - Check VCF has genotypes (not just allele frequencies)
   * - Chromosome naming mismatch
     - Convert: ``sed 's/^chr//' file.vcf`` or ``sed 's/^/chr/' file.vcf``
   * - Different reference builds
     - Ensure both use same genome (hg38 vs GRCh38)

"Sample not found in VCF" error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Error**:

.. code-block:: text

   ValueError: Sample 'sample1' not found in VCF file

**Solution**:

.. code-block:: bash

   # List available samples
   bcftools query -l variants.vcf.gz

   # Use exact name (case-sensitive)
   wasp2-count count-variants sample.bam variants.vcf.gz \
       --samples "NA12878"  # Use exact name from list

"BAM file not sorted" error
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Error**:

.. code-block:: text

   RuntimeError: BAM file is not coordinate-sorted

**Solution**:

.. code-block:: bash

   samtools sort -o sample_sorted.bam sample.bam
   samtools index sample_sorted.bam

Missing BAM index
~~~~~~~~~~~~~~~~~

**Error**:

.. code-block:: text

   FileNotFoundError: [Errno 2] No such file or directory: 'sample.bam.bai'

**Solution**:

.. code-block:: bash

   samtools index sample.bam

Missing VCF index
~~~~~~~~~~~~~~~~~

**Error**:

.. code-block:: text

   FileNotFoundError: variants.vcf.gz.tbi not found

**Solution**:

.. code-block:: bash

   # For VCF.gz
   bcftools index -t variants.vcf.gz
   # Or
   tabix -p vcf variants.vcf.gz

   # For BCF
   bcftools index variants.bcf

Low counts everywhere
~~~~~~~~~~~~~~~~~~~~~

**Symptoms**: Most sites have < 10 reads

**Diagnostic**:

.. code-block:: bash

   # Check overall BAM coverage
   samtools depth sample.bam | awk '{sum+=$3; n++} END {print "Mean depth:", sum/n}'

   # Check alignment rate
   samtools flagstat sample.bam

**Solutions**:

- Verify sequencing depth (need 30M+ reads for RNA-seq)
- Check alignment quality
- Lower ``--min`` threshold if appropriate
- Ensure reads aren't being filtered (check mapping quality)

Mapping Module Issues
---------------------

make-reads produces no output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Symptoms**: No FASTQ files generated

**Diagnostic**:

.. code-block:: bash

   # Check for SNP-overlapping reads
   bedtools intersect -a sample.bam -b <(bcftools view -H variants.vcf.gz | \
       awk '{print $1"\t"$2-1"\t"$2}') -wa | head

**Solutions**:

- Verify BAM and VCF coordinates overlap
- Check VCF has variants in sequenced regions
- Ensure sample has heterozygous SNPs

filter-remapped discards all reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Symptoms**: Output BAM is empty

**Causes**:

1. **Remapping with different parameters**:
   Use identical aligner settings for remapping

2. **Read name format changed**:
   Don't modify read names during remapping

3. **Wrong files specified**:
   Verify JSON file matches remapped BAM

**Diagnostic**:

.. code-block:: bash

   # Check remapped BAM has reads
   samtools view -c wasp_output/sample_remapped.bam

   # Check JSON file exists and is correct
   cat sample_wasp_data_files.json

Analysis Module Issues
----------------------

No significant results
~~~~~~~~~~~~~~~~~~~~~~

**Symptoms**: All FDR values > 0.05

**Possible causes and solutions**:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Cause
     - Solution
   * - Insufficient coverage
     - Lower ``--min`` threshold (e.g., 10 instead of 20)
   * - Few informative SNPs
     - Use more complete VCF (e.g., imputed genotypes)
   * - Biological: no true imbalance
     - This may be a valid result
   * - Conservative statistics
     - Check if binomial test gives different results

Unexpected imbalance patterns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Symptoms**: Results don't match expectations

**Diagnostic steps**:

1. **Check known controls** (imprinted genes should show AI)
2. **Verify strand information** (some tools require specific strand handling)
3. **Check for batch effects** (compare replicates)
4. **Visualize raw counts** before statistical testing

Performance Issues
------------------

Processing is very slow
~~~~~~~~~~~~~~~~~~~~~~~

**Solutions** (in order of impact):

1. **Use PGEN format** (25x faster):

   .. code-block:: bash

      plink2 --vcf variants.vcf.gz --make-pgen --out variants
      wasp2-count count-variants sample.bam variants.pgen ...

2. **Install cyvcf2** (7x faster VCF):

   .. code-block:: bash

      pip install wasp2[cyvcf2]

3. **Process by chromosome**:

   .. code-block:: bash

      for chr in {1..22}; do
          wasp2-count count-variants sample.bam variants.vcf.gz \
              --region chr${chr}.bed --out_file counts_chr${chr}.tsv &
      done
      wait

Out of memory
~~~~~~~~~~~~~

**Error**:

.. code-block:: text

   MemoryError
   # or
   Killed (signal 9)

**Solutions**:

1. **Use PGEN format** (lower memory)
2. **Process by chromosome**
3. **Reduce threads**
4. **Pre-filter VCF**:

   .. code-block:: bash

      bcftools view -s sample1 -g het variants.vcf.gz -O z -o het_only.vcf.gz

"No space left on device"
~~~~~~~~~~~~~~~~~~~~~~~~~

**Error**:

.. code-block:: text

   OSError: [Errno 28] No space left on device

**Solutions**:

1. **Use different temp directory**:

   .. code-block:: bash

      wasp2-count count-variants sample.bam variants.vcf.gz \
          --temp_loc /scratch/large_disk/

2. **Clean up existing temp files**:

   .. code-block:: bash

      rm -rf /tmp/wasp2_*

Single-Cell Issues
------------------

Very sparse matrix
~~~~~~~~~~~~~~~~~~

**Symptoms**: Most cells have 0 counts

This is **normal** for single-cell data. Solutions:

1. **Aggregate by cell type** (recommended)
2. **Focus on high-coverage regions**
3. **Use pseudobulk analysis**

Memory issues with large datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Solutions**:

1. **Process chromosomes separately**
2. **Subsample cells for initial analysis**
3. **Use sparse matrix operations**

.. code-block:: bash

   # Subsample barcodes
   shuf cell_barcodes.txt | head -1000 > subset_barcodes.txt

Error Messages Reference
------------------------

Quick reference for common errors:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Error
     - Module
     - Solution
   * - ``FileNotFoundError: variants.vcf.gz.tbi``
     - count
     - ``bcftools index variants.vcf.gz``
   * - ``ValueError: Sample not found``
     - count
     - Check ``bcftools query -l``
   * - ``RuntimeError: BAM not sorted``
     - count
     - ``samtools sort``
   * - ``OSError: No space left``
     - All
     - Use ``--temp_loc``
   * - ``ImportError: wasp2_rust``
     - All
     - Rebuild or set ``WASP2_DISABLE_RUST=1``

Getting More Help
-----------------

If you can't resolve your issue:

1. **Search existing issues**: `GitHub Issues <https://github.com/Jaureguy760/WASP2-exp/issues>`_

2. **Open a new issue** with:

   - WASP2 version: ``wasp2-count --version``
   - Python version: ``python --version``
   - Operating system
   - Complete error message
   - Command that caused the error
   - Minimal reproducible example
