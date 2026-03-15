# WASP2 Container Usage

## Verified Docker Flow

The Docker image validated for this update is:

```bash
ghcr.io/mcvickerlab/wasp2:1.4.0
```

Pull and inspect the available CLI tools:

```bash
docker pull ghcr.io/mcvickerlab/wasp2:1.4.0

docker run --rm ghcr.io/mcvickerlab/wasp2:1.4.0 wasp2-count --help
docker run --rm ghcr.io/mcvickerlab/wasp2:1.4.0 wasp2-map --help
docker run --rm ghcr.io/mcvickerlab/wasp2:1.4.0 wasp2-analyze --help
docker run --rm ghcr.io/mcvickerlab/wasp2:1.4.0 wasp2-ipscore --help
```

Mount local data when running workflows:

```bash
docker run --rm -v "$PWD":/data ghcr.io/mcvickerlab/wasp2:1.4.0 \
  wasp2-count count-variants /data/sample.bam /data/variants.vcf.gz -o /data/counts.tsv
```

## Intended Singularity / Apptainer Flow

For HPC environments using SIF images:

```bash
singularity pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:1.4.0
singularity exec wasp2.sif wasp2-count --help
```

or:

```bash
apptainer pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:1.4.0
apptainer exec wasp2.sif wasp2-count --help
```

These commands are the intended container path, but they were not executed in
this development environment because neither `singularity` nor `apptainer` was
installed locally.

## Mapping Example

```bash
docker run --rm -v "$PWD":/data ghcr.io/mcvickerlab/wasp2:1.4.0 \
  wasp2-map make-reads /data/sample.bam /data/variants.vcf.gz \
  --samples sample1 \
  --out_dir /data/remap_dir
```

After realigning the swapped FASTQ reads with your aligner of choice:

```bash
docker run --rm -v "$PWD":/data ghcr.io/mcvickerlab/wasp2:1.4.0 \
  wasp2-map filter-remapped /data/remapped.bam \
  --wasp_data_json /data/remap_dir/sample_wasp_data_files.json \
  --out_bam /data/filtered.bam
```

## Counting Example

```bash
docker run --rm -v "$PWD":/data ghcr.io/mcvickerlab/wasp2:1.4.0 \
  wasp2-count count-variants /data/filtered.bam /data/variants.vcf.gz \
  --samples sample1 \
  --region /data/genes.gtf \
  --out_file /data/counts.tsv
```

## Notes

- The image contains the WASP2 package plus `samtools`, `bcftools`, and `bedtools`.
- The documented public mapping workflow is `make-reads -> realign -> filter-remapped`.
- `wasp2-ipscore` is present in the container alongside the main analysis tools.
