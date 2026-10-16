# WASP2 Container Usage

## Docker

Choose an existing tag from the [container registry](https://github.com/mcvickerlab/WASP2/pkgs/container/wasp2)
and replace `PUBLISHED_TAG` below. Development source versions do not imply
that a matching image has been published.

Pull and inspect the available CLI tools:

```bash
docker pull ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG

docker run --rm ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG wasp2-count --help
docker run --rm ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG wasp2-map --help
docker run --rm ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG wasp2-analyze --help
```

Mount local data when running workflows:

```bash
docker run --rm -v "$PWD":/data ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG \
  wasp2-count count-variants /data/sample.bam /data/variants.vcf.gz -o /data/counts.tsv
```

## Singularity / Apptainer

For HPC environments using SIF images:

```bash
singularity pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG
singularity exec wasp2.sif wasp2-count --help
```

or:

```bash
apptainer pull wasp2.sif docker://ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG
apptainer exec wasp2.sif wasp2-count --help
```

## Mapping Example

```bash
docker run --rm -v "$PWD":/data ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG \
  wasp2-map make-reads /data/sample.bam /data/variants.vcf.gz \
  --samples sample1 \
  --out_dir /data/remap_dir
```

After realigning the swapped FASTQ reads with the same aligner, version,
reference, and parameters as the original BAM:

```bash
docker run --rm -v "$PWD":/data ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG \
  wasp2-map filter-remapped /data/remapped.bam \
  --wasp_data_json /data/remap_dir/sample_wasp_data_files.json \
  --out_bam /data/filtered.bam
```

## Counting Example

```bash
docker run --rm -v "$PWD":/data ghcr.io/mcvickerlab/wasp2:PUBLISHED_TAG \
  wasp2-count count-variants /data/filtered.bam /data/variants.vcf.gz \
  --samples sample1 \
  --region /data/genes.gtf \
  --out_file /data/counts.tsv
```

## Notes

- The image contains the WASP2 package plus `samtools`, `bcftools`, and `bedtools`.
- The documented public mapping workflow is `make-reads -> realign -> filter-remapped`.
