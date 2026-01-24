# nf-atacseq: Citations

## Pipeline

If you use nf-atacseq for your analysis, please cite:

> **WASP: Allele-specific software for robust molecular quantitative trait locus discovery**
>
> Bryce van de Geijn, Graham McVicker, Yoav Gilad, Jonathan K Pritchard
>
> _Nature Methods_ 2015 Nov;12(11):1061-3
> doi: [10.1038/nmeth.3582](https://doi.org/10.1038/nmeth.3582)

## Nextflow

> **Nextflow enables reproducible computational workflows**
>
> Paolo Di Tommaso, Maria Chatzou, Evan W. Floden, Pablo Prieto Barja, Emilio Palumbo & Cedric Notredame
>
> _Nature Biotechnology_ 2017 Apr 11;35(4):316-319
> doi: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

## Pipeline components

### Alignment

- **BWA-MEM**

  > Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2
  >
  > [arXiv](https://arxiv.org/abs/1303.3997)

- **Bowtie2**

  > Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012 Mar 4;9(4):357-9.
  >
  > doi: [10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)

### Read Processing

- **fastp**

  > Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018 Sep 1;34(17):i884-i890.
  >
  > doi: [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)

- **Samtools**

  > Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9.
  >
  > doi: [10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352)

- **Picard**

  > Broad Institute. Picard toolkit. [https://broadinstitute.github.io/picard/](https://broadinstitute.github.io/picard/)

### Peak Calling

- **MACS2**

  > Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. Model-based analysis of ChIP-Seq (MACS). Genome Biol. 2008;9(9):R137.
  >
  > doi: [10.1186/gb-2008-9-9-r137](https://doi.org/10.1186/gb-2008-9-9-r137)

### Quality Control

- **FastQC**

  > Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data.
  >
  > [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- **MultiQC**

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8.
  >
  > doi: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)

## BibTeX

```bibtex
@article{vandegeijn2015wasp,
  title={WASP: allele-specific software for robust molecular quantitative trait locus discovery},
  author={van de Geijn, Bryce and McVicker, Graham and Gilad, Yoav and Pritchard, Jonathan K},
  journal={Nature methods},
  volume={12},
  number={11},
  pages={1061--1063},
  year={2015},
  publisher={Nature Publishing Group}
}

@article{ditommaso2017nextflow,
  title={Nextflow enables reproducible computational workflows},
  author={Di Tommaso, Paolo and Chatzou, Maria and Floden, Evan W and Barja, Pablo Prieto and Palumbo, Emilio and Notredame, Cedric},
  journal={Nature biotechnology},
  volume={35},
  number={4},
  pages={316--319},
  year={2017},
  publisher={Nature Publishing Group}
}

@article{li2013bwa,
  title={Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM},
  author={Li, Heng},
  journal={arXiv preprint arXiv:1303.3997},
  year={2013}
}

@article{langmead2012bowtie2,
  title={Fast gapped-read alignment with Bowtie 2},
  author={Langmead, Ben and Salzberg, Steven L},
  journal={Nature methods},
  volume={9},
  number={4},
  pages={357--359},
  year={2012},
  publisher={Nature Publishing Group}
}

@article{chen2018fastp,
  title={fastp: an ultra-fast all-in-one FASTQ preprocessor},
  author={Chen, Shifu and Zhou, Yanqing and Chen, Yaru and Gu, Jia},
  journal={Bioinformatics},
  volume={34},
  number={17},
  pages={i884--i890},
  year={2018},
  publisher={Oxford University Press}
}

@article{li2009samtools,
  title={The sequence alignment/map format and SAMtools},
  author={Li, Heng and Handsaker, Bob and Wysoker, Alec and Fennell, Tim and Ruan, Jue and Homer, Nils and Marth, Gabor and Abecasis, Goncalo and Durbin, Richard},
  journal={Bioinformatics},
  volume={25},
  number={16},
  pages={2078--2079},
  year={2009},
  publisher={Oxford University Press}
}

@article{zhang2008macs,
  title={Model-based analysis of ChIP-Seq (MACS)},
  author={Zhang, Yong and Liu, Tao and Meyer, Clifford A and Eeckhoute, J{\'e}r{\^o}me and Johnson, David S and Bernstein, Bradley E and Nusbaum, Chad and Myers, Richard M and Brown, Myles and Li, Wei and others},
  journal={Genome biology},
  volume={9},
  number={9},
  pages={R137},
  year={2008},
  publisher={BioMed Central}
}

@article{ewels2016multiqc,
  title={MultiQC: summarize analysis results for multiple tools and samples in a single report},
  author={Ewels, Philip and Magnusson, M{\aa}ns and Lundin, Sverker and K{\"a}ller, Max},
  journal={Bioinformatics},
  volume={32},
  number={19},
  pages={3047--3048},
  year={2016},
  publisher={Oxford University Press}
}
```

## Software packaging

- [Bioconda](https://bioconda.github.io/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476.
  >
  > doi: [10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7)

- [BioContainers](https://biocontainers.pro/)

  > da Veiga Leprevost F, Grüning BA, Alber SM, Pireddu L, Bittremieux W, Moreno P, Clements D, Martinez D, Gontier N, Reiter J, Goecks J, Audain E, Perez-Riverol Y, Bowers R, Röst HL. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582.
  >
  > doi: [10.1093/bioinformatics/btx192](https://doi.org/10.1093/bioinformatics/btx192)
