
# WASP2 Pipeline Schematic - BioRender Specification

## Overview
Create a publication-quality pipeline diagram in the style of BiaPy Figure 1.
Three main sections flowing top-to-bottom with grouped boxes.

## Color Palette (Bang Wong colorblind-safe)
- Input boxes: Light blue (#E8F4FD) with blue border (#0072B2)
- Process boxes: Light orange (#FFF3E0) with orange border (#E69F00)
- Output boxes: Light green (#E8F5E9) with green border (#009E73)
- Arrows: Dark gray (#555555)
- Section backgrounds: Light gray (#FAFAFA) with dashed border

## Layout (Top to Bottom)

### Section 1: INPUT (top)
┌─────────────────────────────────────────────────────────┐
│                        INPUT                             │
│  ┌─────────────┐                    ┌─────────────┐     │
│  │    BAM      │                    │    VCF      │     │
│  │  (aligned   │                    │   (het      │     │
│  │   reads)    │                    │  variants)  │     │
│  └─────────────┘                    └─────────────┘     │
└─────────────────────────────────────────────────────────┘
         │                                   │
         └──────────────┬───────────────────┘
                        ▼

### Section 2: WASP REMAPPING (middle, largest)
┌─────────────────────────────────────────────────────────┐
│                   WASP REMAPPING                         │
│                                                          │
│              ┌───────────────────────┐                  │
│              │  Find reads at        │                  │
│              │  variant sites        │                  │
│              └───────────────────────┘                  │
│                        │                                 │
│                        ▼                                 │
│              ┌───────────────────────┐                  │
│              │  Generate alternate   │                  │
│              │  haplotypes           │  ← SNVs + INDELs │
│              └───────────────────────┘                  │
│                        │                                 │
│                        ▼                                 │
│              ┌───────────────────────┐                  │
│              │  Remap (BWA/STAR)     │                  │
│              │  + WASP filter        │                  │
│              └───────────────────────┘                  │
│                                                          │
└─────────────────────────────────────────────────────────┘
                        │
         ┌──────────────┴───────────────┐
         ▼                              ▼

### Section 3: OUTPUT (bottom)
┌─────────────────────────────────────────────────────────┐
│                       OUTPUT                             │
│  ┌─────────────┐                    ┌─────────────┐     │
│  │ Filtered    │                    │  Allele     │     │
│  │    BAM      │                    │  counts     │     │
│  │ (bias-free) │                    │ (ASE-ready) │     │
│  └─────────────┘                    └─────────────┘     │
└─────────────────────────────────────────────────────────┘

## Icons to Include (BioRender library)
- BAM box: DNA helix or sequencing reads icon
- VCF box: Variant/SNP icon (colored dots on line)
- Remap box: Alignment/mapping icon
- Filter: Funnel icon (small, near WASP filter step)
- Output BAM: Clean reads icon
- Allele counts: Bar chart or histogram icon

## Typography
- Section titles: 8pt bold, dark gray
- Box labels: 6pt medium weight
- Annotations: 5pt regular, light gray

## Dimensions
- Total width: ~88mm (single column) or fit within 60mm for panel
- Total height: ~70mm
- Rounded corners on all boxes (radius: 3-4pt)

## Special Elements
1. Dashed section borders to group related steps
2. Merge arrows from BAM+VCF into first process step
3. Split arrows from last process step to two outputs
4. Small "WASP filter" annotation with arrow pointing to remap step
5. Optional: Small "SNVs + INDELs" label near haplotype generation

## Style Notes
- Clean, minimal design (no gradients or shadows)
- Consistent spacing between elements
- Professional scientific aesthetic
- Should complement benchmark plots in panels B and C
