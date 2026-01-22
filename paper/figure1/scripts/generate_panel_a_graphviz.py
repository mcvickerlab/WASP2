#!/usr/bin/env python3
"""
Figure 1 Panel A: WASP2 Architecture Diagram (Graphviz Version)

Alternative implementation using Graphviz for publication-quality pipeline diagrams.
Features:
- Left-to-right flow (rankdir='LR')
- Three main modules: Read Filtering, Variant Counting, AI Detection
- Input nodes: BAM, VCF
- Output nodes: Filtered BAM, Allele Counts, AI Results
- Sub-steps within each module as nested clusters
- Colorblind-safe colors matching existing palette
- 300 DPI output
- Professional rounded boxes with sans-serif fonts

Nature Methods compliant: 180mm width, 300dpi, 7-8pt fonts.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from graphviz import Digraph
from config import get_plot_path, REPO_ROOT

# =============================================================================
# COLOR SCHEME (Colorblind-safe - matching existing Panel A)
# =============================================================================
COLORS = {
    # Main components (Paul Tol bright palette)
    'comp1': '#4477AA',         # Blue - Read Filtering
    'comp1_bg': '#DCE9F4',
    'comp2': '#228833',         # Green - Counting
    'comp2_bg': '#D9EDDB',
    'comp3': '#CCBB44',         # Yellow/Gold - AI Detection
    'comp3_bg': '#F5F2DC',

    # Input/Output
    'input': '#555555',
    'input_bg': '#E8E8E8',
    'output_bg': '#FFFFFF',

    # Container
    'container_bg': '#F5F8FC',
    'container_border': '#B8C9DC',

    # Text and arrows
    'text': '#222222',
    'text_light': '#555555',
    'arrow': '#666666',
    'white': '#FFFFFF',
}

# User-specified colors (requested in task)
REQUESTED_COLORS = {
    'blue': '#0072B2',
    'green': '#009E73',
    'orange': '#E69F00',
}


def create_wasp2_pipeline_diagram():
    """Create the WASP2 architecture diagram using Graphviz."""

    # Create directed graph with left-to-right layout
    dot = Digraph(
        'WASP2_Architecture',
        format='pdf',
        engine='dot',
    )

    # Global graph attributes for publication quality
    dot.attr(
        rankdir='LR',
        splines='ortho',  # Orthogonal edges for clean look
        nodesep='0.4',
        ranksep='0.6',
        compound='true',
        fontname='Helvetica',
        fontsize='8',
        bgcolor='white',
        dpi='300',
        size='7.09,4',  # Nature Methods double column width
        ratio='compress',
        pad='0.3',
    )

    # Default node attributes
    dot.attr(
        'node',
        shape='box',
        style='rounded,filled',
        fontname='Helvetica',
        fontsize='7',
        penwidth='1.2',
        margin='0.15,0.08',
    )

    # Default edge attributes
    dot.attr(
        'edge',
        color=COLORS['arrow'],
        arrowsize='0.7',
        penwidth='1.5',
        fontname='Helvetica',
        fontsize='6',
    )

    # =========================================================================
    # INPUT NODES
    # =========================================================================
    with dot.subgraph(name='cluster_inputs') as inputs:
        inputs.attr(
            label='Inputs',
            labelloc='t',
            labeljust='l',
            fontsize='8',
            fontcolor=COLORS['text'],
            style='rounded,filled',
            fillcolor=COLORS['container_bg'],
            color=COLORS['container_border'],
            penwidth='1.5',
        )

        inputs.node(
            'bam_input',
            'BAM',
            fillcolor=COLORS['input_bg'],
            color=COLORS['input'],
            fontcolor=COLORS['text'],
            fontsize='8',
            fontname='Helvetica-Bold',
        )

        inputs.node(
            'vcf_input',
            'VCF/BCF',
            fillcolor=COLORS['input_bg'],
            color=COLORS['input'],
            fontcolor=COLORS['text'],
            fontsize='8',
            fontname='Helvetica-Bold',
        )

    # =========================================================================
    # MODULE 1: READ FILTERING (Blue)
    # =========================================================================
    with dot.subgraph(name='cluster_filtering') as filtering:
        filtering.attr(
            label='1  Read Filtering',
            labelloc='t',
            labeljust='l',
            fontsize='8',
            fontcolor=COLORS['comp1'],
            fontname='Helvetica-Bold',
            style='rounded,filled',
            fillcolor=COLORS['comp1_bg'],
            color=COLORS['comp1'],
            penwidth='1.5',
        )

        # Sub-steps within Read Filtering
        substeps_1 = [
            ('f1', 'Identify Variant Reads'),
            ('f2', 'Swap Alleles'),
            ('f3', 'Remap Sequences'),
            ('f4', 'Concordance Filter'),
        ]

        for node_id, label in substeps_1:
            filtering.node(
                node_id,
                label,
                fillcolor=COLORS['comp1'],
                color=COLORS['comp1'],
                fontcolor=COLORS['white'],
                fontsize='6',
                fontname='Helvetica-Bold',
            )

        # Feature note
        filtering.node(
            'f_note',
            'Any aligner supported',
            shape='plaintext',
            fillcolor='transparent',
            fontcolor=COLORS['text_light'],
            fontsize='5',
            fontname='Helvetica-Oblique',
        )

        # Internal edges (vertical flow within module)
        filtering.edge('f1', 'f2', style='invis')
        filtering.edge('f2', 'f3', style='invis')
        filtering.edge('f3', 'f4', style='invis')
        filtering.edge('f4', 'f_note', style='invis')

    # =========================================================================
    # MODULE 2: VARIANT COUNTING (Green)
    # =========================================================================
    with dot.subgraph(name='cluster_counting') as counting:
        counting.attr(
            label='2  Variant Counting',
            labelloc='t',
            labeljust='l',
            fontsize='8',
            fontcolor=COLORS['comp2'],
            fontname='Helvetica-Bold',
            style='rounded,filled',
            fillcolor=COLORS['comp2_bg'],
            color=COLORS['comp2'],
            penwidth='1.5',
        )

        # Sub-steps within Variant Counting
        substeps_2 = [
            ('c1', 'Filter Het Variants'),
            ('c2', 'Feature Intersection'),
            ('c3', 'Allele Counting'),
        ]

        for node_id, label in substeps_2:
            counting.node(
                node_id,
                label,
                fillcolor=COLORS['comp2'],
                color=COLORS['comp2'],
                fontcolor=COLORS['white'],
                fontsize='6',
                fontname='Helvetica-Bold',
            )

        # Feature notes
        counting.node(
            'c_note1',
            'Per-SNV & feature-level',
            shape='plaintext',
            fillcolor='transparent',
            fontcolor=COLORS['text_light'],
            fontsize='5',
            fontname='Helvetica-Oblique',
        )
        counting.node(
            'c_note2',
            'Phased haplotype support',
            shape='plaintext',
            fillcolor='transparent',
            fontcolor=COLORS['text_light'],
            fontsize='5',
            fontname='Helvetica-Oblique',
        )

        # Internal edges
        counting.edge('c1', 'c2', style='invis')
        counting.edge('c2', 'c3', style='invis')
        counting.edge('c3', 'c_note1', style='invis')
        counting.edge('c_note1', 'c_note2', style='invis')

    # =========================================================================
    # MODULE 3: AI DETECTION (Yellow/Gold)
    # =========================================================================
    with dot.subgraph(name='cluster_ai') as ai:
        ai.attr(
            label='3  AI Detection',
            labelloc='t',
            labeljust='l',
            fontsize='8',
            fontcolor=COLORS['comp3'],
            fontname='Helvetica-Bold',
            style='rounded,filled',
            fillcolor=COLORS['comp3_bg'],
            color=COLORS['comp3'],
            penwidth='1.5',
        )

        # Sub-steps within AI Detection
        substeps_3 = [
            ('a1', 'Aggregate Counts'),
            ('a2', 'Beta-Binomial Model'),
            ('a3', 'FDR Correction'),
        ]

        for node_id, label in substeps_3:
            ai.node(
                node_id,
                label,
                fillcolor=COLORS['comp3'],
                color=COLORS['comp3'],
                fontcolor=COLORS['white'],
                fontsize='6',
                fontname='Helvetica-Bold',
            )

        # Feature notes
        ai.node(
            'a_note1',
            'Multi-SNV aggregation',
            shape='plaintext',
            fillcolor='transparent',
            fontcolor=COLORS['text_light'],
            fontsize='5',
            fontname='Helvetica-Oblique',
        )
        ai.node(
            'a_note2',
            'Likelihood ratio test',
            shape='plaintext',
            fillcolor='transparent',
            fontcolor=COLORS['text_light'],
            fontsize='5',
            fontname='Helvetica-Oblique',
        )

        # Internal edges
        ai.edge('a1', 'a2', style='invis')
        ai.edge('a2', 'a3', style='invis')
        ai.edge('a3', 'a_note1', style='invis')
        ai.edge('a_note1', 'a_note2', style='invis')

    # =========================================================================
    # OUTPUT NODES
    # =========================================================================
    with dot.subgraph(name='cluster_outputs') as outputs:
        outputs.attr(
            label='Outputs',
            labelloc='t',
            labeljust='l',
            fontsize='8',
            fontcolor=COLORS['text'],
            style='rounded,filled',
            fillcolor=COLORS['container_bg'],
            color=COLORS['container_border'],
            penwidth='1.5',
        )

        outputs.node(
            'filtered_bam',
            'Filtered\\nBAM',
            fillcolor=COLORS['comp1_bg'],
            color=COLORS['comp1'],
            fontcolor=COLORS['text'],
            fontsize='7',
            fontname='Helvetica-Bold',
        )

        outputs.node(
            'allele_counts',
            'Allele\\nCounts',
            fillcolor=COLORS['comp2_bg'],
            color=COLORS['comp2'],
            fontcolor=COLORS['text'],
            fontsize='7',
            fontname='Helvetica-Bold',
        )

        outputs.node(
            'ai_results',
            'AI\\nResults',
            fillcolor=COLORS['comp3_bg'],
            color=COLORS['comp3'],
            fontcolor=COLORS['text'],
            fontsize='7',
            fontname='Helvetica-Bold',
        )

    # =========================================================================
    # MAIN FLOW EDGES (between modules)
    # =========================================================================

    # Inputs to Module 1
    dot.edge('bam_input', 'f1', lhead='cluster_filtering')
    dot.edge('vcf_input', 'f1', lhead='cluster_filtering')

    # Module 1 to Module 2
    dot.edge('f4', 'c1', ltail='cluster_filtering', lhead='cluster_counting')

    # Module 2 to Module 3
    dot.edge('c3', 'a1', ltail='cluster_counting', lhead='cluster_ai')

    # Module 1 to Output (Filtered BAM)
    dot.edge('f4', 'filtered_bam', ltail='cluster_filtering')

    # Module 2 to Output (Allele Counts)
    dot.edge('c3', 'allele_counts', ltail='cluster_counting')

    # Module 3 to Output (AI Results)
    dot.edge('a3', 'ai_results', ltail='cluster_ai')

    return dot


def create_wasp2_pipeline_simple():
    """Create a simpler, cleaner WASP2 architecture diagram."""

    dot = Digraph(
        'WASP2_Architecture',
        format='pdf',
        engine='dot',
    )

    # Global graph attributes
    dot.attr(
        rankdir='LR',
        splines='spline',
        nodesep='0.5',
        ranksep='0.8',
        compound='true',
        fontname='Helvetica',
        fontsize='8',
        bgcolor='white',
        dpi='300',
        size='7.09,3.5',
        pad='0.2',
    )

    # Default node attributes
    dot.attr(
        'node',
        shape='box',
        style='rounded,filled',
        fontname='Helvetica',
        fontsize='7',
        penwidth='1.2',
        margin='0.12,0.06',
    )

    # Default edge attributes
    dot.attr(
        'edge',
        color=COLORS['arrow'],
        arrowsize='0.8',
        penwidth='1.5',
    )

    # =========================================================================
    # INPUT NODES (rank = same to align vertically)
    # =========================================================================
    with dot.subgraph() as s:
        s.attr(rank='same')
        s.node(
            'bam',
            'BAM',
            fillcolor=COLORS['input_bg'],
            color=COLORS['input'],
            fontcolor=COLORS['text'],
            fontsize='8',
            fontname='Helvetica-Bold',
        )
        s.node(
            'vcf',
            'VCF/BCF',
            fillcolor=COLORS['input_bg'],
            color=COLORS['input'],
            fontcolor=COLORS['text'],
            fontsize='8',
            fontname='Helvetica-Bold',
        )

    # =========================================================================
    # MODULE 1: READ FILTERING
    # =========================================================================
    with dot.subgraph(name='cluster_m1') as m1:
        m1.attr(
            label='<<B>1</B>  <B>Read Filtering</B>>',
            labelloc='t',
            fontsize='8',
            fontcolor=COLORS['comp1'],
            style='rounded,filled',
            fillcolor=COLORS['comp1_bg'],
            color=COLORS['comp1'],
            penwidth='2',
            margin='12',
        )

        m1.node('m1_s1', 'Identify Variant Reads', fillcolor=COLORS['comp1'], fontcolor='white', fontsize='6', fontname='Helvetica-Bold')
        m1.node('m1_s2', 'Swap Alleles', fillcolor=COLORS['comp1'], fontcolor='white', fontsize='6', fontname='Helvetica-Bold')
        m1.node('m1_s3', 'Remap Sequences', fillcolor=COLORS['comp1'], fontcolor='white', fontsize='6', fontname='Helvetica-Bold')
        m1.node('m1_s4', 'Concordance Filter', fillcolor=COLORS['comp1'], fontcolor='white', fontsize='6', fontname='Helvetica-Bold')
        m1.node('m1_note', 'Any aligner supported', shape='plaintext', fontcolor=COLORS['text_light'], fontsize='5', fontname='Helvetica-Oblique')

        # Vertical internal flow
        m1.edge('m1_s1', 'm1_s2', color=COLORS['comp1'], arrowsize='0.5', penwidth='1.0')
        m1.edge('m1_s2', 'm1_s3', color=COLORS['comp1'], arrowsize='0.5', penwidth='1.0')
        m1.edge('m1_s3', 'm1_s4', color=COLORS['comp1'], arrowsize='0.5', penwidth='1.0')
        m1.edge('m1_s4', 'm1_note', style='invis')

    # =========================================================================
    # MODULE 2: VARIANT COUNTING
    # =========================================================================
    with dot.subgraph(name='cluster_m2') as m2:
        m2.attr(
            label='<<B>2</B>  <B>Variant Counting</B>>',
            labelloc='t',
            fontsize='8',
            fontcolor=COLORS['comp2'],
            style='rounded,filled',
            fillcolor=COLORS['comp2_bg'],
            color=COLORS['comp2'],
            penwidth='2',
            margin='12',
        )

        m2.node('m2_s1', 'Filter Het Variants', fillcolor=COLORS['comp2'], fontcolor='white', fontsize='6', fontname='Helvetica-Bold')
        m2.node('m2_s2', 'Feature Intersection', fillcolor=COLORS['comp2'], fontcolor='white', fontsize='6', fontname='Helvetica-Bold')
        m2.node('m2_s3', 'Allele Counting', fillcolor=COLORS['comp2'], fontcolor='white', fontsize='6', fontname='Helvetica-Bold')
        m2.node('m2_note1', 'Per-SNV & feature-level', shape='plaintext', fontcolor=COLORS['text_light'], fontsize='5', fontname='Helvetica-Oblique')
        m2.node('m2_note2', 'Phased haplotype support', shape='plaintext', fontcolor=COLORS['text_light'], fontsize='5', fontname='Helvetica-Oblique')

        m2.edge('m2_s1', 'm2_s2', color=COLORS['comp2'], arrowsize='0.5', penwidth='1.0')
        m2.edge('m2_s2', 'm2_s3', color=COLORS['comp2'], arrowsize='0.5', penwidth='1.0')
        m2.edge('m2_s3', 'm2_note1', style='invis')
        m2.edge('m2_note1', 'm2_note2', style='invis')

    # =========================================================================
    # MODULE 3: AI DETECTION
    # =========================================================================
    with dot.subgraph(name='cluster_m3') as m3:
        m3.attr(
            label='<<B>3</B>  <B>AI Detection</B>>',
            labelloc='t',
            fontsize='8',
            fontcolor=COLORS['comp3'],
            style='rounded,filled',
            fillcolor=COLORS['comp3_bg'],
            color=COLORS['comp3'],
            penwidth='2',
            margin='12',
        )

        m3.node('m3_s1', 'Aggregate Counts', fillcolor=COLORS['comp3'], fontcolor='white', fontsize='6', fontname='Helvetica-Bold')
        m3.node('m3_s2', 'Beta-Binomial Model', fillcolor=COLORS['comp3'], fontcolor='white', fontsize='6', fontname='Helvetica-Bold')
        m3.node('m3_s3', 'FDR Correction', fillcolor=COLORS['comp3'], fontcolor='white', fontsize='6', fontname='Helvetica-Bold')
        m3.node('m3_note1', 'Multi-SNV aggregation', shape='plaintext', fontcolor=COLORS['text_light'], fontsize='5', fontname='Helvetica-Oblique')
        m3.node('m3_note2', 'Likelihood ratio test', shape='plaintext', fontcolor=COLORS['text_light'], fontsize='5', fontname='Helvetica-Oblique')

        m3.edge('m3_s1', 'm3_s2', color=COLORS['comp3'], arrowsize='0.5', penwidth='1.0')
        m3.edge('m3_s2', 'm3_s3', color=COLORS['comp3'], arrowsize='0.5', penwidth='1.0')
        m3.edge('m3_s3', 'm3_note1', style='invis')
        m3.edge('m3_note1', 'm3_note2', style='invis')

    # =========================================================================
    # OUTPUT NODES (rank = same to align vertically)
    # =========================================================================
    with dot.subgraph() as s:
        s.attr(rank='same')
        s.node(
            'out_bam',
            'Filtered\nBAM',
            fillcolor=COLORS['comp1_bg'],
            color=COLORS['comp1'],
            fontcolor=COLORS['text'],
            fontsize='7',
            fontname='Helvetica-Bold',
        )
        s.node(
            'out_counts',
            'Allele\nCounts',
            fillcolor=COLORS['comp2_bg'],
            color=COLORS['comp2'],
            fontcolor=COLORS['text'],
            fontsize='7',
            fontname='Helvetica-Bold',
        )
        s.node(
            'out_ai',
            'AI\nResults',
            fillcolor=COLORS['comp3_bg'],
            color=COLORS['comp3'],
            fontcolor=COLORS['text'],
            fontsize='7',
            fontname='Helvetica-Bold',
        )

    # =========================================================================
    # MAIN PIPELINE EDGES
    # =========================================================================

    # Inputs to Module 1
    dot.edge('bam', 'm1_s1', lhead='cluster_m1')
    dot.edge('vcf', 'm1_s1', lhead='cluster_m1')

    # Module 1 to Module 2 (main flow)
    dot.edge('m1_s4', 'm2_s1', ltail='cluster_m1', lhead='cluster_m2', penwidth='2.0')

    # Module 2 to Module 3 (main flow)
    dot.edge('m2_s3', 'm3_s1', ltail='cluster_m2', lhead='cluster_m3', penwidth='2.0')

    # Outputs from each module
    dot.edge('m1_s4', 'out_bam', color=COLORS['comp1'], style='dashed', penwidth='1.2')
    dot.edge('m2_s3', 'out_counts', color=COLORS['comp2'], style='dashed', penwidth='1.2')
    dot.edge('m3_s3', 'out_ai', color=COLORS['comp3'], style='dashed', penwidth='1.2')

    return dot


def create_wasp2_pipeline_horizontal():
    """Create a horizontal-flow WASP2 architecture diagram optimized for LR layout."""

    dot = Digraph(
        'WASP2_Architecture',
        format='pdf',
        engine='dot',
    )

    # Global graph attributes for horizontal layout
    dot.attr(
        rankdir='LR',
        splines='polyline',
        nodesep='0.3',
        ranksep='0.7',
        compound='true',
        fontname='Helvetica',
        fontsize='8',
        bgcolor='white',
        dpi='300',
        size='7.09,4.5',
        pad='0.25',
        newrank='true',
    )

    # Default node attributes
    dot.attr(
        'node',
        shape='box',
        style='rounded,filled',
        fontname='Helvetica',
        fontsize='7',
        penwidth='1.2',
        margin='0.1,0.05',
        height='0.35',
    )

    # Default edge attributes
    dot.attr(
        'edge',
        color=COLORS['arrow'],
        arrowsize='0.7',
        penwidth='1.5',
    )

    # =========================================================================
    # INPUT CLUSTER
    # =========================================================================
    with dot.subgraph(name='cluster_input') as inp:
        inp.attr(
            label='<<FONT POINT-SIZE="8"><B>Inputs</B></FONT>>',
            labelloc='t',
            style='rounded,filled',
            fillcolor=COLORS['container_bg'],
            color=COLORS['container_border'],
            penwidth='1.5',
            margin='15',
        )

        inp.node('bam', '<<B>BAM</B>>', fillcolor=COLORS['input_bg'], color=COLORS['input'], fontsize='8')
        inp.node('vcf', '<<B>VCF/BCF</B>>', fillcolor=COLORS['input_bg'], color=COLORS['input'], fontsize='8')

        # Stack inputs vertically
        inp.edge('bam', 'vcf', style='invis')

    # =========================================================================
    # MODULE 1: READ FILTERING (Mapping)
    # =========================================================================
    with dot.subgraph(name='cluster_mapping') as m1:
        m1.attr(
            label='''<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0">
                <TR><TD><FONT COLOR="{0}" POINT-SIZE="9"><B>1</B></FONT></TD>
                <TD><FONT COLOR="{0}" POINT-SIZE="8"><B>Read Filtering</B></FONT></TD></TR>
            </TABLE>>'''.format(COLORS['comp1']),
            labelloc='t',
            style='rounded,filled',
            fillcolor=COLORS['comp1_bg'],
            color=COLORS['comp1'],
            penwidth='2',
            margin='15',
        )

        # Sub-steps as vertical chain
        m1.node('m1_1', '<<B>Identify Variant Reads</B>>', fillcolor=COLORS['comp1'], fontcolor='white', fontsize='6')
        m1.node('m1_2', '<<B>Swap Alleles</B>>', fillcolor=COLORS['comp1'], fontcolor='white', fontsize='6')
        m1.node('m1_3', '<<B>Remap Sequences</B>>', fillcolor=COLORS['comp1'], fontcolor='white', fontsize='6')
        m1.node('m1_4', '<<B>Concordance Filter</B>>', fillcolor=COLORS['comp1'], fontcolor='white', fontsize='6')
        m1.node('m1_feat', '<<I><FONT POINT-SIZE="5" COLOR="{0}">Any aligner supported</FONT></I>>'.format(COLORS['text_light']),
                shape='plaintext', margin='0')

        # Vertical arrows within module
        m1.edge('m1_1', 'm1_2', color=COLORS['comp1'], arrowsize='0.5', penwidth='1.0')
        m1.edge('m1_2', 'm1_3', color=COLORS['comp1'], arrowsize='0.5', penwidth='1.0')
        m1.edge('m1_3', 'm1_4', color=COLORS['comp1'], arrowsize='0.5', penwidth='1.0')
        m1.edge('m1_4', 'm1_feat', style='invis')

    # =========================================================================
    # MODULE 2: VARIANT COUNTING
    # =========================================================================
    with dot.subgraph(name='cluster_counting') as m2:
        m2.attr(
            label='''<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0">
                <TR><TD><FONT COLOR="{0}" POINT-SIZE="9"><B>2</B></FONT></TD>
                <TD><FONT COLOR="{0}" POINT-SIZE="8"><B>Variant Counting</B></FONT></TD></TR>
            </TABLE>>'''.format(COLORS['comp2']),
            labelloc='t',
            style='rounded,filled',
            fillcolor=COLORS['comp2_bg'],
            color=COLORS['comp2'],
            penwidth='2',
            margin='15',
        )

        m2.node('m2_1', '<<B>Filter Het Variants</B>>', fillcolor=COLORS['comp2'], fontcolor='white', fontsize='6')
        m2.node('m2_2', '<<B>Feature Intersection</B>>', fillcolor=COLORS['comp2'], fontcolor='white', fontsize='6')
        m2.node('m2_3', '<<B>Allele Counting</B>>', fillcolor=COLORS['comp2'], fontcolor='white', fontsize='6')
        m2.node('m2_feat1', '<<I><FONT POINT-SIZE="5" COLOR="{0}">Per-SNV &amp; feature-level</FONT></I>>'.format(COLORS['text_light']),
                shape='plaintext', margin='0')
        m2.node('m2_feat2', '<<I><FONT POINT-SIZE="5" COLOR="{0}">Phased haplotype support</FONT></I>>'.format(COLORS['text_light']),
                shape='plaintext', margin='0')

        m2.edge('m2_1', 'm2_2', color=COLORS['comp2'], arrowsize='0.5', penwidth='1.0')
        m2.edge('m2_2', 'm2_3', color=COLORS['comp2'], arrowsize='0.5', penwidth='1.0')
        m2.edge('m2_3', 'm2_feat1', style='invis')
        m2.edge('m2_feat1', 'm2_feat2', style='invis')

    # =========================================================================
    # MODULE 3: AI DETECTION (Analysis)
    # =========================================================================
    with dot.subgraph(name='cluster_analysis') as m3:
        m3.attr(
            label='''<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0">
                <TR><TD><FONT COLOR="{0}" POINT-SIZE="9"><B>3</B></FONT></TD>
                <TD><FONT COLOR="{0}" POINT-SIZE="8"><B>AI Detection</B></FONT></TD></TR>
            </TABLE>>'''.format(COLORS['comp3']),
            labelloc='t',
            style='rounded,filled',
            fillcolor=COLORS['comp3_bg'],
            color=COLORS['comp3'],
            penwidth='2',
            margin='15',
        )

        m3.node('m3_1', '<<B>Aggregate Counts</B>>', fillcolor=COLORS['comp3'], fontcolor='white', fontsize='6')
        m3.node('m3_2', '<<B>Beta-Binomial Model</B>>', fillcolor=COLORS['comp3'], fontcolor='white', fontsize='6')
        m3.node('m3_3', '<<B>FDR Correction</B>>', fillcolor=COLORS['comp3'], fontcolor='white', fontsize='6')
        m3.node('m3_feat1', '<<I><FONT POINT-SIZE="5" COLOR="{0}">Multi-SNV aggregation</FONT></I>>'.format(COLORS['text_light']),
                shape='plaintext', margin='0')
        m3.node('m3_feat2', '<<I><FONT POINT-SIZE="5" COLOR="{0}">Likelihood ratio test</FONT></I>>'.format(COLORS['text_light']),
                shape='plaintext', margin='0')

        m3.edge('m3_1', 'm3_2', color=COLORS['comp3'], arrowsize='0.5', penwidth='1.0')
        m3.edge('m3_2', 'm3_3', color=COLORS['comp3'], arrowsize='0.5', penwidth='1.0')
        m3.edge('m3_3', 'm3_feat1', style='invis')
        m3.edge('m3_feat1', 'm3_feat2', style='invis')

    # =========================================================================
    # OUTPUT CLUSTER
    # =========================================================================
    with dot.subgraph(name='cluster_output') as out:
        out.attr(
            label='<<FONT POINT-SIZE="8"><B>Outputs</B></FONT>>',
            labelloc='t',
            style='rounded,filled',
            fillcolor=COLORS['container_bg'],
            color=COLORS['container_border'],
            penwidth='1.5',
            margin='15',
        )

        out.node('out_bam', '<<B>Filtered<BR/>BAM</B>>', fillcolor=COLORS['comp1_bg'], color=COLORS['comp1'], fontsize='7')
        out.node('out_counts', '<<B>Allele<BR/>Counts</B>>', fillcolor=COLORS['comp2_bg'], color=COLORS['comp2'], fontsize='7')
        out.node('out_ai', '<<B>AI<BR/>Results</B>>', fillcolor=COLORS['comp3_bg'], color=COLORS['comp3'], fontsize='7')

        # Stack outputs vertically
        out.edge('out_bam', 'out_counts', style='invis')
        out.edge('out_counts', 'out_ai', style='invis')

    # =========================================================================
    # MAIN FLOW EDGES
    # =========================================================================

    # Inputs to Module 1
    dot.edge('bam', 'm1_1', lhead='cluster_mapping', penwidth='2.0')
    dot.edge('vcf', 'm1_1', lhead='cluster_mapping', penwidth='2.0')

    # Module 1 to Module 2
    dot.edge('m1_4', 'm2_1', ltail='cluster_mapping', lhead='cluster_counting', penwidth='2.0')

    # Module 2 to Module 3
    dot.edge('m2_3', 'm3_1', ltail='cluster_counting', lhead='cluster_analysis', penwidth='2.0')

    # Per-stage outputs (dashed to show optional/parallel)
    dot.edge('m1_4', 'out_bam', color=COLORS['comp1'], style='dashed', penwidth='1.2', arrowsize='0.6')
    dot.edge('m2_3', 'out_counts', color=COLORS['comp2'], style='dashed', penwidth='1.2', arrowsize='0.6')
    dot.edge('m3_3', 'out_ai', color=COLORS['comp3'], style='dashed', penwidth='1.2', arrowsize='0.6')

    return dot


def generate_panel_a_graphviz():
    """Generate and save the Graphviz-based Panel A diagram."""

    # Create the diagram
    dot = create_wasp2_pipeline_horizontal()

    # Output paths
    out_dir = get_plot_path(1, 'panel_a_graphviz').parent
    out_dir.mkdir(parents=True, exist_ok=True)

    base_path = out_dir / 'panel_a_graphviz'

    # Render PDF
    dot.format = 'pdf'
    pdf_path = dot.render(str(base_path), cleanup=True)
    print(f"Saved: {pdf_path}")

    # Render PNG (300 DPI is set in graph attributes)
    dot.format = 'png'
    png_path = dot.render(str(base_path), cleanup=True)
    print(f"Saved: {png_path}")

    # Also save the DOT source for reference
    dot_source_path = base_path.with_suffix('.dot')
    with open(dot_source_path, 'w') as f:
        f.write(dot.source)
    print(f"Saved: {dot_source_path}")

    return dot


if __name__ == '__main__':
    generate_panel_a_graphviz()
