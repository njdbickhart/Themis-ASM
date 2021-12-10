#!/usr/bin/env python3

import os
import errno
import sys
import shutil
import subprocess
import argparse


def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A tool to convert a PAF file of assembly alignments into a Circos Jupiter-style plot"
            )
    parser.add_argument('-p', '--paf',
                        help="Input paf alignment file",
                        required=True, type=str,
                        )
    parser.add_argument('-q', '--query',
                        help="Query fasta file (should be samtools faidx indexed)",
                        required=True, type=str,
                        )
    parser.add_argument('-r', '--reference',
                        help="Reference fasta file (should be samtools faidx indexed)",
                        required=True, type=str,
                        )
    parser.add_argument('-b', '--bedfile',
                        help="Optional bed file containing regions to plot on both assemblies",
                        default="None", type=str,
                        )
    parser.add_argument('-o', '--output',
                        help="Output directory path",
                        required=True, type=str,
                        )
    parser.add_argument('-m', '--minimum',
                        help="Optional integer parameter to filter minimum alignment lengths [0.01 of total assembly size]",
                        default=-1, type=int,
                        )

    return parser.parse_args(), parser

class KaryotypeFile:

    def __init__(self, algnDict, targetList, corrKey, ctgLenDict, outDir):
        """ """
        self.algnDict = algnDict
        self.targetList = targetList
        self.corrKey = corrKey
        self.ctgLenDict = ctgLenDict
        self.outDir = outDir

        self.karyotypes = []
        self.ideogramList = []

    def writeKaryotype(self, min_align_length):
        skips = 0
        with open(os.path.join(self.outDir, "karyotype.txt"), 'w') as output:
            for target in self.targetList:
                if target not in self.karyotypes:
                    self.karyotypes.append(target)
                    output.write(f'chr - {target.tag} {target.name} 0 {target.length} {target.color}\n')
                else:
                    continue

                queries_aligned = []
                if target.name not in self.algnDict.keys():
                    continue
                for aln in self.algnDict[target.name]:
                    if aln.length < min_align_length :
                        skips += 1
                        continue
                    #
                    #   =========    target region
                    # -------------- query region
                    #
                    elif (aln.T_start <= target.start and aln.T_end >= target.end) :
                        queries_aligned.append(aln.query)

                    #
                    #   =========    target region
                    # ----- or ----- query region
                    #
                    elif (aln.T_start <= target.start and aln.T_end > target.start and aln.T_end <= target.end) or (aln.T_start >= target.start and aln.T_start <= target.end and aln.T_end >= target.end) :
                        queries_aligned.append(aln.query)
                        #print(aln.T_start, aln.T_end) # DEBUG

                    #
                    #   ========= target region
                    #     ----    query region
                    #
                    elif (aln.T_start >= target.start and aln.T_end <= target.end) :
                        queries_aligned.append(aln.query)
                    else :
                        #print(aln.T_start, aln.T_end) # DEBUG
                        continue
                self.ideogramList.append(f'{target.tag}:{target.start}-{target.end}')
                print(f'Skipped {skips} lines in alignment of target {target.name}')


class Target :
    def __init__(self, name, tag, start, end, color, length) :
        """ """
        self.name = name
        self.tag = tag
        self.start = start
        self.end = end
        self.color = color
        self.length = length

    def __str__(self) :
        return "Target(name={}, tag={}, start={}, end={}, color={}, length={})".format(self.name,self.tag,self.start,self.end,self.color,self.length)

class Coords :
    def __init__(self, name_1, name_2, start_1, start_2, end_1, end_2) :
        """ """
        self.n1 = name_1
        self.s1 = start_1
        self.e1 = end_1
        self.n2 = name_2
        self.s2 = start_2
        self.e2 = end_2

    def __eq__(self, co) :
        return (self.n1 == co.n1 or self.n1 == co.n2) and (self.s1 == co.s1 or self.s1 == co.s2) and (self.e1 == co.e1 or self.e1 == co.e2)

class Alignment :
    """ """
    def __init__(self, query, Q_start, Q_end, T_start, T_end, length) :
        self.query = query
        self.Q_start = Q_start
        self.Q_end = Q_end
        self.T_start = T_start
        self.T_end = T_end
        self.length = length


IDEOGRAM = """
<ideogram>

<spacing>
default = 0.005r
</spacing>

# Ideogram position, fill and outline
radius           = 0.725r
thickness        = 20p
fill             = yes
stroke_color     = dgrey
stroke_thickness = 2p

# Minimum definition for ideogram labels.

show_label       = yes
# see etc/fonts.conf for list of font names
label_font       = default
label_radius     = 1.02r
label_size       = 46
label_parallel   = no

</ideogram>
"""

TICKS = """
show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
color            = black
thickness        = 2p

multiplier       = 1e-6

format           = %.1f

<tick>
spacing        = 100000u
size           = 10p
show_label     = no
label_size     = 15p
label_offset   = 5p
</tick>

<tick>
spacing        = 10000u
size           = 5p
show_label     = no
</tick>

</ticks>
"""

CIRCOS = """
karyotype = ./karyotype.txt

# The chromosomes_unit value is used as a unit (suffix "u") to shorten
# values in other parts of the configuration file. Some parameters,
# such as ideogram and tick spacing, accept "u" suffixes, so instead of

chromosomes_units = 1
chromosomes_display_default = no
chromosomes = <CHROMOSOMES_WILL_GO_HERE>

<links>

<link>
file          = ./links.txt
radius        = 0.995r
bezier_radius = 0r
color         = purple_a2
stroke_color  = black
stroke_thickness = 1
thickness     = 2
ribbon	      = yes

<rules>
<RULES_GO_HERE>
</rules>

</link>

</links>

<<include ideogram.conf>>

<<include ticks.conf>>

<image>
<<include etc/image.conf>>
</image>

<<include etc/colors_fonts_patterns.conf>>

<<include etc/housekeeping.conf>>
"""
