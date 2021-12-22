#!/usr/bin/env python3

import os
import errno
import sys
import shutil
import subprocess as sp
import argparse
from collections import defaultdict, deque


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
                        help="Optional integer parameter to filter minimum alignment lengths [0.05 of total assembly size]",
                        default=-1, type=int,
                        )
    parser.add_argument('-c', '--maxchr',
                        help="Maximum number of chromosomes to plot per assembly.",
                        required=True, type=int,
                        )

    return parser.parse_args(), parser

def main(args, parser):
    print("Not done!")

    # Create main configuration workhorse and generate correspondance tags
    cConf = CircosConf(args.output)
    cConf.createCorrespondance(args.reference, args.query, args.minimum, args.bedfile)

    if args.bedfile != "None":
        cConf.readTargets(args.bedfile)

    # Now create Karyotype file
    kFile = KaryotypeFile(cConf.targetList, cConf.corrKey, cConf.refLenDict, args.output)
    kFile.filterTargetList(args.maxchr)

    kFile.readPAF(args.paf, args.minimum)
    kFile.writeKaryotype(args.minimum)

    # Now create the Link file
    lFile = LinkFile(kFile.algnDict, cConf.targetList, cConf.corrKey, args.output)
    lFile.writeLinks(args.minimum)

    # Finally, create the remaining configuration files
    cConf.createTicksFile()
    cConf.createIdeogramFile()
    cConf.createConfFile(kFile.getIdeogramList(), kFile.generateRuleText())

    # Assuming everything worked, try to run circos
    cConf.run('circos')


class CircosConf:

    def __init__(self, outDir):
        self.outDir = outDir
        self.colors = deque([f'chr{x}' for x in range(1, 25)])

        self.targetList = []
        self.corrKey = defaultdict(dict)
        self.refLenDict = {}
        self.qLenDict = {}

    def replacePattern(self, input_text, patterns, replacements):
        for p, r in zip(patterns, replacements):
            input_text.replace(p, r)

        return input_text

    def createConfFile(self, ideogramtext, rulestext):
        with open(os.path.join(self.outDir, "circos.conf"), 'w') as output:
            output.write(self.replacePattern(CIRCOS,
            ["<CHROMOSOMES_WILL_GO_HERE>", "<RULES_GO_HERE>"],
            [ideogramtext, rulestext]))

    def createTicksFile(self):
        with open(os.path.join(self.outDir, "ticks.conf"), 'w') as output:
            output.write(TICKS)

    def createIdeogramFile(self):
        with open(os.path.join(self.outDir, "ideogram.conf"), 'w') as output:
            output.write(IDEOGRAM)

    def readTargets(self, bedfile):
        if bedfile != "None":
            with open(bedfile, 'r') as input:
                for l in input:
                    s = l.rstrip().split()
                    if len(s) == 0 or s[0] == "" :
                        continue
                    ctg = s[0]
                    start = int(s[1])
                    end = int(s[2])
                    color = self.colors
                    self.colors.rotate(1)
                    try :
                        color = s[3]
                    except :
                        pass
                    self.targetList.append(Target(ctg, self.corrKey["QUERY"][ctg], start, end, color, self.qLenDict[s[0]]))


    def createCorrespondance(self, ref, query, min_align_length, bedfile):
        n = 0
        for f, t in zip([ref, query], ["REF", "QUERY"]):
            if not os.path.exists(f + '.fai'):
                print(f'Building fasta index for {t} file: {f}')
                sp.Popen(f'samtools faidx {f}', shell=True)

            with open(f + '.fai', 'r') as input:
                for l in input:
                    s = l.rstrip().split()
                    if int(s[1]) < min_align_length:
                        continue
                    self.corrKey[t][s[0]] = f'av{n}'
                    if t == "REF":
                        self.refLenDict[s[0]] = int(s[1])
                    else:
                        self.qLenDict[s[0]] = int(s[1])
                    if bedfile == "None":
                        color = self.colors
                        self.colors.rotate(1)
                        self.targetList.append(Target(s[0], self.corrKey[t][s[0]], 1, int(s[1]), color, int(s[1])))
                    n += 1

    def run(self, cmd) :
        cwd = os.getcwd()
        os.chdir(os.path.join(self.outdir))
        print("Running circos in {}".format(os.getcwd()))
        proc = sp.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        proc.communicate()
        print("Done!")
        os.chdir(cwd)

class LinkFile:

    def __init__(self, algnDict, targetList, corrKey, outDir):
        """ """
        self.algnDict = algnDict
        self.targetList = targetList
        self.corrKey = corrKey
        self.outDir = outDir

        self.pairs = []

    def writeLinks(self, min_align_length):
        with open(os.path.join(self.outDir, "links.txt"), 'w') as output:
            for p, target in enumerate(self.targetList) :
                if target.name not in self.algnDict.keys() :
                    continue
                for n, aln in enumerate(self.algnDict[target.name]) :
                    if aln.length < min_align_length :
                        continue
                    link1 = "link{} {} {} {}".format(str(p)+"_"+str(n), target.tag, aln.T_start, aln.T_end)
                    output.write(link1+"\n")
                    link2 = "link{} {} {} {}".format(str(p)+"_"+str(n), self.corrKey["QUERY"][aln.query], aln.Q_start, aln.Q_end)
                    output.write(link2+"\n")
                    self.pairs.append(Coords(target.name, aln.T_start, aln.T_end, aln.query, aln.Q_start, aln.Q_end))


class KaryotypeFile:

    def __init__(self, targetList, corrKey, refLenDict, outDir):
        """ """
        self.targetList = targetList
        self.corrKey = corrKey
        self.refLenDict = refLenDict
        self.outDir = outDir

        self.targets_to_plot = set()
        self.algnDict = {}
        self.karyotypes = []
        self.ideogramList = []

    def getTargetList(self):
        return self.targetList

    def getIdeogramList(self):
        return ';'.join(self.ideogramList)

    def filterTargetList(self, maxChr):
        sorted_chromosomes = [x for x, j in sorted(self.refLenDict.items(), key=lambda item : item[1], reverse=True)]
        self.targets_to_plot = set(sorted_chromosomes[:maxChr])

    def generateRuleText(self):
        text = ''
        for t in self.targetList:
            if t.name in self.targets_to_plot:
                text += f'<rule>\ncondition = to({t.tag})\ncolor={t.color}\nstroke_color={t.color}\n</rule>\n\n'

        return text

    def readPAF(self, paf, min_align_length):
        with open(paf, 'r') as input:
            for l in input:
                s = l.rstrip().split()
                if s[5] not in self.targets_to_plot:
                    continue
                t_len = int(s[6])
                query = s[0]
                q_len = int(s[1])
                q_start = int(s[2])
                q_end = int(s[3])
                t_start = int(s[7])
                t_end = int(s[8])
                length = int(s[10])

                if length <= min_align_length:
                    continue
                if target not in self.algnDict.keys() :
                    self.algnDict[target] = [Alignment(query, q_start, q_end, t_start, t_end, length)]
                else :
                    self.algnDict[target].append(Alignment(query, q_start, q_end, t_start, t_end, length))

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

if __name__ == "__main__":
    args, parser = parse_user_input()
    main(args, parser)
