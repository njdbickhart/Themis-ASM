import subprocess as sp
import numpy as np
import os
import sys

usage = f'python {sys.argv[0]} <BED of gap regions> <Fasta index file> <Output tab file name>'

if len(sys.argv) != 4:
    print(usage)
    sys.exit(-1)



# create chromosome length bed file
with open(sys.argv[2], 'r') as fai, open(sys.argv[2] + '.bed', 'w') as out:
    for l in fai:
        s = l.rstrip().split()
        out.write(f'{s[0]}\t1\t{s[1]}\n')

# define some Variables
sizes = []

# Run bedtools subtraction and process output
cmd = ['bedtools', 'subtract', '-b', sys.argv[1], '-a', f'{sys.argv[2]}.bed']
result = sp.run(cmd, stdout=sp.PIPE ).stdout.decode('utf-8')
for l in result.split('\n'):
    s = l.rstrip().split()
    if len(s) < 3:
        break
    sizes.append(int(s[2]) - int(s[1]))

# calculate final values and print output
sizes.sort(reverse=True)

count = len(sizes)
sum = np.sum(sizes) if count > 0 else 0
mean = np.mean(sizes) if count > 0 else 0
median = np.median(sizes) if count > 0 else 0
stdev = np.std(sizes) if count > 0 else 0

l50v = int(sum / 2)
csum = np.cumsum(sizes)
l50idx = min(csum[csum >= l50v])
c50idx = np.where(csum == l50idx)

n50 = sizes[int(c50idx[0])]
with open(sys.argv[3], 'w') as out:
    out.write("Num\tSumBp\tMeanBp\tMedianBp\tStdevBp\tContigN50\n")
    out.write(f'{count}\t{sum}\t{mean}\t{median}\t{stdev}\t{n50}\n')
