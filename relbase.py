import sys
import random
import math
import matplotlib.pyplot as plt
import matplotlib
import re
import numpy as np

# Input files are in FASTA format. The first sequence is used as a consensus
# while the following sequences are variants of this consensus sequence
# containing mutations. All sequences have the same length and there are no
# indels.
# There are two kinds of sequences
# - control sequence  (Ctrl*.txt)
# - experimental data (Brca2*.txt)
# We count
# - the total number of mutations compared to consensus
# - the kind of mutations
# The result will be graphics which show the mutation frequency for a specific
# kind of mutation after a specific "pseudo-time", which is equivalent to the
# total number of mutations in a sequence.

# Usage:
#  python3 relbase.py --ctrl in/Ctrl_* --exp in/Brca2_*
#  Output in out/: 
#    RelBase_AC.png, RelBase_AT.png ... - specific kind of mutation
#    RelBase.png - all kinds of mutations combined in one image



# order we use to display
cgat = ['G','C','A','T']

# read fasta file
# returns array with all sequences in file (seq as string)
def read_fasta(f):
    f = open(f,'r')
    seq = []
    while True:
        line = f.readline()
        if not line:
            break # done
        if line.startswith('>'):
            continue   # comment line

        # line with sequence
        line = line.rstrip()
        seq.append(line)
    return seq


# Count and accumulate types of mutations of sequences against a consensus
# sequence. dst['mut'] will be a list of buckets with one bucket containing
# a list of information about each sequence which has this specific number of
# mutations. These information are a dictionary k,v where key being the type of
# mutation ('AG','TC'...) and v the number of mutations for this type in this
# sequence. dst['base'] will be a count of the bases in the consensus sequence
# for each base type.
def count_mutations(consensus,seq,dst):
    for b in cgat: 
        dst['base'][b] = 0
    for b in consensus: 
        dst['base'][b] += 1

    mut = dst['mut']
    for s in seq:
        if re.search('[^ACGT]',s):
            print("ignore line " + s)
            continue

        n = 0
        ctr = {}
        for i in range(len(consensus)):
            if consensus[i] != s[i]:
                n += 1
                id = consensus[i] + s[i]
                if id in ctr:
                    ctr[id] += 1
                else:
                    ctr[id] = 1

        if n>=len(mut):
            for i in range(len(mut)-1,n):
                mut.append([])
        mut[n].append(ctr)


# This is used instead of sampling, i.e. each sequence in the bucket is returned
# instead of more or less randomly drawn and averaged sample of sequences.
def select_all(data, percent = None, nmuts = None):
    mut = data['mut']
    if not nmuts:
        nmuts = range(1,len(mut))
    xys = {}
    for n in nmuts:
        if not mut[n]:
            continue
        tot = {}
        t = 0
        for s in mut[n]:
            for k,v in s.items():
                t += v
                if k in tot:
                    tot[k] += v
                else:
                    tot[k] = v
        for k,v in tot.items():
            if not percent:
                True
            elif percent == 'mut' :
                v = v*100/(n*len(mut[n]))
            elif percent == 'base' :
                v = v*100/(len(mut[n])*data['base'][k[0]])
            else:
                raise ValueError("invalid value for percent: " + percent)
            if k in xys:
                xys[k].append([ n, v, 2 ])
            else:
                xys[k] = [ [ n, v, 2 ] ]
    return xys


# scatter plotting
def plot_scatter(ax, k, xys, xlim = [1,7], ylim = [0,25]):
    cmap = [ [ '#b0b0b0','o', xys[0] ], [ 'black','*', xys[1] ] ]
    # cmap = { 'blue': xys[0], 'red': xys[1] }
    x = []
    y = []
    s = []
    c = []
    m = []
    for cmx in cmap:
        color,marker,xysc = cmx
        if not k in xysc: continue
        for xi,yi,si,*_ in xysc[k]:
            if xlim and (xlim[0]>xi or xlim[1]<xi):
                continue
            x.append(xi)
            y.append(yi)
            s.append(si*si*50+2)
            c.append(color)
            m.append(marker)
    if not len(x):
        return
    if ylim: 
        check_bounds(y,ylim)
        ax.set_ylim(ylim)
    if xlim: 
        ax.xaxis.set_ticks(np.arange(xlim[0], xlim[1]+0.1, 1))
        ax.set_xlim(xlim)
    ax.scatter(x,y,s,c,alpha=0.5) # , alpha = 10/len(x) if len(x)>100 else 0.5)
    ax.set_title(k[0] + '>' + k[1])
            

# used to make sure that ylim settings do not cut off values from display
def check_bounds(a,lim):
    if min(a)<lim[0]:
        raise ValueError("minimum={} less than lim[0]={}".
            format(min(a),lim[0]))
    if max(a)>lim[1]:
        raise ValueError("maximum={} greater than lim[1]={}".
            format(max(a),lim[1]))


# create the image files for xys
def create_images(title, name, ylim, xlim, xys):
    print(title)
    matplotlib.rcParams.update({'font.size': 20})

    # one image with all sub-images included
    f, ax = plt.subplots(4,4, sharex=True, sharey=True, figsize=(20,20), dpi=100)
    for l in range(len(cgat)):
        for r in range(len(cgat)):
            if l == r: continue
            k = cgat[l]+cgat[r]
            print('   ' + k)
            plot_scatter(ax[l][r], k, xys, xlim = xlim, ylim = ylim)
    f.suptitle(title)
    plt.savefig('out/' + name + '.png')

    # one file per sub-image
    for l in range(len(cgat)):
        for r in range(len(cgat)):
            if l == r: continue
            k = cgat[l]+cgat[r]
            print('   ' + k)
            f,ax = plt.subplots(1,1, figsize=(5,5), dpi=100)
            plot_scatter(ax, k, xys, xlim = xlim, ylim = ylim)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title('')
            plt.savefig('out/' + name + '_' + k + '.png')


# ctrl and exp
# base: count of each base in consensus sequence (assumes that same consensus is
#   used for all sequences in all files of same type)
# mut: mutations per type and bucket
# see description of count_mutations for more details
ctrl = { 'base': {}, 'mut': [] }
exp  = { 'base': {}, 'mut': [] }

# read FASTA files into ctrl and exp
for f in sys.argv[1:]:
    if f == '--ctrl':
        dst = ctrl
    elif f ==  '--exp':
        dst = exp
    else:
        consensus,*seq = read_fasta(f)
        count_mutations(consensus,seq,dst)

create_images(
    title = "(RelBase) Mutations per Bucket relative to Count of Origin Base",
    name  = 'RelBase',
    ylim  = [0, 1.1],
    xlim  = [1, 5],
    xys   = [ select_all(ctrl, 'base'), select_all(exp, 'base') ],
)


