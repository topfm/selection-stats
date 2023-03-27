#!/usr/bin/env python

import sys
import os
import getopt
import egglib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# This script reads in an alignment and calculates diversity and selection
# statistics based on the window width and window step given by the user.
# Will calculate Fay and Wu's H if an outgroup is provided.

#####
# Edited by M. Youngblom 07/2019 to work with EggLib v3.
# Edited by M. Topf 07/2021 to include work around if theta or pi is NA
#####

# get arguments
def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    alignment = None
    winWidth = 1000
    winStep = 300
    outgroup = None
    try:
        opts, args = getopt.getopt(argv, "a:w:s:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-a':
            alignment = arg
        elif opt == '-w':
            winWidth = int(arg) 
        elif opt == '-s':
            winStep = int(arg)
        elif opt == '-o':
            outgroup = arg
    return (alignment, winWidth, winStep, outgroup)

# usage
def usage():
    print("slidingWindowStats.py\n \
        -a <fasta alignment>\n \
        -w <window width default = 1000>\n \
        -s <window step default = 300>\n \
        -o <outgroup>")

# function to calculate popgen stats using egglib
    # use add_stats to add statistics to be calculated
    # have to divide theta & pi by window size to get per-site value
def calc_stats(a, win):
    statDict = {}
    cs = egglib.stats.ComputeStats()
    cs.add_stats('thetaW', 'Pi', 'D', 'Hsd')
    polyDict = cs.process_align(a, max_missing= .2)
    try:  
        statDict['theta'] = float(polyDict['thetaW'])/int(win)
    except TypeError:
        statDict['theta'] = "NA"
    try:
        statDict['pi'] = float(polyDict['Pi'])/int(win)
    except TypeError:
        statDict['pi'] = "NA"
    statDict['tajimaD'] = polyDict['D']
    statDict['H'] = polyDict['Hsd']
    return statDict

# assign argument values
alignment, winWidth, winStep, outgroup = get_arguments(sys.argv[1:])

# alignment file is required
if alignment is None:
    usage()
    sys.exit()

# write outfile containing header 
outfile = open('windowStats_' + os.path.splitext(alignment)[0] + '.txt', 'w')
outfile.write("Start\tStop\tTheta\tPi\tTajimasD\tFay&WuH\n")

# read in alignment with egglib
align = egglib.io.from_fasta(alignment, alphabet=egglib.alphabets.DNA)

# ? something to do with the outgroup, must test to see this part is still working
for i in range(align.ns):
    align.get_sequence(i)
if outgroup is not None:
    align.group(align.find(outgroup, strict = False), group = 999)

# first window starts at zero, stops at winWidth
start = 0
stop = winWidth
# used for plotting D across genome
location = []
TD = []

# calculating stats for each window and writing to output 
for window in align.slider(winWidth, winStep):
    stats = calc_stats(window, winWidth)
    outfile.write("%i\t%i\t%s\t%s\t%s\t%s\n" % (start, stop, stats['theta'],
    stats['pi'], stats['tajimaD'], stats['H']))
    location.append((start + stop)/2)
    TD.append(stats['tajimaD'])
    start += winStep
    stop += winStep
outfile.close()

plt.plot(location, TD)
plt.xlabel('Location')
plt.ylabel('Tajima\'s D')
plt.savefig("slidingWindowTajimasD.png")
plt.close()
