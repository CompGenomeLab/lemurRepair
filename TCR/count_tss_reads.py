import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-P', type=str, help='plus XR counts bed file')
parser.add_argument('-M', type=str, help='minus XR counts bed file')
parser.add_argument('--tss', type=str,
                    help='biomart export containing TSS coordinates')
parser.add_argument('-O', type=str, help='output directory')
parser.add_argument('-S', type=str, help='sample name')
parser.add_argument('-L', type=str, help='csv file with chromosome names and lengths')

args = parser.parse_args()

samplename = args.S
biomart = args.tss
outdir = args.O
plusreads = args.P
minusreads = args.M
chr_lengths = args.L

chr = {}
tssplus = {}
tssminus = {}
chrnames = []

with open(chr_lengths, "r") as chrlist:
    for line in chrlist:
        ls = line.split(",")
        chr[str(ls[0])[3:]] = int(ls[1])
        tssplus[str(ls[0])[3:]] = []
        tssminus[str(ls[0])[3:]] = []
        chrnames.append(str(ls[0])[3:])

with open(biomart, "r") as mart:
    for line in mart:
        ls = line.split("\t")
        c = str(ls[4])
        if c in chrnames:
            if int(ls[0]) > 9999:
                if str(ls[3]) == "1":
                    if chr[c] > int(ls[0])+10000:
                        tssplus[c].append(
                            [int(ls[0]) - 10000, int(ls[0]) + 10000])
                elif str(ls[3]) == "-1":
                    if chr[c] > int(ls[0])+10000:
                        tssminus[c].append(
                            [int(ls[0]) - 10000, int(ls[0]) + 10000])

readsplus = []
readsminus = []
lastchr = "1"

chrarray = np.zeros(chr[lastchr], np.int8)

with open(plusreads, "r") as plus:
    for line in plus:
        ls = line.split("\t")
        c = str(ls[0])
        mid = int((int(ls[2])-int(ls[1]))/2) + (int(ls[1]))
        if c in chrnames:
            if str(c) != str(lastchr):
                for i in tssplus[lastchr]:
                    readsplus.append(chrarray[i[0]:i[1]])
                for i in tssminus[lastchr]:
                    readsminus.append(chrarray[i[0]:i[1]])
                lastchr = c
                if chr[c] > mid:
                    chrarray = np.zeros(chr[lastchr], np.int8)
                    chrarray[mid] = chrarray[mid] + 1
            else:
                if chr[c] > mid:
                    chrarray[mid] = chrarray[mid] + 1

for i in tssplus[lastchr]:
    readsplus.append(chrarray[i[0]:i[1]])
for i in tssminus[lastchr]:
    readsminus.append(chrarray[i[0]:i[1]])

plusxrplusgenes = [sum(x) for x in zip(*readsplus)]
plusxrminusgenes = [sum(x) for x in zip(*readsminus)]

readsplus2 = []
readsminus2 = []

lastchr = "1"
chrarray = np.zeros(chr[lastchr], np.int8)

with open(minusreads, "r") as plus:
    for line in plus:
        ls = line.split("\t")
        c = str(ls[0])
        mid = int((int(ls[2])-int(ls[1]))/2) + (int(ls[1]))
        if c in chrnames:
            if str(c) != str(lastchr):
                for i in tssplus[lastchr]:
                    readsplus2.append(chrarray[i[0]:i[1]])
                for i in tssminus[lastchr]:
                    readsminus2.append(chrarray[i[0]:i[1]])
                lastchr = c
                if chr[c] > mid:
                    chrarray = np.zeros(chr[lastchr], np.int8)
                    chrarray[mid] = chrarray[mid] + 1
            else:
                if chr[c] > mid:
                    chrarray[mid] = chrarray[mid] + 1
for i in tssplus[lastchr]:
    readsplus2.append(chrarray[i[0]:i[1]])
for i in tssminus[lastchr]:
    readsminus2.append(chrarray[i[0]:i[1]])

minusxrplusgenes = [sum(x) for x in zip(*readsplus2)]
minusxrminusgenes = [sum(x) for x in zip(*readsminus2)]
with open(outdir + "/" + samplename + ".txt", "w") as outf:
    for i in range(20000):
        outf.write(str(plusxrplusgenes[i]) + "\t" + str(plusxrminusgenes[i]) + "\t" + str(
            minusxrplusgenes[i]) + "\t" + str(minusxrminusgenes[i]) + "\n")
