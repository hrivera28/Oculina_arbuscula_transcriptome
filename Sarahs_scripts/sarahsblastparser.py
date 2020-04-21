#! /usr/bin/env python
#usage: sarahsblastparser.py lengththreshold percentIDthreshold blastfile.txt
import sys

lenthresh=int(sys.argv[1]) # the length of the blast overlap threshold
perIDthresh=float(sys.argv[2]) # the percent match threshold
infile = open(sys.argv[3], 'r')
outfile = open('%s_goodmatches.txt' %(sys.argv[3][:-4]), 'w')	

gooddict={}
evaluelist=[]
matches=[]
linenum=0
for line in infile:
	linenum+=1
	cols=line.rstrip().split('\t')
	if int(cols[6])>= lenthresh and float(cols[7])/int(cols[6])*100 >= perIDthresh:
		try:
			if float(gooddict[cols[0]][3]) > float(cols[3]):
				gooddict[cols[0]]=cols
				evaluelist.append(float(cols[3]))
				matches.append(cols[1])
			else:
				continue
		except KeyError:
			gooddict[cols[0]]=cols
			evaluelist.append(float(cols[3]))
			matches.append(cols[1])
for key in sorted(gooddict.keys()):
	outfile.write('%s\n' %('\t'.join(gooddict[key])))

print('Number of blast matches parsed: %d\nNumber of Good matches: %d\nHighest EValue of a Match: %.2e' %(linenum, len(gooddict.keys()), max(evaluelist)))