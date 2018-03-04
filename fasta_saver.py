#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 3 22:32:59 2018

@authors: Elena.Pazhenkova & RagenIV
"""

import sys, re, os

if (len(sys.argv) != 7)and(len(sys.argv) != 8):
	print('Try ./fasta_saver <blastn_res> <chmit.gff3> <ch.fasta> <chmit.fasta> <TYPE> <flag_blast> <min_length>')
	print()
	print('<blastn_res> -- blastn output')
	print('<chmit.gff3> -- mitohondrial gff')
	print('<ch.fasta> -- genome in fasta format')
	print('<chmit.fasta> -- mitohondria in fasta format')
	print('<TYPE> -- name of species')
	print('<flag_blast> -- 0 if genome on mitohodria and else 1')
	print('<min_length> -- minimal acceptable length of results')
	exit()

def peresek(p1s,p1e,p2s,p2e):
	return max(p1s,p2s), min(p1e,p2e)


try:
	tabla = open(sys.argv[1],'r')
	the_tab = tabla.read()
	tabla.close()
except:
	print('Error in <blastn_res>')
	exit()

try:
	chmit_gff3 = open(sys.argv[2],'r')
	gff3 = chmit_gff3.read()
	chmit_gff3.close()
except:
	print('Error in <chmit.gff3>')
	exit()

anno = []

for line in gff3.split('\n'):
	if len(line) == 0:
		continue
	if line[0] == '#':
		continue
	line_con = line.split('\t')
	if line_con[2] != 'gene':
		continue
	# start, stop, name
	naaame = re.findall('gene=(.*);',line_con[8])
	if len(naaame) == 0:
		continue
	anno.append([int(line_con[3]), int(line_con[4]), naaame[0]])

if not((sys.argv[6] == '0')or(sys.argv[6] == '1')):
	print('Error in <flag_blast>')
	exit()

try:
	minlen = int(sys.argv)
except:
	minlen = 50

print('[+] Minimal length: 50')

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
	handle = open(sys.argv[4])
except:
	print('Error in <chmit.fasta>')

try:
	handle2 = open(sys.argv[3])
except:
	print('Error in <ch.fasta>')
for record in SeqIO.parse(handle, "fasta"):
	mitochondr = str(record.seq)

records2 = SeqIO.to_dict(SeqIO.parse(handle2, "fasta"))

if not os.path.exists('RESULTS'):
    os.makedirs('RESULTS')

for line in the_tab.split('\n'):
	if len(line) > 0:
		if line[0] != '#':
			line_con = line.split('\t')

			if sys.argv[6] == '0':
				ident = line_con[0]
				percentage = line_con[2]
				qstart = int(line_con[8])
				qend = int(line_con[9])
				sstart = int(line_con[6])
				send = int(line_con[7])
			else:
				ident = line_con[1]
				percentage = line_con[2]
				qstart = int(line_con[6])
				qend = int(line_con[7])
				sstart = int(line_con[8])
				send = int(line_con[9])

			if percentage == '100.000':
				continue
			if qstart > qend:
				qstart, qend = qend, qstart
			if sstart > send:
				sstart, send = send, sstart

			for gens in anno:
				if (qstart < gens[0] < qend)or(qstart < gens[1] < qend):
					x1, x2 = peresek(gens[0],gens[1],qstart,qend)
					if x2 - x1 < minlen:
						continue
					x1, x2 = x1 -1, x2 -1
					y1 = x1 + send - qend
					y2 = x2 + send - qend

					seq2 = str(records2[ident].seq[y1:y2])
					nuc = SeqRecord(Seq(seq2), id = sys.argv[5]+'_nuclear')


					seq1 = mitochondr[x1:x2]
					mit = SeqRecord(Seq(seq1), id = sys.argv[5]+'_mt')
					fin = [nuc,mit]
					output_handle = open("./RESULTS/"+gens[2]+".fasta", "a")
					for seq in fin:
						SeqIO.write(seq, output_handle, "fasta"), 
					output_handle.close()

print('[!] Check your results in directory RESULTS')
