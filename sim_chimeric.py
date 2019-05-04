#!/usr/bin/env python
import sys
import os
import argparse
import random
import re
import multiprocessing


def GetOpts():
	group = argparse.ArgumentParser()
	group.add_argument('-a', '--a_fasta', help='fasta file', required=True)
	group.add_argument('-b', '--b_fasta', help='fasta file with same chromosome order with fasta file a', required=True)
	group.add_argument('-o', '--output', help='filename of simulated data', required=True)
	group.add_argument('-c', '--chimeric', type=float, help='percentage of chimeric region size, like 5 means 5%%, default: 10', default=10)
	group.add_argument('-t', '--thread', type=int, help='thread number, default: 1', default=1)
	return group.parse_args()


def ReadFasta(inFasta):
	fastaDB = {}
	id_list = []
	with open(inFasta, 'r') as fIn:
		id = ''
		seq = ''
		totalLen = 0
		for line in fIn:
			if line[0] == '>':
				if seq != '':
					fastaDB[id] = seq
				data = line.strip()[1:].split()
				id = data[0]
				id_list.append(id)
				seq = ''
			else:
				seq += line.strip()
		fastaDB[id] = seq
	return fastaDB, id_list


def isOvlp(posList, pos):
	newPosList = sorted(posList)
	s = 0
	e = len(newPosList)-1
	if e < 0:
		return False
	while s<=e:
		m = int((s+e)/2)
		if posList[m][0] > pos[0]:
			e = m-1
		elif posList[m][0] < pos[0]:
			s = m+1
		else:
			return True
	if min(posList[e][1], pos[1])-max(posList[e][0], pos[0]) >= 0:
		return True
	else:
		return False


def subProcess(idListA, fastaA, idListB, fastaB, chimeric):
	for i in range(0, len(idListA)):
		idA = idListA[i]
		idB = idListB[i]
		totalLen = len(fastaA[idA])+len(fastaB[idB])
		totalChiLen = int(totalLen*chimeric)
		posList = []
		while totalChiLen > 0:
			curChiLen = random.randint(5000, 1e6)
			if curChiLen > totalChiLen:
				curChiLen = totalChiLen
			curChiPos = random.randint(0, totalLen-curChiLen-1)
			while isOvlp(posList, [curChiPos, curChiPos+curChiLen]) or curChiPos < len(fastaA[idA])-1 and curChiPos+curChiLen >= len(fastaA[idA]):
				curChiPos = random.randint(0, totalLen-curChiLen-1)
			posList.append([curChiPos, curChiPos+curChiLen])
			totalChiLen -= curChiLen
		posList = sorted(posList)
		seqA = list(fastaA[idA])
		seqB = list(fastaB[idB])
		lenA = len(fastaA[idA])
		lenB = len(fastaB[idB])
		for pos in posList:
			chiLen = pos[1]-pos[0]+1
			if pos[1] < lenA:
				chis = random.randint(0, lenB-chiLen-1)
				chie = chis+chiLen
				seqA[pos[0]: pos[1]] = fastaB[idB][chis: chie]
			else:
				chis = random.randint(0, lenA-chiLen-1)
				chie = chis+chiLen
				seqB[pos[0]-lenA: pos[1]-lenA] = fastaA[idA][chis: chie]
		with open(idA+"-"+idB+".tmp", 'w') as fout:
			fout.write(">%s\n%s\n>%s\n%s\n"%(idA, ''.join(seqA), idB, ''.join(seqB)))



def SimChimeric(aFasta, bFasta, outFa, chimeric, tn):
	random.seed()
	print("Reading fasta file: %s"%(aFasta))
	fastaA, idA = ReadFasta(aFasta)
	print("Reading fasta file: %s"%(bFasta))
	fastaB, idB = ReadFasta(bFasta)
	
	chrCnt = min(len(idA), len(idB))
	chrPerT = int(chrCnt/tn)

	ts = []
	for i in range(0, tn):
		if i < tn-1:
			t = multiprocessing.Process(target=subProcess, args=(idA[i*chrPerT: (i+1)*chrPerT], fastaA, idB[i*chrPerT: (i+1)*chrPerT], fastaB, chimeric))
		else:
			t = multiprocessing.Process(target=subProcess, args=(idA[i*chrPerT:], fastaA, idB[i*chrPerT:], fastaB, chimeric))
		ts.append(t)
	print("Simulating chimeric")
	for t in ts:
		t.start()
	for t in ts:
		t.join()
	print("Merging result")
	with open(outFa, 'w') as fout:
		for i in range(0, chrCnt):
			with open(idA[i]+'-'+idB[i]+'.tmp', 'r') as fin:
				for line in fin:
					fout.write(line)
			os.remove(idA[i]+'-'+idB[i]+'.tmp')
	print("Finished")


if __name__ == "__main__":
	opts = GetOpts()
	aFasta = opts.a_fasta
	bFasta = opts.b_fasta
	outFa = opts.output
	chimeric = opts.chimeric/100.0
	tn = opts.thread
	SimChimeric(aFasta, bFasta, outFa, chimeric, tn)
