#!/usr/bin/env python
import sys
import argparse
import random
import time


def GetOpts():
	group = argparse.ArgumentParser()
	group.add_argument('-a', '--a_contigs', help='first contig fasta file', required=True)
	group.add_argument('-b', '--b_contigs', help='second contig fasta file', required=True)
	group.add_argument('-o', '--output', help='filename of simulated data', required=True)
	group.add_argument('-s', '--blast', help='blast file with format 6, must use first file of input as query and second file as database', required=True)
	group.add_argument('-p', '--prefix', help='prefix of contig file a and contig file b, divided by comma, like: HA, HB, required when contig file contain same IDs', default="")
	group.add_argument('-c', '--collapse', type=float, help='persentage of collapse region size, like 5 means 5%%, default: 10', default=10)

	return group.parse_args()


def ReadFasta(inFasta):
	fastaDB = {}
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
				seq = ''
			else:
				seq += line.strip()
				totalLen += len(line.strip())
		fastaDB[id] = seq
	return fastaDB, totalLen


def ReadBlast(inBlast, prefixList):
	blastDB = {}
	with open(inBlast, 'r') as fBlast:
		for line in fBlast:
			data = line.strip().split()
			queryID = data[0]
			targetID = data[1]
			identity = float(data[2])
			queryRegion = list(map(int, [data[6], data[7]]))
			targetRegion = list(map(int, [data[8], data[9]]))
			if queryID not in blastDB:
				blastDB[queryID] = [targetID, identity, queryRegion, targetRegion]
			else:
				if identity > blastDB[queryID][1]:
					blastDB[queryID] = [targetID, identity, queryRegion, targetRegion]
	
	mapping = {}
	allContigList = []
	for queryID in blastDB:
		targetID = blastDB[queryID][0]
		mapping[prefixList[0]+'-'+queryID] = prefixList[1]+"-"+targetID
		mapping[prefixList[1]+"-"+targetID] = prefixList[0]+'-'+queryID
		allContigList.append(prefixList[0]+'-'+queryID)
		allContigList.append(prefixList[1]+"-"+targetID)
	
	return mapping, allContigList


def SimCollapse(aFasta, bFasta, outFa, blastFile, prefixList, collapse):
	print("\033[32m%s\033[0m Reading first contig file"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	fastaDBA, lenA = ReadFasta(aFasta)
	
	print("\033[32m%s\033[0m Reading species list"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	print("           Reading second contig file")
	fastaDBB, lenB = ReadFasta(bFasta)
	
	lenDB = {}
	for id in fastaDBA:
		lenDB[prefixList[0]+"-"+id] = len(fastaDBA[id])
	for id in fastaDBB:
		lenDB[prefixList[1]+"-"+id] = len(fastaDBB[id])
	
	print("\033[32m%s\033[0m Reading blast file"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	mapping, allContigList = ReadBlast(blastFile, prefixList)
	collapseLen = int((lenA+lenB)*collapse)
	
	print("\033[32m%s\033[0m Removing collapses"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	print("           Total collapse size expected: %d"%(collapseLen))
	removeList = {}
	removeLen = 0
	while collapseLen > 0:
		index = random.randint(0, len(allContigList)-1)
		name = allContigList[index]

		while mapping[name] not in allContigList:
			index = random.randint(0, len(allContigList)-1)
			name = allContigList[index]

		repeatCnt = 0
		isLast = False
		while mapping[name] not in allContigList or lenDB[name] > collapseLen:
			if repeatCnt > 50:
				isLast = True
				break
			if lenDB[name] > collapseLen:
				repeatCnt += 1
			index = random.randint(0, len(allContigList)-1)
			name = allContigList[index]
		
		if isLast:
			break
		pre, ctg = name.split('-')
		collapseLen -= lenDB[name]
		removeLen += lenDB[name]
		
		if pre not in removeList:
			removeList[pre] = []
		removeList[pre].append(ctg)
		
		allContigList.remove(name)
		allContigList.remove(mapping[name])
	
	print("           Total collapse size removed: %d"%(removeLen))
	
	print("\033[32m%s\033[0m Writing result"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	with open(outFa, 'w') as fOut:
		for pre in removeList:
			if pre == prefixList[0]:
				for id in fastaDBA:
					if id not in removeList[pre]:
						fOut.write(">%s_%s\n%s\n"%(pre, id, fastaDBA[id]))
			else:
				for id in fastaDBB:
					if id not in removeList[pre]:
						fOut.write(">%s_%s\n%s\n"%(pre, id, fastaDBB[id]))
	
	print("\033[32m%s\033[0m Finished"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))


if __name__ == "__main__":
	opts = GetOpts()
	aFasta = opts.a_contigs
	bFasta = opts.b_contigs
	outFa = opts.output
	collapse = opts.collapse/100.0
	blastFile = opts.blast
	prefixList = opts.prefix.split(',')
	SimCollapse(aFasta, bFasta, outFa, blastFile, prefixList, collapse)
