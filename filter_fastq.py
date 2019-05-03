#!/usr/bin/env python
import sys
import gc


def is_in_list(filter_list, id):
	s = 0
	e = len(filter_list)-1
	while s<=e:
		mid = (s+e)/2
		if filter_list[mid] > id:
			e = mid-1
		elif filter_list[mid] < id:
			s = mid+1
		else:
			return True
	return False


def filter_fastq(in_fq, in_list, out_fq):
	filter_list = []
	print("Reading filter list")
	with open(in_list, 'r') as fin:
		for line in fin:
			filter_list.append(line.strip())
	
	sorted_filter_list = sorted(filter_list)
	fq_db = {}
	print("Reading fastq")
	with open(in_fq, 'r') as fin:
		for line in fin:
			if line[0] == '@' and len(fq_db[id])==4:
				id = line[1:].strip().split()[0]
				fq_db[id] = []
			fq_db[id].append(line)
	
	print("Writing fastq")
	with open(out_fq, 'w') as fout:
		for id in sorted(fq_db):
			if is_in_list(sorted_filter_list, id) == False:
				fout.write(''.join(fq_db[id]))
	del fq_db
	gc.collect()
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_fastq> <in_list> <out_fastq>")
	else:
		in_fq, in_list, out_fq = sys.argv[1:]
		filter_fastq(in_fq, in_list, out_fq)
