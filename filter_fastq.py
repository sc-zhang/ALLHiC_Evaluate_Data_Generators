#!/usr/bin/env python
import sys
import gc


def filter_fastq(in_fq, in_list, out_fq):
	filter_set = set()
	print("Reading filter list")
	with open(in_list, 'r') as fin:
		for line in fin:
			filter_set.add(line.strip())
	
	fq_db = {}
	print("Reading fastq")
	with open(in_fq, 'r') as fin:
        id = ""
        cnt = 0
		for line in fin:
			if line[0] == '@' and cnt%4==0:
				id = line[1:].strip().split()[0]
				fq_db[id] = []
			fq_db[id].append(line)
            cnt += 1
	
	print("Writing fastq")
	with open(out_fq, 'w') as fout:
		for id in sorted(fq_db):
			if not id in filter_set:
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
