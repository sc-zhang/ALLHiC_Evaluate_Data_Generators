#!/usr/bin/env python
import sys


def split_cross_fastq(in_fq, out_pre):
	f_in = open(in_fq, 'r')
	f_R1 = open(out_pre+"_R1.fastq", 'w')
	f_R2 = open(out_pre+"_R2.fastq", 'w')
	
	tmp_list = []
	for line in f_in:
		tmp_list.append(line)
		if len(tmp_list) == 8:
			f_R1.write(''.join(tmp_list[:4]))
			f_R2.write(''.join(tmp_list[4: 8]))
			tmp_list = []
	
	f_in.close()
	f_R1.close()
	f_R2.close()


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_fq> <out_pre>")
	else:
		in_fq, out_pre = sys.argv[1:]
		split_cross_fastq(in_fq, out_pre)
