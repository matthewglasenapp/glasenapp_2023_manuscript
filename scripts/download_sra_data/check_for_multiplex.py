"""
Check that sequencing reads in fastq files come from a single flowcell lane. 
"""

from joblib import Parallel, delayed
import gzip

# ERR5621404 and SRR6281818 were skipped because they do not inlcude any metadata, only read number and length. 
file_list = [
"/hb/home/mglasena/short_read_data/DRR107784_1.fastq.gz", 
"/hb/home/mglasena/short_read_data/DRR107784_2.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767279_1.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767279_2.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767280_1.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767280_2.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767281_1.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767281_2.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767282_1.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767282_2.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767283_1.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767283_2.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767284_1.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767284_2.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767285_1.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767285_2.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767286_1.fastq.gz", 
#"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR5767286_2.fastq.gz", 
"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR7211988_H2YWCBCXX:1_1.fastq.gz", 
"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR7211988_H2YWCBCXX:1_2.fastq.gz",
"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR7211988_H2YWCBCXX:2_1.fastq.gz",
"/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR7211988_H2YWCBCXX:2_2.fastq.gz"
]

def check_fastq_file(file):
	print("Checking " + file)
	
	read_group_string_list = []

	split_char = ":"
	
	# Extract read group string (<instrument>:<run number>:<flowcell ID>:<lane>)
	# the normal case is three fields separated by the space character (@SRR5767279.1 HS3:147:d0gnlacxx:3:1101:1371:2064 length=100)
	
	if len(gzip.open(file,"rt").readline().split()) == 3:
		
		# Normal case (@SRR5767279.1 HS3:147:d0gnlacxx:3:1101:1371:2064 length=100)
		# Extract HS3:147:d0gnlacxx:3 as read group string
		tmp = gzip.open(file,"rt").readline().split()[1].split(split_char)
		if len(tmp) == 7:
			K = 4
			read_group_string = split_char.join(tmp[:K])
			read_group_string_list.append(read_group_string)
			print("The read group string is " + read_group_string)

		# Edge case for pallidus (SRR5767285): @SRR5767285.1 HS1_0066:8:1101:1237:2149 length=100
		# Extract HS1_0066:8 as read group string 
		elif len(tmp) == 5:
			K = 2
			read_group_string = split_char.join(tmp[:K])
			read_group_string_list.append(read_group_string)
			print("The read group string is " + read_group_string)
	
	# Edge case for pulcherrimus (SRR5767283) and intermedius (SRR5767280) fastq files.
	# There are only two fields in header line: @HS2:148:C0EN2ACXX:3:1101:1327:2227 2:Y:0:
	elif len(gzip.open(file,"r").readline().split()) == 2:
		K = 4
		tmp = gzip.open(file,"r").readline().split()[0][1:].split(split_char)
		read_group_string = split_char.join(tmp[:K])
		read_group_string_list.append(read_group_string)
		print("The read group string is " + read_group_string)

	else:
		print("Something went wrong!")

	# Check that the header for each read has the same read group string as the first header in the file
	for line in gzip.open(file,"rt"):
		if line[0] == "@":
			
			# Normal case of three fields (@SRR5767285.1 HS1_0066:8:1101:1237:2149 length=100)
			if len(line.split()) == 3:
				if line.split()[1][0:len(read_group_string)] not in read_group_string_list:
					new_read_group_string = line.split()[1][0:len(read_group_string)]
					print("Alternate read group string found: {}".format(new_read_group_string))
					read_group_string_list.append(new_read_group_string)

			# Edge case for pulcherrimus (SRR5767283) and intermedius (SRR5767280): @HS2:148:C0EN2ACXX:3:1101:1327:2227 2:Y:0:
			elif len(line.split()) == 2:
				if line.split()[0][1:len(read_group_string)+1] not in read_group_string_list:
					new_read_group_string = line.split()[0][1:len(read_group_string)+1]
					print("Alternate read group string found: {}".format(new_read_group_string))
					read_group_string_list.append(new_read_group_string)

	print(file + " done!")
	return read_group_string_list

def main():
	Parallel(n_jobs=len(file_list))(delayed(check_fastq_file)(file) for file in file_list)

if __name__ == "__main__":
	main()