"""
Check that sequencing reads in fastq files come from a single flowcell lane. 
"""

from joblib import Parallel, delayed

fastq_dir = "/hb/scratch/mglasena/short_read_data/"

# ERR5621404_1 and SRR6281818 were skipped because they do not inlcude any metadata, only read number and length. 
skipped_file_list = ["ERR5621404_1.fastq", "ERR5621404_1_2.fastq", "SRR6281818_1.fastq", "SRR6281818_2.fastq"]

#file_list = ["DRR107784_1.fastq", "DRR107784_2.fastq", "SRR2583947.fastq", "SRR5767279_1.fastq", "SRR5767279_2.fastq", "SRR5767280_1.fastq", "SRR5767280_2.fastq", "SRR5767281_1.fastq", "SRR5767281_2.fastq", "SRR5767282_1.fastq", "SRR5767282_2.fastq", "SRR5767283_1.fastq", "SRR5767283_2.fastq", "SRR5767284_1.fastq", "SRR5767284_2.fastq", "SRR5767285_1.fastq", "SRR5767285_2.fastq", "SRR5767286_1.fastq", "SRR5767286_2.fastq", "SRR7211988_1.fastq", "SRR7211988_2.fastq"]
file_list = ["DRR107784_1.fastq", "DRR107784_2.fastq"]

def check_fastq_file(file):
	file = fastq_dir + file
	print("Checking " + file.split("/")[7])
	split_char = ":"
	
	# Extract read group string (<instrument>:<run number>:<flowcell ID>:<lane>)
	# Normal case. Three fields separated by space (@SRR5767279.1 HS3:147:d0gnlacxx:3:1101:1371:2064 length=100)
	if len(open(file,"r").readline().split()) == 3:
		
		# Normal case (@SRR5767279.1 HS3:147:d0gnlacxx:3:1101:1371:2064 length=100)
		# Extract HS3:147:d0gnlacxx:3 as read group string
		tmp = open(file,"r").readline().split()[1].split(split_char)
		if len(tmp) == 7:
			K = 4
			read_group_string = split_char.join(tmp[:K])
			print("The read group string is " + read_group_string)

		# Edge case for pallidus (SRR5767285): @SRR5767285.1 HS1_0066:8:1101:1237:2149 length=100
		# Extract HS1_0066:8 as read group string 
		elif len(tmp) == 5:
			K = 2
			read_group_string = split_char.join(tmp[:K])
			print("The read group string is " + read_group_string)
	
	# Edge case for pulcherrimus (SRR5767283) and intermedius (SRR5767280) fastq files.
	# There are only two fields in header line: @HS2:148:C0EN2ACXX:3:1101:1327:2227 2:Y:0:
	elif len(open(file,"r").readline().split()) == 2:
		K = 4
		tmp = open(file,"r").readline().split()[0][1:].split(split_char)
		read_group_string = split_char.join(tmp[:K])
		print("The read group string is " + read_group_string)

	else:
		print("Something went wrong!")

	# Check that the header for each read has the same read group string as the first header in the file
	for line in open(file,"r"):
		if line[0] == "@":
			
			# Normal case of three fields (@SRR5767285.1 HS1_0066:8:1101:1237:2149 length=100)
			if len(line.split()) == 3:
				if line.split()[1][0:len(read_group_string)] != read_group_string:
					print ("error!")
					print(line)
					break

			# Edge case for pulcherrimus (SRR5767283) and intermedius (SRR5767280): @HS2:148:C0EN2ACXX:3:1101:1327:2227 2:Y:0:
			elif len(line.split()) == 2:
				if line.split()[0][1:len(read_group_string)+1] != read_group_string:
					print("error!")
					print(line)
					break

	print(file.split("/")[7] + " done!")

def main():
	Parallel(n_jobs=len(file_list))(delayed(check_fastq_file)(file) for file in file_list)

if __name__ == "__main__":
	main()