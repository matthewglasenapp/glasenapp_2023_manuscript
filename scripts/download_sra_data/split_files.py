from joblib import Parallel, delayed
import gzip

short_read_directory = "/hb/home/mglasena/short_read_data/"

forward_file = "/hb/home/mglasena/short_read_data/DRR107784_1.fastq.gz"

working_dir = "/hb/home/mglasena/short_read_data/"

output_dir = "/hb/scratch/mglasena/"

file_list = ["DRR107784_1.fastq.gz","DRR107784_1.fastq.gz"]

read_group_1_forward = open(output_dir + "DRR107784_C1J2AACXX:5_1.fastq", "a")
read_group_1_reverse = open(output_dir + "DRR107784_C1J2AACXX:5_2.fastq", "a")
read_group_2_forward = open(output_dir + "DRR107784_C1JJYACXX:4_1.fastq", "a")
read_group_2_reverse = open(output_dir + "DRR107784_C1JJYACXX:4_2.fastq", "a")
read_group_3_forward = open(output_dir + "DRR107784_C1J2AACXX:4_1.fastq", "a")
read_group_3_reverse = open(output_dir + "DRR107784_C1J2AACXX:4_2.fastq", "a")
read_group_4_forward = open(output_dir + "DRR107784_C1MA0ACXX:8_1.fastq", "a")
read_group_4_reverse = open(output_dir + "DRR107784_C1MA0ACXX:8_2.fastq", "a")

def find_all_read_group_strings(file):
	print("Sorting " + file)

	read_group_string_set = set()

	split_char = ":"
	K = 4
	
	tmp = gzip.open(file,"rt").readline().split()[1].split(split_char)
	read_group_string = split_char.join(tmp[:K])

	read_group_string_set.add(read_group_string)

	for line in gzip.open(file,"rt"):
		if line[0:4] == "@SRR" or line[0:4] == "@DRR":
			for item in read_group_string_set:
				if item in line:
					continue
			else:
				tmp = line.split()[1].split(split_char)
				new_read_group_string = split_char.join(tmp[:K])
				read_group_string_set.add(new_read_group_string)

	print(read_group_string_set)
	return read_group_string_set

def write_new_files(file):
	
	print("Demultiplexing " + file)

	file = working_dir + file

	read_group_string_1 = 'HWI-ST462R:262:C1J2AACXX:5'
	read_group_string_2 = 'HWI-1KL134:197:C1JJYACXX:4'
	read_group_string_3 = 'HWI-ST462R:262:C1J2AACXX:4'
	read_group_string_4 = 'HWI-1KL134:198:C1MA0ACXX:8'

	if file == "/hb/home/mglasena/short_read_data/DRR107784_1.fastq.gz":
		with gzip.open(file,"rt") as f:
			for line in f:
				if line[0:4] == "@DRR":
					if read_group_string_1 in line:
						read_group_1 = "True"
						read_group_2 = "False"
						read_group_3 = "False"
						read_group_4 = "False"
					elif read_group_string_2 in line:
						read_group_2 = "True"
						read_group_1 = "False"
						read_group_3 = "False"
						read_group_4 = "False"
					elif read_group_string_3 in line:
						read_group_3 = "True"
						read_group_1 = "False"
						read_group_2 = "False"
						read_group_4 = "False"
					elif read_group_string_4 in line:
						read_group_4 = "True"
						read_group_1 = "False"
						read_group_2 = "False"
						read_group_3 = "False"

				if read_group_1 == "True":
					read_group_1_forward.write(line)
				elif read_group_2 == "True":
					read_group_2_forward.write(line)
				elif read_group_3 == "True":
					read_group_3_forward.write(line)
				elif read_group_4 == "True":
					read_group_4_forward.write(line)

	elif file == "/hb/home/mglasena/short_read_data/DRR107784_2.fastq.gz":
		with gzip.open(file,"rt") as f:
			for line in f:
				if line[0:4] == "@SRR":
					if read_group_string_1 in line:
						read_group_1 = "True"
						read_group_2 = "False"
						read_group_3 = "False"
						read_group_4 = "False"
					elif read_group_string_2 in line:
						read_group_2 = "True"
						read_group_1 = "False"
						read_group_3 = "False"
						read_group_4 = "False"
					elif read_group_string_3 in line:
						read_group_3 = "True"
						read_group_1 = "False"
						read_group_2 = "False"
						read_group_4 = "False"
					elif read_group_string_4 in line:
						read_group_4 = "True"
						read_group_1 = "False"
						read_group_2 = "False"
						read_group_3 = "False"

				if read_group_1 == "True":
					read_group_1_reverse.write(line)
				elif read_group_2 == "True":
					read_group_2_reverse.write(line)
				elif read_group_3 == "True":
					read_group_3_reverse.write(line)
				elif read_group_4 == "True":
					read_group_4_reverse.write(line)
				

def main():
	#Parallel(n_jobs=2)(delayed(write_new_files)(file) for file in file_list)
	write_new_files("DRR107784_1.fastq.gz")

if __name__ == "__main__":
	main()
