#from joblib import Parallel, delayed

#SRR7211988

#@SRR7211988.1 HISEQ:348:H2YWCBCXX:1:1101:1057:2031 length=150

#@SRR7211988.155831747 HISEQ:348:H2YWCBCXX:2:1101:1212:2021 length=150

#forward_file = "/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR7211988_1.fastq"
#reverse_file = "/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR7211988_2.fastq"

#read_group_strings = ['HISEQ:348:H2YWCBCXX:1', 'HISEQ:348:H2YWCBCXX:2']

working_dir = "/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/"

output_dir = "/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/test/"

file_list = ["SRR7211988_1.fastq","SRR7211988_2.fastq"]

lane_1_forward = open(output_dir + "SRR7211988_lane1_1.fastq","a")
lane_1_reverse = open(output_dir + "SRR7211988_lane1_2.fastq","a")
lane_2_forward = open(output_dir + "SRR7211988_lane2_1.fastq","a")
lane_2_reverse = open(output_dir + "SRR7211988_lane2_2.fastq","a")

def find_all_read_group_strings(file):
	print("Sorting " + file.split("/")[7])

	read_group_string_set = set()

	split_char = ":"
	K = 4
	
	tmp = open(file,"r").readline().split()[1].split(split_char)
	read_group_string = split_char.join(tmp[:K])

	read_group_string_set.add(read_group_string)

	for line in open(file,"r"):
		if line[0:4] == "@SRR":
			for item in read_group_string_set:
				if item in line:
					continue
			else:
				tmp = line.split()[1].split(split_char)
				new_read_group_string = split_char.join(tmp[:K])
				read_group_string_set.add(new_read_group_string)

	print(read_group_string_set)

def write_new_files(file):
	
	print("Demultiplexing " + file)

	file = working_dir + file
	
	lane_1_read_group_string =  "HISEQ:348:H2YWCBCXX:1"
	lane_2_read_group_string = "HISEQ:348:H2YWCBCXX:2"
	
	if file == "/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR7211988_1.fastq":
		with open(file,"r") as f:
			for line in f:
				if line[0:4] == "@SRR":
					if lane_1_read_group_string in line:
						lane_1 = "TRUE"
					elif lane_2_read_group_string in line:
						lane_1 = "FALSE"

				if lane_1 == "TRUE":
					lane_1_forward.write(line)
				else:
					lane_2_forward.write(line)

	elif file == "/hb/groups/pogson_group/dissertation/data/raw_sequencing_reads/SRR7211988_2.fastq":
		with open(file,"r") as f:
			for line in f:
				if line[0:4] == "@SRR":
					if lane_1_read_group_string in line:
						lane_1 = "TRUE"
					elif lane_2_read_group_string in line:
						lane_1 = "FALSE"

				if lane_1 == "TRUE":
					lane_1_reverse.write(line)
				else:
					lane_2_reverse.write(line)
				

def main():
	#find_all_read_group_strings(forward_file)
	for file in file_list:
		write_new_files(file)
	#Parallel(n_jobs=2)(delayed(write_new_files)(file) for file in file_list)

if __name__ == "__main__":
	main()
