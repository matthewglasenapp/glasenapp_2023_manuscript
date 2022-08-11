# Create a partition file for window based tree building with iqtree. 

number_snv_sites_alignment = 100001

window_size = 50

output_file = "partition_file.txt"

def write_partition_file(output_file, length):
	with open(output_file,"a") as f:
		f.write("#nexus" + "\n")
		f.write("begin sets;" + "\n")

		counter = 0
		intervals = length // (window_size + 1)
		window_counter = 0 
		start_int = 1

		for i in range(1,intervals+1):
			counter += 1
			window_counter += 1
			stop_interval = start_int + window_size
			string = "charset window_{} = {}-{};".format(counter, start_int, stop_interval)
			f.write("\t" + string + "\n") 
			start_int += (window_size + 1)

		f.write("end;" + "\n")

write_partition_file(output_file, number_snv_sites_alignment)

