output_file = "/Users/matt/desktop/test.txt"

def write_partition_file(output_file, length):
	with open(output_file,"a") as f:
		f.write("#nexus" + "\n")
		f.write("begin sets;" + "\n")

		counter = 0
		intervals = length // 51
		window_counter = 0 
		start_int = 1

		for i in range(1,intervals+1):
			counter += 1
			window_counter += 1
			stop_interval = start_int + 50 
			string = "charset window_{} = {}-{};".format(counter, start_int, stop_interval, stop_interval-start_int)
			f.write("\t" + string + "\n") 
			start_int += 51

		f.write("end;" + "\n")

write_partition_file(output_file, 10001)

