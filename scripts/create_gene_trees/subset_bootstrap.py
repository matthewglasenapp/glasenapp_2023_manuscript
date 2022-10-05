boot_file = "test.txt"
output_file = "output"

# Create test file
#with open("test.txt","a") as f:
	#for n in range(0,999999):
		#f.write(str(n) + "\n")

def subset_boot_file():
	line_counter = 0 

	with open(output_file,"a") as f:
		for line in open(boot_file,"r"):
			if ((line_counter + 1000) % 1000) == 0:
				count = 0
			if count <= 99:
				f.write(line)
				count += 1 
			line_counter += 1

def main():
	subset_boot_file()

if __name__ == "__main__":
	main()





