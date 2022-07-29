# My main program does:

1) Get a list of lists of [output_file,coordinate_file] using get_file_paths_pairs_list()

2) Run process_single_scaffold on for each list of output file, coordinate_file

Process_single_scaffold appends results to a bunch of lists defined at the top of the script

# Chunk this out further 
process_single_scaffold calls:

	find_tracts() - returns index_tract_list and coordinate_tract_list - lists of lists of tracts in format of  [[start_index,stop_index,length in base pairs], []] and [[start_coordinate,stop_coordinate,length], []]

	sort_tract_list() - Sorts a tract list in descending order of length in base pairs

	get_tract_length_dist() - get list of all tract lengths for a tract list
