# My main program does:

1) Get a list of lists of [output_file,coordinate_file] using get_file_paths_pairs_list()

2) Run process_single_scaffold on for each list of output file, coordinate_file

Process_single_scaffold appends results to a bunch of lists defined at the top of the script

# Chunk this out further 
process_single_scaffold calls:

	find_tracts() - returns index_tract_list - list of lists of tracts in format of  [[start_index,stop_index,length in base pairs], []]

	get_coordinate_tract_list_and_tract_lengths() - returns coordinate_tract_list - list of lists of tracts in format of [[start_coordinate,stop_coordinate,length], []]

	sort_index_tract_list() - Sorts index_tract_list in descending order of length in base pairs

	sorted_coordinate_tract_list() - sorts coordinate_tract_list in descending order of length in base pairs 

	get_tract_length_dist() - get list of all tract lengths for




Notes: I think process single scaffold should only call on helper functions. Chunk stuff out into helper functions. 