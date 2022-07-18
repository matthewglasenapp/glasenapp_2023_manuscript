"""Pre-process raw fastq files. Steps inlcluded are:
1. Picard FastqToSam
2. Picard MarkIlluminaAdapters
3. Picard SamToFastq | bwa mem | gatk MergeBamAlignment
4. gatk MarkDuplicatesSpark
5. gatk HaplotypeCaller
6. bcftools norm
7. gatk IndexFeatureFile

fasterq-dump is not included. This script assumes there is a directiory with the raw fastq files for samples. Intermediate files are deleted. The final result is a vcf file.
"""

import os 
from joblib import Parallel, delayed

# Root directory for output files 
root_dir = "/hb/groups/pogson_group/dissertation/data/"

# Path to S. purpuratus reference genome file
reference_genome = "/hb/groups/pogson_group/dissertation/data/purpuratus_reference/GCF_000002235.5_Spur_5.0_genomic.fna"

# Temporary directory for intermediate files
tmp_dir = "/hb/scratch/mglasena/"

# Directory containing raw fastq read files
raw_fastq_dir = root_dir + "do_not_delete/raw_sequencing_reads/"

# Directory for unmapped bam files
ubam_dir = root_dir + "unmapped_bam_files/"
make_ubam_dir = "mkdir -p {}".format(ubam_dir)
os.system(make_ubam_dir)

# Directory for unampped bam files with adapters marked
uBAM_XT_dir = root_dir + "uBAM_XT/"
make_uBAM_XT_dir = "mkdir -p {}".format(uBAM_XT_dir)
os.system(make_uBAM_XT_dir)

# Directory for mapped bam files
mapped_bam_dir = root_dir + "mapped_bam_files/"
make_mapped_bam_dir = "mkdir -p {}".format(mapped_bam_dir)
os.system(make_mapped_bam_dir)

# Directory for mapped bam files with duplicates marked
dedup_bam_dir = root_dir + "dedup_mapped_bam_files/"
make_dedup_bam_dir = "mkdir -p {}".format(dedup_bam_dir)
os.system(make_dedup_bam_dir)

# Directory for vcf files
vcf_dir = root_dir + "vcf_files/"
make_vcf_dir = "mkdir -p {}".format(vcf_dir)
os.system(make_vcf_dir)

dict = {
"SRR5767279" : ["SRR5767279","fragilis","QB3KMK013","SAMN07269103","VJCQB3KMK013","HS3:147:d0gnlacxx:3","d0gnlacxx:3"],
"SRR5767281" : ["SRR5767281","nudus","QB3KMK011","SAMN07269101","VJCQB3KMK011","HS2:148:C0EN2ACXX:4","C0EN2ACXX:4"],
"SRR5767282" : ["SRR5767282","franciscanus","QB3KMK010","SAMN07269100","VJCQB3KMK010","HS2:148:C0EN2ACXX:5","C0EN2ACXX:5"],
#"DRR107786" : ["DRR107786","pulcherrimus","SAMD00098133","SAMD00098133","DRR107786","HWI-ST462R:262:C1J2AACXX:3","C1J2AACXX:3"],
#"SRR5767284" : ["SRR5767284","depressus","QB3KMK015","SAMN07269098","VJCQB3KMK015","HS3:171:d0le4acxx:2","d0le4acxx:2"],
#"SRR5767285" : ["SRR5767285","pallidus","QB3KMK002","SAMN07269097","VJCQB3KMK002","HS1_0066:8","HS1_0066:8"],
#"SRR5767286" : ["SRR5767286","droebachiensis","QB3KMK014","SAMN07269096","VJCQB3KMK014","HS3:171:d0le4acxx:1","d0le4acxx:1"],
"SRR6281818" : ["SRR6281818","purpuratus","S.purpuratus#1","SAMN08013506","A630_1","SRR6281818","SRR6281818"],
#"ERR5671699" : ["ERR5671699","lividus","4","ERS2351987","4089_HiC_PI_3","ERR5671699","ERR5671699"],
"SRR5767283" : ["SRR5767283","pulcherrimus","QB3KMK016","SAMN07269099","VJCQB3KMK016","HS3:171:d0le4acxx:3","d0le4acxx:3"],
"SRR5767280" : ["SRR5767280","intermedius","QB3KMK012","SAMN07269102","VJCQB3KMK012","HS2:148:C0EN2ACXX:3","C0EN2ACXX:3"],
"SRR7211988" : ["SRR7211988","purpuratus","SPUR.00","SAMN00829422","CIT_GEC_SP_1",["HISEQ:348:H2YWCBCXX:1","HISEQ:348:H2YWCBCXX:2"],["H2YWCBCXX:1","H2YWCBCXX:2"]]
}

class Accessions:

	def __init__(self, accession, species, sample_name, ncbi_biosample, library_name, read_group_string, platform_unit):
		self.accession = accession
		self.species = species
		self.sample = sample_name
		self.biosample = ncbi_biosample
		self.library = library_name
		self.read_group_string = read_group_string
		self.platform_unit = platform_unit

		print("Processing {}: {}!".format(self.accession, self.species))

	def __repr__(self):
		return self.acceession

	def convert_fastq_to_unmapped_bam(self):
		print("Picard FastqToSam. Converting fastq files to unampped BAM format.")

		# Multiplex
		if len(self.read_group_string) == 2:
			lanes = [item.split(":")[-1] for item in self.read_group_string]
			def convert_lane(lane):
				number_string = "_lane{}".format(lane)
				platform_unit = self.platform_unit[lanes.index(lane)]
				file1 = raw_fastq_dir + self.accession + number_string + "_1.fastq.gz"
				file2 = raw_fastq_dir + self.accession + number_string + "_2.fastq.gz"
				output_file = ubam_dir + self.species + "_" + self.accession + number_string + "_unaligned_reads.bam"
				fastq_to_sam = "gatk FastqToSam -F1 {} -F2 {} -O {} -PL ILLUMINA -RG {} -SM {} -LB {} -PU {} -SO queryname --TMP_DIR {}".format(file1, file2, output_file, platform_unit, self.sample, self.library, platform_unit, tmp_dir)
				os.system(fastq_to_sam)

			Parallel(n_jobs=len(lanes))(delayed(convert_lane)(lane) for lane in lanes)

		# No multiplex
		else:
			file1 = raw_fastq_dir + self.accession + "_1.fastq.gz"
			file2 = raw_fastq_dir + self.accession + "_2.fastq.gz"
			output_file = ubam_dir + self.species + "_" + self.accession + "_unaligned_reads.bam"
			fastq_to_sam = "gatk FastqToSam -F1 {} -F2 {} -O {} -PL ILLUMINA -RG {} -SM {} -LB {} -PU {} -SO queryname --TMP_DIR {}".format(file1, file2, output_file, self.platform_unit, self.sample, self.library, self.platform_unit, tmp_dir)
			os.system(fastq_to_sam)

	def mark_illumina_adapters(self):
		print("Picard MarkIlluminaAdapters. Marking adapter sequences in unmapped BAM file.")

		# Multiplex
		if len(self.read_group_string) == 2:
			lanes = [item.split(":")[-1] for item in self.read_group_string]
			def mark_adapters(lane):
				number_string = "_lane{}".format(lane)
				ubam_file = ubam_dir + self.species + "_" + self.accession + number_string + "_unaligned_reads.bam"
				metrics_file = uBAM_XT_dir + self.species + "_" + self.accession + number_string + "_adapter_metrics.txt"
				output_file = uBAM_XT_dir + self.species + "_" + self.accession + number_string + "_unaligned_reads_XT.bam"
				mark_adapters = "gatk MarkIlluminaAdapters -I {} -M {} -O {} --TMP_DIR {}".format(ubam_file, metrics_file, output_file, tmp_dir)
				os.system(mark_adapters)

			Parallel(n_jobs=len(lanes))(delayed(mark_adapters)(lane) for lane in lanes)

		# No multiplex
		else:
			ubam_file = ubam_dir + self.species + "_" + self.accession + "_unaligned_reads.bam"
			metrics_file = uBAM_XT_dir + self.species + "_" + self.accession + "_adapter_metrics.txt"
			output_file = uBAM_XT_dir + self.species + "_" + self.accession + "_unaligned_reads_XT.bam"
			mark_adapters = "gatk MarkIlluminaAdapters -I {} -M {} -O {} --TMP_DIR {}".format(ubam_file, metrics_file, output_file, tmp_dir)
			os.system(mark_adapters)

	def align_to_reference(self):
		print("gatk SamToFastq | bwa mem | gatk MergeBamAlingment. Aligning reads to reference genome.")
	
		# Multiplex
		if len(self.read_group_string) == 2:
			lanes = [item.split(":")[-1] for item in self.read_group_string]
			def align(lane):
				number_string = "_lane{}".format(lane)
				uBAM_XT = uBAM_XT_dir + self.species + "_" + self.accession + number_string + "_unaligned_reads_XT.bam"
				unmapped_BAM = ubam_dir + self.species + "_" + self.accession + number_string + "_unaligned_reads.bam"
				output_file = mapped_bam_dir + self.species + "_" + self.accession + number_string + "_aligned_reads.bam"
				tmp = "{}{}_lane{}/".format(tmp_dir, self.accession, lane)
				align = "gatk SamToFastq -I {} --FASTQ /dev/stdout --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 __INTERLEAVE true --NON_PF true --TMP_DIR {} | bwa mem -M -t 8 -p {} /dev/stdin | gatk MergeBamAlignment --ALIGNED_BAM /dev/stdin --UNMAPPED_BAM {} -R {} -O {} --ADD_MATE_CIGAR true --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS --SORT_ORDER coordinate --CREATE_INDEX --TMP_DIR {}".format(uBAM_XT, tmp, reference_genome, unmapped_BAM, reference_genome, output_file, tmp)
				os.system(align)

				# Clean up intermediate_files
				#os.system("rm " + uBAM_XT)
				#os.system("rm " + unmapped_BAM)

			Parallel(n_jobs=len(lanes))(delayed(align)(lane) for lane in lanes)

		# No multiplex
		else:
			uBAM_XT = uBAM_XT_dir + self.species + "_" + self.accession + "_unaligned_reads_XT.bam"
			unmapped_BAM = ubam_dir + self.species + "_" + self.accession + "_unaligned_reads.bam"
			output_file = mapped_bam_dir + self.species + "_" + self.accession + "_aligned_reads.bam"
			tmp = "{}{}/".format(tmp_dir, self.accession)
			align = "gatk SamToFastq -I {} --FASTQ /dev/stdout --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true -NON_PF true --TMP_DIR {} | bwa mem -M -t 16 -p {} /dev/stdin | gatk MergeBamAlignment --ALIGNED_BAM /dev/stdin --UNMAPPED_BAM {} -R {} -O {} --ADD_MATE_CIGAR true --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS --SORT_ORDER coordinate --CREATE_INDEX --TMP_DIR {}".format(uBAM_XT, tmp, reference_genome, unmapped_BAM, reference_genome, output_file, tmp)
			os.system(align)

			# Clean up intermediate_files
			#os.system("rm " + uBAM_XT)
			#os.system("rm " + unmapped_BAM)

	def mark_duplicates(self):
		print("Picard MarkDuplicates. Marking Duplicate reads in BAM files.")

		# Multiplex
		if len(self.read_group_string) == 2:
			lanes = [item.split(":")[-1] for item in self.read_group_string]
			number_string_one = "_lane{}".format(lanes[0])
			number_string_two = "_lane{}".format(lanes[1])
			input_bam_1 = mapped_bam_dir + self.species + "_" + self.accession + number_string_one + "_aligned_reads.bam"
			input_bam_2 = mapped_bam_dir + self.species + "_" + self.accession + number_string_two + "_aligned_reads.bam"
			output_file = dedup_bam_dir + self.species + "_" + self.accession + "_dedup_aligned_reads.bam"
			metrics_file = dedup_bam_dir + self.species + "_" + self.accession + "_dedup_metrics.txt"
			mark_duplicates = "gatk MarkDuplicates -I {} -I {} -O {} --TMP_DIR {} -M {}".format(input_bam_1, input_bam_2, output_file, tmp_dir, metrics_file)
			os.system(mark_duplicates)

			# Clean up intermediate_files
			#os.system("rm " + input_bam_1)
			#os.system("rm " + input_bam_2)

		# No multiplex
		else:
			input_bam = mapped_bam_dir + self.species + "_" + self.accession + "_aligned_reads.bam"
			output_file = dedup_bam_dir + self.species + "_" + self.accession + "_dedup_aligned_reads.bam"
			metrics_file = dedup_bam_dir + self.species + "_" + self.accession + "_dedup_metrics.txt"
			mark_duplicates = "gatk MarkDuplicates -I {} -O {} --TMP_DIR {} -M {}".format(input_bam, output_file, tmp_dir, metrics_file)
			os.system(mark_duplicates)

			# Clean up intermediate_files
			#os.system("rm " + input_bam)

	def get_alignment_stats(self):
		print("Samtools flagstat. Getting alignment metrics from clean BAM files.")
	
		input_file = dedup_bam_dir + self.species + "_" + self.accession + "_dedup_aligned_reads.bam"
		output_file = dedup_bam_dir + self.species + "_" + self.accession + "_dedup_aligned_reads_flagstat.tsv"
		flagstat = "samtools flagstat -@ 16 -O tsv {} > {}".format(input_file, output_file)
		os.system(flagstat)

	def call_variants(self):
		print("gatk HaplotypeCaller. Calling variants.")
	
		input_file = dedup_bam_dir + self.species + "_" + self.accession + "_dedup_aligned_reads.bam"
		output_file = vcf_dir + self.species + "_" + self.accession + ".g.vcf.gz"
		haplotype_caller = "gatk HaplotypeCaller -R {} -I {} --native-pair-hmm-threads 10 -O {} -ERC GVCF".format(reference_genome, input_file, output_file)
		os.system(haplotype_caller)

		# Delete unnecesary BAM file. If need to retain BAM file, comment out this line
		#os.system("rm " + input_file)
		#os.system("rm " + input_file + ".bai")
		#os.system("rm " + input_file + ".sbi")

	def normalize_indels(self):
		print("bcftools norm. Left-align and normalize indels.")
	
		input_file = vcf_dir + self.species + "_" + self.accession + ".g.vcf.gz"
		output_file = vcf_dir + self.species + "_" + self.accession + "_norm.g.vcf.gz"
		norm = "bcftools norm -f {} -m- -Oz -o {} {}".format(reference_genome, output_file, input_file)
		os.system(norm)

		#Remove the original vcf file and index
		os.system("rm " + input_file)
		os.system("rm " + input_file + ".tbi")

	def index_vcf(self):
		print("gatk IndexFeatureFile. Create index for normalized vcf files.")
	
		input_file = vcf_dir + self.species + "_" + self.accession + "_norm.g.vcf.gz"
		index = "gatk IndexFeatureFile -I {}".format(input_file)
		os.system(index)

def main():
	array_id = os.environ["array_id"]
	accession_list = list(dict)
	accession_id = accession_list[int(array_id)]
	accession = Accessions(accession_id,dict[accession_id][1],dict[accession_id][2],dict[accession_id][3],dict[accession_id][4],dict[accession_id][5],dict[accession_id][6])

	accession.convert_fastq_to_unmapped_bam()
	
	accession.mark_illumina_adapters()
	accession.align_to_reference()
	accession.mark_duplicates()
	accession.get_alignment_stats()
	accession.call_variants()
	accession.normalize_indels()
	accession.index_vcf()

if __name__ == "__main__":
	main()
