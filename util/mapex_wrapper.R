# Gregory Way 2017
# util/mapex_wrapper.R
#
# Will run MAPEX on input sample VCF files

library(optparse)
library(mapexr)

# Load command arguments

option_list <- list(
  optparse::make_option(c("-b", "--path_to_bam"),
                        type = "character",
                        help = "location of BAM file"),
  optparse::make_option(c("-i", "--path_to_bam_index"),
                        type = "character",
                        help = "location of BAM BAI file"),
  optparse::make_option(c("-v", "--path_to_vcf"),
                        type = "character",
                        help = "location of VCF file"),
  optparse::make_option(c("-o", "--blast_output"),
                        type = "character",
                        help = "location to write output blast file"),
  optparse::make_option(c("-l", "--blast"),
                        type = "character",
                        default = "blastn",
                        help = "location of command 'blastn'"),
  optparse::make_option(c("-f", "--results_output"),
                        type = "character",
                        help = "name of the output results file"),
  optparse::make_option(c("-d", "--blast_db"),
                        type = "character",
                        default = "/lorax/sanchezlab/shared/pdx_exomeseq/reference/combined_g1k_v37_mm9.fasta",
                        help = "location of combined mouse/human reference"),
  optparse::make_option(c("-t", "--num_threads"),
                        type = "numeric",
                        default = 1,
                        help = "the number of threads to run mapexr using"),
  optparse::make_option(c("-q", "--mapq"),
                        type = "numeric",
                        default = 1,
                        help = "minimum mapq score")
)

opt_parser <- optparse::OptionParser(option_list = option_list);
opt <- optparse::parse_args(opt_parser);

# Load arguments
bam <- opt$path_to_bam
bamidx <- opt$path_to_bam_index
variants <- opt$path_to_vcf
blastout <- opt$blast_output # Read level blast results, not required
blastpath <- opt$blast # if blast is in the users path, just "blastn" here
results_output <- opt$results_output
blastdb <- opt$blast_db
threads <- opt$num_threads # number of threads consumed by blastn
mapq <- opt$mapq # this is the default for minimum mapq score

print(bam)
print(bamidx)
print(variants)
print(blastout)
print(blastpath)
print(results_output)
print(blastdb)
print(threads)
print(mapq)


results <- mapexr::run_mapex(bam_file = bam,
                             bam_idx = bamidx,
                             variant_file = variants,
                             variant_file_type = 'vcf',
                             blast_out_file = blastout,
                             blastn_path = blastpath,
                             blast_threads = threads,
                             min_mapq = mapq)

write.table(results, results_output, sep='\t')
