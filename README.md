# mapq_plotter

A script to help you visualize the distribution of MAPQ scores for the reads around a given position in your bam file.

# Example usage at single location:
python mapq_around_locus.py
--bam_file mother.bam
--chromosome 20
--position 16059512
--output_dir outputs

# Example usage at several locations specified in a tab-delimited file with a 'Chromosome' and 'Start' column:
python mapq_around_locus.py
--bam_file mother.bam
--positions_file positions.txt
--output_dir outputs

# Example positions.txt file contents (tab-delimited):
Chromosome	Start

20	16059226

20	16061080

20	16064664

20	16067971

20	16073033

20	16078668

