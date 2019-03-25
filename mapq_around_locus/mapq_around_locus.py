import matplotlib
matplotlib.use('Agg')
import argparse
import pysam
import subprocess
import sys
import math

def mapq_around_locus(bam_file, chromosome, position, window):
    start_position = position-math.floor(window/2)
    end_position = position+math.ceil(window/2)
    pos = '{}:{}-{}'.format(chromosome, start_position, end_position)
    print(pos)
    get_qualities_command = 'samtools view {} {}'.format(bam_file, pos)
    print(get_qualities_command)

    # Run the samtools view command using subprocess, and place the stdout output into the qualities variable as a string
    samtools_view = subprocess.check_output(get_qualities_command.split())

    # Split the sam_view output string into each line comprising it
    sam_lines = samtools_view.decode("utf-8") .split('\n')
    quality_scores = []
    for line in sam_lines:
        try:
            quality = int(line.split('\t')[4])
            quality_scores.append(quality)
        except:
            sys.stdout.write("Couldn't get quality score for {}\n".format(line.split('\t')))
    print(quality_scores)


def main():
    parser = argparse.ArgumentParser(description='Run mapq_around_locus')
    parser.add_argument('--bam_file', metavar='bam_file', type=str)
    parser.add_argument('--chromosome', metavar='chromosome', type=str)
    parser.add_argument('--position', metavar='position', type=int)
    parser.add_argument('--window', metavar='window', type=int, default=50)
    parser.add_argument('--slide', metavar='slide', type=int, default=25)
    parser.add_argument('--max_distance', metavar='max_distance', type=int, default=100)
    args = parser.parse_args()

    bam_file = args.bam_file
    chromosome = args.chromosome
    position = args.position
    window = args.slide
    slide = args.slide
    max_distance = args.max_distance

    for i in range(0, max_distance, slide):a
        mapq_around_locus(bam_file, chromosome, position+i*slide, window)


if __name__ == '__main__':
    main()