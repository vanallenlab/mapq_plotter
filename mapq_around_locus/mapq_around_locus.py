import matplotlib.pyplot as plt
#matplotlib.use('Agg')
import argparse
import pysam
import subprocess
import sys
import math
import pandas as pd
import seaborn as sns
import random


def mapq_around_locus(bam_file, chromosome, position, window, offset):
    """Get the distribution of quality scores for all reads falling within the window of size 'window' centered around
    'position', and return a dataframe column containing these quality scores."""
    start_position = position-math.floor(window/2)
    end_position = position+math.ceil(window/2)
    pos = '{}:{}-{}'.format(chromosome, start_position, end_position)
    get_qualities_command = 'samtools view {} {}'.format(bam_file, pos)

    # Run the samtools view command using subprocess
    samtools_view = subprocess.check_output(get_qualities_command.split())

    # Split the sam_view output string into each line comprising it
    sam_lines = samtools_view.decode("utf-8") .split('\n')
    quality_scores = []
    for line in sam_lines:
        try:
            quality = int(line.split('\t')[4])
            quality_scores.append(quality)
        except:
            sys.stderr.write("Couldn't get quality score for {}\n".format(line.split('\t')))
    #quality_scores = random.sample(range(0, 60), 30)
    return pd.DataFrame(quality_scores, columns=[offset])


def plot_mapq_around_locus(bam_file, chromosome, position, window, slide, max_distance):
    """Plot the distribution of read quality scores"""
    dfs = []
    for offset in range(-max_distance, max_distance+slide, slide):
        dfs.append(mapq_around_locus(bam_file, chromosome, position + offset, window, offset))

    merged_df = pd.DataFrame()
    for df in dfs:
        merged_df = pd.concat([merged_df, df], axis=1)

    quantiles_df = pd.DataFrame()
    for df in dfs:
        quantiles = df.quantile([.2, .4, .6, .8, 1], axis=0)
        quantiles_df = pd.concat([quantiles_df, quantiles], axis=1)

    sns.set_style('whitegrid')
    plt.xticks(size=8)
    plt.ylabel('MapQ Score')
    plt.xlabel('Offset (base pairs)')
    plt.title('Read MAPQ quantiles in {} bp range around {}:{}'.format(max_distance*2, chromosome, position))

    quantiles_values = ['20%', '40%', '60%', '80%', '100%']
    quantiles_df.index = pd.Index(quantiles_values)

    quantiles_df = quantiles_df.T

    plt.plot(quantiles_df)
    plt.savefig('{}:{}.png'.format(chromosome, position))
    plt.figure()


def main():
    parser = argparse.ArgumentParser(description='Run mapq_around_locus')
    parser.add_argument('--bam_file', metavar='bam_file', type=str)
    parser.add_argument('--chromosome', metavar='chromosome', type=str)
    parser.add_argument('--position', metavar='position', type=int)
    parser.add_argument('--positions_file', metavar='positions_file', type=str, default=None)
    parser.add_argument('--window', metavar='window', type=int, default=100)
    parser.add_argument('--slide', metavar='slide', type=int, default=50)
    parser.add_argument('--max_distance', metavar='max_distance', type=int, default=400)
    args = parser.parse_args()

    bam_file = args.bam_file
    chromosome = args.chromosome
    position = args.position
    positions_file = args.positions_file
    window = args.window
    slide = args.slide
    max_distance = args.max_distance

    if positions_file:
        # If a positions file is specified, output plots for each of the positions listed in the file
        positions = pd.read_csv(positions_file, sep='\t', header='infer')
        for row in positions.iterrows():
            row = row[1]
            chromosome = row.Chromosome
            position = row.Start
            plot_mapq_around_locus(bam_file, chromosome, position, window, slide, max_distance)
    else:
        # If a positions file is not specified, but chromosome and position are specified, just do a one-time plotting:
        plot_mapq_around_locus(bam_file, chromosome, position, window, slide, max_distance)


if __name__ == '__main__':
    main()