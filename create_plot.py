import argparse
import numpy as np
import gzip
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')

def get_quality_scores_from_fastq(file_path, limit=None):
    quality_scores = []
    with gzip.open(file_path, 'rt') as f:
        line_count = 0
        for line in f:
            line = line.strip()
            line_count += 1
            
            # Check if the line is a quality scores line
            if line_count % 4 == 0:
                quality_scores.append(line)
            if limit and len(quality_scores) >= limit:
                break
                
    return quality_scores

def qscores_list_to_mat(qscores_list):
    bchar_vals_list = []
    max_len = 0
    for line in qscores_list:
        if len(line) > max_len:
            max_len = len(line)
    for line in qscores_list:
        bchar_vals = np.asarray([ord(item)-33 for item in list(line)])
        bchar_vals = np.pad(bchar_vals, (0, max_len-len(line)), 'constant')
        bchar_vals_list.append(bchar_vals)
    bchar_vals_mat = np.asarray(bchar_vals_list)
    return bchar_vals_mat.mean(0)

def collect_bam_qscores(bam_scores_file, limit=None):
    with open(bam_scores_file) as bfile:
        if limit:
            blines = bfile.readlines()[0:limit]
        else:
            blines = bfile.readlines()
        blines = [item.strip() for item in blines]
        return qscores_list_to_mat(blines)

def collect_fastq_scores(fastq_file, limit):
    qscores = get_quality_scores_from_fastq(fastq_file, limit)
    return qscores_list_to_mat(qscores)

def create_plots(bam_scores_file, read1, read2, limit=None, output_fname='plot.png'):
    print(f"Fetching bam qscores ... {bam_scores_file}")
    mean_bam_scores = collect_bam_qscores(bam_scores_file, limit)

    print(f"Fetching read1 qscores ... ")
    read1_qscores = collect_fastq_scores(read1, limit)
    print(f"Fetching read2 qscores ... ")
    read2_qscores = collect_fastq_scores(read2, limit)

    #print("Plotting bam scores")
    plt.plot(mean_bam_scores, color='blue', label='Alignment inferred quality scores')

    # Create a plot for the second array (green color)
    plt.plot(read1_qscores, color='green', label='Read 1: instrument reported phred scores')

    # Create a plot for the third array (red color)
    #print("Plotting read2 scores")
    plt.plot(read2_qscores, color='red', label='Read 2: instrument reported phred scores')

    # Set x and y axis labels
    plt.xlabel('Read position')
    plt.ylabel('Q scores')

    plt.legend()

    print("Saving the plot now ... ")
    plt.savefig(output_fname, dpi=200)
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam_qc_output_file")
    parser.add_argument("-1", "--read1")
    parser.add_argument("-2", "--read2")
    parser.add_argument("-l",'--limit', default=None)
    parser.add_argument("-o", "--output_file", default="plot.png")
    args = parser.parse_args()
    
    create_plots(args.bam_scores_file, args.read1, args.read2, args.limit, args.output_file)

if __name__ == '__main__':
    main()