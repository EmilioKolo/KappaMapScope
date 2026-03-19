# KappaMapScope

Calculates the mappability of a DNA sequence over different values of k and makes a figure of it.

Includes a separate script to plot the mappability setting limits and log-scales.

# Requirements

Tools: 
* Bowtie2

Python3 libraries: 
* matplotlib
* pysam

# Example script

Calculate the mappability score of a single fasta file on a single k value
```
python cli_single.py \
    -i ./input/sequences.fasta \
    -k 100 \
    -x ./input/sequences.index \
    --output-prefix ./results/sequences_mappability \
    -t 4 --cleanup --create-bam
```

Calculate the mappability score of a single fasta file over a list of k values
```
python cli_sweep.py \
    -i ./input/sequences.fasta \
    -x ./input/sequences.index \
    --k-list 50,100,200,300,400,500,750,1000 \
    --output-prefix ./results/sequences_mappability \
    -p 2 -t 4
```
