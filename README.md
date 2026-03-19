# KappaMapScope

Calculates the mappability of a DNA sequence over different values of k and makes a figure of it.

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
    --output-prefix ./results/sequences_k100 \
    -t 4 --cleanup --create-bam
```

Calculate the mappability score of a single fasta file over a list of k values
```

```
