# Analyzing monkeypox sequencing data

# Clone the repository

```
git clone git@github.com:fredericlemoine/mpx_apobec3f.git
cd mxpv_analysis
```

# Required dependencies

Install java and [singularity](https://sylabs.io/singularity/).

# Install [Nextflow](https://www.nextflow.io/)

If not already done, first install nextflow:
```
mkdir $HOME/bin
pushd $HOME/bin
curl -s https://get.nextflow.io | bash
popd
echo 'export PATH=$PATH:$HOME/bin' >> $HOME/.bashrc
```

# Run the workflow

```
nextflow run main.nf --results results_batch_x/ --data 'data/batch_x/*fas' --reference data/reference.fasta
```

This will create the folder `results_batch_x`, with the following files:
- `align.fasta`: The "raw" MSA
- `align_clean.fasta`: The MSA without 70% N containing sequences
- `mutations.txt`: File containing the number of occurences of each mutations for each samples, with the following columns:
    - Sample: Name of the file (Ex: `PLATE_1_D02_M13uni-21-D02`)
    - Annot1: Long annotation (Ex: `D02`)
    - Annot: Short annotation (Ex: `D`)
    - Type: Type of mutation (C-T, G-A or Others)
    - Position: Position of the mutation on the reference sequence (0-based)
    - Ref: Reference nucleotide (no rev comp)
    - Comp: Mutated nucleotide (no rev comp)
    - TrinucSeq1-TrinucSeq2: Reference context -1,0,+1 - Mutated context -1,0,+1 (no rev comp)
    - Strand: If G-A, then strand=-
    - RefContext: Reference context, if strand=- (G-A), then rev comp
- `reference_cut.fasta`: The sub sequence of interest of the reference sequence
