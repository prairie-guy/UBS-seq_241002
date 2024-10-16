## reference/

- Directory for large references, indexes and genomes

- Except for this file, nothing else in this dir will be stored in github

- To include additional files or dirs, edit `.gitignore` to append '!' to expended pathname

### References used for BAT-seq repos

The dir `~/uofc_data/bat_seq/Projects/reference/` should be used by specific repo versions of BAT-seq

#### fasta/
- Origial primary sequences for BAT-seq are located in: https://github.com/prairie-guy/Genomic_References/tree/main/BAT_seq
   - The Genomic_References repo lives locally at: ~/uofc_data/Genomic_References

#### hisat3n
- histat3n indexes for fasta sequences. Sample build:
```
hisat-3n-build --base-change C,T ../fasta/lambda.fa lambda
```
