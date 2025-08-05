## directives


#### rule all: and expand()

Let's say we have the same manifest file as before:

```
loci/locus1.fasta
loci/locus2.fasta
loci/locus3.fasta
...
loci/locus1453.fasta
```

How will Snakemake know what to insert for `{locus_id}` in our `mafft_align` rule?