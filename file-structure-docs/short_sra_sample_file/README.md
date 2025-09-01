These files were created by starting with `short_sra_test.fastq`

```fastq
@short_sra_test.1/1
ACGT
+
!mn~
```

Then running:

```sh
latf-load --quality PHRED_33 short_sra_test.fastq -o short_sra_test.sra.dir
kar --create short_sra_test.sra --directory short_sra_test.sra.dir/
hexdump -C short_sra_test.sra > short_sra_test.sra.hexdump.txt
tree short_sra_test.sra.dir/ > short_sra_test.sra.dir.tree
```

