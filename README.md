# Vcffilter

*Vcffilter* provides a tool set for fast extraction of genotype data from VCF files. The set includes:

- **vcffilter** for filtering and extracting information from VCF
- **restorevcf** to restore a VCF file from extracted information with *vcffilter* and optionally apply further filters
- **myzcat** for a (probably) faster unpacking of gzipped files to stdout

### Prerequisites:

Vcffilter requires a standard C/C++ build chain and the packages *zlib* and *libboost-program-options-dev*.
For building simply type `make` in the *Release* subfolder.
For pre- and post-processing *bcftools* and *bgzip* might be useful.


## vcffilter

*vcffilter* requires a valid **unpacked** VCF-file stream from stdin. Usually, VCF files are compressed as `.vcf.gz`. So, to generate a valid VCF-file stream, you need to use *zcat* or *myzcat* (see below).
*vcffilter* only extracts the following information from the input file:

- the chromosome (only printed once in the header of the output, taken from the first data line after the header from the `CHR` column)
- `POS` column
- `ID` column
- `REF` column
- `ALT` column
- `QUAL` column
- `FILTER` column
- `INFO` column
- the genotypes `GT` from the genotype columns
- *optional:* `GQ` from the genotype columns (if present and the `--gq` switch is provided)

You might want to compress the information after extraction again, which you could do using *gzip*.

#### Example:

```
zcat input.vcf.gz | vcffilter | gzip -c > compressed_extraction.gz
```

Note, that *vcffilter* prints some additional information to *stderr* during the run.

### Important!! Requirements for input VCFs:

The genotype (`GT`) field **must be the first field** in the genotype columns **and it MUST NOT be the only information** in the genotype columns (as the parser searches for the double-colon ":" character for the end of the GT field).


## restorevcf

If you want to restore your compressed and extracted data, you can use *restorevcf* to restore a valid VCF file.
Per default, *restorevcf* also restores all information from the `INFO` column with the only exception that the `AF`, `AC` and `AN` fields will be replaced by `OrgAF`, `OrgAC` and `OrgAN` containing the original information, and `AF`, `AC` and `AN` will be recreated with the actual re-calculated allele frequency, allele count and allele number during restoration.

#### Important note:

*restorevcf* does not produce a VCF header which is required for a valid VCF. 
You need to provide the header in a seperate file first.
If you want to extract a header from a VCF file, you could use `bcftools view -h input.vcf.gz > header.vcf`.

#### Example:

```
zcat compressed_extraction.gz | restorevcf > uncompressed_extraction
cat header.vcf uncompressed_extraction |  bgzip -c > restored.vcf.gz
rm uncompressed_extraction
```

Note the use of *bgzip* instead of *gzip* as valid VCF files need to be able to be indexed, which is not possible using *gzip*.

Alternatively, you can use `bcftools convert` to generate a file output in your desired format, e.g. to generate a `.bcf` file:

```
zcat compressed_extraction.gz | restorevcf > uncompressed_extraction
cat header.vcf uncompressed_extraction | bcftools convert - -Ob -o restored.bcf
rm uncompressed_extraction
```

#### Optional filter and conversion options:

*restorevcf* provides optional filter and conversion options which will be applied on-the-fly during restoration. For a full list of options type `restorevcf --help`.

Here is an excerpt of the most important filters:

- `--fpass` returns only variants with `PASS` in the `FILTER` column
- `--rminfo` removes all fields from the `INFO` column (but still re-generates `AF`, `AC` and `AN`)
- `keepaa` if present, the `AAScore` field will be kept if `--rminfo` was applied
- `--macfilter` keeps only variants with a minor allele count greater or equal the provided number
- `--maffilter` keeps only variants with a minor allele frequency greater or equal the provided number
- `--aafilter` keeps only variants with an `AAScore` greater or equal the provided number
- `--missfilter` keeps only variants with a missingness rate below the provided number
- `--filterunknown` removes all variants with unknown alleles (named `*`)
- `--splitma` splits multi-allelic variants into several bi-allelic ones, filling up with the reference `0` (implies `--rminfo`).

#### Example:

The following command restores a VCF with on-the-fly filtering for `PASS` in the `FILTER` column, a missingness rate below 0.1, an `AAScore` >= 0.8 and a MAC >= 4, removal of unknown alleles and splitting of multi-allelic variants to bi-allelics (which implies the removal of all `INFO` fields besides `AF`, `AC` and `AN`, but also keeps `AAScore` due to the `--keepaa` switch).

```
zcat compressed_extraction.gz | restorevcf --fpass --missfilter 0.1 --aafilter 0.8 --macfilter 4 --filterunknown --splitma --keepaa > uncompressed_extraction
cat header.vcf uncompressed_extraction |  bgzip -c > restored.vcf.gz
rm uncompressed_extraction
```

## removesamples

*removesamples* requires a file as argument that contains the IDs of samples (one exclusively in each line). *removesamples* reads an uncompressed VCF file from *stdin* and writes uncompressed VCF to *stdout*, the samples in the input file are removed during this process (if they are found). 
Informational, warning and error messages are written to *stderr*.

**Note:** *removesamples* removes all information from the *INFO* column. However, it recalculates and sets the tags for *AC (allele count)* and *AN (allele number)*. Note, that multi-allelics are probably not counted correctly as *AC* reflects the number of known (i.e. not missing) non-zero alleles. Unknown (i.e. missing) alleles are still counted for *AN*.

**Note:** Further note, that *removesamples* requires the *genotype (GT)* to be the first field in each sample column (which is the usual case). The *FORMAT* column will not be checked.

#### Optional filters:

- `--macfilter` keeps only variants with a minor allele count greater or equal the provided number
- `--maffilter` keeps only variants with a minor allele frequency greater or equal the provided number
- `--missfilter` keeps only variants with a missingness rate below the provided number

#### Example:

*removesamples* is much faster than applying *bcftools* with the *-S* option, followed by the *bcftools +fill-tags* plugin and a final *bcftools filter* call. However, compression and decompression will quickly become the bottleneck (which is again slow in *bcftools*). Thus, I recommend to use the *vcf.gz* format for VCF input and output and use the *bgzip* tool with the *--threads* option for compression and decompression:

```
bgzip -d -c --threads 4 input.vcf.gz | \
  removesamples my_exclude_samples_file --macfilter 4 --missfilter 0.1 | \
  bgzip --threads 4 > output.vcf.gz
```

## myzcat

*myzcat* can be used to replace *zcat*. It might be a little bit faster than the original *zcat* as it pre-allocates a large buffer of 1 GB for unpacking at the beginning, what *zcat* usually doesn't do. However, using *bgzip* with *--threads* might still be faster when in a multi-threaded environment.
