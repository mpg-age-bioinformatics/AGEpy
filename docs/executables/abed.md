## Intro

`abed` is annotation tool for bed files.

It annotates bed files with gene names, gene ids, and feature information from a
provided annotations file (GTF).

## Examples

```
$ abed -b K27AC_chip1_peaks.bed -g hg38.83.gtf -s hg38.83.genome \
  -c "chr,start,end,name,signal_value,strand,fold_change,p_value,Benjamini_Hochberg_FDR" \
  -p 1000,200 -o annotated.bed.tsv
```

## Output files

```
chr     start   end     name    signal_value    strand  fold_change     p_value Benjamini_Hochberg_FDR  annotated_gene_features
1       1205710 1205930 chip1_peak_2799 241.68113891    .       12.1486518993   1.18025200195e-06       0.00504618822857        TNFRSF18: promoter
1       1616560 1616780 chip1_peak_5889 71.050614487    .       3.57152066778   5.81256160902e-06       0.0113928182561 RP11-345P4.9: promoter; MIB2: five_prime_utr, exon, promoter, CDS
1       1892440 1892660 chip1_peak_3527 243.582136289   .       12.2442098543   1.93910073064e-06       0.00651852531479        RP1-140A9.1: exon
1       2212540 2212870 chip1_peak_25   81.3040545107   .       4.08693314134   9.22431065693e-12       4.99354651514e-06       RP11-181G12.4: exon, promoter; FAAP20: exon, promoter
1       3712500 3712720 chip1_peak_6234 38.8954679096   .       1.95516912169   6.52541908518e-06       0.0120595193967 RP5-1092A11.2; TP73
1       3772780 3773000 chip1_peak_4768 120.93909338    .       6.0792784787    3.69904369916e-06       0.00905667169324        SMIM1: five_prime_utr, exon, promoter
1       4680280 4680500 chip1_peak_7707 110.246753848   .       5.54180372353   9.96101899764e-06       0.014735318932  AJAP1
1       5652020 5652240 chip1_peak_526  145.04841094    .       7.29118813738   2.00633103721e-08       0.000477258277665       RP11-154H17.1
1       6330720 6330940 chip1_peak_7213 153.651918104   .       7.72366298467   8.71110008982e-06       0.0138290963917 ACOT7
1       6362730 6362950 chip1_peak_7153 87.6949697748   .       4.40818702657   8.54621275409e-06       0.0136908694171 ACOT7
1       6421360 6421580 chip1_peak_5440 279.339262378   .       14.0416230896   4.90344812744e-06       0.0104638228061 HES2
1       6423890 6424220 chip1_peak_3597 398.194019503   .       20.0161276677   2.01692509039e-06       0.00664985446589        HES2; ESPN: promoter
1       6859710 6860040 chip1_peak_53   79.510332304    .       3.99676761667   4.53361213715e-11       1.10018291319e-05       CAMTA1
```

## Help

```
$ abed --help

usage: abed [-h] [-b BED] [-g GTF] [-s SIZES] [-c COLUMNS] [-p PROMOTER]
            [-o OUTPUT]

abed is an annotation tool for bed files.

optional arguments:
  -h, --help            show this help message and exit
  -b BED, --bed BED     /path/to/file.bed (default: None)
  -g GTF, --gtf GTF     /path/to/file.gtf (default: None)
  -s SIZES, --sizes SIZES
                        /path/to/file.genome. Tab separated values of
                        'chromosome name' and 'size' information. (default:
                        None)
  -c COLUMNS, --columns COLUMNS
                        A comma separated string of column headers to use when
                        reading in the bed file. eg.: 'chr,start,end,name'.
                        (default: None)
  -p PROMOTER, --promoter PROMOTER
                        A comma separated list containing the upstream start
                        of the promoter region from the TSS and the downstream
                        end of the promoter region from the TSS. eg.:
                        '1000,200'. (default: None)
  -o OUTPUT, --output OUTPUT
                        /path/to/output.tsv. (default: None)
```
