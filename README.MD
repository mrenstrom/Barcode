The source code to process barcodes. This program only works for 20 base pair barcodes. The input source must be FASTQ.
This program assumes each fastq line containes: fixed start tag, 20 base pair barcode, end tag

The fastq file may contain sequence data from different tests, which are demultiplexed (by the start tag) in the metafile.

The metafile has the following format:   the first line is either "SamplePrep" or "MasterPrep"
-----------------------------------------------------------
SamplePrep
logOutputDirectory
sourcefile.fastq
startTAG1 endTag1 outputFile1
...
