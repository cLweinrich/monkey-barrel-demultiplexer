# Monkey Barrel Demultiplexer
A demultiplexing pipeline designed specifically to demultiplex barcoded sequences that have been ligated together end-to-end.

# Amplicons are sometimes (unintentionally) ligated end-to-end during library prep
During preparation for a sequencing run on an ONT MinION, the adapater ligation step requires the addition of a T4 ligase. 
The T4 ligase is intended to ligate the adapters to sequences which are connected together by a 1bp A-T overhang.

However, read length distributions from the MinION show a peak at the expected product length, but then one or more smaller peaks at perfect multiples of the expected length (e.g. a 300bp product has peaks at 300, 600, 900, etc.) Analysis of the final sequences reveals this is a result of products being occasionally being joined together end-to-end. 

This is presumably due the T4 ligase having low blunt-end ligation activity. Normally, the blunt end activity is negligible because products should be A-tailed. However, if the previous A-tailing step went poorly, there will be products without A-overhangs that can be ligated together by this blunt-end activity.

Although the resultant reads are multimeric amalgamations of various products, the information should still be perfectly salvageable since the constituent monomers with identifying barcodes are intact, albeit strung together like a chain of monkeys.

Hence, I have written this script to search for and extract monomers with a valid barcode pairing from within a larger sequence. Monomers can then be processed as usual. 

(As a sidenote, since the MinION is primarily designed for long read sequencing, *intentionally* joining shorter reads together might not be a bad idea. If one were to try this, I would imagine that *BEFORE* the A-tailing step, one would perform an additional ligation step to chain products together, *then* do A-tailing and finish up with a second ligation step to connect adapters.)

# Requirements
The demultiplexing script requires [regex](https://pypi.org/project/regex/) and [Biopython](https://biopython.org/). 

The full Snakemake pipeline additionally requires [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [VSEARCH](https://github.com/torognes/vsearch). 

(All these dependencies are available via conda/mamba install)
