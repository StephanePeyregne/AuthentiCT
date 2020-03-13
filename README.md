AuthentiCT is a command-line tool to estimate the proportion of present-day DNA contamination in ancient DNA datasets generated from single-stranded libraries.
It estimates contamination using an ancient DNA damage model and assumes that the contaminant is not deaminated.

Installation:
=============
AuthentiCT requires Python 3.6+ and relies on a number of non-standard libraries. Here is the list of these dependencies with the versions that we used:
- numpy (version 1.17.3)
- numdifftools (version 0.9.39)
- pandas (version 0.25.2)
- scipy (version 1.3.1)

After downloading or cloning this repository, you can install AuthentiCT by typing from the terminal whic should install all necessary dependencies:
```
pip3 install [path to the downloaded repository]
```
It should also be installed if you run the setup.py file like so: 
```
python3 setup.py install
```
The installation AuthentiCT should run on both linux and macOS. If you have any issue with the installation, contact me at: stephane\_peyregne@eva.mpg.de

Inputs:
=======
The input file should be in the SAM format. The presence of a MD field for each sequence is necessary to identify substitutions to the reference genome. This MD field should be in the standard format as generated with the samtools calmd command.

You can also provide a configuration file specifying the values of the parameters. This configuration file should list on each line the name and value of a parameter in a tab-separated format. Here is an example:
```
e	0.003170
rss	0.602679
lo	0.261435
lss	0.142886
lds	0.008414
rds	0.035792
contam	0.001000
o	0.637999
o2	0.477368
```

Here is what they correspond to:
```
e = error rate,
rss = rate of C-to-T substitutions in single-stranded regions,
lo = length of single-stranded overhangs,
lss = length of internal single-stranded regions,
lds = length of double-stranded regions,
contam = contamination estiamte,
o = the frequency of 5' single-stranded overhangs,
o2 is proportional to the frequency of 3' single-stranded overhangs
```

Note that AuthentiCT only applies to sequences generated from single-stranded libraries and is not tailored for the more complex deamination patterns of sequences generated from double-stranded libraries. Also, it cannot estimate contamination from samples treated with Uracil-DNA Glycosylase (UDG) or other enzymes that alter the deamination patterns. For the same reason, it is not possible to estimate contamination in a filtered dataset of only deaminated sequences.

Commands:
=========
AuthentiCT has 3 sub-commands: deam2cont to estimate contamination, deamination to print C-to-T substitution frequencies, simulation to simulate sequences under the ancient DNA damage model. We describe the options of each command below.

deam2cont
---------
	-o	name of output file (by default, it uses the standard output)
	-t	estimates contamination from the frequencies of C-to-T substitutions at the first and last positions of sequences
	-m	mapping quality cutoff
	-l	sequence length cutoff
	-b	base quality cutoff
	-s	maximum number of sequences used to fit the deamination model (the default is 100,000 sequences)
	-p	only keep sequences overlapping provided positions
	--decoding	prints the posterior probabilities of each state, line by line for each informative site

Example of a command:
```
AuthentiCT deam2cont -o [output name] -s 10000 [input.sam]
```
or if you have a file in the BAM format:
```
samtools view input.bam | AuthentiCT deam2cont -o [output name] -s 10000 -
```

deamination
-----------
	-o	name of output file (by default, it uses the standard output)
	-d	computes the observed and expected distance between internal C-to-T substitutions and disregard the specified number of positions at both ends of sequences
	-m	mappinq quality cutoff
	-l	sequence length cutoff
	-b	base quality cutoff

simulation
----------
	-o	name of output file (by default, it uses the standard output)
	-N	number of simulated sequences
	-GC	proportion of G and C bases in the sequences
	-L	average length of sequences
	-l	minimum length of seauences

Outputs:
========

deam2cont
---------
deam2cont outputs by default a list of the maximum likelihood estimates of the parameter values (2nd column) with their associated standard error (3rd column). 
It then outputs information about the sequences:

| sequence name | sequence of observations | flag (0 for forward, 16 for reverse) | sequence length | vector of posterior probabilities for the 5' single-stranded state | vector of posterior probabilities for the double-stranded state | vector of posterior probabilities for the single-stranded state | vector of posterior probabilities for the 3' single-stranded state | vector of informative positions in the sequence | probability that the sequence is endogenous | likelihood from the endogenous model | likelihood from the contaminant model | score (likelihood ratio) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

One can also print the posterior probabilities line by line for each site using the option --decoding. This will output the position, the length of the sequence, the base at that position, the posterior probabilities of 5’ single-stranded overhang, double-stranded, single-stranded and 3’ single-stranded overhang states, respectively, and, finally, the sequence name.

simulation
----------
The output for this command is in the SAM format. We add a ST field (e.g. ST:Z:555dddddsssdddd3333) that corresponds to the state for each base with 5, d, s and 3 corresponding to the 5’ single-stranded overhang, double-stranded, single-stranded and 3’ single-stranded overhang states, respectively.

citation
--------
The method is described in a preprint that will be soon available on bioRxiv.

