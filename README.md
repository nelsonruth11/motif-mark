## **Overview**

This script inputs a fasta file and a text file containing motifs and generates a .png figure that showing where each motif prefers to bind on every gene in the file. 
<br>

The motif file should contain one motif per line. Motifs can be ambiguous if they follow [standard IUPAC rules](https://en.wikipedia.org/wiki/Nucleic_acid_notation). One exception to this is that Uracil and Thymine are treated interchangeably here.
<br>

## **Details**

Using PyCairo to draw, the script uses object oriented programming to draw motif objects onto gene objects. 

Argparse is required here to designate your fasta file and motif file, and the output .png graph will inherit the name of your fasta file. For example, Fig_1.fasta will create Fig_1.png.

There is a small element of randomness involved with graphing. Each motif is colored using a series of three random numbers ranging from (0-1). If the user is unhappy with the way the graph looked, simply rerun the scrip. If the user wants to keep the current color scheme, they can set a seed following [these instructions](https://docs.python.org/3/library/random.html). The y values for each motif binding site are also influenced by a random number to allow viewing of overlapping motifs.
<br>

This script can process any number of genes of any length, as well as any number of motifs, though after about 10 different motifs things get messy due to excessive sequence overlap and motif color similarity.
<br>

## **Usage**

The user does not need to change any code in the file, simply add the argparse options for the fasta and motif file then run the script. The output figure will be added to the folder containing the script.
<br>
