#!/usr/bin/env python
import re, cairo, random, argparse

def get_args():
    parser = argparse.ArgumentParser(description="This script inputs a fasta file and a file of motifs and goutputs a graph showing where the motifs can bind in each sequence.")
    parser.add_argument("-f", "--input_file", help="Your fasta file.", required = True)
    parser.add_argument("-m", "--motif_file", help="Your motifs file, must be one motif per line.", required = True)

    return parser.parse_args()

args = get_args()
fasta_file = args.input_file
motif_file = args.motif_file

fasta = "/home/ndr/bgmp/bioinfo/Bi625/OOP/Figure_1.fasta"
#motifs = "/home/ndr/bgmp/bioinfo/Bi625/OOP/Fig_1_motifs.txt"

ambigs_dict_2 = {'A':['[A]'], 'C':['[C]'], 'G':['[G]'], 'T':['[TU]'], 'U':['[UT]'], 'W':['[ATU]'], 'S':['[CG]'],
'M':['[AC]'], 'K':['[GTU]'], 'R':['[AG]'], 'Y':['[CTU]'],'B':['[CGTU]'], 'D':['[AGTU]'], 'H':['[ACTU]'],
'V':['[ACG]'], 'N':['[ACGTU]']}

def create_motif_list(motifs):
    """Creates list from motif text file."""
    motif_list = []
    with open (motif_file, "r") as fh:
        while True:
            line = fh.readline()
            if line =='':
                break
            if line[-1] != "\n":
                line = line + "\n"
            motif_list.append(line)
    return(motif_list)
motif_list = create_motif_list(motif_file)

def convert_to_two_line_fasta_dict(input_file):
    """Converts multilined FASTA files to two-lined dictionary, where keys = gene name, values = sequence."""
    seq = ""     #initialize seq to ""
    first_line=True
    full_dict = {}
    with open(fasta_file, "r") as fh:
        while True:
            line = fh.readline().strip()
            if line =='':
                break
            if line[0] == '>':
                if first_line:
                    header=line
                    first_line=False
                else:
                    full_dict[header] = seq
                    header = line
                    seq = ""
            else:
                seq += line
    full_dict[header] = seq
    return full_dict
fasta_dict = convert_to_two_line_fasta_dict(fasta_file)

#determines max length of gene infasta file
max_gene_len = int()
for gene in fasta_dict.values():
    if max_gene_len < len(gene):
        max_gene_len = len(gene)

#used for coloring the motifs
random_list = []
for x in range(10):
    random_list_2 = []
    for y in range(3):
        random_list_2.append(random.random())
    random_list.append(random_list_2)
    random_list_2 = []

class Motif:
    def __init__(self, motif, seq):
        """
        *Maps and draws motifs: 
        *Different color per motif
        *Accounts for overlap in sequences and graph
        *Input is a motif and sequence.
        """
        self.match = [] #regex pattern of motif
        self.start = [] #start position of motif in sequence
        self.length = int() #length of motif
        self.color = [] #color of motif
        self.name = '' #name of motif

        pos_list = []
        flat_list = []
        motif = motif.strip()

        #create regex pattern for each motif
        self.name = motif
        self.length = len(motif)
        individual_list = []
        for symbol in motif.upper():
            individual_list.append(ambigs_dict_2[symbol])
        flat_list = [x for l in individual_list for x in l]
        match = ''.join(flat_list)
        self.match = match

        #finds where each motif matches in every sequence
        for pos in re.finditer(rf"(?={self.match})", seq, re.I):
            self.start.append(pos.span()[0])

    def draw(self, ctx):
        #drawing the exons
        for start_pos in self.motifs.start:
            start_pos = 10 + start_pos
            end_pos = start_pos + self.motifs.length 
            y_pos = (100 + group[2] * 85) - (30 * random.random())
            ctx.rectangle(start_pos, y_pos, self.motifs.length, -2)
            ctx.set_source_rgba(self.color[0], self.color[1], self.color[2]) #save this
            ctx.fill()  
            
gene_group_list = [] #list of [gene name, sequence, index, and exon start position and length] for each entry in fasta file
class GeneGroup:
    def __init__(self, fasta_dict):
        """
        *Draws genes: 
        *Exons shown as large outline boxes.
        *Input is fasta file dictionary.
        """
        self.name = '' #gene name
        self.seq = '' #gene sequence
        self.motifs = [] #contains start positions of motifs in sequence
        self.gene_number = int() #used to offset drawing code
        self.exon = [] #contains start position and length of each exon 
        for index, (gene, seq) in enumerate(fasta_dict.items()):
            self.gene_number = index
            self.name = gene
            self.seq = seq

            total_length = len(seq)
            exon = "[A-Z]"
            exon_start = re.search(exon, seq)
            exon_end = re.search(exon, seq[::-1])
            exon_start = exon_start.span()[0]
            exon_length = (total_length - exon_end.span()[0] - 1) - exon_start
            self.exon = exon_start, exon_length
            list_entry = [gene, seq, index, (exon_start, exon_length)]
            gene_group_list.append(list_entry)

 
    def draw(self, ctx):
        """
        *Draws gene: 
        *Denotes exons and introns. 
        *Creates legend for motifs.
        *Inputs a PyCairo surface.
        """
        ctx.rectangle(10, 100 + (group[2] * 85), len(group[1]), 1)
        #draw genes
        ctx.rectangle(10 + group[3][0], 100 + (group[2] * 85), group[3][1], -10)
        ctx.set_source_rgba(0.0, 0.0, 0.0)
        ctx.fill_preserve()
        ctx.set_source_rgba(1, 1, 1)
        ctx.set_line_width(0.5)
        ctx.stroke()

        #add gene names
        ctx.set_source_rgba(1, 1, 1)
        ctx.move_to(10, 100 + (group[2] * 85) + 15) 
        ctx.show_text(group[0])

        for index, motif in enumerate(motif_list):
            self.color = random_list[index]
            self.motifs = Motif(motif, group[1])
            
            #create legend
            ctx.set_source_rgba(self.color[0], self.color[1], self.color[2])
            ctx.rectangle((index) * 100 + 10,  25 , 10, 10)
            ctx.fill()
            ctx.move_to((index * 100 + 25), 32)
            ctx.show_text(self.motifs.name)

            #draw motifs 
            Motif.draw(self, ctx)

##########
x = GeneGroup(fasta_dict)

surface = cairo.ImageSurface(cairo.FORMAT_RGB24, max_gene_len + 50, len(gene_group_list) * 95)
ctx = cairo.Context(surface)
for group in gene_group_list:
    x.draw(ctx)

surface.write_to_png(f'{fasta_file}.png')
##########