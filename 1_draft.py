#!/usr/bin/env python
import re, cairo, math, random
#add argparse stuff for later

fasta = "/home/ndr/bgmp/bioinfo/Bi625/OOP/Figure_1.fasta"
motifs = "/home/ndr/bgmp/bioinfo/Bi625/OOP/Fig_1_motifs.txt"

ambigs_dict_2 = {'A':['[A]'], 'C':['[C]'], 'G':['[G]'], 'T':['[TU]'], 'U':['[UT]'], 'W':['[ATU]'], 'S':['[CG]'],
'M':['[AC]'], 'K':['[GTU]'], 'R':['[AG]'], 'Y':['[CTU]'],'B':['[CGTU]'], 'D':['[AGTU]'], 'H':['[ACTU]'],
'V':['[ACG]'], 'N':['[ACGTU]']}

def create_motif_list(motifs):
    """Creates list from motif text file."""
    motif_list = []
    with open (motifs, "r") as fh:
        while True:
            line = fh.readline()
            if line =='':
                break
            if line[-1] != "\n":
                line = line + "\n"
            motif_list.append(line)
    return(motif_list)
motif_list = create_motif_list(motifs)
#print(motif_list)

#converting fasta to two liner
#updated to return a dictionary instead of two lists
def convert_to_two_line_fasta_dict(input_file):
    """Converts multilined FASTA files to two-lined"""
    seq = ""     #initialize seq to ""
    first_line=True
    full_dict = {}
    with open(input_file, "r") as fh:
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
full_dict = convert_to_two_line_fasta_dict(fasta)
#print(len(full_dict))

class Finder:
    def __init__(self, name):
        "Determines sequence ambiguity of motifs."
        self.name = name
        #self.motif_pos = self.find_pos()
        #self.fasta = full_dict
        self.matches = []
        self.color = (random.random(), random.random(), random.random())

    def regex(self, motif): #add stuff from other objects?
        """Inputs ambiguous motif sequence and translates it to regex C/G/A/T/U search query."""
        flat_list = []
        #for motif in motif_list
        #motif = motif.strip()
        individual_list = []
        for symbol in motif.upper():
            individual_list.append(ambigs_dict_2[symbol])
        flat_list = [x for l in individual_list for x in l]
        match = ''.join(flat_list)
        return match

    def find_pos(self, match, gene, seq):
        """Finds positions of binding motifs in a fasta file."""
        #matches_list = []
        self.matches.append(motif)

        #for header, seq in fasta.items():
        self.matches.append(gene)
        pos_list = []

        for pos in re.finditer(match, seq, re.I):
            #pos_list = []
            pos_list.append(pos.span()[0])
        self.matches.append(pos_list)


class Gene:
    def __init__(self, gene):
        """Determines exon locations of a DNA sequence and length of total sequence."""
        self.gene_length = int()
        self.exon_start = int()
        self.exon_length = int()
        #self.gene = gene (maybe)
        #self.matches = Finder(find_pos)

    def traits(self, gene):        
        """Finds start position and length of exon, and length of entire gene."""
        total_length = len(gene)
        exon = "[A-Z]"
        exon_start = re.search(exon, gene)
        exon_end = re.search(exon, gene[::-1])
        self.gene_length = total_length
        self.exon_start = exon_start.span()[0]
        self.exon_length = (total_length - exon_end.span()[0] - 1) - exon_start.span()[0]

#Finds longest gene in fasta file.
max_gene_len = int()
for gene in full_dict.values():
    if max_gene_len < len(gene):
        max_gene_len = len(gene)
#print(max_gene_len)

surface = cairo.ImageSurface(cairo.FORMAT_RGB24, max_gene_len + 50, 500)
ctx = cairo.Context(surface)
#spacer used to space genes out depending on number of genes
spacer = 50
for gene, seq in full_dict.items():
    seq_obj = Gene(seq)
    seq_obj.traits(seq)
    num_genes = len(full_dict)
    spacer += round(400/num_genes)

    #adds genes as rectangles (grey)
    ctx.rectangle(10, spacer, seq_obj.gene_length, 1)
    ctx.rectangle(10 + seq_obj.exon_start, spacer, seq_obj.exon_length, - 10)
    ctx.set_source_rgb(0.6, 0.6, 0.6)
    ctx.fill()
    
    #adds gene names (yellow)
    ctx.set_source_rgb(1, 1, 0)
    ctx.move_to(10, spacer + 15)
    ctx.show_text(gene)

#halp 
    for motif in motif_list:
        motif = motif.strip()
        x = Finder(motif)
        x.regex(motif)
        regex = x.regex(motif)
        x.find_pos(regex, gene, seq)
    #print(x.color)
    #print(x.name)

        # if x.matches[1] == gene:
        #     #print(x.name, x.matches[2])

        #     for start in x.matches[2]:
        #         #print(start)
        #         ctx.rectangle(10 + start, spacer, start + len(x.name), -5)
        #         ctx.set_source_rgb(0, 0, 1)
        #         ctx.fill()


    # for w in range(0, len(x.matches),2):
    #     print(x.name)
    #     print(x.matches[w], x.matches[w+1])
    # for count, info in enumerate(x.matches):
    #     #ctx.rectangle

    #     print(count, info)

        # print(x.name)
        # print(x.matches)

# graphing_dict = {} #keys = gene name, values = [motif, binding sites]
# for gene in full_dict.keys():
#     graphing_dict[gene] = []
# print(graphing_dict)


    #print(x.matches[1])
    #print(x.matches[2])

    #if x.matches
    # for y in x.matches[2]:
    #     print(y)
    # ctx.fill()
        # ctx.rectangle(10, 10, 100, 100)
        # ctx.set_source_rgb(1, 0, 0)
        # ctx.fill()

surface.write_to_png("first_attempt.png")


# for motif in motif_list:
#     motif = motif.strip()
#     x = Finder(motif)
#     x.regex(motif)
#     regex = x.regex(motif)
#     #for header, seq in full_dict.items():
#     x.find_pos(regex, full_dict)

    #print(x.name)


# ctx.rectangle(25, 50, 50, 120)
# ctx.set_source_rgb(1, 0, 0)
# ctx.fill()

# surface.write_to_png('example.png')

# for name in x.name:
#     print(name)
#maybe += regex stuff onto a string
#make find pos take a particular sequence, then feed individual motifs into it
#maybe store in dictionary with sequence: 3 lists of matches for each sequence
#find pos is frigged up: make it do for one motif

# for motif in motif_list:
#     motif = motif.strip()
#     x = Finder(motif)
#     x.regex(motif)
#     y = x.regex(motif)
    #print(x.name)
    # matches_list = []

    # for header, seq in full_dict.items():
    #     matches_list = []
    #     matches_list.append(motif)
    #     matches_list.append(header)
    #     for pos in re.finditer(y, seq, re.I):
    #         matches_list.append(pos.span())
    #     print(matches_list)



    #print(re.finditer(x.find_pos(motif), full_dict.values(), re.I))

    #motif_list_2.append(find_pos(motif))
    #motif.find_pos(motif)
    #motif_list_2.append(motif.yell)
#print(motif_list_2)

# for x in motif_list_2:
#     x = x.motif_name
#     print(x)
    #x.find_pos()

# for x in motif_list_2:
#     print(x.motif_name)

#maybe combine Motif and Gene?


# location_dict = {} #sequence header : intron and exon information
# for header, seq in full_dict.items():
#     location_dict[header] = Gene()
# print(location_dict)

# for gene in location_dict.values():
#     print(gene.gene_length)

# class Draw:
#     def __init__(self):
#         """Inputs information from previous classes and creates the graph"""

# lil_seqs = ["cgGTGTTCCc", "cactgCATCCAAGGTca", "tcctcagcFCAAGGTcacttgtgct"]

#surface = cairo.ImageSurface(cairo.FORMAT_RGB24, 500, 500)
# ctx = cairo.Context(surface)

# ctx.rectangle(25, 50, 50, 120)
# ctx.set_source_rgb(1, 0, 0)
# ctx.fill()

# surface.write_to_png('example.png')

# ctx.rectangle(125, 50, 50, 120)
# ctx.set_source_rgb(0, 1, 1)
# ctx.set_line_width(4)
# ctx.stroke()

# ctx.rectangle(225, 50, 50, 120)
# ctx.set_source_rgb(0, 0, 1)
# ctx.fill_preserve()
# ctx.set_source_rgb(1, 1, 0)
# ctx.set_line_width(4)
# ctx.stroke()