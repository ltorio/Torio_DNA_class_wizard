#!/usr/bin/env python
# coding: utf-8

# # A class called seq 
# ## This class should accept the following attributes:
# ### name, organism, sequence, and type. When making an isntance of seq, the input for name will be a string (e.g. 'RAS_G12D') the input for organism will be a string (e.x. 'human'), the input for sequence will be a DNA, RNA, or protein sequence (e.g. 'ATCGAAATC') and the type will be either 'DNA', 'RNA', or 'Protein'
# 
# ## This class should have three methods: 
# ### 1. info -this should print the name, type, organism, and sequence of the instance
# ### 2. length -this should count the length of the sequence string
# ### 3. fasta_out -this should write the name, organism, type, and sequence as a fasta file. 

class seq:
    # Call instance attributes
    def __init__(self, name, organism, sequence, type):
        self.name=name
        self.organism=organism
        self.sequence=sequence
        self.type=type
     
    # define the info function 
    def info(self):
        print(self.name)
        print(self.type)
        print(self.organism)
        print(self.sequence)
    
    # define the length function
    def length(self):
        print(len(self.sequence))
    
    # the fasta_out method creats a new fa file
    # it writes to a file using the sequence name as part of the file name
    # the first line of the file includes the name of the sequence, and the organism,
    # the second line of the file is the sequence
    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(">" + self.name + "_" + self.organism + "_" + self.type + "\n" + self.sequence)
        f.close()
       


# # A new child class of seq called protein
# ## This should have a new attribute called size (For instances of protein, size values will be in kDa, like '52')
# ## Overwrite the parent class seq function fasta_out to include the protein size in the first line of the fasta file

class protein(seq):
    def __init__(self, name, organism, sequence, type, size):
        self.size=size # new attribute size
        super().__init__(name, organism, sequence, type) # inherit parent methods info() and length()
    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(">" + self.name + "_" + self.organism + "_" + self.type + "_"+str(self.size)+"\n" + self.sequence)
        f.close()
    
    #eturns the total molecular weight of the protein sequence
    def mol_weight(self):
        # a dictionary of amino acids and their molecular weight copied from helpful_variables.txt
        aa_mol_weights={'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
                'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
                'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
                'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19}
        weight = 0.0 # will store the total molecular weight of the sequence
        
        # loops thru each amino acid in the sequence and sums the molecular weight of each amino acid
        for aa in self.sequence:
            weight=weight+aa_mol_weights[aa]
        print(weight)
            
        


# # A child class of seq called nucleotide
# ## Make a new method called gc_content that calculates the percent of letters that are G or C and then prints the gc content percentage
# 
# ## for the gc_content method, there are multiple ways to do it, but you may find this page helpful https://www.geeksforgeeks.org/python-count-occurrences-of-a-character-in-string/ 

# In[46]:


# Write the new nucleotide class here
class nucleotide(seq):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type) # fully inherits parent's methods
        
    # calculates the percent of letters that are G or C and then prints the gc content percentage
    def gc_content(self):
        totalContent = len(self.sequence)
        gcContent = self.sequence.count("G")+self.sequence.count("C")
        gcPercent = gcContent/totalContent*100
        print(gcPercent)


# DNA is a child class of nucleotide that is used for DNA sequences
# a new methods: transcribe, six_frames, reverse_complement
class DNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type) # fully inherits parent's __init__
    
    # the transcribe method transcribes DNA to RNA and prints the transcribed sequence 
    # (aka replace the Ts in the DNA sequence with Us)
    def transcribe(self):
        rna = self.sequence.replace("T","U")
        print(rna)
        return(rna)
    
    # the six_frames method prints all 6 coding frames of the sequence 
    # (3 frames on the forward strand, 3 frames on the reverse complement strand)
    def six_frames(self):
        forwardFrames=["","",""] # will store the 3 frames for the forward strand
        reverseFrames=["","",""] # will store the 3 frames for the reverse complement strand
        reverseStrand = self.sequence[::-1] # the last 5 bases of the forward strand in reverse order
        reverseCompStart="" # will store the first 5 bases of the reverse complement strand
        
        forwardFrames[0]=self.sequence[0:] # first reading frame on the forward strand
        forwardFrames[1]=self.sequence[1:] # second reading frame on the forward strand
        forwardFrames[2]=self.sequence[2:] # third reading frame on the forward strand
        
        # loops thru a strand and creates a complement strand
        for letter in reverseStrand:
            if letter == "A":
                reverseCompStart=reverseCompStart+"T"
            elif letter == "C":
                reverseCompStart=reverseCompStart+"G"
            elif letter =="T":
                reverseCompStart=reverseCompStart+"A"
            elif letter == "G":
                reverseCompStart=reverseCompStart+"C"

        reverseFrames[0]=reverseCompStart[0:] # first reading frame on the reverse complement strand
        reverseFrames[1]=reverseCompStart[1:] # second reading frame on the reverse complement strand
        reverseFrames[2]=reverseCompStart[2:] # third reading frame on the reverse complement strand
        
        print(forwardFrames[0])
        print(forwardFrames[1])
        print(forwardFrames[2])
        print(reverseFrames[0])
        print(reverseFrames[1])
        print(reverseFrames[2])
        
    # returns the reverse complement of the sequence
    def reverse_complement(self):
        reverseCompliment="" # will store the reverse compliment strand
        reverse=self.sequence[::-1] # a reversed string from the sequence
        
        # loops thru a strand and creates a complement strand
        for letter in reverse:
            if letter == "A":
                reverseCompliment=reverseCompliment+"T"
            elif letter == "C":
                reverseCompliment=reverseCompliment+"G"
            elif letter =="T":
                reverseCompliment=reverseCompliment+"A"
            elif letter == "G":
                reverseCompliment=reverseCompliment+"C"
        print(reverseCompliment)
        return(reverseCompliment)
        
        
        


# RNA is a child class of nucleotide that is used for RNA sequences
# a new methods: start, translate

class RNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type) # inheritance

    # returns and prints the index of the start sequence AUG
    def start(self):
        startIndex = self.sequence.find("AUG")
        print(startIndex)
        return(startIndex)
    
    # returns and prints the sequence after translation    
    def translate(self):
        # a dictionary for RNA to amino acid translations from helpful_variables.txt
        standard_code = {
         "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
         "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
         "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
         "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
         "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
         "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
         "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
         "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
         "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
         "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
         "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
         "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
         "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
        
        codonList=[] # will store codons from an rna sequence
        aaSeq = "" # will store the translated amino acid sequence
        startCodonIndex = self.start() # finds the start codon (AUG) position
        codingRegion=self.sequence[startCodonIndex:] # the region of rna to be translated
        
        # loops thru a seqence breaks the sequence into 3 letter codons 
        for i in range(0,len(codingRegion),3):
            codonList.append(codingRegion[i:i+3])
        if len(codonList[-1])!=3: # handles cases where when the last codon in the list is not 3 bases long and thus invalid
            codonList.pop(-1)
        
        # loops thru the codonList and creates the translated amino acid sequence
        for codon in codonList:
            aaSeq=aaSeq+standard_code[codon]
        print(aaSeq)
        return(aaSeq)
        
        


# # Part B: test out your new methods

# Q8 assign the sequence: "CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA" to a variable of class DNA called uidA. This variable attributes should be name uidA, organism bacteria, and type DNA.
uidA = DNA(name = 'uidA', sequence= 'CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA', 
                      organism = 'bacteria', type = 'DNA')

# Q9 Use the fasta_out function to write the sequence and information for uidA DNA to a fasta file
uidA.fasta_out()

# Q10 10. Use the six_frames function and reverse_complement functions to output the six coding frames and reverse complement of the uidA DNA sequence
uidA.six_frames()
uidA.reverse_complement()

# Q11 Using the transcribe function for the DNA class, transcribe the uidA DNA sequence to an RNA sequence.
# Q12 Save this RNA sequence as a RNA class object called uidA_RNA with the same other attributes except the name should be uidA_RNA and the type should be RNA
uidA_RNA = RNA(name = 'uidA_RNA', sequence= uidA.transcribe(), 
                      organism = 'bacteria', type = 'RNA')

# Q13 Use the fasta_out() function to write the RNA sequence and information for uidA to a fasta file
uidA_RNA.fasta_out()

# Q14 Use the translate method on the RNA object uidA_RNA to translate the RNA sequence
# Q15 Save this amino acid sequence as a protein class object called uidA_protein. Set the name as uidA_protein and the type as protein. You can set the size attribute as any value. Use the fasta_out() function to write this protein sequence and information to a new fasta file
uidA_protein = protein(name = 'uidA_protein', sequence= uidA_RNA.translate(), 
                      organism = 'bacteria', type = 'protein', size = 9999)
uidA_protein.fasta_out()

# Q16 Use the method mol_weight to output the molecular weight of the amino acid sequence uidA_protein
uidA_protein.mol_weight()

