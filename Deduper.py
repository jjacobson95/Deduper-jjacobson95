#!/usr/bin/env python3.8
#my python3.9 disapeared. I think this will be okay.

import argparse
from os import close, write
import re
import sys

parser = argparse.ArgumentParser(description = "Argparse")
parser.add_argument('-f', '--file', help="Absolute file path to SAM file [required]", required =False)
parser.add_argument('-paired', help="Use 'True' if file is pair-ended [Unavailable]", required =False)
parser.add_argument('-umi', help="Absolute path to UMI file [required]", required =False)

single_or_paired = False
help_message = False
umi_file = False
sam_input = False

args = parser.parse_args()
sam_input = args.file
single_or_paired = args.paired
umi_file = args.umi

print("TTTTTTTTT", umi_file)

#conditional messages
if single_or_paired != False and single_or_paired != "False":
    print ("\nError! Unable to comprehend pair-ended data.")
    sys.exit("Please use single-end data or indicate '-p False'.")
if sam_input == False:
    sys.exit("\nError! Sam file is required!")

#messages
print("\nDeduper.py is used to de-duplicate identical reads cause by PCR amplification during library prep. \nOptions: \n -f, --file: Input absolute path to the input sam file. \n -p, --paired: Input <True> if file is pair-ended. Note - this script is not yet capable of using pair-end data. \n -u, --umi: Input absolute path to UMI file. Currently this is a required parameter. \n -h, --help: Display argument details.")
print("\nDeduper.py has begun.")

#create umi set
if umi_file == False:
    UMI_present = False
else:
    UMI_present = True
    umi_set = set()
    with open(umi_file, "r") as file_umi:
        for line in file_umi:
            umi_line = line.strip()
            umi_set.add(umi_line)
    print("\nUMIs identified.")


#Function 1: Strand Checker - check which strand the read is aligned to
def strand_check(bit_flag):
    """
    This function checks which strand the read is on from the input of the bitwise flag.
    It returns forward_strand = Bool. It also returns a toss_read function if unable to determine strand.
    """
#this is not necessary for a file with only two options of a bitwise flag. If 16th bit = 0 it is forward, if 1 it is reverse. (double check this)
    #bitwise = '{0:b}'.format(int(bit_flag))
    #bits = [int(x) for x in str(bitwise)]

    #16 indicates that this is on the reverse strand
    if int(bit_flag) == 16:
        forward_strand = False
    #0 indicates that this is on the forward strand
    elif int(bit_flag) == 0:
        forward_strand = True
    else:
        print(bit_flag)
    return(forward_strand)


#Function 2: Clipped Position Start Finder - Positive Strand
def FS_pos_finder(CIGAR, input_pos):
    """
    This function finds the starting postion of a read for the Positive/Forward strand.
    """
    cigar_dict = {}
    bases_skipped = 0

#this splits cigar strings by number and letter - aka: '2S10M3S' -> ['2S','10M','3S']
    cigs = re.split(r"([0-9]+[A-Z,a-z]+)", CIGAR)
    cigs = [string for string in cigs if string != ""]

#This adds all values of identical letters/operators. The above string would now be in dict: {M: 10, S: 5}
    for item in cigs:
        if item[-1] in cigar_dict:
            cigar_dict[item[-1]] += int(item[:-1])
        else:
            cigar_dict[item[-1]] = int(item[:-1])

#This adds softclips to dictionary. The above dict would now be: {M: 10, S: 5, left_S: 2, right_S: 3}
    if cigs[0][-1] ==  "S":
        cigar_dict["left_S"] = int(cigs[0][:-1])
    if cigs[-1][-1] == "S":
        cigar_dict["right_S"] = int(cigs[-1][:-1])
        
    tru_position = int(input_pos)

#this finds the true starting position of the positive/forward strand
#this also includes bases_skipped - this is used for choosing the best duplicate
    for key in cigar_dict:
        if key == "left_S":
            tru_position -= cigar_dict[key]
            bases_skipped += cigar_dict[key]
        elif key == "D":
            bases_skipped += cigar_dict[key]
            continue
        elif key == "N":
            bases_skipped += cigar_dict[key]
            continue
        elif key == "right_S":
            bases_skipped += cigar_dict[key]
            continue
    return(tru_position, bases_skipped)


#Function 3: Clipped Position Start Finder - Negative Strand
def RS_pos_finder(CIGAR, input_pos):
    """
    This function finds the starting postion of a read for the Negative/Reverse strand.
    """
    cigar_dict = {}
    bases_skipped = 0

#this splits cigar strings by number and letter - aka: '2S10M3S' -> ['2S','10M','3S']
    cigs = re.split(r"([0-9]+[A-Z,a-z]+)", CIGAR)
    cigs = [string for string in cigs if string != ""]


#This adds all values of identical letters/operators. The above string would now be in dict: {M: 10, S: 5}
    for item in cigs:
        if item[-1] in cigar_dict:
            cigar_dict[item[-1]] += int(item[:-1])
        else:
            cigar_dict[item[-1]] = int(item[:-1])

#This adds softclips to dictionary. The above dict would now be: {M: 10, S: 5, left_S: 2, right_S: 3}
    if cigs[-1][-1] == "S":
        cigar_dict["right_S"] = int(cigs[-1][:-1])
    if cigs[0][-1] == "S":
        cigar_dict["left_S"] = int(cigs[0][:-1])

        
    tru_position = int(input_pos)

#this finds the true starting position of the negative/reverse strand
#this also includes bases_skipped - this is used for choosing the best duplicate
    for key in cigar_dict:
        if key == "M":
            tru_position += cigar_dict[key]
        elif key == "I":
            continue
        elif key == "D":
            tru_position += cigar_dict[key]
            bases_skipped += cigar_dict[key]
        elif key == "N":
            tru_position += cigar_dict[key]
            bases_skipped += cigar_dict[key]
        elif key == "right_S":
            tru_position += cigar_dict[key]
            bases_skipped += cigar_dict[key]
        elif key == "left_S":
            bases_skipped += cigar_dict[key]
            continue
        elif key == "H":
            continue
        elif key == "S":
            continue
        elif key == "P":
            continue
        elif key == "X":
            tru_position += cigar_dict[key]
        elif key == "=":
            tru_position += cigar_dict[key]
        else:
            print("Improper Cigar Operator: ", key, ". Found in Cigar: ", CIGAR)
    return(tru_position, bases_skipped)



#these are counters for output
duplicate_counter = 0
wrong_umi = 0
header_counter = 0
total_reads = 0

#this will hold all the reads that we want to keep.
PCR_reads_dict = {}

#this will count all output reads per chromosome
chrom_counter_dict = {}

#output file name
output_file_name = (str(sam_input) + "_deduped")
ff = open(output_file_name, "w")

#these are used for update messages
head_check = 0
deduplicate_check = 0
update_check = 0


with open(sam_input, "r") as file_sam:
    for line in file_sam:
        #first write out headers to new file
        component = line.split()
        if ("@") in component[0][0]:
            ff.write(line)
            header_counter +=1
            continue

        head_check +=1
        if head_check == 1:
            print("Headers identified.")

        #Check if UMIs are correct
        the_umi = component[0][-8:]
        if the_umi not in umi_set and UMI_present == True:
            wrong_umi +=1
            continue
       
       #Assign values needed for determining if PCR duplicate
        chrom = component[2]
        input_pos = component[3]
        CIGAR = component[5]
        bit_flag = component[1]

        #check the strand with the strand_check function
        forward_strand = strand_check(bit_flag)

        #determine the true starting position
        #also determine the number of deletions/clips for choosing best duplicate
        if forward_strand == True:
            real_position = FS_pos_finder(CIGAR,input_pos)[0]
            skipped_bases = FS_pos_finder(CIGAR,input_pos)[1]
        elif forward_strand == False:
            real_position = RS_pos_finder(CIGAR,input_pos)[0]
            skipped_bases = RS_pos_finder(CIGAR,input_pos)[1]

        #unique title for PCR_reads_dict - each unique PCR read will have a unique title. Duplicates will overwrite eachother.
        title = ("chromosome: " + str(chrom) + " Position: " + str(real_position) + " UMI: " + str(the_umi) + " Forward Strand: " + str(forward_strand))

        #This read is a duplicate
        if title in PCR_reads_dict:
            #this chooses the best duplicate based off having the fewest number of deletions, Ns,and soft-clips.
            if skipped_bases < previous_skipped_bases:
                PCR_reads_dict[title] = line
                previous_skipped_bases = skipped_bases

            #count duplicates and print out some progress updates.
            duplicate_counter +=1
            if duplicate_counter == 10000:
                print("\nDeduplicator in progress.\n    10,000 duplicates removed.")
            if duplicate_counter == 100000:
                print("   100,000 duplicates removed.")
            if duplicate_counter == 1000000:
                print(" 1,000,000 duplicates removed.")
            if duplicate_counter == 2500000:
                print(" 2,500,000 duplicates removed.")
            if duplicate_counter == 5000000:
                print(" 5,000,000 duplicates removed.")
            if duplicate_counter == 10000000:
                print(" 10,000,000 duplicates removed.")

        #This is a unique read (or the first of a duplicate).
        else:
            PCR_reads_dict[title] = line
            #counters
            if chrom not in chrom_counter_dict:
                chrom_counter_dict[chrom] = 0
            chrom_counter_dict[chrom] += 1
            total_reads += 1
            #this assists the best duplicate feature
            previous_skipped_bases = skipped_bases


#Print all chromosome read counts
print("Writing output to SAM file.")
for v in PCR_reads_dict.values():
    ff.write(v)

#Result messages
print("\nSuccess.\n\nResults:")

print("chromosome dictionary: \n", chrom_counter_dict)
print("\ntotal duplicates removed: ", duplicate_counter)
print("incorrect UMIs found: ", wrong_umi)
print("Header counts: ", header_counter)
print("Total reads after deduplication removal: ", total_reads)
print("Output file is saved as: ", output_file_name, "\n")


ff.close()
 