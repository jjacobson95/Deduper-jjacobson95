#!/usr/bin/env python3.8

import argparse
from os import close, write
import re


def argparse_options():
    parser = argparse.ArgumentParser(description = "Argparse")
    parser.add_argument ("--file", help="Absolute file path", required =False, type = str)
    parser.add_argument ("--paired", help="Use 'True' if file is pair-ended", required =False, type = str)
    parser.add_argument ("--umi", help="Absolute path to UMI file, if exists.", required =False, type = str)
    parser.add_argument ("--help", help="Use 'True' to print a useful help message", required =False, type = str)
    return parser.parse_args()



single_or_paired = False
help_message = False
umi_file = False
sam_input = False

args = argparse_options()
sam_input = args.file
single_or_paired = args.paired
umi_file = args.umi
help_message = args.help


if help_message != False:
    print("Deduper.py is used to de-duplicate identical reads cause by PCR amplification during library prep. \n Options: \n --file: Input absolute path to the input sam file. \n --paired: Input <True> if file is pair-ended. Note - this script is not yet capable of using pair-end data. \n --umi: Input absolute path to UMI file. Currently this is a required parameter. /n --help: Display this message.")

if single_or_paired != False:
    print ("Error! Unable to comprehend pair-ended data.")

if umi_file == False:
    print ("Error! Deduper.py currently requires a Umi file.")

if sam_input == False:
    print ("Error! Sam file is required!")


#create umi set
umi_set = set()
with open(umi_file, "r") as file_umi:
    for line in file_umi:
        umi_set.add(line)

print(umi_set)


#Function 1: Strand Checker - check which strand the read is aligned to
def strand_check(bit_flag):
    """
    This function checks which strand the read is on from the input of the bitwise flag.
    It returns forward_strand = Bool. It also returns a toss_read function if unable to determine strand.
    """
    bitwise = '{0:b}'.format(bit_flag)
    bits = [int(x) for x in str(bitwise)]
    if len(bits) >=5:
        if bits[-5] == 0:
            forward_strand = True
        if bits[-5] == 1:
            forward_strand = False
    else:
        forward_strand = True
    return(forward_strand)


#Function 2: Clipped Position Start Finder - Positive Strand
def FS_pos_finder(CIGAR, input_pos):
    rightmost_softclip = False
    leftmost_softclip = False
    cigar_dict = {}
    cigs = re.split(r"([0-9]+[A-Z,a-z]+)", CIGAR)
    cigs = [string for string in cigs if string != ""]
    
    for item in cigs:
        if item[-1] in cigar_dict:
            cigar_dict[item[-1]] = int(cigar_dict[item[-1]]) + int(item[:-1])
        else:
            cigar_dict[item[-1]] = int(item[:-1])

    if cigs[-1][-1] ==  "S":
        cigar_dict["right_S"] = cigs[-1][:-1]
        rightmost_softclip = True
    if cigs[0][-1] ==  "S":
        cigar_dict["left_S"] = cigs[0][:-1]
        leftmost_softclip = True
        
    tru_position = input_pos

    for key in cigar_dict:
        if key == "left_s" and leftmost_softclip == True:
            tru_position - cigar_dict[key]
    return(tru_position)

#Function 3: Clipped Position Start Finder - Negative Strand
def RS_pos_finder(CIGAR, input_pos):
    rightmost_softclip = False
    leftmost_softclip = False
    cigar_dict = {}
    cigs = re.split(r"([0-9]+[A-Z,a-z]+)", CIGAR)
    cigs = [string for string in cigs if string != ""]

    for item in cigs:
        if item[-1] in cigar_dict:
            cigar_dict[item[-1]] = int(cigar_dict[item[-1]]) + int(item[:-1])
        else:
            cigar_dict[item[-1]] = int(item[:-1])

    if cigs[-1][-1] == "S":
        cigar_dict["right_S"] = cigs[-1][:-1]
        rightmost_softclip = True
    if cigs[0][-1] == "S":
        cigar_dict["left_S"] = cigs[0][:-1]
        leftmost_softclip = True
        
    tru_position = input_pos 

    for key in cigar_dict:
        if key == "M":
            #add the Ms
            tru_position + cigar_dict[key]
        elif key == "I":
            #ignore the Insertions
            continue
        elif key == "D":
            #add the Deletions
            tru_position + cigar_dict[key]
        elif key == "N":
            #add the Ns (skipped regions)
            tru_position + cigar_dict[key]
        elif key == "right_S" and rightmost_softclip == True:
            #add the rightmost softclip
            tru_position + cigar_dict[key]
        elif key == "left_S" and leftmost_softclip == True:
            #ignore the leftmost softclip
            continue
        elif key == "H":
            #ignore hard clipping
            continue
        elif key == "P":
            #ignore padding
            continue
        elif key == "X":
            #add the Xs (seq match)
            tru_position + cigar_dict[key]
        elif key == "=":
            #add the '='s (seq mismatch)
            tru_position + cigar_dict[key]
        else:
            print("Improper Cigar Operator: ", key, ". Found in Cigar: ", CIGAR)
    return(tru_position)




#Function 4: Best_Duplicate - choose the best duplicate based on a scoring system I devise.
#This function will only activate if strand, adjusted position, chrom, UMI are all identical.
#Scoring system -> best has the most Ms. If same amount of Ms, check average quality. 
#If same quality, choose the first one.
#Then compare this to the next one. and so on. 


#Some logic - first I will sort the file. Then I will run through it. Lines will only write to file if 
#position no longer matches calulated fixed postion of the previous run
#I will also maintain a dictionary(?) of previously found Chrom, Pos, Umi, Strand.



wrong_umi = 0
PCR_reads_dict = {}

ff = open("Deduped_file.Sam", "w")

with open(sam_input, "r") as file_sam:
    for line in file_sam:
        #first write out headers to new file
        if "@SQ" in line:
            ff.write(line)
            continue

        #then do all this stuff
        component = line.split()
        the_umi = component[0][-8:]
        if the_umi not in umi_set:
            wrong_umi +=1
        chrom = component[2]
        input_pos = component[3]
        CIGAR = component[5]
        bit_flag = component[1]
        forward_strand = strand_check(bit_flag)

        if forward_strand == True:
            real_position = FS_pos_finder(CIGAR,input_pos)
        elif forward_strand == False:
            real_position = RS_pos_finder(CIGAR,input_pos)

        title = ("chromosome: " + str(chrom) + " Position: " + str(real_position) + " UMI: " + str(the_umi) + " Forward Strand: " + str(forward_strand))
        
        if title in PCR_reads_dict:
            #this is a duplicate
            #run function best_duplicate
            continue
        else:
            PCR_reads_dict[title] = line


for v in PCR_reads_dict.values():
    ff.write(v)


#close output file
 






    #first write header to sam_output
    #then do the real stuff

