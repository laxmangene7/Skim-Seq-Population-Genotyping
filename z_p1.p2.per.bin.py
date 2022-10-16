from csv import reader
import pandas as pd
import csv
import os


#Folder path containing input files to iterate through
path = "/bulk/laxman7/Monococcum-assembly-line/Variant-calling-wgs/laxman-trim-parent-wgs-file/combined-nex0021.nex0025.validated.aegiloref.snpcall.parentsnpposition/Merge-vcf/file_separate2_code/Separate_col/non-missing/"

#Specify the suffix that will be used to look through files with. for instance, just '.txt' will grab all .txt files.
input_file_suffix = "NA.no.H.txt"

#Windows bin file containaing chromosome, the beginning and end of the bin.
#NOTE: If you are using a .bed file, make sure that the delimiter is \t instead of ,
#"Chr1A","0","1000000"
#"Chr1A","1000000","2000000"
#"Chr1A","2000000","3000000"
#..., ..., ...

windows_file = "new.bin.csv"
windows_delimiter = ','

#The output file is the same as the input file, with two added columns: Reference Allele count and Alternative Allele count. Headerless.

####################################################################################################################################################

def determine_ref_alt_per_bin(windows_file,input_file,output_file):

    with open(input_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        s_chrom = []
        s_pos = []
        s_allele = []
        #If there is a header in the input file, uncomment this
        #header = next(csv_reader)
        for lines in csv_reader:
            s_chrom.append(lines[0])
            s_pos.append(int(float(lines[1])))
            s_allele.append(lines[2])

    sequence_counter = 0
    with open(windows_file, 'r') as csvinput:
        with open(output_file, 'w') as csvoutput:
            writer = csv.writer(csvoutput, lineterminator='\n', delimiter='\t')
            reader = csv.reader(csvinput, delimiter=windows_delimiter)

            all=[]
            
            # Iterate over each row after the header in the csv
            rows = list(reader)
            #for row in reader:
            for row in rows:
                # row variable is a list that represents a row in big parents csv
                allele_array = []
                if (sequence_counter != len(s_chrom)):
                    for n in range(sequence_counter,len(s_chrom)):
                        #print("begin row:" + str(n))
                        if (s_chrom[n] == row[0]) and (s_pos[n]) in range(int(float(row[1])),int(float(row[2]))):
                            #USED TO BE REF
                            if (s_allele[n] == "P1"):
                                allele_array.append(0)
                            #USED TO BE ALT
                            elif (s_allele[n] == "P2"):
                                allele_array.append(1)
                            else:
                                print("ERROR: Found positional match, but not Ref or Alt allele:" + s_allele[n])
                            #Found a sequence match! moving onto the next one          
                            sequence_counter = sequence_counter + 1
                        else:
                            #print("Sequence not found in bin, averaging depths and moving onto next bin.")
                            break
                    if len(allele_array) == 0:
                        #print("NOTE: NO sequences in bin, moving to next bin")
                        row.append(0)
                        row.append(0)
                        all.append(row)
                    else:
                        #print("Calculating counts of each allele in bin")
                        #Because alt alleles are stored as a 1, and reference alleles stored as a 0.
                        #So, summing our allele array will give us our alt count, and the length of our array minus that is the reference count.
                        alt_allele_count = sum(allele_array)
                        ref_allele_count = len(allele_array) - alt_allele_count
                        row.append(ref_allele_count)
                        row.append(alt_allele_count)
                        all.append(row)
                else:
                    #print("END OF SEQUENCES AT BIN: ")
                    break

            #All rows are being written at once, but I do not believe that this is incredibly memomry heavy.
            writer.writerows(all)
    #print("end")

os.chdir(path)
print("PATH:")
print(path)
# iterate through all of file
for file in os.listdir():
    # Check whether file is in text format or not, and create the input path and output path (add majorAllele.txt at end of file name)
    if file.endswith(input_file_suffix):
        print(file)
        in_file_path = f"{path}{file}"
        out_file_path = (in_file_path[:-3] + "allelesCounted.txt")
        determine_ref_alt_per_bin(windows_file,in_file_path,out_file_path)
