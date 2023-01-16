#!/usr/bin/env python
# coding: utf-8
import os, subprocess
from subprocess import PIPE
import pandas as pd
import argparse
import matplotlib.pyplot as plt

def replace_char_at_index(old_str, index, replacement):
        new_str=old_str
        if index < len(old_str):
            new_str=old_str[0:index]+replacement+old_str[index+1:]
        return new_str

def listToString(s):
        str1=' '
        return(str1.join(s))

def tRNAmap(modified_data, unmodified_data, reference_fasta_file, position_analyze, out_folder):
    print('Input file (modified): '+modified_data)
    print('Input file (unmodified): '+unmodified_data)  
    print('Analyzing position: '+str(position_analyze))

    proc=subprocess.run('conda activate tRNAmap', shell=True, stdout=PIPE, stderr=PIPE)

    reference_fasta = pd.read_csv(reference_fasta_file, header=None)
    print(reference_fasta)

    seq_name_list=[]
    for i in range(0, len(reference_fasta), 2):
        seq_name=reference_fasta.iloc[i, 0]
        seq_name=seq_name.replace('>', '')
        seq_name_list.append(seq_name)
        
    seq_list=[]
    for i in range(1, len(reference_fasta), 2):
        RNA_sequence=reference_fasta.iloc[i, 0]
        seq_list.append(RNA_sequence)

    Combined_data=[]

    for sequence, name in zip(seq_list, seq_name_list):
        AtoT=replace_char_at_index(sequence, position_analyze, 'T')
        AtoG=replace_char_at_index(sequence, position_analyze, 'G')
        AtoC=replace_char_at_index(sequence, position_analyze, 'C')
        print('target_name: '+name)
        print('target_sequence (WT): '+sequence)
        print('target_sequence_A'+str(position_analyze)+'T: '+AtoT)
        print('target_sequence_A'+str(position_analyze)+'G: '+AtoG)
        print('target_sequence_A'+str(position_analyze)+'C: '+AtoC+'\n')
        data=[]
        data=[name]
        data.append(sequence)
        analyzing_list=[]
        analyzing_list=[sequence, AtoT, AtoG, AtoC]
        analyzing_name=[]
        analyzing_name=['WT', 'AtoT', 'AtoG', 'AtoC']
        for var, name in zip(analyzing_list, analyzing_name):
            #for Modified
            argument='seqkit grep -s -i -p {0} '.format(var)+modified_data+' | grep @NB -c'
            print('argument_modified('+name+'): '+argument)
            proc=subprocess.run(argument, shell=True, stdout=PIPE, stderr=PIPE)
            output_modified=proc.stdout.decode('utf8')
            output_modified.rstrip('\n')
            print('Read number_modified('+name+'):'+output_modified)
            data.append(int(output_modified))
            print()
            #for Unmodified
            argument='seqkit grep -s -i -p {0} '.format(var)+unmodified_data+' | grep @NB -c'
            print('argument_unmodified('+name+'): '+argument)
            proc=subprocess.run(argument, shell=True, stdout=PIPE, stderr=PIPE)
            output_unmodified=proc.stdout.decode('utf8')
            output_unmodified.rstrip('\n')
            print('Read number_unmodified('+name+'):'+output_unmodified)
            data.append(int(output_unmodified))
            print()
        try:
            total_read_modified=int(data[2])+int(data[4])+int(data[6])+int(data[8])
            Mrate_modified=(int(total_read_modified)-int(data[2]))/total_read_modified
        except ZeroDivisionError:
            Mrate_modified=0
        try:
            total_read_unmodified=int(data[3])+int(data[5])+int(data[7])+int(data[9])
            Mrate_unmodified=(int(total_read_unmodified)-int(data[3]))/total_read_unmodified
        except ZeroDivisionError:
            Mrate_unmodified=0
        Mrate=Mrate_modified-Mrate_unmodified
        print('total read (modified): ' + str(total_read_modified))
        print('Mrate (modified): ' + str(Mrate_modified))
        print('total read (unmodified): ' + str(total_read_unmodified))
        print('Mrate (unmodified): ' + str(Mrate_unmodified))
        print('Mrate: ' + str(Mrate))
        data.append(total_read_modified)
        data.append(Mrate_modified)
        data.append(total_read_unmodified)
        data.append(Mrate_unmodified)
        data.append(Mrate)
        Combined_data.append(data)
        print('-------------------Done!-------------------\n')  
        

            
    Compiled_df=pd.DataFrame(Combined_data, columns=["sequence_name", "sequence", "Read_modified_WT", "Read_unmodified_WT",
                                                    "Read_modified_AtoT","Read_unmodified_AtoT", "Read_modified_AtoG","Read_unmodified_AtoG",
                                                    "Read_modified_AtoC","Read_unmodified_AtoC", "total_read_modified",
                                                    "Mrate_modified", "total_read_unmodified", "Mrate_unmodified", "Mrate"])
    print(Compiled_df)        

    #create a new directry
    dirname = out_folder+'/'
    os.makedirs(dirname, exist_ok=True)
    csv_name = dirname+'tRNAmap_data.csv'
    Compiled_df.to_csv(csv_name, index = False, header=True)

    #graph
    fig_size = len(reference_fasta)/14
    png_name = dirname+'tRNAmap_fig.png'
    print("Drawing figure...\n")
    x = Compiled_df['sequence_name']
    y = Compiled_df['Mrate']
    f = plt.figure()
    f.set_figwidth(fig_size)
    plt.bar(x, y, color="#4f85a6", alpha=1)
    plt.xlabel('tRNA variants')
    plt.xticks(rotation=90)
    plt.ylabel('Mutation rate')
    plt.savefig(png_name, dpi = 1000)
    #plt.show() 
    plt.clf() 
    print('-------------------Done!-------------------\n')  

def main():
    parser = argparse.ArgumentParser(description = "This script performs mutational profiling at a specified m1A position. The script counts mutations using Seqkit (https://bioinf.shenwei.me/seqkit/usage/). Example code: tRNAmap -m Data/tRNAmap_plus_trim_filt.fastq.gz -u Data/tRNAmap_minus_trim_filt.fastq.gz -r Data/reference_sequence.fa -p 21 -o Output. The example code is for G. stearothermophilus m1A22 methyltransferase TrmK which modifies A22 in G. stearothermophilus tRNALeu.") 

    parser.add_argument('-m', '--modified_data', type=str, required=True, help='Specify a path to modified data.')
    parser.add_argument('-u', '--unmodified_data', type=str, required=True, help='Specify a path to unmodified data.')
    parser.add_argument('-r', '--reference_fasta', type=str, required=True, help='Specify a path to a reference fasta file.')
    parser.add_argument('-p', '--position', type=int, required=True, help='Specify an m1A position that you want to analyze.')
    parser.add_argument('-o', '--output_folder', type=str, required=True, help='Specify an output folder name.')

    args = parser.parse_args()

    tRNAmap(args.modified_data, args.unmodified_data, args.reference_fasta, args.position, args.output_folder)

if __name__ == '__main__':
        main()