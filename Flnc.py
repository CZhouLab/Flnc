#!/usr/bin/env python
import copy
import argparse
import textwrap
import os
import re
import subprocess
import logging
import datetime
import json
import Mode1
import Mode2
import Share
import logging

def process_bed(input_file, library, output_dir, model, gtf_file, user_json_file):
        LOG_FILENAME = output_dir + "/logfile.log"
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO)
        logging.info("-----------------------------------------")
        logging.info(str(datetime.datetime.now()))
        logging.info("Input file (bed format): " + input_file)
        logging.info("Library: " + library)
        logging.info("Output directory: " + output_dir)
        logging.info("Machine learning model: " + model)
        logging.info("GTF file: " + gtf_file)

        os.system("awk '{print $0,$4}' FS=\"\t\" OFS=\"\t\" " + input_file + " >" + output_dir + "/input.bed")
        Mode2.TrueLncRNA_PredictFrom_PutativeLncRNABED(output_dir+"/input.bed", library, output_dir, model, user_json_file, step_index = 1)
        
        logging.info("Flnc pipeline is done!\noutput files are under " + output_dir + "/output/ ")
        logging.info("-----------------------------------------")

def process_fasta(input_file, library, output_dir, model, gtf_file, user_json_file):
        LOG_FILENAME = output_dir + "/logfile.log"
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO)
        logging.info("-----------------------------------------")
        logging.info(str(datetime.datetime.now()))
        logging.info("Input file (fasta format): " + input_file)
        logging.info("Library: " + library)
        logging.info("Output directory: " + output_dir)
        logging.info("Machine learning model: " + model)
        logging.info("GTF file: " + gtf_file)
        
        putative_lncRNA_bed = Mode2.fasta_2_bed(input_file, output_dir, library, user_json_file)
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 1: fasta to bed " + TimeNow)
        Mode2.TrueLncRNA_PredictFrom_PutativeLncRNABED(putative_lncRNA_bed, library, output_dir, model, user_json_file, step_index = 2)
        
        logging.info("Flnc pipeline is done!\noutput files are under " + output_dir + "/output/ ")
        logging.info("-----------------------------------------")

def process_singleend_fastq(input_file, library, output_dir, strand, model, gtf_file, user_json_file):
        LOG_FILENAME = output_dir + "/logfile.log"
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO)
        logging.info("-----------------------------------------")
        logging.info(str(datetime.datetime.now()))
        logging.info("Input file (fastq format): " + input_file)
        logging.info("Library: " + library)
        logging.info("Output directory: " + output_dir)
        logging.info("Strand: " + strand)
        logging.info("Machine learning model: " + model)
        logging.info("GTF file: " + gtf_file)

        Replicate_list = input_file.split(",")
        
        GTF_List=[]
        for i in range(0,len(Replicate_list)):
                Rep_index = i+1
                Rep_fastq = Replicate_list[i].strip() 
                Mode1.HISAT2_Alignment_SingleEnd(str(Rep_index), Rep_fastq, strand, output_dir, library, user_json_file) 
                StringTie_GTF = Mode1.StringTie_Assembly(str(Rep_index), strand, output_dir, library, user_json_file)
                Strawberry_GTF = Mode1.Strawberry_Assembly(str(Rep_index), strand, output_dir, library, user_json_file)
                GTF_List.append(StringTie_GTF)
                GTF_List.append(Strawberry_GTF)
        GTF_List_Value = "\n".join(GTF_List)
        with open(output_dir + '/GTF_MergeList.txt', 'w') as rsh:
                rsh.write(GTF_List_Value)

        GTF_input = Mode1.StringTie_Merge(output_dir+'/GTF_MergeList.txt', output_dir, library, user_json_file)
        Share.MergeGTF_2_NoncodingLongTranscripts_BED(GTF_input, output_dir, library, user_json_file)
       
        os.system("rm -f " + output_dir + "/expressed_transcripts.id")
        for i in range(0,len(Replicate_list)):
                Rep_index = i+1
                Rep_bam = output_dir + "/Rep" + str(Rep_index) + ".sortedbycoord.bam"
                Mode1.FPKM(str(Rep_index), Rep_bam, GTF_input, strand, output_dir, library)
                Mode1.Cat_FPKM_ReadCount(str(Rep_index), GTF_input, output_dir, library)
                Expressed_Transcripts_file = Mode1.SelectExpress(str(Rep_index), output_dir, library)
                
                cmd_1 = '{Expressed_Transcripts_file}'
                cmd_1_R = cmd_1.format(Expressed_Transcripts_file = Expressed_Transcripts_file)
                f1 = open(output_dir + "/expressed_transcripts.id", "a")
                cmd_1_Run = subprocess.Popen(['cat'] + cmd_1_R.split(), stdout = f1)
                cmd_1_Run.communicate()
                f1.close()

        putative_lncRNA_bed = Mode1.PutativeLncRNA_Generation(output_dir + "/expressed_transcripts.id", output_dir, library)
        Mode1.TrueLncRNA_PredictFrom_PutativeLncRNABED(len(Replicate_list), putative_lncRNA_bed, library, output_dir, model, user_json_file, step_index = 15)
        os.system("cp " + output_dir +"/*HISAT2 " + output_dir + "/output")
        logging.info("Flnc pipeline is done!\noutput files are under " + output_dir + "/output/ ")
        logging.info("-----------------------------------------")

def process_pairend_fastq(input_R1_file, input_R2_file, library, output_dir, strand, model, gtf_file, user_json_file):
        Replicate_list_R1 = input_R1_file.split(",")
        Replicate_list_R2 = input_R2_file.split(",")
        assert len(Replicate_list_R1) == len(Replicate_list_R2)
        
        LOG_FILENAME = output_dir + "/logfile.log"
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO)
        logging.info("-----------------------------------------")
        logging.info(str(datetime.datetime.now()))
        logging.info("Input R1 file (fastq format): " + input_R1_file)
        logging.info("Input R2 file (fastq format): " + input_R2_file)
        logging.info("Library: " + library)
        logging.info("Output directory: " + output_dir)
        logging.info("Strand: " + strand)
        logging.info("Machine learning model: " + model)
        logging.info("GTF file: " + gtf_file)

        GTF_List=[]
        for i in range(0,len(Replicate_list_R1)):
                Rep_index = i+1
                Rep_fastq_R1 = Replicate_list_R1[i].strip()
                Rep_fastq_R2 = Replicate_list_R2[i].strip()
                Mode1.HISAT2_Alignment_PairedEnd(str(Rep_index), Rep_fastq_R1, Rep_fastq_R2, strand, output_dir, library, user_json_file) 
                StringTie_GTF = Mode1.StringTie_Assembly(str(Rep_index), strand, output_dir, library, user_json_file)
                Strawberry_GTF = Mode1.Strawberry_Assembly(str(Rep_index), strand, output_dir, library, user_json_file)
                GTF_List.append(StringTie_GTF)
                GTF_List.append(Strawberry_GTF)
        GTF_List_Value = "\n".join(GTF_List)
        with open(output_dir + '/GTF_MergeList.txt', 'w') as rsh:
                rsh.write(GTF_List_Value)

        GTF_input = Mode1.StringTie_Merge(output_dir+'/GTF_MergeList.txt', output_dir, library, user_json_file)
        GTF_input = output_dir + "/stringtie_merged.changeID.gtf"
        Share.MergeGTF_2_NoncodingLongTranscripts_BED(GTF_input, output_dir, library, user_json_file)
       
        os.system("rm -f " + output_dir + "/expressed_transcripts.id")
        for i in range(0,len(Replicate_list_R1)):
                Rep_index = i+1
                Rep_bam = output_dir + "/Rep" + str(Rep_index) + ".sortedbycoord.bam"
                Mode1.FPKM(str(Rep_index), Rep_bam, GTF_input, strand, output_dir, library)
                Mode1.Cat_FPKM_ReadCount(str(Rep_index), GTF_input, output_dir, library)
                Expressed_Transcripts_file = Mode1.SelectExpress(str(Rep_index), output_dir, library)
                
                cmd_1 = '{Expressed_Transcripts_file}'
                cmd_1_R = cmd_1.format(Expressed_Transcripts_file = Expressed_Transcripts_file)
                f1 = open(output_dir + "/expressed_transcripts.id", "a")
                cmd_1_Run = subprocess.Popen(['cat'] + cmd_1_R.split(), stdout = f1)
                cmd_1_Run.communicate()
                f1.close()

        putative_lncRNA_bed = Mode1.PutativeLncRNA_Generation(output_dir + "/expressed_transcripts.id", output_dir, library)
        Mode1.TrueLncRNA_PredictFrom_PutativeLncRNABED(len(Replicate_list_R1), putative_lncRNA_bed, library, output_dir, model, user_json_file, step_index = 15)
        os.system("cp " + output_dir +"/*HISAT2 " + output_dir + "/output")
        logging.info("Flnc pipeline is done!\noutput files are under " + output_dir + "/output/ ")
        logging.info("-----------------------------------------")

shared_parser = argparse.ArgumentParser(add_help=True, formatter_class=argparse.RawTextHelpFormatter)
shared_parser.add_argument('-v', '--version', action = 'version', version = 'v1.0')
shared_parser.add_argument('-g', '--gtf_file', type = str, default = '', help = textwrap.dedent('''\
Full path of the reference gene annotation file in GTF format.
Default: gencode.v29.annotation.gtf in the LIB folder.
'''))
shared_parser.add_argument('-m', '--model', type = str, default = "rf", help = textwrap.dedent('''\
Choose the abbreviation of one of the following models:
rf: random forest.
lr: logistic regression.
nb: naive Bayes.
dt: decision tree.
knn: k-nearest neighbors.
rbfsvm: support vector machines with RBF kernel.
lsvm: support vector machines with linear kernel.
ensemble: the common result predicted by all models 
Default: rf.
'''), choices = ["rf", "lr", "nb", "dt", "knn", "rbfsvm", "lsvm", "ensemble"])
shared_parser.add_argument('-s', '--strand', type = str, choices = ["first", "second", "unstrand"], default = "first", help = textwrap.dedent('''\
This option is required only if "-f fastq", otherwise this argument is not needed.
Specify strand-specific information with the following three options: 
first: corresponds to fr-firststrand of the -library-type option in the TopHat tool for stranded RNA-seq data.
second: corresponds to fr-secondstrand of the -library-type option in the TopHat tool for stranded RNA-seq data.
unstrand: specific for unstranded RNA-seq data.
Default: first.
'''))
shared_parser.add_argument('-l', '--library', type = str, required = True, help = textwrap.dedent('''\
Mandatory argument.
Full path of the LIB folder, which can be downloaded from Zenodo: https://zenodo.org/record/5711975/files/LIB.zip?download=1 
'''))
shared_parser.add_argument('-o', '--output_dir', type = str, required = True, help = textwrap.dedent('''\
Mandatory argument.
Please specify the name of the output folder. This must be specified as a full path. For example, "-o /home/username/Flnc_sample1_output".
'''))
shared_parser.add_argument('-f','--format', required = True, help = textwrap.dedent('''\
Mandatory argument.
The format of the input file: fastq, or fasta or bed.
If using the pair subcommand, the format must be "fastq".
If using single subcommand, the format can be fastq, or fasta, or bed.
'''), choices = ["fastq", "fasta", "bed"])
parser = argparse.ArgumentParser(usage="python Flnc.py [options]", formatter_class=argparse.RawTextHelpFormatter)
subparser = parser.add_subparsers(help='help for subcommand', dest='subcommand')

parser_a = subparser.add_parser("pair", add_help=False, parents=[shared_parser], formatter_class=argparse.RawTextHelpFormatter, description='When running Flnc with paired RNA-seq data, it is critical that the *_1 files and the *_2 files of replicates appear in separate comma-delimited lists, and that the order of the files in the two lists is the same')
parser_a.add_argument('-1', type = str, required = True, help = textwrap.dedent('''\
This argument is mandatory if using the pair subcommand.
Full path of the mate 1 file of paired FASTQ files, paired with the mate 2 file specified with "-2 " option.
The mate 1 of replicates can be input through comma delimitation, e.g., "<path>/Rep1_1.fastq,<path>/Rep2_1.fastq".
'''), dest='file1')
parser_a.add_argument('-2', type = str, required = True, help = textwrap.dedent('''\
This argument is mandatory if using the pair subcommand.
Full path of the mate 2 file of paired FASTQ files, paired with the mate 1 file specified with "-1 " option.
The mate 2 of replicates can be input through comma delimitation, e.g., "<path>/Rep1_2.fastq,<path>/Rep2_2.fastq".
'''), dest='file2')

parser_b = subparser.add_parser("single", add_help=False, parents=[shared_parser], formatter_class=argparse.RawTextHelpFormatter)
parser_b.add_argument("-u", required = True, help=textwrap.dedent('''
This argument is mandatory if using the single subcommand
Full path of the single input file.
If "-f fastq", please input the full path of FASTQ file of single-end RNA-seq data. FASTQ files for replicates can be input through comma delimitation, For example, "<path>/Rep1.fastq,<path>/Rep2.fastq". 
If "-f fasta", please input the full path of files with transcripts in FASTA format.
If "-f bed", please input the full path of files with transcripts in BED format.
'''), dest='file')

args = parser.parse_args()
if args.model != "":
        model = args.model
else:
        model = "rf"
        
if args.gtf_file != "":
        gtf_file = args.gtf_file
else:
        gtf_file = args.library + "/gencode.v29.annotation.gtf"

if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

user_json_file = Share.load_reference(args)

if args.format in ["fastq"]:
        assert args.strand in ["first", "second", "unstrand"]
        if args.subcommand == "single":
                process_singleend_fastq(args.file, args.library, args.output_dir, args.strand, model, gtf_file, user_json_file)
        elif args.subcommand == "pair":
                process_pairend_fastq(args.file1, args.file2, args.library, args.output_dir, args.strand, model, gtf_file, user_json_file)
elif args.format in ["bed"]:
        assert args.subcommand == "single"
        process_bed(args.file, args.library, args.output_dir, model, gtf_file, user_json_file)
        
elif args.format in ["fasta"]:
        assert args.subcommand == "single"
        process_fasta(args.file, args.library, args.output_dir, model, gtf_file, user_json_file)
