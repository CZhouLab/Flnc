import os
import re
import subprocess
import argparse
import logging
import datetime
import json
import textwrap
import Share

def HISAT2_Strandness_option(Strandness, SeqType):
        if Strandness == "first":
                if SeqType == "PE":
                        HISAT2_Strandness = "RF"
                elif SeqType == "SE":
                        HISAT2_Strandness = "R"
        elif Strandness == "second":
                if SeqType == "PE":
                        HISAT2_Strandness = "FR"
                elif SeqType == "SE":
                        HISAT2_Strandness = "F"
        
        return HISAT2_Strandness

def StringTie_Strandness_option(Strandness):
        if Strandness == "first":
                StringTie_Strandness = "rf"
        elif Strandness == "second":
                StringTie_Strandness = "fr"
        
        return StringTie_Strandness

def HTSeq_Strandness_option(Strandness):
        if Strandness == "first":
                HTSeq_Strandness = "reverse"
        elif Strandness == "second":
                HTSeq_Strandness = "yes"
        elif Strandness == "unstrand":
                HTSeq_Strandness = "no"

        return HTSeq_Strandness

def HISAT2_Alignment_SingleEnd(Replicate_index, R1_file, Strandness, OutputDic, LIB, json_file):
        with open(json_file) as input_config:
            config = json.load(input_config)
        
        HISAT2Index = config['file_config']['hisat2index']
        Known_splice_site = config['file_config']['known_splice_site']
       
        sed_1_path = LIB + "/sed_1.sh"
        sed_2_path = LIB + "/sed_2.sh"
        sed_3_path = LIB + "/sed_3.sh"

        HISAT2_novel_splice_site_path = OutputDic + "/Rep" + Replicate_index + "_HISAT2_novel_splice_site.txt"
        SeqType="SE"

        if (Strandness == "first") or (Strandness == "second"):
                HISAT2_Strandness = HISAT2_Strandness_option(Strandness,SeqType)
                run_hisat2_1 = subprocess.Popen([LIB+'/centos_0.1.sif', 'hisat2', '-p', '10', '--dta', '-x', HISAT2Index,
                                                '-U', R1_file,
                                                '--add-chrname', '--rna-strandness', HISAT2_Strandness, '--fr', '--known-splicesite-infile',
                                                Known_splice_site, '--novel-splicesite-outfile',HISAT2_novel_splice_site_path,
                                                '--novel-splicesite-infile', HISAT2_novel_splice_site_path, '--seed', '168', '--phred33',
                                                '--min-intronlen', '20', '--max-intronlen', '500000'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elif Strandness == "unstrand":
                run_hisat2_1 = subprocess.Popen([LIB+'/centos_0.1.sif', 'hisat2', '-p', '10', '--dta', '-x', HISAT2Index,
                                                '-U', R1_file,
                                                '--add-chrname', '--fr', '--known-splicesite-infile',
                                                Known_splice_site, '--novel-splicesite-outfile',HISAT2_novel_splice_site_path,
                                                '--novel-splicesite-infile', HISAT2_novel_splice_site_path, '--seed', '168', '--phred33',
                                                '--min-intronlen', '20', '--max-intronlen', '500000'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                                
        run_hisat2_2 = subprocess.Popen([LIB+'/centos_0.1.sif', 'sambamba', 'view', '-h', '-S', '-f', 'bam', '-t', '10', '/dev/stdin'], stdin=run_hisat2_1.stdout, stdout=subprocess.PIPE, )
        run_hisat2_3_str = '--sort-by-name -t 10 --tmpdir ' + OutputDic + "/Rep" + Replicate_index + ' -o ' + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.bam /dev/stdin'
        run_hisat2_3 = subprocess.Popen([LIB+'/centos_0.1.sif', 'sambamba', 'sort'] + run_hisat2_3_str.split(), stdin=run_hisat2_2.stdout, )
        run_hisat2_3.communicate()
                
        stderr_value = run_hisat2_1.communicate()[1]
        HISAT2_log = OutputDic + "/Rep" + Replicate_index + '.HISAT2'
        file(HISAT2_log, 'w').write(stderr_value.encode('utf-8'))

        run_sambamba_1_str = '-t 10 -h ' + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.bam'
        run_sambamba_1 = subprocess.Popen([LIB+'/centos_0.1.sif', 'sambamba', 'view'] + run_sambamba_1_str.split(), stdout=subprocess.PIPE, )
        run_sambamba_2 = subprocess.Popen([sed_1_path], stdin=run_sambamba_1.stdout, stdout=subprocess.PIPE, )
        run_sambamba_3 = subprocess.Popen([sed_2_path], stdin=run_sambamba_2.stdout, stdout=subprocess.PIPE, )
        run_sambamba_4 = subprocess.Popen([sed_3_path], stdin=run_sambamba_3.stdout, stdout=subprocess.PIPE, )
        run_sambamba_5_str = '-t 10 -f bam -S -o ' + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.renamed.bam /dev/stdin'
        run_sambamba_5 = subprocess.Popen([LIB+'/centos_0.1.sif', 'sambamba', 'view'] + run_sambamba_5_str.split(), stdin=run_sambamba_4.stdout)
        run_sambamba_5.communicate()
        os.system("rm -f " + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.bam')

        run_sambamba_6_str = '-t 10 --tmpdir=' + OutputDic + "/Rep" + Replicate_index + ' -o ' + OutputDic + "/Rep" + Replicate_index + '.sortedbycoord.bam ' + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.renamed.bam'
        run_sambamba_6 = subprocess.Popen([LIB+'/centos_0.1.sif', 'sambamba', 'sort'] + run_sambamba_6_str.split())
        run_sambamba_6.communicate()
        os.system("rm -f " + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.renamed.bam')
        
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 1: Rep" + str(Replicate_index) + " reads alignment (HISAT2) " + TimeNow)
        
def HISAT2_Alignment_PairedEnd(Replicate_index, R1_file, R2_file, Strandness, OutputDic, LIB, json_file):
        with open(json_file) as input_config:
            config = json.load(input_config)

        HISAT2Index = config['file_config']['hisat2index']
        Known_splice_site = config['file_config']['known_splice_site']

        sed_1_path = LIB + "/sed_1.sh"
        sed_2_path = LIB + "/sed_2.sh"
        sed_3_path = LIB + "/sed_3.sh"

        HISAT2_novel_splice_site_path = OutputDic + "/Rep" + Replicate_index + "_HISAT2_novel_splice_site.txt"
        SeqType="PE"

        if (Strandness == "first") or (Strandness == "second"):
                HISAT2_Strandness = HISAT2_Strandness_option(Strandness,SeqType)
                run_hisat2_1 = subprocess.Popen([LIB+'/centos_0.1.sif', 'hisat2', '-p', '10', '--dta', '-x', HISAT2Index,
                                                '-1', R1_file,'-2', R2_file,
                                                '--add-chrname', '--rna-strandness', HISAT2_Strandness, '--fr', '--known-splicesite-infile',
                                                Known_splice_site, '--novel-splicesite-outfile',HISAT2_novel_splice_site_path,
                                                '--novel-splicesite-infile', HISAT2_novel_splice_site_path, '--seed', '168', '--phred33',
                                                '--min-intronlen', '20', '--max-intronlen', '500000'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elif Strandness == "unstrand":
                run_hisat2_1 = subprocess.Popen([LIB+'/centos_0.1.sif', 'hisat2', '-p', '10', '--dta', '-x', HISAT2Index,
                                                '-1', R1_file,'-2', R2_file,
                                                '--add-chrname', '--fr', '--known-splicesite-infile',
                                                Known_splice_site, '--novel-splicesite-outfile',HISAT2_novel_splice_site_path,
                                                '--novel-splicesite-infile', HISAT2_novel_splice_site_path, '--seed', '168', '--phred33',
                                                '--min-intronlen', '20', '--max-intronlen', '500000'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        run_hisat2_2 = subprocess.Popen([LIB+'/centos_0.1.sif', 'sambamba', 'view', '-h', '-S', '-f', 'bam', '-t', '10', '/dev/stdin'], stdin=run_hisat2_1.stdout, stdout=subprocess.PIPE, )
        run_hisat2_3_str = '--sort-by-name -t 10 --tmpdir ' + OutputDic + "/Rep" + Replicate_index + ' -o ' + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.bam /dev/stdin'
        run_hisat2_3 = subprocess.Popen([LIB+'/centos_0.1.sif', 'sambamba', 'sort'] + run_hisat2_3_str.split(), stdin=run_hisat2_2.stdout, )
        run_hisat2_3.communicate()

        stderr_value = run_hisat2_1.communicate()[1]
        HISAT2_log = OutputDic + "/Rep" + Replicate_index + '.HISAT2'
        file(HISAT2_log, 'w').write(stderr_value.encode('utf-8'))

        run_sambamba_1_str = '-t 10 -h ' + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.bam'
        run_sambamba_1 = subprocess.Popen([LIB+'/centos_0.1.sif', 'sambamba', 'view'] + run_sambamba_1_str.split(), stdout=subprocess.PIPE, )
        run_sambamba_2 = subprocess.Popen([sed_1_path], stdin=run_sambamba_1.stdout, stdout=subprocess.PIPE, )
        run_sambamba_3 = subprocess.Popen([sed_2_path], stdin=run_sambamba_2.stdout, stdout=subprocess.PIPE, )
        run_sambamba_4 = subprocess.Popen([sed_3_path], stdin=run_sambamba_3.stdout, stdout=subprocess.PIPE, )
        run_sambamba_5_str = '-t 10 -f bam -S -o ' + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.renamed.bam /dev/stdin'
        run_sambamba_5 = subprocess.Popen([LIB+'/centos_0.1.sif', 'sambamba', 'view'] + run_sambamba_5_str.split(), stdin=run_sambamba_4.stdout)
        run_sambamba_5.communicate()
        os.system("rm -f " + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.bam')

        run_sambamba_6_str = '-t 10 --tmpdir=' + OutputDic + "/Rep" + Replicate_index + ' -o ' + OutputDic + "/Rep" + Replicate_index + '.sortedbycoord.bam ' + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.renamed.bam'
        run_sambamba_6 = subprocess.Popen([LIB+'/centos_0.1.sif', 'sambamba', 'sort'] + run_sambamba_6_str.split())
        run_sambamba_6.communicate()
        os.system("rm -f " + OutputDic + "/Rep" + Replicate_index + '.sortedbyname.renamed.bam')

        TimeNow = str(datetime.datetime.now())
        logging.info("Step 1: Rep" + str(Replicate_index) + " reads alignment (HISAT2) " + TimeNow)

def StringTie_Assembly(Replicate_index, Strandness, OutputDic, LIB, json_file):
        with open(json_file) as input_config:
            config = json.load(input_config)
        GTF_file = config['file_config']['StringTie_gtf']

        StringTie_GTF = OutputDic + "/Rep" + Replicate_index + '.stringtie.gtf'
        BamFile = OutputDic + "/Rep" + Replicate_index + '.sortedbycoord.bam'

        if (Strandness == "first") or (Strandness == "second"):
                StringTie_Strandness = StringTie_Strandness_option(Strandness)
                cmd = '--{StringTie_Strandness} -p 20 -G {GTF_file} -o {StringTie_GTF} -l {SampleName} -f 0 -m 200 -a 10 -j 1 -M 1 -g 50 {BamFile}'
                cmd_R = cmd.format(StringTie_Strandness = StringTie_Strandness, GTF_file = GTF_file, StringTie_GTF = StringTie_GTF, SampleName = "Rep"+Replicate_index, BamFile = BamFile)
        elif Strandness == "unstrand":
                cmd = '-p 20 -G {GTF_file} -o {StringTie_GTF} -l {SampleName} -f 0 -m 200 -a 10 -j 1 -M 1 -g 50 {BamFile}'
                cmd_R = cmd.format(GTF_file = GTF_file, StringTie_GTF = StringTie_GTF, SampleName = "Rep"+Replicate_index, BamFile = BamFile)
        cmd_Run = subprocess.Popen([LIB+'/centos_0.4.sif', 'stringtie'] + cmd_R.split())
        cmd_Run.communicate()

        TimeNow = str(datetime.datetime.now())
        logging.info("Step 2: Rep" + str(Replicate_index) + " transcripts assembly (StringTie) " + TimeNow)

        return StringTie_GTF

def Strawberry_Assembly(Replicate_index, Strandness,  OutputDic, LIB, json_file):
        with open(json_file) as input_config:
            config = json.load(input_config)
        GTF_file = config['file_config']['Strawberry_gtf']

        Strawberry_GTF = OutputDic + "/Rep" + Replicate_index + '.strawberry.gtf'
        BamFile = OutputDic + "/Rep" + Replicate_index + '.sortedbycoord.bam'

        if os.path.exists(Strawberry_GTF):
                os.remove(Strawberry_GTF)
                
        if (Strandness == "first") or (Strandness == "second"):
                StringTie_Strandness = StringTie_Strandness_option(Strandness)
                cmd = '--{StringTie_Strandness} -g {GTF_file} -o {Strawberry_GTF} -p 20 -m 0 -t 200 -s 10 -d 50 --no-quant --min-depth-4-transcript 0.1 {BamFile}'
                cmd_R = cmd.format(StringTie_Strandness = StringTie_Strandness, GTF_file = GTF_file, Strawberry_GTF = Strawberry_GTF, BamFile = BamFile)
        elif Strandness == "unstrand":
                cmd = '-g {GTF_file} -o {Strawberry_GTF} -p 20 -m 0 -t 200 -s 10 -d 50 --no-quant --min-depth-4-transcript 0.1 {BamFile}'
                cmd_R = cmd.format(GTF_file = GTF_file, Strawberry_GTF = Strawberry_GTF, BamFile = BamFile)
        cmd_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'strawberry'] + cmd_R.split())
        cmd_Run.communicate()

        TimeNow = str(datetime.datetime.now())
        logging.info("Step 3: Rep" + str(Replicate_index) + " transcripts assembly (Strawberry) " + TimeNow)
        
        return Strawberry_GTF

def GTF_fromStringTieMerge_ChangeID(GTF_ori, GTF_new, LIB):
                cmd_1 = LIB+'/src/gff3sort-master/gff3sort.pl {GTF_ori}'
                cmd_1_R = cmd_1.format(GTF_ori = GTF_ori)
                
                fo = open(GTF_ori + '.sort', "w")
                cmd_1_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'perl'] + cmd_1_R.split(), stdout = fo)
                cmd_1_Run.communicate()
                fo.close()

                cmd_2 = LIB+'/src/GTF_fromStringTieMerge_ChangeID.py {GTF_ori}.sort {GTF_new}'
                cmd_2_R = cmd_2.format(GTF_ori = GTF_ori, GTF_new = GTF_new)
                cmd_2_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_2_R.split())
                cmd_2_Run.communicate()


def StringTie_Merge(GTF_MergeList_File, OutputDic, LIB, json_file):
        with open(json_file) as input_config:
            config = json.load(input_config)
        GTF_file = config['file_config']['StringTie_gtf']

        Stringtie_merged_gtf = OutputDic + '/stringtie_merged.gtf'

        cmd = '--merge -p 20 -G {GTF_file} -F 0 -g 0 -f 0 -i -m 0 -T 0 -c 0 -o {Stringtie_merged_gtf} {GTF_MergeList_File}'
        cmd_R = cmd.format(GTF_file = GTF_file, Stringtie_merged_gtf = Stringtie_merged_gtf, GTF_MergeList_File = GTF_MergeList_File)
        cmd_Run = subprocess.Popen([LIB+'/centos_0.4.sif', 'stringtie'] + cmd_R.split())
        cmd_Run.communicate()

        Stringtie_merged_gtf_tmp = OutputDic + '/stringtie_merged.gtf.tmp'
        cmd_1 = "awk '$1!~/chrKI/ && $1!~/^GL/ && $7!=\".\" {print}' FS=\"\t\" OFS=\"\t\" " + Stringtie_merged_gtf + " >" + Stringtie_merged_gtf_tmp
        os.system(cmd_1)
        os.system("mv " + Stringtie_merged_gtf_tmp + " " + Stringtie_merged_gtf)

        Stringtie_merged_gtf_changeID = OutputDic + '/stringtie_merged.changeID.gtf'
        GTF_fromStringTieMerge_ChangeID(Stringtie_merged_gtf, Stringtie_merged_gtf_changeID, LIB)
        
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 4: Transcript assembly (merged) " + TimeNow)
    
        return Stringtie_merged_gtf_changeID

def GTF_RelatedFile_Generation(gtf_file, LIB, json_file):
        with open(json_file) as input_config:
            config = json.load(input_config)
        hg38_fasta = config['file_config']['genome_sequence']

        cmd_1 = LIB+'/src/Gene_Transcript_from_genegtf.py {gtf_file}'
        cmd_1_R = cmd_1.format(gtf_file = gtf_file)
        f1 = open(gtf_file + '.Gene_Transcript', "w")
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()
        
        cmd_2 = LIB+'/src/gtf2bed.pl {gtf_file}'
        cmd_2_R = cmd_2.format(gtf_file = gtf_file)
        f2 = open(gtf_file + '.bed', "w")
        cmd_2_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'perl'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        cmd_5 = LIB+'/src/ExonLength_forGene_fromgtf.R {gtf_file} {gtf_file}.ExonLength'
        cmd_5_R = cmd_5.format(gtf_file = gtf_file)
        cmd_5_Run = subprocess.Popen([LIB+'/centos_0.2.sif', 'Rscript'] + cmd_5_R.split())
        cmd_5_Run.communicate()
       
def PureLoci_RelatedTo_FilterRNA(BED_file, LIB, json_file, Pseudogenes_included = "T"):
        with open(json_file) as input_config:
            config = json.load(input_config)
        config = config['file_config']

        if Pseudogenes_included == "T":
                Filtered_RNA = config['Six_types']
        elif Pseudogenes_included == "F":
                Filtered_RNA = config['Five_types']

        cmd_1 = 'intersect -a {BED_file} -b {Filtered_RNA} -s -v -wa'
        cmd_1_R = cmd_1.format(BED_file = BED_file, Filtered_RNA = Filtered_RNA)
        f1 = open(BED_file + '.FilterRNA.NotOverlap', "w")
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.4.sif', 'bedtools'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()

        cmd_2 = 'intersect -a {BED_file} -b {Filtered_RNA} -s -wa'
        cmd_2_R = cmd_2.format(BED_file = BED_file, Filtered_RNA = Filtered_RNA)
        f2 = open(BED_file + '.FilterRNA.Overlap', "w")
        cmd_2_Run = subprocess.Popen([LIB+'/centos_0.4.sif', 'bedtools'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        cmd_3 = 'intersect -a {BED_file}.FilterRNA.NotOverlap -b {BED_file}.FilterRNA.Overlap -s -v -wa'
        cmd_3_R = cmd_3.format(BED_file = BED_file)
        f3 = open(BED_file + '.FilterRNA.NotOverlap_Double', "w")
        cmd_3_Run = subprocess.Popen([LIB+'/centos_0.4.sif', 'bedtools'] + cmd_3_R.split(), stdout = f3)
        cmd_3_Run.communicate()
        f3.close()
    
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 5: Transcripts selection (not overlapped with six/five types of annotated RNA on the same strand) " + TimeNow)
        return BED_file + '.FilterRNA.NotOverlap_Double'

def BED_2_fasta(BED_file, LIB, json_file):
        with open(json_file) as input_config:
            config = json.load(input_config)
        hg38_fasta = config['file_config']['genome_sequence']
        
        cmd_1 = 'getfasta -fi {hg38_fasta} -bed {BED_file} -name -split -s -fo {BED_file}.fa.tmp' # if "Segmentation fault" is reported, submit the job using bsub 
        cmd_1_R = cmd_1.format(BED_file = BED_file, hg38_fasta = hg38_fasta)
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.4.sif', 'bedtools'] + cmd_1_R.split())
        cmd_1_Run.communicate()

        cmd_2 = '-w 60 {BED_file}.fa.tmp'
        cmd_2_R = cmd_2.format(BED_file = BED_file)
        f2 = open(BED_file + '.fa', "w")
        cmd_2_Run = subprocess.Popen(['fold'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        return BED_file + '.fa'
        
def CPAT(FASTA_file, OutputDic, LIB, json_file):
        with open(json_file) as input_config:
            config = json.load(input_config)
        CPAT_hex = config['file_config']['CPAT_hex']
        CPAT_logitModel = config['file_config']['CPAT_logitmodel']
        
        cmd_1 = '--gene={FASTA_file} -o {OutputDic}/cpat.out --hex={CPAT_hex} --logitModel={CPAT_logitModel}'
        cmd_1_R = cmd_1.format(FASTA_file = FASTA_file, OutputDic = OutputDic, CPAT_hex = CPAT_hex, CPAT_logitModel = CPAT_logitModel)
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'cpat.py'] + cmd_1_R.split())
        cmd_1_Run.communicate()
        
        cmd_2 = '{OutputDic}/cpat.out.r'
        cmd_2_R = cmd_2.format(OutputDic = OutputDic)
        cmd_2_Run = subprocess.Popen([LIB+'/centos_0.2.sif', 'Rscript'] + cmd_2_R.split())
        cmd_2_Run.communicate()
       
        cmd_3 = "awk 'NR!=1 && $6>=0.364 {print $1}' FS=\"\t\" OFS=\"\t\" " + OutputDic + "/cpat.out | cut -d\"(\" -f 1 | sort | uniq >" + OutputDic + "/cpat.out.CodingTranscript.id"
        os.system(cmd_3)
        
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 6: Coding potential evaluation (CPAT) " + TimeNow)
        return  OutputDic + "/cpat.out.CodingTranscript.id"

def LGC(FASTA_file, OutputDic, LIB):
        cmd_1 = '{FASTA_file} {OutputDic}/LGC.out'
        cmd_1_R = cmd_1.format(FASTA_file = FASTA_file, OutputDic = OutputDic)
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.3.sif', 'lgc-1.0.0.py'] + cmd_1_R.split())
        cmd_1_Run.communicate()
        
        cmd_2 = "awk '$5==\"Coding\" {print $1}' FS=\"\t\" OFS=\"\t\" " + OutputDic + "/LGC.out | cut -d\"(\" -f 1 | sort | uniq >" + OutputDic + "/LGC.out.CodingTranscript.id"
        os.system(cmd_2)

        TimeNow = str(datetime.datetime.now())
        logging.info("Step 7: Coding potential evaluation (LGC) " + TimeNow)
        return OutputDic + "/LGC.out.CodingTranscript.id"

def PLEK(FASTA_file, OutputDic, LIB):
        cmd_1 = 'PLEK.py -fasta {FASTA_file} -out {OutputDic}/PLEK.out -thread 10'
        cmd_1_R = cmd_1.format(FASTA_file = FASTA_file, OutputDic = OutputDic)
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.1.sif'] + cmd_1_R.split())
        cmd_1_Run.communicate()

        cmd_2 = "awk '$1==\"Coding\" {print $3}' FS=\"\t\" OFS=\"\t\" " + OutputDic + "/PLEK.out | cut -d\">\" -f 2 | cut -d\"(\" -f 1 | sort | uniq >" + OutputDic + "/PLEK.out.CodingTranscript.id"
        os.system(cmd_2)
        
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 8: Coding potential evaluation (PLEK) " + TimeNow)
        return OutputDic + "/PLEK.out.CodingTranscript.id"

def CPPred(FASTA_file, OutputDic, LIB, json_file):
        with open(json_file) as input_config:
            config = json.load(input_config)

        Human_Hexamer = config['file_config']['CPPred_human_hexamer']
        Human_Range = config['file_config']['CPPred_human_range']
        Human_Model = config['file_config']['CPPred_human_model']
        
        CPPred_path = config['software_config']['CPPred_path']
        
        Current_path = os.getcwd()
        os.chdir(CPPred_path)
        
        cmd_1 = 'CPPred.py -i {FASTA_file} -hex {Human_Hexamer} -r {Human_Range} -mol {Human_Model} -spe Human -o {OutputDic}/CPPred.out'
        cmd_1_R = cmd_1.format(FASTA_file = FASTA_file, Human_Hexamer = Human_Hexamer, Human_Range = Human_Range, Human_Model = Human_Model, OutputDic = OutputDic)
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.3.sif', 'python'] + cmd_1_R.split())
        cmd_1_Run.communicate()
       
        os.chdir(Current_path)

        cmd_2 = "awk '$(NF-1)==\"coding\" {print $1}' FS=\"\t\" OFS=\"\t\" " + OutputDic + "/CPPred.out | cut -d\"(\" -f 1 | sort | uniq >" + OutputDic + "/CPPred.out.CodingTranscript.id"
        os.system(cmd_2)

        TimeNow = str(datetime.datetime.now())
        logging.info("Step 9: Coding potential evaluation (CPPred) " + TimeNow)
        return OutputDic + "/CPPred.out.CodingTranscript.id"

def Get_CodingTranscript_BED(gtf_file, BED_file, CodingTranscript_list, LIB):
        cmd_1 = LIB+'/src/Transcript_Gene_Transcript_inGTF.py {CodingTranscript_list} {gtf_file}'
        cmd_1_R = cmd_1.format(CodingTranscript_list = CodingTranscript_list, gtf_file = gtf_file)
        f1 = open(CodingTranscript_list + '.out', "w")
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()
        
        cmd_2 = LIB+'/src/my_join.pl -a {BED_file} -b {CodingTranscript_list}.out -F 4 -f 1'
        cmd_2_R = cmd_2.format(BED_file = BED_file, CodingTranscript_list = CodingTranscript_list)
        f2 = open(CodingTranscript_list + '.out.bed.tmp', "w")
        cmd_2_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'perl'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        cmd_3 = "awk '$13!=\"\"' FS=\"\t\" OFS=\"\t\" " + CodingTranscript_list + ".out.bed.tmp | cut -f 1-12 >" + CodingTranscript_list + ".out.bed"
        os.system(cmd_3)

        os.system("rm -f " + CodingTranscript_list + '.out.bed.tmp')
        return CodingTranscript_list + ".out.bed"

def PureLoci_RelatedTo_CodingRNA(gtf_file, BED_file, FASTA_file, OutputDic, LIB, json_file):
        CPAT_CodingTranscript_list = CPAT(FASTA_file, OutputDic, LIB, json_file)
        CPAT_CodingTranscript_BED = Get_CodingTranscript_BED(gtf_file, BED_file, CPAT_CodingTranscript_list, LIB)

        LGC_CodingTranscript_list = LGC(FASTA_file, OutputDic, LIB)
        LGC_CodingTranscript_BED = Get_CodingTranscript_BED(gtf_file, BED_file, LGC_CodingTranscript_list, LIB)

        PLEK_CodingTranscript_list = PLEK(FASTA_file, OutputDic, LIB)
        PLEK_CodingTranscript_BED = Get_CodingTranscript_BED(gtf_file, BED_file, PLEK_CodingTranscript_list, LIB)
        
        CPPred_CodingTranscript_list = CPPred(FASTA_file, OutputDic, LIB, json_file)
        CPPred_CodingTranscript_BED = Get_CodingTranscript_BED(gtf_file, BED_file, CPPred_CodingTranscript_list, LIB)
        
        cmd_0 = '{CPAT_CodingTranscript_BED} {LGC_CodingTranscript_BED} {PLEK_CodingTranscript_BED} {CPPred_CodingTranscript_BED}'
        cmd_0_R = cmd_0.format(CPAT_CodingTranscript_BED = CPAT_CodingTranscript_BED, LGC_CodingTranscript_BED = LGC_CodingTranscript_BED, PLEK_CodingTranscript_BED = PLEK_CodingTranscript_BED, CPPred_CodingTranscript_BED = CPPred_CodingTranscript_BED)
        output_file = OutputDic + '/codingpotential_transcript.bed'
        f0 = open(output_file, "w")
        cmd_0_Run = subprocess.Popen(['cat'] + cmd_0_R.split(), stdout = f0)
        cmd_0_Run.communicate()
        f0.close()
        
        cmd_1 = 'intersect -a {BED_file} -b {CodingRNA} -s -v -wa'
        cmd_1_R = cmd_1.format(BED_file = BED_file, CodingRNA = output_file)
        f1 = open(BED_file + '.CodingRNA.NotOverlap', "w")
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.4.sif', 'bedtools'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()

        cmd_2 = 'intersect -a {BED_file} -b {CodingRNA} -s -wa'
        cmd_2_R = cmd_2.format(BED_file = BED_file, CodingRNA = output_file)
        f2 = open(BED_file + '.CodingRNA.Overlap', "w")
        cmd_2_Run = subprocess.Popen([LIB+'/centos_0.4.sif', 'bedtools'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        cmd_3 = 'intersect -a {BED_file}.CodingRNA.NotOverlap -b {BED_file}.CodingRNA.Overlap -s -v -wa'
        cmd_3_R = cmd_3.format(BED_file = BED_file)
        f3 = open(BED_file + '.CodingRNA.NotOverlap_Double', "w")
        cmd_3_Run = subprocess.Popen([LIB+'/centos_0.4.sif', 'bedtools'] + cmd_3_R.split(), stdout = f3)
        cmd_3_Run.communicate()
        f3.close()
      
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 10: Transcripts selection (lacking coding potential) " + TimeNow)

        return BED_file + '.CodingRNA.NotOverlap_Double'

def SelectLongRNA(BED_file, LIB):
        cmd = LIB+'/src/LengthSelect_from_BED.py {BED_file} 200'
        cmd_R = cmd.format(BED_file = BED_file)
        f0 = open(BED_file + '.Long', "w")
        cmd_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_R.split(), stdout = f0)
        cmd_Run.communicate()
        f0.close()
        
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 11: Transcripts selection (size >= 200nt) " + TimeNow)
        return BED_file + '.Long'

def ReadCount(HISAT2_AlignSummary_file, condition = "Total"):
        number = 1        
        
        fi = open(HISAT2_AlignSummary_file, "r")
        for content in fi:
                if re.search("^Warning", content):
                        pass
                        number = number + 1
        fi.close()

        fi = open(HISAT2_AlignSummary_file, "r")
        for i in range(1,number):
                fi.readline()

        Total_ReadCount = fi.readline().strip().split(" ")[0]
        fi.readline()
        fi.readline()
        UniqueMapped_ReadCount = fi.readline().strip().split("(")[0].replace(" ","")
        MultipleMapped_ReadCount = fi.readline().strip().split("(")[0].replace(" ","")
        fi.close()

        if condition == "Total":
                return Total_ReadCount
        elif condition == "Unique":
                return UniqueMapped_ReadCount
        elif condition == "TotalMapped":
                return int(UniqueMapped_ReadCount)+int(MultipleMapped_ReadCount)
        else:
                pass
        

def FPKM(Replicate_index, BamFile, gtf_file, Strandness, OutputDic, LIB, ReadsNumber = "TotalMapped"):
        HTSeq_Strandness = HTSeq_Strandness_option(Strandness)
        cmd_1 = '-f bam -r name --stranded={HTSeq_Strandness} -m union {BamFile} {gtf_file}'
        cmd_1_R = cmd_1.format(HTSeq_Strandness = HTSeq_Strandness, BamFile = BamFile, gtf_file = gtf_file)
        f1 = open(OutputDic + "/Rep" + Replicate_index + '_HTSEQ_count.txt', "w")
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'htseq-count'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()
                
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 12: Rep" + str(Replicate_index) + " read count calculation " + TimeNow)

        count = ReadCount(OutputDic + "/Rep" + Replicate_index + ".HISAT2", condition = ReadsNumber) 
                
        cmd_2 = LIB+'/src/rpkm.pl {OutputDic}/Rep{Replicate_index}_HTSEQ_count.txt {gtf_file}.ExonLength {count}'
        cmd_2_R = cmd_2.format(OutputDic = OutputDic, Replicate_index = Replicate_index, gtf_file = gtf_file, count = count)
        f2 = open(OutputDic + "/Rep" + Replicate_index + '_Gene_FPKM.txt', "w")
        cmd_2_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'perl'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()
                
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 13: Rep" + str(Replicate_index) + " FPKM/RPKM calculation " + TimeNow)                

def Cat_FPKM_ReadCount(Replicate_index, gtf_file, OutputDic, LIB):
        FPKM_file = OutputDic + "/Rep" + Replicate_index + '_Gene_FPKM.txt'
        HTSeq_output = OutputDic + "/Rep" + Replicate_index + '_HTSEQ_count.txt'
                
        cmd_1 = LIB+'/src/Transcript_Gene_FPKM_ReadCounts.py {gtf_file} {FPKM_file} {HTSeq_output}'
        cmd_1_R = cmd_1.format(gtf_file = gtf_file, FPKM_file = FPKM_file, HTSeq_output = HTSeq_output)
        f1 = open(OutputDic + "/Rep" + Replicate_index + "_Transcript_FPKM_ReadCounts.txt", "w")
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()

        return OutputDic + "/Rep" + Replicate_index + "_Transcript_FPKM_ReadCounts.txt"

def SelectExpress(Replicate_index, OutputDic, LIB):
        BED_file = OutputDic + "/stringtie_merged.changeID.gtf.bed.FilterRNA.NotOverlap_Double.CodingRNA.NotOverlap_Double.Long"
        Expression_file = OutputDic + "/Rep" + Replicate_index + "_Transcript_FPKM_ReadCounts.txt"
    
        cutoff = os.popen(LIB+"/centos_0.2.sif Rscript "+LIB+"/src/FPKM_Cutoff.R {Expression_file} {BED_file}".format(Expression_file = Expression_file, BED_file = BED_file)).read()
        cutoff = cutoff.split('\n')[-2].split()[-1]

        cmd_1 = ("awk 'NR!=1 && ($3>%s && $4>10) {print $1}' FS=\"\t\" OFS=\"\t\" " % cutoff) + Expression_file + " >" + OutputDic + "/Rep" + Replicate_index + "_expressed_transcripts.id"
        os.system(cmd_1)           
        
        return OutputDic + "/Rep" + Replicate_index + "_expressed_transcripts.id"

def PutativeLncRNA_Generation(Expressed_transcripts_file, OutputDic, LIB):
        BED_file = OutputDic + "/stringtie_merged.changeID.gtf.bed.FilterRNA.NotOverlap_Double.CodingRNA.NotOverlap_Double.Long"

        cmd_3 = LIB+'/src/my_join.pl -a {BED_file} -b {Expressed_transcripts_file} -F 4 -f 1'
        cmd_3_R = cmd_3.format(BED_file = BED_file, Expressed_transcripts_file = Expressed_transcripts_file)
        f3 = open(BED_file + '.Expressed_tmp1', "w")
        cmd_3_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'perl'] + cmd_3_R.split(), stdout = f3)
        cmd_3_Run.communicate()
        f3.close()

        cmd_4 = "awk '$13!=\"\"' FS=\"\t\" OFS=\"\t\" " + BED_file + ".Expressed_tmp1 | cut -f 1-12 >" + BED_file + ".Expressed_tmp2"
        os.system(cmd_4)
        
        cmd_5 = LIB+'/src/my_join.pl -a {BED_file}.Expressed_tmp2 -b {OutputDic}/stringtie_merged.changeID.gtf.Gene_Transcript -F 4 -f 2'
        cmd_5_R = cmd_5.format(BED_file = BED_file, OutputDic = OutputDic)
        f5 = open(BED_file + '.Expressed_tmp3', "w")
        cmd_5_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'perl'] + cmd_5_R.split(), stdout = f5)
        cmd_5_Run.communicate()
        f5.close()

        cmd_6 = "awk '$1!~/chrKI/ && $1!~/^GL/ && $14!=\"\"' FS=\"\t\" OFS=\"\t\" " + BED_file + ".Expressed_tmp3 | cut -f 1-13 | sort | uniq >" + OutputDic + "/lncRNA.bed"
        os.system(cmd_6)

        os.system("rm -f " + BED_file + ".Expressed_tmp1")
        os.system("rm -f " + BED_file + ".Expressed_tmp2")
        os.system("rm -f " + BED_file + ".Expressed_tmp3")

        TimeNow = str(datetime.datetime.now())
        logging.info("Step 14: Transcripts selection (expressed) " + TimeNow)

        return OutputDic + "/lncRNA.bed"

def True_LncRNA_Infor(Replicate_count, OutputDic, LIB, user_json_file):
        Putaive_LncRNA_BED = OutputDic + "/lncRNA.bed"
        Putaive_LncRNA_Feature = OutputDic + "/lncRNA.features"
        
        with open(user_json_file) as input_config:
                config = json.load(input_config)
        
        gtf_file = config['file_config']["StringTie_gtf"]

        FPKM_file_list = []
        ReadCount_file_list = []
        FPKM_Dict = {}
        ReadCount_Dict = {}
        for i in range(0,int(Replicate_count)):
                Rep_index = i+1
                FPKM_file_list.append("Rep"+str(Rep_index)+"_FPKM")
                ReadCount_file_list.append("Rep"+str(Rep_index)+"_ReadCount")
        
                fj = open(OutputDic + "/Rep" + str(Rep_index) + "_Transcript_FPKM_ReadCounts.txt")
                fj.readline()
                for line_j in fj:
                        line_j = line_j.strip()
                        element_j = line_j.split("\t")
                        transcript_id = element_j[0]
                        fpkm = element_j[2]
                        read_count = element_j[3]

                        if transcript_id in FPKM_Dict:
                                FPKM_Dict[transcript_id].append(fpkm)
                        else:
                                FPKM_Dict[transcript_id] = list()
                                FPKM_Dict[transcript_id].append(fpkm)

                        if transcript_id in ReadCount_Dict:
                                ReadCount_Dict[transcript_id].append(read_count)
                        else:
                                ReadCount_Dict[transcript_id] = list()
                                ReadCount_Dict[transcript_id].append(read_count)
                fj.close()

        
        PutativeLncRNA_BED = OutputDic + "/lncRNA.bed"
        PutativeLncRNA_GTF = OutputDic + "/lncRNA.gtf"
        os.system(LIB+"/centos_0.1.sif python " + LIB + "/src/bed2gtf.py " + PutativeLncRNA_BED + " >" + PutativeLncRNA_GTF)
        os.system(LIB+"/centos_0.1.sif gffcompare -r " + gtf_file + " -o " + PutativeLncRNA_GTF + "_gffcompare " + PutativeLncRNA_GTF)

        cmd_2 = "cut -f 2,5 " + PutativeLncRNA_GTF + "_gffcompare.tracking | " + "awk '{print $3,$1}' FS=\"[\\t|]\" OFS=\"\\t\" >" + PutativeLncRNA_GTF + ".Transcript_Locus"
        os.system(cmd_2)

        Transcript_Locus_Dict = {}
        fy = open(PutativeLncRNA_GTF + ".Transcript_Locus", "r")
        for line_y in fy:
                line_y = line_y.strip()
                element_y = line_y.split("\t")
                tran = element_y[0]
                locus = element_y[1]
                Transcript_Locus_Dict[tran] = locus
        fy.close()

        Features = {}
        fk = open(OutputDic +"/lncRNA.features")
        fk.readline()
        for line_k in fk:
                line_k = line_k.strip()
                element_k = line_k.split("\t")
                tran_id = element_k[0]
                SE_ME = element_k[1]
                divergent = element_k[2]
                intergenic = element_k[3]
                antisense = element_k[4]
                tran_length = element_k[5]
                promoter = element_k[6]

                if tran_id in Features:
                        Features[tran_id].append(SE_ME + "\t" + divergent + "\t" + antisense + "\t" + intergenic + "\t" + promoter + "\t" + tran_length)
                else:
                        Features[tran_id] = list()
                        Features[tran_id].append(SE_ME + "\t" + divergent + "\t" + antisense + "\t" + intergenic + "\t" + promoter + "\t" + tran_length)
        fk.close()

        fy = open(PutativeLncRNA_BED, "r")
        fp = open(OutputDic + "/output/putative_lncRNAs.bed", "w")
        fp.write("chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts\n")
        fo = open(OutputDic + "/output/putative_lncRNA_infor.txt", "w")
        fo.write("TranscriptID\tLociID\tMulti_Exon\tDivergent\tAntisense\tIntergenic\tPromoter\tTranscriptLength\t"+"\t".join(FPKM_file_list)+"\t"+"\t".join(ReadCount_file_list)+"\n")
        for line_y in fy:
                line_y = line_y.strip()
                element_y = line_y.split("\t")
                tran_id = element_y[3]
                element_y.pop()

                if tran_id in Features:
                        feature = "".join(Features[tran_id])
                        fpkm = "\t".join(FPKM_Dict[tran_id])
                        read_count = "\t".join(ReadCount_Dict[tran_id])
                        locus_id = Transcript_Locus_Dict[tran_id]
                        fo.write(tran_id + "\t" + locus_id + "\t" + feature + "\t" + fpkm + "\t" + read_count + "\n")
                        fp.write("\t".join(element_y) + "\n")

        fy.close()
        fo.close()
        fp.close()

        for model_name in ["rf","lr","nb","dt","knn","rbfsvm","lsvm", "ensemble"]:
                True_LncRNA_Predict_file = OutputDic + "/" + model_name + "_lncRNA.predict"
                True_LncRNA_BED = OutputDic + "/true_lncRNA." + model_name + ".bed"
                True_LncRNA_infor = OutputDic + "/true_lncRNA_infor." + model_name + ".txt" 

                cmd = LIB+'/src/my_join.pl -a {OutputDic}/output/putative_lncRNAs.bed -b {True_LncRNA_Predict_file} -F 4 -f 1'
                cmd_R = cmd.format(OutputDic = OutputDic, True_LncRNA_Predict_file = True_LncRNA_Predict_file)
                f = open(True_LncRNA_BED+".tmp", "w")
                cmd_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'perl'] + cmd_R.split(), stdout = f)
                cmd_Run.communicate()
                f.close()

                cmd_1 = "awk 'NR==1 || $14==1' FS=\"\t\" OFS=\"\t\" " + True_LncRNA_BED + ".tmp | cut -f 1-12 >" + True_LncRNA_BED
                os.system(cmd_1)
                os.system("rm -f " + True_LncRNA_BED + ".tmp")

                cmd = LIB+'/src/my_join.pl -b {True_LncRNA_Predict_file} -a {OutputDic}/output/putative_lncRNA_infor.txt -F 1 -f 1'
                cmd_R = cmd.format(True_LncRNA_Predict_file = True_LncRNA_Predict_file, OutputDic = OutputDic)
                f = open(True_LncRNA_infor+".tmp", "w")
                cmd_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'perl'] + cmd_R.split(), stdout = f)
                cmd_Run.communicate()
                f.close()

                cmd_1 = "awk 'NR==1 || $NF==\"1\" {print}' FS=\"\\t\" OFS=\"\\t\" " + True_LncRNA_infor + ".tmp | rev | cut -f 2- | rev >" + True_LncRNA_infor
                os.system(cmd_1)
                os.system("rm -f " + True_LncRNA_infor + ".tmp")

def TrueLncRNA_PredictFrom_PutativeLncRNABED(Replicate_count, PutativeLncRNABED, library, OutputDic, model, user_json_file, step_index = 1):
        Share.Features(PutativeLncRNABED, OutputDic, library, user_json_file)
        TimeNow = str(datetime.datetime.now())
        logging.info("Step " + str(step_index) + ": Features identification of putative lncRNAs " + TimeNow)
        step_index = step_index + 1

        Share.True_LncRNA_Predict(OutputDic, library, user_json_file)
        True_LncRNA_Infor(Replicate_count, OutputDic, library, user_json_file)
        TimeNow = str(datetime.datetime.now())
        logging.info("Step " + str(step_index) + ": Identification of true lncRNAs " + TimeNow)

        if model != '':
                True_LncRNA_BED = OutputDic + "/true_lncRNA." + model + ".bed"
                True_LncRNA_infor = OutputDic + "/true_lncRNA_infor." + model + ".txt"
                if model in ["rf","lr","nb","dt","knn","rbfsvm","lsvm", "ensemble"]:
                        os.system("cp " + True_LncRNA_BED + " " + OutputDic + "/output/")
                        os.system("cp " + True_LncRNA_infor + " " + OutputDic + "/output/")
                else:
                        return
        else:
                True_LncRNA_BED = OutputDic + "/true_lncRNA.rf.bed"
                True_LncRNA_infor = OutputDic + "/true_lncRNA_infor.rf.txt"
                os.system("cp " + True_LncRNA_BED + " " + OutputDic + "/output/")
                os.system("cp " + True_LncRNA_infor + " " + OutputDic + "/output/")
