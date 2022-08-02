import os
import re
import subprocess
import argparse
import logging
import datetime
import json
import textwrap
import Share

def fasta_2_bed(Fasta_file, OutputDic, LIB, json_file):
        with open(json_file) as input_config:
            config = json.load(input_config)
        hg38_fasta = config['file_config']['genome_sequence']

        cmd_0 = '{hg38_fasta} {Fasta_file} {OutputDic}/blat.psl'
        cmd_0_R = cmd_0.format(hg38_fasta = hg38_fasta, Fasta_file = Fasta_file, OutputDic = OutputDic)
        cmd_0_Run = subprocess.Popen([LIB+'/centos_0.5.sif', 'blat'] + cmd_0_R.split())
        cmd_0_Run.communicate()

        cmd = LIB+'/src/psl_2_bed.py {OutputDic}/blat.psl {json_file}'
        cmd_R = cmd.format(OutputDic = OutputDic, json_file = json_file)
        f = open(OutputDic+"/input.bed", "w")
        cmd_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_R.split(), stdout = f)
        cmd_Run.communicate()
        f.close()
        
        return OutputDic + "/input.bed"

def True_LncRNA_Infor(OutputDic, LIB):
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

                cmd_1 = "awk 'NR==1 || $14==1' FS=\"\\t\" OFS=\"\\t\" " + True_LncRNA_BED + ".tmp | cut -f 1-12 >" + True_LncRNA_BED
                os.system(cmd_1)
                os.system("rm -f " + True_LncRNA_BED + ".tmp")

                cmd = LIB+'/src/my_join.pl -b {True_LncRNA_Predict_file} -a {OutputDic}/output/putative_lncRNA_infor.txt -F 1 -f 1'
                cmd_R = cmd.format(True_LncRNA_Predict_file = True_LncRNA_Predict_file, OutputDic = OutputDic)
                f = open(True_LncRNA_infor+".tmp", "w")
                cmd_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'perl'] + cmd_R.split(), stdout = f)
                cmd_Run.communicate()
                f.close()

                cmd_1 = "awk 'NR==1 || $NF==\"1\" {print}' FS=\"\\t\" OFS=\"\\t\" " + True_LncRNA_infor + ".tmp | cut -f 1-10 >" + True_LncRNA_infor
                os.system(cmd_1)
                os.system("rm -f " + True_LncRNA_infor + ".tmp")

def TrueLncRNA_PredictFrom_PutativeLncRNABED(PutativeLncRNABED, library, OutputDic, model, user_json_file, step_index = 1):
        Share.Features(PutativeLncRNABED, OutputDic, library, user_json_file)
        TimeNow = str(datetime.datetime.now())
        logging.info("Step " + str(step_index) + ": Features identification of putative lncRNAs " + TimeNow)
        step_index = step_index + 1

        os.system("mkdir " + OutputDic + "/output")
       
        os.system('echo -e TranscriptID"\t"LocusID"\t"Multi_Exon"\t"Divergent"\t"Antisense"\t"Intergenic"\t"Promoter"\t"TranscriptLength"\t"FPKM"\t"ReadCount >' + OutputDic + '/output/putative_lncRNA_infor.txt')
        cmd = "awk 'NR!=1 {print $1,\"NA\",$2,$3,$5,$4,$7,$6,\"NA\",\"NA\"}' FS=\"\\t\" OFS=\"\\t\" " + OutputDic + "/lncRNA.features >>" + OutputDic + "/output/putative_lncRNA_infor.txt"
        os.system(cmd)
        os.system("cut -f 1-12 " + OutputDic + "/input.bed >" + OutputDic + "/output/putative_lncRNAs.bed")
        os.system("sed -i '1i chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts' " + OutputDic + "/output/putative_lncRNAs.bed") 

        Share.True_LncRNA_Predict(OutputDic, library, user_json_file)
        True_LncRNA_Infor(OutputDic, library)
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
