import os
import re
import subprocess
import argparse
import logging
import datetime
import json
import textwrap
import Mode1
import Mode2
from collections import OrderedDict


def GTF_SortByAlphabet(gtf_file, LIB, output_dir):
        GTF_SortByAlphabet_file = output_dir + "/SortByAlphabet.user.gtf"
        cmd = LIB+"/centos_0.1.sif perl "+LIB+"/src/gff3sort-master/gff3sort.pl " + gtf_file +" >" + GTF_SortByAlphabet_file
        os.system(cmd)
        return GTF_SortByAlphabet_file

def ProteinCoding_GTF_and_TSS(gtf_file, LIB, output_dir):
        ProteinCoding_gtf = output_dir + "/protein_coding.user.gtf"
        ProteinCoding_tss = output_dir + "/protein_coding.user.TSS"
        ProteinCoding_bed = output_dir + "/protein_coding.user.bed"

        cmd = 'egrep \"transcript_type \\"protein_coding\" ' + gtf_file + " >" + ProteinCoding_gtf
        os.system(cmd)

        cmd_1 = LIB+"/centos_0.1.sif perl "+LIB+"/src/gtf2bed.pl " + ProteinCoding_gtf + " >" + ProteinCoding_bed
        os.system(cmd_1)

        cmd_2 = LIB+"/centos_0.1.sif python "+LIB+"/src/TSS_fromBED.py " + ProteinCoding_bed + " >" + ProteinCoding_tss
        os.system(cmd_2)

        os.system("rm -f " + ProteinCoding_bed)
        return ProteinCoding_gtf, ProteinCoding_tss

def DifferentType_RNA(gtf_file, LIB, output_dir):
        Different_type_rna = output_dir + "/Different_Type_RNA.user.bed"
        cmd = LIB+"/centos_0.1.sif python "+LIB+"/src/Different_Type_RNA.py " + gtf_file +" >" + Different_type_rna
        os.system(cmd)
        return Different_type_rna

def Splice_site(gtf_file, LIB, output_dir):
        Splice_site_file = output_dir + "/splice_site.user.txt"
        cmd = LIB+"/centos_0.2.sif Rscript "+LIB+"/src/SpliceSite_fromGTF.R " + gtf_file + " " + Splice_site_file
        os.system(cmd)
        return Splice_site_file

def load_reference(args):
        GTF_SortByAlphabet_file = ''
        ProteinCoding_gtf = ''
        ProteinCoding_tss = ''
        if args.gtf_file != '':
                if os.path.exists(args.gtf_file):
                        GTF_SortByAlphabet_file = args.output_dir + "/SortByAlphabet.user.gtf"
                        if not os.path.exists(GTF_SortByAlphabet_file):
                                GTF_SortByAlphabet_file = GTF_SortByAlphabet(args.gtf_file, args.library, args.output_dir)

                        ProteinCoding_gtf = args.output_dir + "/protein_coding.user.gtf"
                        ProteinCoding_tss = args.output_dir + "/protein_coding.user.TSS"
                        if (not os.path.exists(ProteinCoding_gtf)) or (not os.path.exists(ProteinCoding_tss)):
                                ProteinCoding_infor = ProteinCoding_GTF_and_TSS(args.gtf_file, args.library, args.output_dir)
                                ProteinCoding_gtf = ProteinCoding_infor[0]
                                ProteinCoding_tss = ProteinCoding_infor[1]

                        Different_type_rna = args.output_dir + "/Different_Type_RNA.user.bed"
                        if not os.path.exists(Different_type_rna):
                                Different_type_rna = DifferentType_RNA(args.gtf_file, args.library, args.output_dir)

                        Splice_site_file = args.output_dir + "/splice_site.user.txt"
                        if not os.path.exists(Splice_site_file):
                                Splice_site_file = Splice_site(args.gtf_file, args.library, args.output_dir)

        ref_dict = {
            "file_config": {
                "hisat2index": "{Library}/HISAT2_index/genome".format(Library = args.library),
                "known_splice_site": "{Library}/gencode.v29.annotation.splice_site.txt".format(Library = args.library) if args.gtf_file=='' else Splice_site_file,
                "StringTie_gtf": "{Library}/gencode.v29.annotation.SortByAlphabet.gtf".format(Library = args.library) if args.gtf_file == '' else GTF_SortByAlphabet_file,
                "Strawberry_gtf": "{Library}/gencode.v29.annotation.SortByAlphabet.gtf".format(Library = args.library) if args.gtf_file == '' else GTF_SortByAlphabet_file,
                "genome_sequence": "{Library}/hg38.fa".format(Library = args.library),
                "Six_types": "{Library}/FiveTypeRNA/Six_types.RNA.bed".format(Library = args.library),
                "Five_types": "{Library}/FiveTypeRNA/Five_types.RNA.bed".format(Library = args.library),
                "CPAT_hex": "{Library}/CPAT/Human_Hexamer.tsv".format(Library = args.library),
                "CPAT_logitmodel": "{Library}/CPAT/Human_logitModel.RData".format(Library = args.library),
                "CPPred_human_hexamer": "{Library}/CPPred_20180516/Hexamer/Human_Hexamer.tsv".format(Library = args.library),
                "CPPred_human_range": "{Library}/CPPred_20180516/Human_Model/Human.range".format(Library = args.library),
                "CPPred_human_model": "{Library}/CPPred_20180516/Human_Model/Human.model".format(Library = args.library),
                "protein_coding_tss": "{Library}/protein_coding_gencodev29.TSS".format(Library = args.library) if args.gtf_file == '' else ProteinCoding_tss,
                "protein_coding_gtf": "{Library}/protein_coding_gencodev29.gtf".format(Library = args.library) if args.gtf_file == '' else ProteinCoding_gtf,
                "transcript_type": "{Library}/Different_Type_RNA.bed".format(Library = args.library) if args.gtf_file == '' else Different_type_rna,
                "chr_length": "{Library}/chrNameLength.txt".format(Library = args.library),
                "transcript_length_background": "{Library}/TranscriptLength_Background.txt".format(Library = args.library),
                "LogisticRegression_model": "{Library}/LogisticRegression.sav".format(Library = args.library),
                "DecisionTree_model": "{Library}/DecisionTree.sav".format(Library = args.library),
                "NaiveBayes_model": "{Library}/NaiveBayes.sav".format(Library = args.library),
                "KNN_model": "{Library}/KNN.sav".format(Library = args.library),
                "RBFSVM_model": "{Library}/RBFSVM.sav".format(Library = args.library),
                "LinearSVM_model": "{Library}/Linear_SVM.sav".format(Library = args.library),
                "RandomForest_model": "{Library}/RandomForest.sav".format(Library = args.library)
                },
            "software_config": {
                "CPPred_path": "{Library}/CPPred_20180516/bin".format(Library = args.library)
            }
        }
        with open(args.output_dir + '/user.json', "w") as outjson:
                json.dump(ref_dict, outjson, indent=4)

        TimeNow = str(datetime.datetime.now())
        logging.info("jason file (" + args.output_dir + "/user.json) generation " + TimeNow)

        return args.output_dir + '/user.json'

def Genomic_Structure(BED_file, LIB, json_file):
        with open(json_file) as input_config:
            config = json.load(input_config)

        Protein_coding_TSS = config['file_config']['protein_coding_tss']
        GTF = config['file_config']['protein_coding_gtf']

        cmd_1 = LIB+'/src/TSS_fromBED.py {BED_file}'
        cmd_1_R = cmd_1.format(BED_file = BED_file)
        f1 = open(BED_file + ".TSS", "w")
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()

        cmd_2 = '-a {Protein_coding_TSS} -b {BED_file}.TSS -w 2000 -Sm'
        cmd_2_R = cmd_2.format(BED_file = BED_file, Protein_coding_TSS = Protein_coding_TSS)
        f2 = open(BED_file + ".divergent.tmp", "w")
        cmd_2_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'windowBed'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        cmd_3 = "cut -f 10 " + BED_file + ".divergent.tmp |  sort | uniq >" + BED_file + ".divergent"
        os.system(cmd_3)

        cmd_4 = LIB+'/src/my_join.pl -a {BED_file} -b {BED_file}.divergent -F 4 -f 1'
        cmd_4_R = cmd_4.format(BED_file = BED_file)
        f4 = open(BED_file + '.non-divergent.tmp1', "w")
        cmd_4_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'perl'] + cmd_4_R.split(), stdout = f4)
        cmd_4_Run.communicate()
        f4.close()

        cmd_5 = "awk '$14==\"\"' FS=\"\t\" OFS=\"\t\" " + BED_file + ".non-divergent.tmp1 | cut -f 1-13 >" + BED_file + ".non-divergent.tmp2"
        os.system(cmd_5)

        os.system(LIB+"/centos_0.1.sif python "+LIB+"/src/bed2gtf.py " + BED_file + ".non-divergent.tmp2 >" + BED_file + ".non-divergent.gtf")
        os.system("echo " + BED_file + ".non-divergent.gtf >" + BED_file + ".non-divergent_gtf")

        cmd_6 = '-r {GTF} -i {BED_file}.non-divergent_gtf -o {BED_file}.non-divergent_gffcompare'
        cmd_6_R = cmd_6.format(GTF = GTF, BED_file = BED_file)
        cmd_6_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'gffcompare'] + cmd_6_R.split())
        cmd_6_Run.communicate()

        cmd_7 = LIB+'/src/gffcompare_tracking.py {BED_file}.non-divergent_gffcompare.tracking {BED_file} {BED_file}.divergent {json_file}'
        cmd_7_R = cmd_7.format(BED_file = BED_file, json_file = json_file)
        f7 = open(BED_file + '.Genomic_Structure', "w")
        cmd_7_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_7_R.split(), stdout = f7)
        cmd_7_Run.communicate()
        f7.close()

        os.system("rm -f " + BED_file + '.divergent')
        os.system("rm -f " + BED_file + '.divergent.tmp')
        os.system("rm -f " + BED_file + '.non-divergent*')

        return BED_file + '.Genomic_Structure'

def Single_Multi_Exon(BED_file):
        cmd_1 = "awk '$10==1 {print $4,\"SE\"}' FS=\"\t\" OFS=\"\t\" " + BED_file + " >" + BED_file + ".SE"
        cmd_2 = "awk '$10>1 {print $4,\"ME\"}' FS=\"\t\" OFS=\"\t\" " + BED_file + " >" + BED_file + ".ME"
        cmd_3 = "cat " + BED_file + ".SE " + BED_file + ".ME >" + BED_file + ".Single_Multi_Exon"
        os.system(cmd_1)
        os.system(cmd_2)
        os.system(cmd_3)

        os.system("rm -f " + BED_file + '.SE')
        os.system("rm -f " + BED_file + '.ME')

        return BED_file + '.Single_Multi_Exon'

def PromoterSignature(BED_file, TSSG_output, LIB, json_file, TSS_Extend_Length = 1000):
        with open(json_file) as input_config:
            config = json.load(input_config)
        hg38_fasta = config['file_config']['genome_sequence']

        if not os.path.isdir(TSSG_output):
                os.mkdir(TSSG_output)
        else:
                os.system("rm -rf " + TSSG_output)
                os.mkdir(TSSG_output)

        cmd_1 = LIB+'/src/Promoter_fromBED.py {BED_file} {TSS_Extend_Length} {json_file}'
        cmd_1_R = cmd_1.format(BED_file = BED_file, TSS_Extend_Length = TSS_Extend_Length, json_file = json_file)
        f1 = open(BED_file + ".promoter", "w")
        cmd_1_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()

        cmd_2 = 'getfasta -fi {hg38_fasta} -bed {BED_file}.promoter -name -split -s -fo {BED_file}.promoter.fa.tmp'
        cmd_2_R = cmd_2.format(hg38_fasta = hg38_fasta, BED_file = BED_file)
        cmd_2_Run = subprocess.Popen([LIB+'/centos_0.4.sif', 'bedtools'] + cmd_2_R.split())
        cmd_2_Run.communicate()

        cmd_3 = 'fold -w 60 ' + BED_file +".promoter.fa.tmp >" + BED_file + ".promoter.fa"
        os.system(cmd_3)

        cmd_4 = LIB+'/src/Split_FastaFile.py {BED_file}.promoter.fa {TSSG_output}'
        cmd_4_R = cmd_4.format(BED_file = BED_file, TSSG_output = TSSG_output)
        cmd_4_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_4_R.split())
        cmd_4_Run.communicate()

        fa_list = os.listdir(TSSG_output)
        current_directory = os.getcwd()
        os.chdir(LIB + "/TSSG")
        RunTSSG_path = "./run_tssg.sh"
        for File in fa_list:
                cmd_5 = '{TSSG_output}/{File} {TSSG_output}/{File}.out'
                cmd_5_R = cmd_5.format(TSSG_output = TSSG_output, File = File)
                cmd_5_Run = subprocess.Popen([RunTSSG_path] + cmd_5_R.split())
                cmd_5_Run.communicate()

        os.chdir(current_directory)
        os.system("rm -f " + BED_file + ".PromoterSignature")
        for File in fa_list:
                if not re.search("out$", File):
                        cmd_6 = LIB+'/src/TSSG.py {TSSG_output}/{File}.out'
                        cmd_6_R = cmd_6.format(TSSG_output = TSSG_output, File = File)
                        f6 = open(BED_file + ".PromoterSignature", "a")
                        cmd_6_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_6_R.split(), stdout = f6)
                        cmd_6_Run.communicate()
                        f6.close()

        os.system("rm -f " + BED_file + '.promoter*')
        os.system("rm -rf " + TSSG_output)

        return BED_file + ".PromoterSignature"

def Transcript_Length(BED_file):
        fo =open(BED_file + ".Transcript_Length", "w")
        fi = open(BED_file, "r")
        for line in fi:
                line=line.strip()
                element=line.split("\t")
                T=element[3]
                ExonSize=element[10]
                ExonSize_List=ExonSize.split(",")
                ExonSize_Value=[]
                for i in range(0,len(ExonSize_List)-1):
                        ExonSize_Value.append(int(ExonSize_List[i]))

                TranscriptLength=sum(ExonSize_Value)
                fo.write(T+"\t"+str(TranscriptLength)+"\n")
        fi.close()
        fo.close()

        return BED_file + ".Transcript_Length"

def Features(PutativeLncRNA_BED, output_dir, LIB, json_file):
        TSSG_name = os.path.basename(output_dir)
        
        Genomic_Structure_output = Genomic_Structure(PutativeLncRNA_BED, LIB, json_file)
        Single_Multi_Exon_output = Single_Multi_Exon(PutativeLncRNA_BED)
        PromoterSignature_output = PromoterSignature(PutativeLncRNA_BED, LIB+"/TSSG/input/" + TSSG_name, LIB, json_file)
        Transcript_Length_output = Transcript_Length(PutativeLncRNA_BED)

        f1 = open(Single_Multi_Exon_output, "r")
        Single_Multi_Exon_Dict = {}
        for line in f1:
                line = line.strip()
                element = line.split("\t")
                Transcript = element[0]
                Type = element[1]

                Type_new = "NA"
                if Type == "SE":
                        Type_new = "0"
                elif Type  == "ME":
                        Type_new = "1"
                Single_Multi_Exon_Dict[Transcript] = Type_new
        f1.close()

        f2 = open(Genomic_Structure_output, "r")
        Genomic_Struture_Dict = {}
        for line in f2:
                line = line.strip()
                element = line.split("\t")
                Transcript = element[0]
                Type = element[1]

                if Type != "NO_USE":
                        Type_new = "NA\tNA\tNA"
                        if Type == "divergent":
                                Type_new = "1\t0\t0"
                        elif Type == "intergenic":
                                Type_new = "0\t1\t0"
                        elif Type == "antisense":
                                Type_new = "0\t0\t1"

                        Genomic_Struture_Dict[Transcript] = Type_new
        f2.close()

        f4 = open(Transcript_Length_output, "r")
        Length_Dict = {}
        for line in f4:
                line = line.strip()
                element = line.split("\t")
                Transcript = element[0]
                Length = element[1]

                Length_Dict[Transcript] = Length
        f4.close()

        f5 = open(PromoterSignature_output, "r")
        PromoterSignature_Dict = {}
        for line in f5:
                line = line.strip()
                element = line.split("\t")
                Transcript = element[0]
                PromoterSig = element[1]

                if int(PromoterSig) >= 1:
                        PromoterSignature_Dict[Transcript] = "1"
                else:
                        PromoterSignature_Dict[Transcript] = "0"
        f5.close()

        fo = open(output_dir + "/lncRNA.features", "w")
        fo.write("ID\tSE_ME\tdivergent\tintergenic\tantisense\tTranscriptLength\tPromoter\n")
        for T in Genomic_Struture_Dict:
                SE_ME = Single_Multi_Exon_Dict[T]
                Genomic_Struture = Genomic_Struture_Dict[T]
                Size = Length_Dict[T]
                Promoter = PromoterSignature_Dict[T]
                output = T+"\t"+SE_ME+"\t"+str(Genomic_Struture)+"\t"+str(Size)+"\t"+str(Promoter)
                fo.write(output+"\n")

def MergeGTF_2_NoncodingLongTranscripts_BED(GTF_input, output_dir, library, user_json_file):
        Mode1.GTF_RelatedFile_Generation(GTF_input, library, user_json_file)
        PureLoci_RelatedTo_FilterRNA_BED = Mode1.PureLoci_RelatedTo_FilterRNA(GTF_input + '.bed', library, user_json_file)
        PureLoci_RelatedTo_FilterRNA_fasta = Mode1.BED_2_fasta(PureLoci_RelatedTo_FilterRNA_BED, library, user_json_file)
        PureLoci_RelatedTo_CodingRNA_BED = Mode1.PureLoci_RelatedTo_CodingRNA(GTF_input, PureLoci_RelatedTo_FilterRNA_BED, PureLoci_RelatedTo_FilterRNA_fasta, output_dir, library, user_json_file)
        Mode1.SelectLongRNA(PureLoci_RelatedTo_CodingRNA_BED, library)

def EnsembleModel(output_dir):
        rf = output_dir + "/rf_lncRNA.predict"
        lr = output_dir + "/lr_lncRNA.predict"
        nb = output_dir + "/nb_lncRNA.predict"
        dt = output_dir + "/dt_lncRNA.predict"
        knn = output_dir + "/knn_lncRNA.predict"
        rbfsvm = output_dir + "/rbfsvm_lncRNA.predict"
        lsvm = output_dir + "/lsvm_lncRNA.predict"

        Label_Dict = OrderedDict()
        input_file_list = [rf, dt, knn, lr, rbfsvm, nb, lsvm]
        for sub_file in input_file_list:
                fi = open(sub_file, "r")
                fi.readline()
                for line in fi:
                        line = line.strip()
                        element = line.split("\t")
                        ID = element[0]
                        Label = element[1]
                        
                        if ID in Label_Dict:
                                Label_Dict[ID].append(Label)
                        else:
                                Label_Dict[ID] = list()
                                Label_Dict[ID].append(Label)
                fi.close

        ensemble_output = output_dir + "/ensemble_lncRNA.predict"
        fo = open(ensemble_output, "w")
        fo.write("ID\tPredict_Label\n")
        for ID in Label_Dict:
                value = "\t".join(Label_Dict[ID])
                if value == "1\t1\t1\t1\t1\t1\t1":
                        fo.write(ID+"\t1\n")
                else:
                        fo.write(ID+"\t0\n")
        fo.close()

        return ensemble_output

def True_LncRNA_Predict(output_dir, LIB, json_file):
        Putaive_LncRNA_Feature = output_dir + "/lncRNA.features"

        with open(json_file) as input_config:
            config = json.load(input_config)

        gtf_file = config['file_config']["StringTie_gtf"]

        ModelName_in_short = ["rf","lr","nb","dt","knn","rbfsvm","lsvm"]
        ModelName_in_jason = ["RandomForest_model","LogisticRegression_model","NaiveBayes_model","DecisionTree_model","KNN_model","RBFSVM_model","LinearSVM_model"]
        for model_index in range(0,len(ModelName_in_jason)):
                model_name = ModelName_in_short[model_index]
                FullLength_LncRNA_Predict = output_dir + "/" + model_name + "_lncRNA.predict"
                
                model_file = config['file_config'][ModelName_in_jason[model_index]]
                
                cmd_0 = LIB+'/src/FullLength_lncRNA_predict.py {model_file} {Putaive_LncRNA_Feature} {FullLength_LncRNA_Predict} {json_file}'
                cmd_0_R = cmd_0.format(model_file = model_file, Putaive_LncRNA_Feature = Putaive_LncRNA_Feature, FullLength_LncRNA_Predict = FullLength_LncRNA_Predict, json_file = json_file)
                cmd_0_Run = subprocess.Popen([LIB+'/centos_0.1.sif', 'python'] + cmd_0_R.split())
                cmd_0_Run.communicate()

        ensemble_predict_file = EnsembleModel(output_dir)
