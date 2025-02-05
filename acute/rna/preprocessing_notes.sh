###################################
####GENERATING THE SAMPLE SHEET####
###################################

#cd for file in 2428__*; do echo "$file" | awk -F'2428__|-read' '{print $2}'; done | sort | uniq | grep "Chrom*" > /scratch/tmp/ed7982/chromium/01_snpcalling/raw_longevity_list.txt
cd /Genomics/ayroleslab2/scott/git/chromium/data/acute/rna/
for file in 2427__*fastq.gz; do echo "$file" | awk -F'2427__|-read' '{print $2}'; done | sort | uniq | grep "Chrom*" > /Genomics/argo/users/ed7982/crvi-analysis/acute/rna/raw_rna_list.txt
cd /Genomics/argo/users/ed7982/crvi-analysis/acute/rna

# based on #https://docs.google.com/spreadsheets/d/1m-OEZSrV5PRa-gBCsA7YdVIYIL5dI4wK/edit?gid=1330189754#gid=1330189754
# Chrom_40 is acute and Chrom_50 is crvi

# we want sample name, treatment, plate, head/body
awk -F'_' '{
    # Determine treatment based on the second field
    if ($2 == 40) treatment = "control";
    else if ($2 == 50) treatment = "3mM";
    else { print "Error: Invalid second field value in line: " $0 > "/dev/stderr"; next; }

    # Determine tissue based on the fourth field
    if ($4 == "b") tissue = "body";
    else if ($4 == "h") tissue = "head";
    else { print "Error: Invalid fourth field value in line: " $0 > "/dev/stderr"; next; }

    # Create plate as fields 2, 3, and 4 joined by underscores
    plate = $2 "_" $3 "_" $4;

    # Determine date based on plate
    if (plate == "40_1_b" || plate == "40_2_b" || plate == "50_1_h" || plate == "50_2_h") date = "8_4_22";
    else if (plate == "40_3_b" || plate == "40_4_b" || plate == "50_3_h" || plate == "50_4_h") date = "8_5_22";
    else if (plate == "40_1_h" || plate == "40_2_h" || plate == "50_1_b" || plate == "50_2_b") date = "8_8_22";
    else if (plate == "40_3_h" || plate == "40_4_h" || plate == "50_3_b" || plate == "50_4_b") date = "8_9_22";
    else { print "Error: Invalid plate value in line: " $0 > "/dev/stderr"; next; }

    # Print the output
    print $0 "\t" treatment "\t" tissue "\t" plate "\t" date;
}' raw_rna_list.txt > rna_sample_sheet_noheader.txt

echo -e "sample_name\ttreatment\ttissue\tplate\tdate" > rna_sample_sheet_headeronly.txt

cat rna_sample_sheet_headeronly.txt rna_sample_sheet_noheader.txt > rna_sample_sheet.tsv



###################################
#####GENERATING THE UNIT SHEET#####
###################################


#need to convert luisa's trimmomatic arguments to cutadapt arguments
#luisa says list of trimmers (see manual)
    # trimmers:
    #   - "ILLUMINACLIP:AdaptersTrim.fasta:1:30:7"
    #   - "SLIDINGWINDOW:4:15"
    #   - "MINLEN:20"
#these 3 numbers after ILLUMINACLIP are
#seedMismatches: specifies the maximum mismatch count which will still allow a full
# match to be performed
# palindromeClipThreshold: specifies how accurate the match between the two 'adapter
# ligated' reads must be for PE palindrome read alignment.
# simpleClipThreshold: specifies how accurate the match between any adapter etc.
# sequence must be against a read.

#cutadapt says list of trimmers (see manua

# --poly-a
# -e 1
# --nextseq-trim=20
# --trim-n

cd /Genomics/ayroleslab2/scott/git/chromium/data/acute/rna/
for file in 2427__*fastq.gz; do echo "$file" | awk -F'2427__|-read' '{print $2}'; done | sort | uniq | grep "Chrom*" > /Genomics/argo/users/ed7982/crvi-analysis/acute/rna/raw_rna_list.txt
cd /Genomics/argo/users/ed7982/crvi-analysis/acute/rna


#-a file:/Genomics/argo/users/ed7982/crvi-analysis/acute/rna/AdaptersTrim.fasta --poly-a --quality-cutoff 10,0 -e 1 --trim-n --nextseq-trim=20 -m 20 --overlap 2

awk '{
    sample_name=$1
    fastq1="/Genomics/ayroleslab2/scott/git/chromium/data/acute/rna/2427__" sample_name "-read-1.fastq.gz"
    #fastq4="/Genomics/ayroleslab2/scott/git/chromium/data/acute/rna/2427__" sample_name "-read-4.fastq.gz"
    adapters="-a file:/Genomics/argo/users/ed7982/crvi-analysis/acute/rna/AdaptersTrim.fasta --poly-a --quality-cutoff 10,0 -e 1 --trim-n --nextseq-trim=20 -m 20 --overlap 2"
    print sample_name "\t" sample_name "\t" fastq1 "\t" "\t" "\t" adapters "\t" "none" #removed fq4
}' raw_rna_list.txt > rna_unit_sheet_noheader.txt

echo -e "sample_name\tunit_name\tfq1\tfq2\tsra\tadapters\tstrandedness" > rna_unit_sheet_headeronly.txt

cat rna_unit_sheet_headeronly.txt rna_unit_sheet_noheader.txt > rna_unit_sheet.tsv
rm rna_unit_sheet_noheader.txt rna_unit_sheet_headeronly.txt rna_sample_sheet_noheader.txt rna_sample_sheet_headeronly.txt



###################################
#####GET BODY/HEAD SHEETS ONLY#####
###################################
awk 'NR==1 || $1 ~ /_b_/' rna_sample_sheet.tsv > rna_sample_sheet_body.tsv
awk 'NR==1 || $1 ~ /_b_/' rna_unit_sheet.tsv > rna_unit_sheet_body.tsv

awk 'NR==1 || $1 ~ /_h_/' rna_sample_sheet.tsv > rna_sample_sheet_head.tsv
awk 'NR==1 || $1 ~ /_h_/' rna_unit_sheet.tsv > rna_unit_sheet_head.tsv


#moved over manually

head -n 3 rna_sample_sheet_head.tsv > ./head/rna_sample_sheet_head_top3.tsv
head -n 3 rna_unit_sheet_head.tsv > ./head/rna_unit_sheet_head_top3.tsv 
