#!/bin/bash


#rename sampleIDs
##############################################################################################################

#
echo "run script for rna.seq.analysis"

####### rename dir ########


mv NG-14833_Pas_E_1_lib233711_5747_1 Pas.15dpa.1
mv NG-14833_Pas_E_2_lib233712_5747_1 Pas.15dpa.2
mv NG-14833_Pas_E_3_lib233713_5767_1 Pas.15dpa.3
mv NG-14833_Pas_M_1_lib233714_5767_1 Pas.25dpa.1
mv NG-14833_Pas_M_2_lib233715_5767_1 Pas.25dpa.2
mv NG-14833_Pas_M_3_lib233716_5767_1 Pas.25dpa.3

mv NG-14833_Reg_15dpa_1_lib233717_5767_1 Reg.15dpa.1
mv NG-14833_Reg_15dpa_2_lib233718_Merge Reg.15dpa.2
mv NG-14833_Reg_15dpa_3_lib233719_Merge Reg.15dpa.3
mv NG-14833_Reg_25dpa_1_lib233720_Merge Reg.25dpa.1
mv NG-14833_Reg_25dpa_2_lib233721_5747_6 Reg.25dpa.2
mv NG-14833_Reg_25dpa_3_lib233722_5747_6 Reg.25dpa.3

mv NG-14833_Tit_15dpa_1_lib233705_5747_1 Tit.15dpa.1
mv NG-14833_Tit_15dpa_2_lib233706_Merge Tit.15dpa.2
mv NG-14833_Tit_15dpa_3_lib233707_5747_1 Tit.15dpa.3
mv NG-14833_Tit_25dpa_1_lib233708_5747_1 Tit.25dpa.1
mv NG-14833_Tit_25dpa_2_lib233709_5747_1 Tit.25dpa.2
mv NG-14833_Tit_25dpa_3_lib233710_Merge Tit.25dpa.3


####### rename files ########

find . -iname "*NG-14833_Pas_E_1_lib233711_5747_1*" -exec rename 's/NG-14833_Pas_E_1_lib233711_5747_1/Pas.15dpa.1/' '{}' \;
find . -iname "*NG-14833_Pas_E_2_lib233712_5747_1*" -exec rename 's/NG-14833_Pas_E_2_lib233712_5747_1/Pas.15dpa.2/' '{}' \;
find . -iname "*NG-14833_Pas_E_3_lib233713_5767_1*" -exec rename 's/NG-14833_Pas_E_3_lib233713_5767_1/Pas.15dpa.3/' '{}' \;
find . -iname "*NG-14833_Pas_M_1_lib233714_5767_1*" -exec rename 's/NG-14833_Pas_M_1_lib233714_5767_1/Pas.25dpa.1/' '{}' \;
find . -iname "*NG-14833_Pas_M_2_lib233715_5767_1*" -exec rename 's/NG-14833_Pas_M_2_lib233715_5767_1/Pas.25dpa.2/' '{}' \;
find . -iname "*NG-14833_Pas_M_3_lib233716_5767_1*" -exec rename 's/NG-14833_Pas_M_3_lib233716_5767_1/Pas.25dpa.3/' '{}' \;

find . -iname "*NG-14833_Reg_15dpa_1_lib233717_5767_1*" -exec rename 's/NG-14833_Reg_15dpa_1_lib233717_5767_1/Reg.15dpa.1/' '{}' \;

find . -iname "*NG-14833_Reg_25dpa_2_lib233721_5747_6*" -exec rename 's/NG-14833_Reg_25dpa_2_lib233721_5747_6/Reg.25dpa.2/' '{}' \;
find . -iname "*NG-14833_Reg_25dpa_3_lib233722_5747_6*" -exec rename 's/NG-14833_Reg_25dpa_3_lib233722_5747_6/Reg.25dpa.3/' '{}' \;

find . -iname "*NG-14833_Tit_15dpa_1_lib233705_5747_1*" -exec rename 's/NG-14833_Tit_15dpa_1_lib233705_5747_1/Tit.15dpa.1/' '{}' \;

find . -iname "*NG-14833_Tit_15dpa_3_lib233707_5747_1*" -exec rename 's/NG-14833_Tit_15dpa_3_lib233707_5747_1/Tit.15dpa.3/' '{}' \;
find . -iname "*NG-14833_Tit_25dpa_1_lib233708_5747_1*" -exec rename 's/NG-14833_Tit_25dpa_1_lib233708_5747_1/Tit.25dpa.1/' '{}' \;
find . -iname "*NG-14833_Tit_25dpa_2_lib233709_5747_1*" -exec rename 's/NG-14833_Tit_25dpa_2_lib233709_5747_1/Tit.25dpa.2/' '{}' \;



#

mv Reg.15dpa.2/Merge_1.fq.gz Reg.15dpa.2/Reg.15dpa.2-sortmerna-trimmomatic_1.fq.gz
mv Reg.15dpa.2/Merge_2.fq.gz Reg.15dpa.2/Reg.15dpa.2-sortmerna-trimmomatic_2.fq.gz

mv Reg.15dpa.3/Merge_1.fq.gz Reg.15dpa.3/Reg.15dpa.3-sortmerna-trimmomatic_1.fq.gz
mv Reg.15dpa.3/Merge_2.fq.gz Reg.15dpa.3/Reg.15dpa.3-sortmerna-trimmomatic_2.fq.gz

mv Reg.25dpa.1/Merge_1.fq.gz Reg.25dpa.1/Reg.25dpa.1-sortmerna-trimmomatic_1.fq.gz
mv Reg.25dpa.1/Merge_2.fq.gz Reg.25dpa.1/Reg.25dpa.1-sortmerna-trimmomatic_2.fq.gz

mv Tit.15dpa.2/Merge_1.fq.gz Tit.15dpa.2/Tit.15dpa.2-sortmerna-trimmomatic_1.fq.gz
mv Tit.15dpa.2/Merge_2.fq.gz Tit.15dpa.2/Tit.15dpa.2-sortmerna-trimmomatic_2.fq.gz

mv Tit.25dpa.3/Merge_1.fq.gz Tit.25dpa.3/Tit.25dpa.3-sortmerna-trimmomatic_1.fq.gz
mv Tit.25dpa.3/Merge_2.fq.gz Tit.25dpa.3/Tit.25dpa.3-sortmerna-trimmomatic_2.fq.gz


#

find . -iname "*-sortmerna-trimmomatic*" -exec rename 's/-sortmerna-trimmomatic/-sort-trim/' '{}' \;
##########################

echo "Script done all...."
