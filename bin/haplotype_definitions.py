#!/usr/bin/env python3
"""
MHC Haplotype definitions for non-human primates.
Adapted from the original miseq_genotyping_without_labkey.py Colab notebook.
"""

# Indian Rhesus (Mamu) haplotype definitions
indian_rhesus = {'PREFIX': 'Mamu'}

indian_rhesus['MHC_A_HAPLOTYPES'] = {
    'A001.01': ['A1_001'],
    'A002.01': ['A1_002_01'],
    'A003.01': ['A1_003'],
    'A004.01': ['A1_004'],
    'A005.01': ['A1_005_01'],
    'A006.01': ['A1_006_01'],
    'A007.01': ['A1_007_01'],
    'A008.01': ['A1_008_01'],
    'A009.01': ['A1_009_01'],
    'A010.01': ['A1_010_01'],
    'A011.01': ['A1_011'],
    'A012.01': ['A1_012_01'],
    'A013.01': ['A1_013'],
    'A014.01': ['A1_014_01'],
    'A015.01': ['A1_015_01'],
    'A016.01': ['A1_016_01'],
    'A017.01': ['A1_017_01'],
    'A018.01': ['A1_018_01'],
    'A019.01': ['A1_019'],
    'A020.01': ['A1_020'],
    'A021.01': ['A1_021'],
    'A022.01': ['A1_022'],
    'A023.01': ['A1_023'],
    'A024.01': ['A1_024'],
    'A025.01': ['A1_025'],
    'A026.01': ['A1_026'],
    'A027.01': ['A1_027'],
    'A028.01': ['A1_028'],
    'A029.01': ['A1_029'],
    'A030.01': ['A1_030'],
    'A031.01': ['A1_031_01'],
    'A032.01': ['A1_032_01'],
    'A033.01': ['A1_033_01'],
    'A034.01': ['A1_034_01'],
    'A035.01': ['A1_035_01'],
    'A036.01': ['A1_036_01'],
    'A037.01': ['A1_037_01'],
    'A038.01': ['A1_038_01'],
    'A039.01': ['A1_039_01'],
    'A040.01': ['A1_040_01'],
    'A041.01': ['A1_041_01'],
    'A042.01': ['A1_042_01'],
    'A043.01': ['A1_043_01'],
    'A044.01': ['A1_044_01'],
    'A045.01': ['A1_045_01'],
    'A046.01': ['A1_046_01'],
    'A047.01': ['A1_047_01'],
    'A048.01': ['A1_048_01'],
    'A049.01': ['A1_049_01'],
    'A050.01': ['A1_050_01'],
    'A051.01': ['A1_051_01'],
    'A052.01': ['A1_052_01'],
    'A053.01': ['A1_053_01'],
    'A054.01': ['A1_054_01'],
    'A055.01': ['A1_055_01'],
    'A056.01': ['A1_056_01'],
    'A057.01': ['A1_057_01'],
    'A058.01': ['A1_058_01'],
    'A059.01': ['A1_059_01'],
    'A060.01': ['A1_060_01'],
    'A061.01': ['A1_061_01'],
    'A062.01': ['A1_062_01'],
    'A063.01': ['A1_063_01'],
    'A064.01': ['A1_064_01'],
    'A065.01': ['A1_065_01'],
    'A066.01': ['A1_066_01'],
    'A067.01': ['A1_067_01'],
    'A068.01': ['A1_068_01'],
    'A069.01': ['A1_069_01'],
    'A070.01': ['A1_070_01']
}

indian_rhesus['MHC_B_HAPLOTYPES'] = {
    'B001.01': ['B_001', 'B_007', 'B_030'],
    'B001.03': ['B_001_02', 'B_094', 'B_095'],
    'B002.01': ['B_002'],
    'B003.01': ['B_003', 'B_029'],
    'B004.01': ['B_004', 'B_006'],
    'B005.01': ['B_005', 'B_031'],
    'B008.01': ['B_008', 'B_032'],
    'B009.01': ['B_009', 'B_036'],
    'B010.01': ['B_010', 'B_046'],
    'B011.01': ['B_011', 'B_045'],
    'B012.01': ['B_012', 'B_049'],
    'B013.01': ['B_013', 'B_028'],
    'B014.01': ['B_014', 'B_023'],
    'B015.01': ['B_015'],
    'B016.01': ['B_016'],
    'B017.01': ['B_017', 'B_047'],
    'B018.01': ['B_018'],
    'B019.01': ['B_019'],
    'B020.01': ['B_020'],
    'B021.01': ['B_021'],
    'B022.01': ['B_022'],
    'B024.01': ['B_024'],
    'B025.01': ['B_025'],
    'B026.01': ['B_026'],
    'B027.01': ['B_027'],
    'B033.01': ['B_033'],
    'B034.01': ['B_034'],
    'B035.01': ['B_035'],
    'B037.01': ['B_037'],
    'B038.01': ['B_038'],
    'B039.01': ['B_039'],
    'B040.01': ['B_040'],
    'B041.01': ['B_041'],
    'B042.01': ['B_042'],
    'B043.01': ['B_043'],
    'B044.01': ['B_044'],
    'B048.01': ['B_048'],
    'B050.01': ['B_050'],
    'B051.01': ['B_051'],
    'B052.01': ['B_052'],
    'B053.01': ['B_053'],
    'B054.01': ['B_054'],
    'B055.01': ['B_055'],
    'B056.01': ['B_056'],
    'B057.01': ['B_057'],
    'B058.01': ['B_058'],
    'B059.01': ['B_059'],
    'B060.01': ['B_060'],
    'B061.01': ['B_061'],
    'B062.01': ['B_062'],
    'B063.01': ['B_063'],
    'B064.01': ['B_064'],
    'B065.01': ['B_065'],
    'B066.01': ['B_066'],
    'B067.01': ['B_067'],
    'B068.01': ['B_068'],
    'B069.01': ['B_069'],
    'B070.01': ['B_070'],
    'B071.01': ['B_071'],
    'B072.01': ['B_072'],
    'B073.01': ['B_073'],
    'B074.01': ['B_074'],
    'B075.01': ['B_075'],
    'B076.01': ['B_076'],
    'B077.01': ['B_077'],
    'B078.01': ['B_078'],
    'B079.01': ['B_079'],
    'B080.01': ['B_080'],
    'B081.01': ['B_081'],
    'B082.01': ['B_082'],
    'B083.01': ['B_083'],
    'B084.01': ['B_084'],
    'B085.01': ['B_085'],
    'B086.01': ['B_086'],
    'B087.01': ['B_087'],
    'B088.01': ['B_088'],
    'B089.01': ['B_089'],
    'B090.01': ['B_090'],
    'B091.01': ['B_091'],
    'B092.01': ['B_092'],
    'B093.01': ['B_093'],
    'B096.01': ['B_096'],
    'B097.01': ['B_097'],
    'B098.01': ['B_098'],
    'B099.01': ['B_099'],
    'B100.01': ['B_100'],
    'B101.01': ['B_101'],
    'B102.01': ['B_102'],
    'B103.01': ['B_103'],
    'B104.01': ['B_104'],
    'B105.01': ['B_105'],
    'B106.01': ['B_106'],
    'B107.01': ['B_107'],
    'B108.01': ['B_108'],
    'B109.01': ['B_109'],
    'B110.01': ['B_110'],
    'B111.01': ['B_111'],
    'B112.01': ['B_112'],
    'B113.01': ['B_113'],
    'B114.01': ['B_114'],
    'B115.01': ['B_115'],
    'B116.01': ['B_116'],
    'B117.01': ['B_117'],
    'B118.01': ['B_118'],
    'B119.01': ['B_119'],
    'B120.01': ['B_120'],
    'B121.01': ['B_121'],
    'B122.01': ['B_122'],
    'B123.01': ['B_123'],
    'B124.01': ['B_124'],
    'B125.01': ['B_125'],
    'B126.01': ['B_126'],
    'B127.01': ['B_127'],
    'B128.01': ['B_128'],
    'B129.01': ['B_129'],
    'B130.01': ['B_130'],
    'B131.01': ['B_131'],
    'B132.01': ['B_132'],
    'B133.01': ['B_133'],
    'B134.01': ['B_134'],
    'B135.01': ['B_135'],
    'B136.01': ['B_136'],
    'B137.01': ['B_137'],
    'B138.01': ['B_138'],
    'B139.01': ['B_139'],
    'B140.01': ['B_140'],
    'B141.01': ['B_141'],
    'B142.01': ['B_142'],
    'B143.01': ['B_143'],
    'B144.01': ['B_144'],
    'B145.01': ['B_145'],
    'B146.01': ['B_146'],
    'B147.01': ['B_147'],
    'B148.01': ['B_148'],
    'B149.01': ['B_149'],
    'B150.01': ['B_150'],
    'B151.01': ['B_151'],
    'B152.01': ['B_152'],
    'B153.01': ['B_153'],
    'B154.01': ['B_154'],
    'B155.01': ['B_155'],
    'B156.01': ['B_156'],
    'B157.01': ['B_157'],
    'B158.01': ['B_158'],
    'B159.01': ['B_159'],
    'B160.01': ['B_160'],
    'B161.01': ['B_161'],
    'B162.01': ['B_162'],
    'B163.01': ['B_163'],
    'B164.01': ['B_164'],
    'B165.01': ['B_165'],
    'B166.01': ['B_166'],
    'B167.01': ['B_167'],
    'B168.01': ['B_168'],
    'B169.01': ['B_169']
}

indian_rhesus['MHC_DRB_HAPLOTYPES'] = {
    'DR01.01': ['DRB1_04_06_01', 'DRB5_03_01'],
    'DR01.03': ['DRB1_04_11', 'DRB5_03_09'],
    'DR02.01': ['DRB1_03_01', 'DRB5_01_01'],
    'DR02.02': ['DRB1_03_04', 'DRB5_01_03'],
    'DR03.01': ['DRB1_01_01', 'DRB3_02_01'],
    'DR03.02': ['DRB1_01_02', 'DRB3_02_01'],
    'DR04.01': ['DRB1_02_01', 'DRB4_01_01'],
    'DR04.02': ['DRB1_02_02', 'DRB4_01_01'],
    'DR04.03': ['DRB1_02_04', 'DRB4_01_01'],
    'DR05.01': ['DRB1_11_01', 'DRB3_02_02'],
    'DR06.01': ['DRB1_13_01', 'DRB3_03_01'],
    'DR06.02': ['DRB1_13_02', 'DRB3_03_01'],
    'DR07.01': ['DRB1_07_01', 'DRB4_02_01'],
    'DR08.01': ['DRB1_12_01', 'DRB3_02_03'],
    'DR09.01': ['DRB1_04_01', 'DRB5_02_01'],
    'DR10.01': ['DRB1_10_01', 'DRB4_03_01'],
    'DR11.01': ['DRB1_08_01', 'DRB4_01_02'],
    'DR12.01': ['DRB1_09_01', 'DRB4_02_02'],
    'DR13.01': ['DRB1_14_01', 'DRB3_02_04'],
    'DR14.01': ['DRB1_04_07', 'DRB5_02_02'],
    'DR15.01': ['DRB1_15_01', 'DRB5_01_04'],
    'DR16.01': ['DRB1_16_01', 'DRB5_02_03'],
    'DR17.01': ['DRB1_17_01', 'DRB3_02_05'],
    'DR18.01': ['DRB1_18_01', 'DRB3_02_06'],
    'DR19.01': ['DRB1_04_02', 'DRB5_02_04'],
    'DR20.01': ['DRB1_04_03', 'DRB5_01_02'],
    'DR21.01': ['DRB1_04_08', 'DRB5_02_05'],
    'DR22.01': ['DRB1_04_09', 'DRB5_03_02'],
    'DR23.01': ['DRB1_04_10', 'DRB5_03_03'],
    'DR24.01': ['DRB1_04_04', 'DRB5_03_04'],
    'DR25.01': ['DRB1_04_05', 'DRB5_03_05'],
    'DR26.01': ['DRB1_03_02', 'DRB5_01_05'],
    'DR27.01': ['DRB1_03_03', 'DRB5_01_06'],
    'DR28.01': ['DRB1_01_03', 'DRB3_02_07'],
    'DR29.01': ['DRB1_01_04', 'DRB3_02_08'],
    'DR30.01': ['DRB1_01_05', 'DRB3_02_09']
}

indian_rhesus['MHC_DQA_HAPLOTYPES'] = {
    '01_02': ['DQA1_01_02'],
    '01_07': ['DQA1_01_07'],
    '02_01': ['DQA1_02_01'],
    '02_02': ['DQA1_02_02'],
    '02_03': ['DQA1_02_03'],
    '02_04': ['DQA1_02_04'],
    '02_05': ['DQA1_02_05'],
    '02_06': ['DQA1_02_06'],
    '02_07': ['DQA1_02_07'],
    '03_01': ['DQA1_03_01'],
    '03_02': ['DQA1_03_02'],
    '04_01': ['DQA1_04_01'],
    '05_01': ['DQA1_05_01'],
    '06_01': ['DQA1_06_01'],
    '06_02': ['DQA1_06_02'],
    '06_03': ['DQA1_06_03'],
    '06_04': ['DQA1_06_04'],
    '07_01': ['DQA1_07_01'],
    '07_02': ['DQA1_07_02'],
    '08_01': ['DQA1_08_01'],
    '09_01': ['DQA1_09_01'],
    '10_01': ['DQA1_10_01'],
    '11_01': ['DQA1_11_01'],
    '12_01': ['DQA1_12_01'],
    '13_01': ['DQA1_13_01'],
    '14_01': ['DQA1_14_01'],
    '15_01': ['DQA1_15_01'],
    '16_01': ['DQA1_16_01'],
    '17_01': ['DQA1_17_01'],
    '18_01': ['DQA1_18_01'],
    '19_01': ['DQA1_19_01'],
    '20_01': ['DQA1_20_01'],
    '21_01': ['DQA1_21_01'],
    '22_01': ['DQA1_22_01'],
    '23_01': ['DQA1_23_01'],
    '24_01': ['DQA1_24_01'],
    '25_01': ['DQA1_25_01'],
    '26_01': ['DQA1_26_01']
}

indian_rhesus['MHC_DQB_HAPLOTYPES'] = {
    '06_01': ['DQB1_06_01'],
    '06_07': ['DQB1_06_07'],
    '06_08': ['DQB1_06_08'],
    '06_09': ['DQB1_06_09'],
    '06_10': ['DQB1_06_10'],
    '06_11': ['DQB1_06_11'],
    '06_12': ['DQB1_06_12'],
    '06_13': ['DQB1_06_13'],
    '06_14': ['DQB1_06_14'],
    '06_15': ['DQB1_06_15'],
    '06_16': ['DQB1_06_16'],
    '06_17': ['DQB1_06_17'],
    '06_18': ['DQB1_06_18'],
    '18_01': ['DQB1_18_01'],
    '18_02': ['DQB1_18_02'],
    '18_03': ['DQB1_18_03'],
    '18_04': ['DQB1_18_04'],
    '18_05': ['DQB1_18_05'],
    '18_06': ['DQB1_18_06'],
    '18_07': ['DQB1_18_07'],
    '18_08': ['DQB1_18_08'],
    '18_09': ['DQB1_18_09'],
    '18_10': ['DQB1_18_10'],
    '18_11': ['DQB1_18_11'],
    '18_12': ['DQB1_18_12'],
    '18_13': ['DQB1_18_13'],
    '18_14': ['DQB1_18_14'],
    '18_15': ['DQB1_18_15'],
    '18_16': ['DQB1_18_16'],
    '18_17': ['DQB1_18_17'],
    '18_18': ['DQB1_18_18'],
    '15_01': ['DQB1_15_01'],
    '15_02': ['DQB1_15_02'],
    '15_03': ['DQB1_15_03'],
    '15_04': ['DQB1_15_04'],
    '15_05': ['DQB1_15_05'],
    '15_06': ['DQB1_15_06'],
    '15_07': ['DQB1_15_07'],
    '15_08': ['DQB1_15_08'],
    '15_09': ['DQB1_15_09'],
    '15_10': ['DQB1_15_10'],
    '15_11': ['DQB1_15_11'],
    '15_12': ['DQB1_15_12'],
    '15_13': ['DQB1_15_13'],
    '15_14': ['DQB1_15_14'],
    '15_15': ['DQB1_15_15'],
    '08_01': ['DQB1_08_01'],
    '08_02': ['DQB1_08_02'],
    '08_03': ['DQB1_08_03'],
    '08_04': ['DQB1_08_04'],
    '08_05': ['DQB1_08_05'],
    '08_06': ['DQB1_08_06'],
    '08_07': ['DQB1_08_07'],
    '08_08': ['DQB1_08_08'],
    '08_09': ['DQB1_08_09'],
    '08_10': ['DQB1_08_10'],
    '08_11': ['DQB1_08_11'],
    '08_12': ['DQB1_08_12'],
    '08_13': ['DQB1_08_13'],
    '08_14': ['DQB1_08_14'],
    '17_01': ['DQB1_17_01'],
    '17_02': ['DQB1_17_02'],
    '17_03': ['DQB1_17_03'],
    '17_04': ['DQB1_17_04'],
    '17_05': ['DQB1_17_05'],
    '17_06': ['DQB1_17_06'],
    '17_07': ['DQB1_17_07'],
    '17_08': ['DQB1_17_08'],
    '17_09': ['DQB1_17_09'],
    '17_10': ['DQB1_17_10'],
    '17_11': ['DQB1_17_11'],
    '17_12': ['DQB1_17_12'],
    '05_01': ['DQB1_05_01'],
    '05_02': ['DQB1_05_02'],
    '05_03': ['DQB1_05_03'],
    '05_04': ['DQB1_05_04'],
    '05_05': ['DQB1_05_05'],
    '05_06': ['DQB1_05_06'],
    '05_07': ['DQB1_05_07'],
    '05_08': ['DQB1_05_08'],
    '05_09': ['DQB1_05_09'],
    '05_10': ['DQB1_05_10'],
    '05_11': ['DQB1_05_11'],
    '05_12': ['DQB1_05_12'],
    '05_13': ['DQB1_05_13'],
    '05_14': ['DQB1_05_14']
}

indian_rhesus['MHC_DPA_HAPLOTYPES'] = {
    '02_03': ['DPA1_02_03'],
    '02_08': ['DPA1_02_08'],
    '02_09': ['DPA1_02_09'],
    '02_10': ['DPA1_02_10'],
    '02_11': ['DPA1_02_11'],
    '02_12': ['DPA1_02_12'],
    '02_13': ['DPA1_02_13'],
    '02_14': ['DPA1_02_14'],
    '02_15': ['DPA1_02_15'],
    '02_16': ['DPA1_02_16'],
    '02_17': ['DPA1_02_17'],
    '02_18': ['DPA1_02_18'],
    '02_19': ['DPA1_02_19'],
    '02_20': ['DPA1_02_20'],
    '03_01': ['DPA1_03_01'],
    '03_02': ['DPA1_03_02'],
    '03_03': ['DPA1_03_03'],
    '03_04': ['DPA1_03_04'],
    '03_05': ['DPA1_03_05'],
    '04_01': ['DPA1_04_01'],
    '05_01': ['DPA1_05_01'],
    '06_01': ['DPA1_06_01'],
    '07_01': ['DPA1_07_01'],
    '08_01': ['DPA1_08_01'],
    '09_01': ['DPA1_09_01']
}

indian_rhesus['MHC_DPB_HAPLOTYPES'] = {
    '01g1': ['DPB1_01g1'],
    '01g2': ['DPB1_01g2'],
    '01g3': ['DPB1_01g3'],
    '02_01': ['DPB1_02_01'],
    '03_01': ['DPB1_03_01'],
    '04_01': ['DPB1_04_01'],
    '05_01': ['DPB1_05_01'],
    '06_01': ['DPB1_06_01'],
    '07_01': ['DPB1_07_01'],
    '08_01': ['DPB1_08_01'],
    '09_01': ['DPB1_09_01'],
    '10_01': ['DPB1_10_01'],
    '11_01': ['DPB1_11_01'],
    '12_01': ['DPB1_12_01'],
    '13_01': ['DPB1_13_01'],
    '14_01': ['DPB1_14_01'],
    '15_01': ['DPB1_15_01'],
    '16_01': ['DPB1_16_01'],
    '17_01': ['DPB1_17_01'],
    '18_01': ['DPB1_18_01'],
    '19_01': ['DPB1_19_01'],
    '20_01': ['DPB1_20_01'],
    '21_01': ['DPB1_21_01'],
    '22_01': ['DPB1_22_01'],
    '23_01': ['DPB1_23_01'],
    '24_01': ['DPB1_24_01'],
    '25_01': ['DPB1_25_01'],
    '26_01': ['DPB1_26_01'],
    '27_01': ['DPB1_27_01'],
    '28_01': ['DPB1_28_01'],
    '29_01': ['DPB1_29_01'],
    '30_01': ['DPB1_30_01']
}

# MCM (Mafa) haplotype definitions
# Updated with full FASTA headers from 31588_MCM_MHC_miseq_deduplicated_annotated.fasta
mcm = {'PREFIX': 'Mafa'}

mcm['MHC_A_HAPLOTYPES'] = {
    'M1A': [
        '01_Mafa_A1_063:01:01:01;Mafa_A1_063:02:01:01_M1M2M3',
        '02_Mafa_A2_05:46:01:01;Mafa_A2_05:04:01:01;Mafa_A2_05:11:01:01;Mafa_A2_05:01:01:01_M1M2M3M4M5',
        '02_Mafa_A4_01:09:01:01N;Mafa_A4_01:01:01:01_M1M2M3M6',
        '16_Mafa_AG1_05:36:01:01N_M1',
        '16_Mafa_AG2_01:05:01:01;Mafa_AG2_01:03:01:01;Mafa_AG2_01:07:02:01;Mafa_AG2_01:14:01:01;Mafa_AG2_01:03:03:01_M1M2M3M4',
        '16_Mafa_AG3_02:16:01:01;Mafa_AG3_02:03:02:02;Mafa_AG3_02:03:05:01_M1M2M4',
        '16_Mafa_AG5_06:05:01:01;Mafa_AG5_07_M4nov;Mafa_AG5_07:03:02:01;Mafa_AG5_07:18:01:01_M1M2M7',
        '16_Mafa_AG6_04:05:03:03;Mafa_AG6_04:05:03:02;Mafa_AG6_04:05:03:01_M1M3M5M6',
        '16_Mafa_AG6_04:08:01:01_M1M3'
    ],
    'M2A': [
        '01_Mafa_A1_063:01:01:01;Mafa_A1_063:02:01:01_M1M2M3',
        '02_Mafa_A2_05:46:01:01;Mafa_A2_05:04:01:01;Mafa_A2_05:11:01:01;Mafa_A2_05:01:01:01_M1M2M3M4M5',
        '02_Mafa_A4_01:09:01:01N;Mafa_A4_01:01:01:01_M1M2M3M6',
        '16_Mafa_AG2_01:05:01:01;Mafa_AG2_01:03:01:01;Mafa_AG2_01:07:02:01;Mafa_AG2_01:14:01:01;Mafa_AG2_01:03:03:01_M1M2M3M4',
        '16_Mafa_AG3_02:16:01:01;Mafa_AG3_02:03:02:02;Mafa_AG3_02:03:05:01_M1M2M4',
        '16_Mafa_AG5_06:04:01:01;Mafa_AG1_05:02:03:01_M2M3',
        '16_Mafa_AG5_06:05:01:01;Mafa_AG5_07_M4nov;Mafa_AG5_07:03:02:01;Mafa_AG5_07:18:01:01_M1M2M7'
    ],
    'M3A': [
        '01_Mafa_A1_063:01:01:01;Mafa_A1_063:02:01:01_M1M2M3',
        '02_Mafa_A2_05:46:01:01;Mafa_A2_05:04:01:01;Mafa_A2_05:11:01:01;Mafa_A2_05:01:01:01_M1M2M3M4M5',
        '02_Mafa_A4_01:09:01:01N;Mafa_A4_01:01:01:01_M1M2M3M6',
        '16_Mafa_AG1_05:09:02:01_M3',
        '16_Mafa_AG2_01:05:01:01;Mafa_AG2_01:03:01:01;Mafa_AG2_01:07:02:01;Mafa_AG2_01:14:01:01;Mafa_AG2_01:03:03:01_M1M2M3M4',
        '16_Mafa_AG3_03:03:02:01;Mafa_AG3_03:01:04:01_M3M4',
        '16_Mafa_AG5_06:04:01:01;Mafa_AG1_05:02:03:01_M2M3',
        '16_Mafa_AG6_04:05:03:03;Mafa_AG6_04:05:03:02;Mafa_AG6_04:05:03:01_M1M3M5M6',
        '16_Mafa_AG6_04:08:01:01_M1M3'
    ],
    'M4A': [
        '01_Mafa_A1_031:01:01:01_M4',
        '02_Mafa_A2_05:46:01:01;Mafa_A2_05:04:01:01;Mafa_A2_05:11:01:01;Mafa_A2_05:01:01:01_M1M2M3M4M5',
        '02_Mafa_A5_30:01:01_nov;Mafa_A5_30:01:01:01_M4',
        '16_Mafa_AG1_05:35:01:01;Mafa_AG1_05:11:01:02_M4',
        '16_Mafa_AG2_01:05:01:01;Mafa_AG2_01:03:01:01;Mafa_AG2_01:07:02:01;Mafa_AG2_01:14:01:01;Mafa_AG2_01:03:03:01_M1M2M3M4',
        '16_Mafa_AG3_02:16:01:01;Mafa_AG3_02:03:02:02;Mafa_AG3_02:03:05:01_M1M2M4',
        '16_Mafa_AG3_03:03:02:01;Mafa_AG3_03:01:04:01_M3M4',
        '16_Mafa_AG5_06:05:01:01;Mafa_AG5_07_M4nov;Mafa_AG5_07:03:02:01;Mafa_AG5_07:18:01:01_M1M2M7',
        '16_Mafa_AG6_04:27:01:01N;Mafa_AG6_04:26:01:01_M4M7'
    ],
    'M5A': [
        '01_Mafa_A1_033:01:01:01_M5',
        '02_Mafa_A2_05:46:01:01;Mafa_A2_05:04:01:01;Mafa_A2_05:11:01:01;Mafa_A2_05:01:01:01_M1M2M3M4M5',
        '16_Mafa_AG1_05:23:02:01;Mafa_AG1_05:10:02:02_M5M6',
        '16_Mafa_AG3_03_M5nov',
        '16_Mafa_AG5_06_M5nov',
        '16_Mafa_AG6_04:05:03:03;Mafa_AG6_04:05:03:02;Mafa_AG6_04:05:03:01_M1M3M5M6'
    ],
    'M6A': [
        '01_Mafa_A1_032:01:01:01_M6',
        '01_Mafa_A1_047:02:01:01_M6',
        '02_Mafa_A4_01:09:01:01N;Mafa_A4_01:01:01:01_M1M2M3M6',
        '16_Mafa_AG1_05:23:02:01;Mafa_AG1_05:10:02:02_M5M6',
        '16_Mafa_AG3_02:06:01:04_M6',
        '16_Mafa_AG5_06:02:03:03;Mafa_AG2_01:15:01:01N_M6M7',
        '16_Mafa_AG6_04:05:03:03;Mafa_AG6_04:05:03:02;Mafa_AG6_04:05:03:01_M1M3M5M6'
    ],
    'M7A': [
        '01_Mafa_A1_060:05:01:01_M7',
        '16_Mafa_AG1_05:34:01:01_M7',
        '16_Mafa_AG3_02:17:01:01_M7',
        '16_Mafa_AG3_02:18:01:01_M7',
        '16_Mafa_AG5_06:02:03:03;Mafa_AG2_01:15:01:01N_M6M7',
        '16_Mafa_AG5_06:05:01:01;Mafa_AG5_07_M4nov;Mafa_AG5_07:03:02:01;Mafa_AG5_07:18:01:01_M1M2M7',
        '16_Mafa_AG6_04:25:01:01_M7',
        '16_Mafa_AG6_04:27:01:01N;Mafa_AG6_04:26:01:01_M4M7'
    ]
}

mcm['MHC_B_HAPLOTYPES'] = {
    'M1B': [
        '03_Mafa_B11L_01_M1nov',
        '03_Mafa_B17_01_M1nov;Mafa_B17_01_M5nov;Mafa_B17_01_M7nov;Mafa_B17_01_M4nov;Mafa_B17_01:05:01:01_M3M4',
        '03_Mafa_B19Ps_01_M1nov'
    ],
    'M2B': [
        '03_Mafa_B_046:01:01:01_M2',
        '03_Mafa_B_057:02:01:01_M2',
        '03_Mafa_B_060:05:01:01;Mafa_B_060:31:01:01;Mafa_B_060:06:01:01_M2M5',
        '03_Mafa_B_064:01:01:01_M2',
        '03_Mafa_B_104:01:01:01_M2',
        '03_Mafa_B_109:13:01:01N_M2',
        '03_Mafa_B_134:02:01:01_M2',
        '03_Mafa_B_144:02:01:01_M2',
        '03_Mafa_B11L_01_M5nov;Mafa_B11L_01_M2nov',
        '03_Mafa_B16_01:13:01:01;Mafa_B16_01:15:01:01_ext_M2'
    ],
    'M3B': [
        '03_Mafa_B_019:03:01:01_M3',
        '03_Mafa_B_079:03:01:01;Mafa_B_079:01:01:01_M3M4',
        '03_Mafa_B_098:01:01:01;Mafa_B_098:04:01:03_M3M6',
        '03_Mafa_B_109:04:01:01_M3',
        '03_Mafa_B_109:31:01:01_M3',
        '03_Mafa_B_148:01:01:01_M3',
        '03_Mafa_B_150:01:01:01_M3',
        '03_Mafa_B_162:03:01:01_M3',
        '03_Mafa_B17_01_M1nov;Mafa_B17_01_M5nov;Mafa_B17_01_M7nov;Mafa_B17_01_M4nov;Mafa_B17_01:05:01:01_M3M4',
        '03_Mafa_B22_01:01:01:01N;Mafa_B22_01:01:02:01N_M3M4'
    ],
    'M4B': [
        '03_Mafa_B_011:01:01:01_M4',
        '03_Mafa_B_070:02:01:01_M4',
        '03_Mafa_B_075:01:01:01_M4',
        '03_Mafa_B_079:03:01:01;Mafa_B_079:01:01:01_M3M4',
        '03_Mafa_B_098:05:01:01_M4',
        '03_Mafa_B_109:02:01:01_M4',
        '03_Mafa_B_207:01:01:01N_M4',
        '03_Mafa_B02Ps_01:03:01:01_M4',
        '03_Mafa_B11L_01:06:01:01N;Mafa_B11L_01_M6nov03_M4',
        '03_Mafa_B11L_01_M4nov',
        '03_Mafa_B17_01_M1nov;Mafa_B17_01_M5nov;Mafa_B17_01_M7nov;Mafa_B17_01_M4nov;Mafa_B17_01:05:01:01_M3M4',
        '03_Mafa_B19Ps_01_M4nov',
        '03_Mafa_B22_01:01:01:01N;Mafa_B22_01:01:02:01N_M3M4'
    ],
    'M5B': [
        '03_Mafa_B_027:02:01:01_ext_M5',
        '03_Mafa_B_051:03:01:01_M5',
        '03_Mafa_B_060:05:01:01;Mafa_B_060:31:01:01;Mafa_B_060:06:01:01_M2M5',
        '03_Mafa_B_088:01:01:01_M5',
        '03_Mafa_B_147:01:01:01_M5',
        '03_Mafa_B02Ps_01:11:02:01_M5',
        '03_Mafa_B11L_01_M5nov;Mafa_B11L_01_M2nov',
        '03_Mafa_B16_01:04:01:01_ext;Mafa_B16_01:18:01:01N_M5M6',
        '03_Mafa_B17_01_M1nov;Mafa_B17_01_M5nov;Mafa_B17_01_M7nov;Mafa_B17_01_M4nov;Mafa_B17_01:05:01:01_M3M4',
        '03_Mafa_B19Ps_01_M5nov',
        '03_Mafa_B21Ps_01:05:02:01_M5'
    ],
    'M6B': [
        '03_Mafa_B_033:01:01:01_M6',
        '03_Mafa_B_036:01:02:01_M6',
        '03_Mafa_B_037:01:01:01_M6',
        '03_Mafa_B_045:01:01:01_M6',
        '03_Mafa_B_045:03:01:01_M6',
        '03_Mafa_B_046:09:01:01_M6',
        '03_Mafa_B_050:04:01:01_M6',
        '03_Mafa_B_051:33:01:01_M6',
        '03_Mafa_B_060:01:01:01;Mafa_B_060:04:01:01_M6',
        '03_Mafa_B_095:01:01:01_M6',
        '03_Mafa_B_097:07:01:01_M6',
        '03_Mafa_B_098:01:01:01;Mafa_B_098:04:01:03_M3M6',
        '03_Mafa_B_098:06:01:01_M6',
        '03_Mafa_B_149:01:01:01_M6',
        '03_Mafa_B_151:01:01:01_M6',
        '03_Mafa_B_167:01:01:01N_M6',
        '03_Mafa_B11L_01:06:01:01N;Mafa_B11L_01_M6nov03_M4',
        '03_Mafa_B11L_01_M6nov01',
        '03_Mafa_B11L_01_M6nov02',
        '03_Mafa_B14Ps_175:01:01:01_M6',
        '03_Mafa_B16_01:04:01:01_ext;Mafa_B16_01:18:01:01N_M5M6',
        '03_Mafa_B16_01:14:01:01_M6',
        '03_Mafa_B16_01:16:01:01_ext_M6',
        '03_Mafa_B17_01_M6nov',
        '03_Mafa_B19Ps_01_M6nov',
        '03_Mafa_B22_01:02:01:01_M6'
    ],
    'M7B': [
        '03_Mafa_B11L_01_M7nov',
        '03_Mafa_B17_01_M1nov;Mafa_B17_01_M5nov;Mafa_B17_01_M7nov;Mafa_B17_01_M4nov;Mafa_B17_01:05:01:01_M3M4'
    ]
}

mcm['MHC_DRB_HAPLOTYPES'] = {
    'M1DR': ['09_Mafa_DRB9_01_M1nov'],
    'M2DR': ['05_Mafa_DRB1_10:01:01:01_M2', '09_Mafa_DRB9_01_M2nov'],
    'M3DR': ['05_Mafa_DRB1_10:02:01:01_M3', '09_Mafa_DRB6_01:09:01:01_M3', '09_Mafa_DRB9_01:01_M3'],
    'M4DR': ['05_Mafa_DRB1_04N_M4nov;Mafa_DRB1_04N_M5nov', '08_Mafa_DRB5_03:01:01:02;Mafa_DRB5_03:01:01:01_M4M5', '09_Mafa_DRB9_01_M4nov'],
    'M5DR': ['05_Mafa_DRB1_04N_M4nov;Mafa_DRB1_04N_M5nov', '08_Mafa_DRB5_03:01:01:02;Mafa_DRB5_03:01:01:01_M4M5', '09_Mafa_DRB6_01:13:02_M5nov;Mafa_DRB6_01:13:02_ext', '09_Mafa_DRB9_01_M5nov'],
    'M6DR': [],
    'M7DR': ['05_Mafa_DRB1_03:02:01:01_M7']
}

mcm['MHC_DQA_HAPLOTYPES'] = {
    'M1DQA': ['10_Mafa_DQA1_24:03:01:02_M1'],
    'M2DQA': ['10_Mafa_DQA1_01:04:01:01_M2'],
    'M3DQA': ['10_Mafa_DQA1_05:03:01:01_M3'],
    'M4DQA': [],
    'M5DQA': [],
    'M6DQA': [],
    'M7DQA': ['10_Mafa_DQA1_23:01:01:01_M7']
}

mcm['MHC_DQB_HAPLOTYPES'] = {
    'M1DQB': ['11_Mafa_DQB1_18:01:01:01_M1'],
    'M2DQB': ['11_Mafa_DQB1_06:01:01:01;Mafa_DQB1_06:01:02_ext_M2'],
    'M3DQB': ['11_Mafa_DQB1_16:01:01:01_M3'],
    'M4DQB': ['11_Mafa_DQB1_06:08:01:01_M4'],
    'M5DQB': ['11_Mafa_DQB1_06:11:01:01_M5'],
    'M6DQB': [],
    'M7DQB': ['11_Mafa_DQB1_18:14:01:01_M7']
}

mcm['MHC_DPA_HAPLOTYPES'] = {
    'M1DPA': ['12_Mafa_DPA1_07:02:01:01_M1'],
    'M2DPA': ['12_Mafa_DPA1_07:01:01:01_M2'],
    'M3DPA': ['12_Mafa_DPA1_13:01:01:01;Mafa_DPA1_13:03_ext_M3'],
    'M4DPA': ['12_Mafa_DPA1_04:01:01:01;Mafa_DPA1_04:01_M7nov_M4'],
    'M5DPA': [],
    'M6DPA': [],
    'M7DPA': ['12_Mafa_DPA1_04:01:01:01;Mafa_DPA1_04:01_M7nov_M4']
}

mcm['MHC_DPB_HAPLOTYPES'] = {
    'M1DPB': ['13_Mafa_DPB2_01_M1nov'],
    'M2DPB': ['13_Mafa_DPB2_01_M2nov'],
    'M3DPB': ['13_Mafa_DPB1_09:02:01:01_M3', '13_Mafa_DPB2_01:01:01:01N_M3'],
    'M4DPB': ['13_Mafa_DPB1_12:02N_M4nov;Mafa_DPB1_12:02N_M7nov', '13_Mafa_DPB2_01_M4nov;Mafa_DPB2_01_M7nov'],
    'M5DPB': ['13_Mafa_DPB2_01_M5nov;Mafa_DPB2_01_M6nov'],
    'M6DPB': ['13_Mafa_DPB1_04:01_ext;Mafa_DPB1_04_M6nov', '13_Mafa_DPB2_01_M5nov;Mafa_DPB2_01_M6nov'],
    'M7DPB': ['13_Mafa_DPB1_03:03_M7nov;Mafa_DPB1_03:03_ext', '13_Mafa_DPB1_12:02N_M4nov;Mafa_DPB1_12:02N_M7nov', '13_Mafa_DPB2_01_M4nov;Mafa_DPB2_01_M7nov']
}

# Create master haplotype dictionary
haplotype_dict = {
    'MAMU': indian_rhesus,
    'MCM': mcm,
    'MANE': indian_rhesus,  # Pig-tailed uses same as Indian Rhesus
    'Rhesus': indian_rhesus,
    'Cynomolgus': mcm,
    'Pig-tailed': indian_rhesus
}

def get_species_haplotypes(species):
    """
    Get haplotype definitions for a given species.
    
    Args:
        species (str): Species name ('Rhesus', 'Cynomolgus', 'Pig-tailed')
        
    Returns:
        dict: Haplotype definitions for the species
    """
    return haplotype_dict.get(species, indian_rhesus)  # Default to Rhesus

def call_haplotypes(locus, locus_haplotype_definitions, df):
    """
    Specify locus (e.g., Mamu-A) to haplotype and provide dictionary of haplotype definitions.
    Update specified pandas dataframe.
    
    Args:
        locus (str): MHC locus (e.g., 'Mamu-A', 'Mafa-B')
        locus_haplotype_definitions (dict): Haplotype definitions for the locus
        df (pd.DataFrame): Dataframe containing allele information
        
    Returns:
        tuple: (haplotype1, haplotype2)
    """
    # Convert alleles in dataframe to string
    # Makes it possible to search for values more easily
    allele_str = df['allele'].to_string(header=False, index=False)
    
    # Create list to store haplotypes for a sample
    sample_haplotypes = []
    
    # Loop through haplotypes
    for haplotype, alleles in locus_haplotype_definitions.items():
        # If all diagnostic alleles for a haplotype are present
        # Save haplotype name to list
        if all((x) in allele_str for x in alleles):
            sample_haplotypes.append(haplotype)
    
    # Evaluate haplotypes
    if len(sample_haplotypes) == 0:  # If there are no haplotypes, that isn't possible
        # Special exception for MCM A1*063 which is often undercalled but is very important
        if locus == 'Mafa-A' and '05_M1M2M3_A1_063g' in allele_str:
            h1 = 'A1_063'
            h2 = '-'
        else:
            h1 = 'ERR: NO HAP'
            h2 = 'ERR: NO HAP'
    elif len(sample_haplotypes) == 1:
        # Special exception for MCM A1*063 which is often undercalled but is very important
        # If A1_063 is present in genotypes but M1, M2, and M3 are not present in single called haplotype, add A1_063 to haplotype2
        if locus == 'Mafa-A' and '05_M1M2M3_A1_063g' in allele_str and not any(y in sample_haplotypes[0] for y in ('M1A', 'M2A', 'M3A')):
            h1 = sample_haplotypes[0]
            h2 = 'A1_063'
        else:
            h1 = sample_haplotypes[0]
            h2 = '-'
    elif len(sample_haplotypes) == 2:
        h1 = sample_haplotypes[0]
        h2 = sample_haplotypes[1]
    elif len(sample_haplotypes) > 2:
        h1 = 'ERR: TMH (' + ', '.join(sample_haplotypes) + ')'
        h2 = 'ERR: TMH (' + ', '.join(sample_haplotypes) + ')'

    # DPA, DPB, DQA, DQB can only have two genotypes though other loci can have more
    # For these loci, error if more than two genotypes are reported

    # Test number of rows that match locus, if > 2 set h1 and h2 to error
    # Need to explicitly cast allele as string
    diploid_loci = ['DPA', 'DPB', 'DQA', 'DQB']

    for i in diploid_loci:
        if (i in locus):
            if df.allele.astype(str).str.contains(i).sum() > 2:
                h1 = 'ERR: TMG'
                h2 = 'ERR: TMG'

    return (h1, h2)