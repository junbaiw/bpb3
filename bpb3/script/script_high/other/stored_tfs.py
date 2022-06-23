#this is a python version of check accuracy for 67 snps.
#that were used in paper https://pubmed.ncbi.nlm.nih.gov/31001324/
#and paper https://pubmed.ncbi.nlm.nih.gov/26202972/
#TF to its name alianse
tfNameMap={}
tfNameMap['NFKB']=('NFKB', "NF-kappaB",'NFKB1','NFKB2')
tfNameMap['NRF2']=('NRF2', 'NFE2L2', 'NFE2', 'NRF-2','NFE2L1')
tfNameMap['NFY']=('NFY', "NF-Y",'NFYA','NFYB','NFYC')
tfNameMap['YY1']=('YY1',)
tfNameMap['USF1']=('USF','USF2')
tfNameMap['USF-1']=('USF','USF1','USF2')
tfNameMap['USF2']=('USF','USF1')
tfNameMap['c-MYB']=('MYB','MYBL1','MYBL2')
tfNameMap['c-MYB']=('MYB','MYBL1','MYBL2')
tfNameMap['C-MYB']=('MYB','MYBL1','MYBL2')
tfNameMap['IRF1']=('IRF1','IRF-1','IRF2','IRF3','IRF4','IRF5','IRF6','IRF7','IRF8','IRF9')
tfNameMap['HLF']=('HLF',)
tfNameMap['P53']=('P53','TP53')
tfNameMap['p53']=('P53','TP53')
tfNameMap['HNF-1']=('HNF1','HNF1B','HNF1A')
tfNameMap['HNF.1']=('HNF1','HNF1A','HNF1B')
tfNameMap['USF']=('USF','USF1','USF2')
tfNameMap['SPI1']=('SPI1','SPIB','SPIC')
tfNameMap['PU.1']=('SPI1','SPIB','SPIC')
tfNameMap['beta-catenin']=('CTNNB1','CTNNB','EVR7')

tfNameMap['AP-2alpha']=('AP-2','AP2','TFAP2','AP2TF','TFAP2A','TFAP2B','TFAP2C','TFAP2D','TFAP2E')
tfNameMap['AP-2-alpha']=('AP-2','AP2','TFAP2','AP2TF','TFAP2A','TFAP2B','TFAP2C','TFAP2D','TFAP2E')
tfNameMap['AP.2alpha']=('AP-2','AP2','TFAP2','AP2TF','TFAP2A','TFAP2B','TFAP2C','TFAP2D','TFAP2E')
tfNameMap['AP-2gamma']=('AP-2','AP2','TFAP2','AP2TF','TFAP2C','TFAP2B','TFAP2C','TFAP2D','TFAP2E')
tfNameMap['AP-2-ALPHA']=('AP-2','AP2','TFAP2','AP2TF','TFAP2A','TFAP2B','TFAP2C','TFAP2D','TFAP2E')
tfNameMap['AP-2ALPHA']=('AP-2','AP2','TFAP2','AP2TF','TFAP2A','TFAP2B','TFAP2C','TFAP2D','TFAP2E')
tfNameMap['CCAAT']=('CRP1', 'CCAAT', 'CEBPZ', 'CEBPA', 'CEBPG', 'CEBPB', 'CEBP', "C/EBP",'CEBPE')
tfNameMap['CACCC']=('CACCC','EKLF','KLF1','KLF-1','KLF','KLF4','KLF3','KLF8','KLF6','KLF5','KLF9','KLF13','KLF12','KLF14','KLF15','KLF16')
tfNameMap['SP1']=('SP1','SP-1','TFSP1')
tfNameMap['SP-1']=('SP1','SP-1','TFSP1')
tfNameMap['Sp1']=('SP1','SP-1','TFSP1')
tfNameMap['MSX1']=('MSX1','MSX-1','HOX7','ECTD3','HYD1','STHAG1','MSX2')
tfNameMap['GATA-1']=('GATA','GATA1','GATA2','GATA3','GATA4','GATA5','GATA6')
tfNameMap['GATA-A']=('GATA','GATA1','GATA2','GATA3','GATA4','GATA5','GATA6')
tfNameMap['GATA.1']=('GATA','GATA1','GATA2','GATA3','GATA4','GATA5','GATA6')
tfNameMap['GATA.A']=('GATA','GATA1','GATA2','GATA3','GATA4','GATA5','GATA6')
tfNameMap['HNF-4']=('HNF4a','HNF4g','HNF4','TCF14','TCF-14','MODY','MODY1','HNF4A','HNF4G')
tfNameMap['HNF.4']=('HNF4a','HNF4g','HNF4','TCF14','TCF-14','MODY','MODY1','HNF4A','HNF4G')
tfNameMap['C/EBP']=('CRP1', 'CEBP', 'CCAAT', "C/EBP",'CEBPA','CEBPB','CEBPD','CEBPE','CEBPG')
tfNameMap['C.EBP']=('CRP1', 'CEBP', 'CCAAT', "C/EBP",'CEBPA','CEBPB','CEBPD','CEBPE','CEBPG')
tfNameMap['CEBP']=('CRP1', 'CEBP', 'CCAAT', "C/EBP",'CEBPA','CEBPB','CEBPD','CEBPE','CEBPG')
tfNameMap['CEBP-delta']=('CRP1', 'CEBP', 'CCAAT', "C/EBP",'CEBPD','CEBPA','CEBPB','CEBPE','CEBPG')
tfNameMap['C/EBP-delta']=('CRP1', 'CEBP', 'CCAAT', "C/EBP",'CEBPA','CEBPB','CEBPD','CEBPE','CEBPG')
tfNameMap['CEBP-DELTA']=('CRP1', 'CEBP', 'CCAAT', "C/EBP",'CEBPD','CEBPA','CEBPB','CEBPE','CEBPG')
tfNameMap['SIX3']=('SIX3','HPE2','SIX','SIX-3','HPE-2','SIX1','SIX2','SIX4','SIX5','SIX6')
tfNameMap['Six3']=('SIX3','HPE2','SIX','SIX-3','HPE-2','SIX1','SIX2','SIX4','SIX5','SIX6')

tfNameMap['AR']=('_AR_', "^ar$",'AR')
tfNameMap['BCL11A']=('BCL11A',)
tfNameMap['c-Fos']=('AP1_', 'AP-1','FOS','JUND','JUN','JUNB')
tfNameMap['c-FOS']=('AP1_', 'AP-1','FOS','JUND','JUN','JUNB')
tfNameMap['C-FOS']=('AP1_', 'AP-1','FOS','JUND','JUN','JUNB')
tfNameMap['JunD']=('AP1_','JUND','JUN','JUNB')
tfNameMap['c-Myc']=('MYC','MAX')
tfNameMap['Max']=('MYC',)
tfNameMap['c-Rel']=('REL_', 'c-Rel','REL','RELA','RELB')
tfNameMap['C-REL']=('REL_', 'c-Rel','REL','RELA','RELB')
#not here ccnt2
tfNameMap['CCNT2']=('CCNT2',)
tfNameMap['E2F6']=('E2F','E2F2','E2F3','E2F4','E2F5','E2F6','E2F7','E2F8','E2F9')
tfNameMap['HA-E2F1']=('E2F','E2F1','E2F2','E2F3','E2F4','E2F5','E2F6','E2F7','E2F8','E2F9')
tfNameMap['EBF']=('EBF','EBF1')
tfNameMap['EBF1']=('EBF',)
tfNameMap['Egr1']=('EGR1', "Egr-1",'EGR3')
tfNameMap['EGR1']=('EGR1', "Egr-1",'EGR3')
tfNameMap['ELF1']=('ELF1','ELF2','ELF3','ELF4','ELF5')
tfNameMap['ETS']=('ETS','ETS1','ETS2','ETS3','ETV1','ELK1','ERG','ETV4','ETV7','ETV6','ELK3','ELK4')
tfNameMap['FOXA1']=('FOXA','FOXA1','FOXA2','FOXA3')
tfNameMap['GABP']=('GABP','GABPA')
tfNameMap['HEY1']=('HEY1','HEY2')
tfNameMap['MEF2A']=('MEF2','MEF2B','MEF2C','MEF2D')
#not here myf5
tfNameMap['MYF-5']=('MYF','MYF6')
tfNameMap['Oct-1']=('POU2F2', "OCT",'POU2F1','POU2F3')
tfNameMap['Oct-2']=('POU2F2', "OCT",'POU2F1','POU2F3')
tfNameMap['OCT-1']=('POU2F2', "OCT",'POU2F1','POU2F3')
tfNameMap['POU2F2']=('POU2F2', "OCT",'POU2F1','POU2F3')
#not here p300
tfNameMap['p300']=('P300','EP300')
tfNameMap['PAX5-C20']=('PAX5',)
tfNameMap['RXRA']=('RXRA',)
tfNameMap['SOX10']=('SOX10',)
tfNameMap['SRY']=('SRY',)
tfNameMap['STAT1']=('STAT','STAT2','STAT3','STAT4','STAT5A','STAT5B','STAT6')
tfNameMap['STAT3']=('STAT','STAT2','STAT3','STAT4','STAT5A','STAT5B','STAT6')
#not here TATA
tfNameMap['TAF1']=('TATA','TBP')
tfNameMap['TBP']=('TATA', "TBP")
tfNameMap['TCF4']=('TCF4', "TCF-4",'TCF3','TCF7')
tfNameMap['TEL2']=('ETV7','ETV1','ETV2','ETV3','ETV4','ETV5','ETV6','ETS','ETV7')
tfNameMap['TR4']=('NR2C2',)
tfNameMap['ZBTB7A']=('ZBTB7A','ZBTB7B','ZBTB7C')


#67 snps records
fredrikssonMuts = ("TERT_ETS|GAIN_1", "TERT_ETS|GAIN_2", "TERT_ETS|GAIN_3", "TERT_ETS|GAIN_4")

hgmdMuts = ("ADAMTS18_TEL2|GAIN", "CDKN2B-AS1_STAT1|LOSS", "CHGB_SRY|LOSS_YY1|LOSS", 
             "CHGB_c-FOS|LOSS", "CHRNA3_Oct-1|LOSS", "FABP4_C/EBP|LOSS", "FAS_c-Rel|GAIN", 
             "FMO1_YY1|LOSS", "FOXE1_MYF-5|LOSS", "FOXE1_USF1|GAIN_USF2|GAIN", 
             "FZD1_Egr1|GAIN", "GFI1B_Oct-1|LOSS", "GFI1B_GATA-1|LOSS", "GSTM1_AP-2-alpha|LOSS", 
             "HBG2_GATA-1|LOSS", "IGF1_C/EBP-delta|LOSS", 
             "IL18_OCT-1|GAIN", "IL4_Oct-1|LOSS", "LEP_Sp1|GAIN", "MPZ_SOX10|LOSS", 
             "MYC_TCF4|GAIN_beta-catenin|GAIN", "PVT1_YY1|LOSS", 
             "SIRT1_p53|LOSS", "SLC47A1_Sp1|LOSS", "SLC7A1_Sp1|LOSS", "SORT1_C/EBP|GAIN", 
             "TMPRSS2_AR|LOSS", "UGT1A7_TBP|LOSS", "UGT2B17_FOXA1|GAIN")

andersonEpsteinMuts = ("AFP_HNF-1_1|GAIN", "AFP_HNF-1_2|GAIN", "AGTRL1_SP1|GAIN", 
                        "ALOX15_SPI1|GAIN", "CETP_SP1|UNKNOWN", "COL1A1_Sp1|UNKNOWN", 
                        "F9_C/EBP|LOSS", "F9_HLF|UNKNOWN", "F9_HNF-4|UNKNOWN", "FECH_SP1|LOSS", 
                        "FLT1_P53|GAIN", "FTH1_NFY|UNKNOWN", "GP1BB_GATA-1|LOSS", "GPD2_NRF2|UNKNOWN", 
                        "HBB_CACCC_1|UNKNOWN", "HBB_CACCC_2|UNKNOWN", "HBB_CACCC_3|UNKNOWN", 
                        "HBB_CACCC_4|UNKNOWN", "HBB_CACCC_5|UNKNOWN", "HBB_CCAAT|UNKNOWN", 
                        "HBM_GATA-1|GAIN", "HOXB7_USF1|UNKNOWN", "IRF2_IRF1|LOSS", "IRF6_AP-2alpha|LOSS", 
                        "ITGA2_SP1|LOSS", "LIPC_USF|LOSS", "NFKBIL1_USF1|LOSS", "OPRM1_NFKB|UNKNOWN", 
                        "PKLR_GATA-A|LOSS", "PTGS2_c-MYB|GAIN", "SFTPB_SP-1|GAIN", "SOX9_MSX1|UNKNOWN", 
                        "SP1_NFY|UNKNOWN", "TCOf1_YY1|LOSS")


