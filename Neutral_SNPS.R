##############################
## Determining SNP overlap between outlier analyses
##
## Matt Brachmann (PhDMattyB)
##
## 2019-07-05
##
##############################

setwd('~/PhD/SNP Demographic modelling/Outliers_directory/')

library(tidyverse)
library(wesanderson)
library(patchwork)
library(janitor)
library(devtools)
library(skimr)
library(tvthemes)
library(rsed)
library(data.table)

theme_set(theme_bw())

# Other packages to load

library(vegan)
library(psych)
library(adegenet)
library(ggman)
library(qvalue)
library(pcadapt)
library(OutFLANK)
library(LEA)
library(diveRsity)

## MAP file ####
## load in the map file so we know where outliers are within the 
## genome
#MAP = read_tsv('icelandic_pops_plink.map')
MAP = read_tsv('Feb202019_Poly_Plink_input.map')
MAP$`#Chromosome`
chr_names = c('AC01', 'AC02', 'AC03', 'AC04p', 'AC04q.1:29', 'AC04q.2', 'AC05', 'AC06.1', 'AC06.2', 'AC07', 'AC08', 'AC09', 'AC10', 'AC11', 'AC12', 'AC13', 'AC14', 'AC15', 'AC16', 'AC17', 'AC18', 'AC19', 'AC20', 'AC21', 'AC22', 'AC23', 'AC24', 'AC25', 'AC26', 'AC27', 'AC28', 'AC30', 'AC31', 'AC32', 'AC33', 'AC34', 'AC35', 'AC36', 'AC37', 'contigs')

t(chr_names)
length(chr_names)
MAP = mutate(.data = MAP,
             CHROME = as.factor(case_when(
               `#Chromosome` == '1' ~ 'AC01',
               `#Chromosome` == '2' ~ 'AC02',
               `#Chromosome` == '3' ~ 'AC03',
               `#Chromosome` == '4' ~ 'AC04p',
               `#Chromosome` == '5' ~ 'AC04.1:29',
               `#Chromosome` == '6' ~ 'AC04q.2',
               `#Chromosome` == '7' ~ 'AC05',
               `#Chromosome` == '8' ~ 'AC06.1',
               `#Chromosome` == '9' ~ 'AC06.2',
               `#Chromosome` == '10' ~ 'AC07',
               `#Chromosome` == '11' ~ 'AC08',
               `#Chromosome` == '12' ~ 'AC09',
               `#Chromosome` == '13' ~ 'AC10',
               `#Chromosome` == '14' ~ 'AC11',
               `#Chromosome` == '15' ~ 'AC12',
               `#Chromosome` == '16' ~ 'AC13',
               `#Chromosome` == '17' ~ 'AC14',
               `#Chromosome` == '18' ~ 'AC15',
               `#Chromosome` == '19' ~ 'AC16',
               `#Chromosome` == '20' ~ 'AC17',
               `#Chromosome` == '21' ~ 'AC18',
               `#Chromosome` == '22' ~ 'AC19',
               `#Chromosome` == '23' ~ 'AC20',
               `#Chromosome` == '24' ~ 'AC21',
               `#Chromosome` == '25' ~ 'AC22',
               `#Chromosome` == '26' ~ 'AC23',
               `#Chromosome` == '27' ~ 'AC24',
               `#Chromosome` == '28' ~ 'AC25',
               `#Chromosome` == '29' ~ 'AC26',
               `#Chromosome` == '30' ~ 'AC27',
               `#Chromosome` == '31' ~ 'AC28',
               `#Chromosome` == '32' ~ 'AC30',
               `#Chromosome` == '33' ~ 'AC31',
               `#Chromosome` == '34' ~ 'AC32',
               `#Chromosome` == '35' ~ 'AC33',
               `#Chromosome` == '36' ~ 'AC34',
               `#Chromosome` == '37' ~ 'AC35',
               `#Chromosome` == '38' ~ 'AC36',
               `#Chromosome` == '39' ~ 'AC37',
               `#Chromosome` > '39' ~ 'Contigs')))

is.na(MAP$CHROME)
MAP$CHROME[is.na(MAP$CHROME)] = 'Contigs'
MAP$CHROME
#View(MAP3$CHROME3)

## Need to load the ped file with only polymorphic populations 
## Once loaded we can split the genotype file up into each populations
PED = read_tsv('Feb202019_Poly_Plink_input.ped')

## POPN PED FILES #####
## This allows us to create ped files for each population!!
# Pop_PED = PED %>% 
#   filter(`#FamilyID` %in% c('8', '4', '5'))
# write_tsv(Svin_PED, 'Svinavatn_Genotype.ped', col_names = T)

## POPS ####
## We will use this ped file to delinate populations

POPS = PED %>% select(`#FamilyID`) %>% 
filter(`#FamilyID` %in% c())

POPS = mutate(.data = POPS,
              POP_name = as.factor(case_when(
                `#FamilyID` == "1" ~ "T.LGB",
                `#FamilyID` == '2' ~ 'V.BR',
                `#FamilyID` == '3' ~ 'V.SIL',
                `#FamilyID` == '8' ~ 'S.LGB',
                `#FamilyID` == '4'~ 'S.PL',
                `#FamilyID` == '5' ~ 'S.PI',
                `#FamilyID` == '6' ~ 'T.PL',
                `#FamilyID` == '7' ~ 'T.SB',
                `#FamilyID` == '9' ~ 'G.SB',
                `#FamilyID` == '10' ~ 'G.PI'))) 
## Fst outlier overlap #####
## testing to see how many outliers overlap between the methods
## That we have applied to determine general outliers within the arctic 
## charr genome. 

## Load in the outlier loci across multipl different tests
## outlier loci common to multipl populations
PolyPopn_outliers = read_csv('.csv')

## outlier loci based on PCAdapt analysis
PCAdapt_out = read_csv('.csv') %>% 
  arrange(MARKER_ID)

## outlier loci based on sNMF analysis
snmf_out = read_csv('.csv') %>%
  rename(MARKER_ID = Marker.ID,
         CHROME = X.Chromosome,
         GDIST = Genetic.distance,
         PDIST = Physical.position,
         PVAL = snmf_pval.pvalues) %>%
  #select(-qval) %>%
  arrange(MARKER_ID)

## outlier loci based on BayeScan analysis
Baye_out = read_csv('.csv') %>% 
  rename(MARKER_ID = 2, 
         CHROME = `#Chromosome`, 
         GDIST = `Genetic distance`) %>% 
  #select(-qval) %>% 
  arrange(MARKER_ID)


## Determining the amount of SNP overlap between the three methods
Two_methods = PCAdapt_out$MARKER_ID[PCAdapt_out$MARKER_ID %in% Baye_out$MARKER_ID]
length(Two_method_OL)

Three_methods = Two_method_OL[Two_method_OL %in% snmf_out$MARKER_ID]
length(Three_method_OL)

## neutral SNPS ####
PCAdapt_out$MARKER_ID[PCAdapt_out$MARKER_ID %in% PolyPopn_outliers$MARKER_ID]
Baye_poly_overlap = BayeScan_out$MARKER_ID[BayeScan_out$MARKER_ID %in% PolyPopn_outliers$MARKER_ID]
length(Baye_poly_overlap)

## !!! Make Sure to remove any outlier loci that overlap between
## PolyPop_outliers and the other outlier analyses !!!!
PolyPopn_outliers = PolyPopn_outliers %>% 
  dplyr::filter(MARKER_ID %!in% c())


Outliers = bind_rows(PCAdapt_out, BayeScan_out, PolyPopn_outliers)
## not in function
## Function was adapted from Dr. Tony Kess
'%!in%' <- function(x,y)!('%in%'(x,y)) # not in

Neutral_SNP_list = data.frame(MAP$`Marker ID` %!in% Outliers$MARKER_ID,
                              stringsAsFactors = F) %>% 
  as_tibble() %>% 
  rename(NEUTRAL = 'MAP..Marker.ID....in..Outliers.MARKER_ID') 

bind_cols(MAP, Neutral_SNP_list) %>% 
  filter(NEUTRAL == TRUE) %>% 
  dplyr::select(`Marker ID`) %>% 
  write_tsv('.txt', 
            col_names = F)

## use --extract in plink to extract all of the neutral snps from the list