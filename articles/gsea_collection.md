# Build GSEA Gene Set Collections

``` r
suppressPackageStartupMessages(library(deseq2pip))
```

    ## Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
    ## 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'

``` r
## create a temporary directory for the GSEA database files
temp_dir <- tempdir()
#temp_dir <- "/Users/hungm/Documents/development/deseq2pip/inst/extdata/"
```

## Create gene set collections for `run_gsea()`

``` r
# create gene set collections for human and mouse
collection_name <- c("HALLMARK", "GOBP", "KEGG", "REACTOME", "BIOCARTA", "TFT")
for(i in c("hs_gs", "mm_gs")){
    if(i == "hs_gs"){
        org = "Homo sapiens"}
    else if(i == "mm_gs"){
        org = "Mus musculus"}
    gs <- list()

    # c1 gene set
    HM<-msigdbr(species = org, category ="H")
    gs[[1]] <- HM %>% dplyr::select(., gs_name, gene_symbol) %>% mutate(collection = collection_name[1])

    # c5 gene set
    GO<-msigdbr(species = org, category ="C5")
    gs[[2]] <-GO %>% dplyr::select(., gs_name, gene_symbol) %>% filter(stringr::str_detect(gs_name, 'GOBP')) %>% mutate(collection = collection_name[2])

    # c2 gene set
    C2<-msigdbr(species = org, category ="C2")
    gs[[3]] <-C2 %>% dplyr::select(., gs_name, gene_symbol) %>% filter(stringr::str_detect(gs_name, 'KEGG')) %>% mutate(collection = collection_name[3])
    gs[[4]] <-C2 %>% dplyr::select(., gs_name, gene_symbol) %>% filter(stringr::str_detect(gs_name, 'REACTOME')) %>% mutate(collection = collection_name[4])
    gs[[5]] <-C2 %>% dplyr::select(., gs_name, gene_symbol) %>% filter(stringr::str_detect(gs_name, 'BIOCARTA')) %>% mutate(collection = collection_name[5])

    # c3 gene set
    TFT <-msigdbr(species = org, category ="C3")
    gs[[6]] <-TFT %>% dplyr::select(., gs_name, gene_symbol) %>% mutate(collection = collection_name[6])

    gs <- bind_rows(gs)
    assign(i, gs)}
```

    ## Warning: The `category` argument of `msigdbr()` is deprecated as of msigdbr 10.0.0.
    ## ℹ Please use the `collection` argument instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Using human MSigDB with ortholog mapping to mouse. Use `db_species = "MM"` for mouse-native gene sets.
    ## This message is displayed once per session.

``` r
# Save TSV files to temporary directory
write.table(mm_gs, file.path(temp_dir, "20250319_mouse_msigdbr.tsv"), sep = "\t", quote = F, row.names = F)
write.table(hs_gs, file.path(temp_dir, "20250319_human_msigdbr.tsv"), sep = "\t", quote = F, row.names = F)
```

``` r
# show the first few rows of the human and mouse gene set collections
head(hs_gs)
```

    ## # A tibble: 6 × 3
    ##   gs_name               gene_symbol collection
    ##   <chr>                 <chr>       <chr>     
    ## 1 HALLMARK_ADIPOGENESIS ABCA1       HALLMARK  
    ## 2 HALLMARK_ADIPOGENESIS ABCB8       HALLMARK  
    ## 3 HALLMARK_ADIPOGENESIS ACAA2       HALLMARK  
    ## 4 HALLMARK_ADIPOGENESIS ACADL       HALLMARK  
    ## 5 HALLMARK_ADIPOGENESIS ACADM       HALLMARK  
    ## 6 HALLMARK_ADIPOGENESIS ACADS       HALLMARK

``` r
head(mm_gs)
```

    ## # A tibble: 6 × 3
    ##   gs_name               gene_symbol collection
    ##   <chr>                 <chr>       <chr>     
    ## 1 HALLMARK_ADIPOGENESIS Abca1       HALLMARK  
    ## 2 HALLMARK_ADIPOGENESIS Abcb8       HALLMARK  
    ## 3 HALLMARK_ADIPOGENESIS Acaa2       HALLMARK  
    ## 4 HALLMARK_ADIPOGENESIS Acadl       HALLMARK  
    ## 5 HALLMARK_ADIPOGENESIS Acadm       HALLMARK  
    ## 6 HALLMARK_ADIPOGENESIS Acads       HALLMARK

## Create GMT files for EnrichmentMap

``` r
# convert mouse gene set collection to a list of gene symbols
gs_name <- mm_gs$gs_name
mm_gs_list <- mm_gs %>%
    select(-c(collection, gs_name)) %>%
    .$gene_symbol %>%
    split(., gs_name)
mm_gs_list %>% head(.)
```

    ## $AAACCAC_MIR140
    ##  [1] "Abcc4"     "Actn4"     "Acvr1"     "Adam9"     "Adamts5"   "Bach1"    
    ##  [7] "Bcl2l2"    "Cebpa"     "Cit"       "Csk"       "Ctcf"      "Cul3"     
    ## [13] "Dpp4"      "Dtna"      "E2f3"      "Ehd1"      "Gys1"      "Hspa13"   
    ## [19] "Igfbp5"    "Map3k8"    "Mgat1"     "Mmd"       "Nr3c1"     "Pfn1"     
    ## [25] "Pitx2"     "Reep1"     "Rnf19a"    "Siah1a"    "Slc25a13"  "Slc38a2"  
    ## [31] "Snx2"      "Stag1"     "Trp53inp2" "Tssk2"     "Vegfa"     "Wnt1"     
    ## [37] "Zbtb10"   
    ## 
    ## $AAAGACA_MIR511
    ##  [1] "Abcg8"    "Ace"      "Ago2"     "Ago4"     "Alcam"    "Ankzf1"  
    ##  [7] "Atl2"     "Atrx"     "Btg1"     "Bub3"     "Calm1"    "Ccnd1"   
    ## [13] "Cdh2"     "Cltc"     "Cnot4"    "Dctn4"    "Dsc1"     "Dusp6"   
    ## [19] "E2f3"     "Enpp1"    "Fgf13"    "Fn1"      "Foxn3"    "Gad2"    
    ## [25] "Gja1"     "Has2"     "Kcnma1"   "Klf9"     "Klhl24"   "Lats1"   
    ## [31] "Luc7l3"   "Nhlh2"    "Nr4a2"    "Peli1"    "Pik3r3"   "Pmepa1"  
    ## [37] "Ppargc1a" "Prpf4b"   "Psma1"    "Rab22a"   "Rab2a"    "Rev3l"   
    ## [43] "Rhot1"    "Rnf19a"   "Selenop"  "Slc6a6"   "Socs2"    "Ss18"    
    ## [49] "Syt11"    "Tial1"    "Tnrc6b"   "Tob1"     "Vav3"     "Wnt16"   
    ## 
    ## $AAAGGAT_MIR501
    ##  [1] "Acaca"    "Acadsb"   "Adcyap1"  "Adipor2"  "Als2"     "Ammecr1" 
    ##  [7] "Atp6v1h"  "Bcl6"     "Celf2"    "Celsr2"   "Cul1"     "Dcx"     
    ## [13] "Errfi1"   "H2ax"     "Has2"     "Jun"      "Kif2a"    "Lpin1"   
    ## [19] "Map3k8"   "Mef2c"    "Myb"      "Mylk"     "Nfasc"    "Nfil3"   
    ## [25] "Nr4a3"    "Pdk1"     "Pitx2"    "Plxnb1"   "Pnn"      "Prkce"   
    ## [31] "Pura"     "Rab22a"   "Ret"      "Rnf11"    "Slc25a3"  "Smc1a"   
    ## [37] "Srsf2"    "Synj1"    "Tnni2"    "Togaram1" "Tomm70a"  "Usp12"   
    ## [43] "Vdac2"   
    ## 
    ## $AAAGGGA_MIR204_MIR211
    ##  [1] "Adcy6"     "Ago4"      "Akap1"     "Angpt1"    "Aph1a"     "Arcn1"    
    ##  [7] "Arhgap29"  "Atp2b1"    "Baz2a"     "Bcl11b"    "Bcl2"      "Brpf3"    
    ## [13] "Cdc25b"    "Cdc42"     "Cdh2"      "Clip1"     "Cox5a"     "Dennd5a"  
    ## [19] "Dlg5"      "Dtx1"      "Dyrk1a"    "Edem1"     "Eef1e1"    "Ell2"     
    ## [25] "Elmod3"    "Elovl6"    "Esr1"      "Ezr"       "Farp1"     "Fbn2"     
    ## [31] "Fbxw7"     "Fjx1"      "Gapvd1"    "Gga2"      "Hook3"     "Hs2st1"   
    ## [37] "Igf2r"     "Kcna3"     "Kmt5a"     "Lats1"     "Lpar1"     "Mall"     
    ## [43] "Med13l"    "Metap1"    "Mon2"      "Mrpl35"    "Ncoa7"     "Nr3c1"    
    ## [49] "Nr4a2"     "Plag1"     "Ppargc1a"  "Ppp3r1"    "Prdm2"     "Prrx1"    
    ## [55] "Rab14"     "Rab1a"     "Rbsn"      "Reep1"     "Rhobtb3"   "Rhot1"    
    ## [61] "Rictor"    "Rps6ka3"   "Rps6ka5"   "Sec24d"    "Shc1"      "Sirt1"    
    ## [67] "Spop"      "Ssrp1"     "Tcf7l1"    "Tgfbr2"    "Tnrc6b"    "Trp53inp1"
    ## [73] "Triap1"    "Wee1"      "Wsb1"      "Ywhag"    
    ## 
    ## $AAANWWTGC_UNKNOWN
    ##  [1] "Actb"     "App"      "Atp1b1"   "Atp2b4"   "Bcl6"     "Calm1"   
    ##  [7] "Cd14"     "Cdh13"    "Cited2"   "Cntfr"    "Dcaf11"   "Ddit3"   
    ## [13] "Dmd"      "Dusp1"    "Efna5"    "Elf4"     "Ephb2"    "Esr1"    
    ## [19] "Fbxw7"    "Fn1"      "Foxn3"    "Fzd7"     "Gata3"    "Gpc3"    
    ## [25] "Hhex"     "Inhba"    "Leap2"    "Lox"      "Luc7l3"   "Mef2c"   
    ## [31] "Meis1"    "Mgll"     "Mid1"     "Mmp3"     "Mrap"     "Mrps18b" 
    ## [37] "Myh2"     "Mylk"     "Nfe2l2"   "Ntrk3"    "Olfm1"    "Parp9"   
    ## [43] "Pax6"     "Pcsk1"    "Pdgfrb"   "Pik3r3"   "Pou2f1"   "Ppargc1a"
    ## [49] "Prpf4b"   "Rora"     "Sgcd"     "Skil"     "Skp2"     "Sparcl1" 
    ## [55] "Tfdp2"    "Tgif1"    "Thbs2"    "Tle4"     "Tlk1"     "Zw10"    
    ## 
    ## $AAAYRNCTG_UNKNOWN
    ##   [1] "Acvr1"     "Adam12"    "Add3"      "Aggf1"     "Apbb2"     "Arid4a"   
    ##   [7] "Arpc2"     "Atoh8"     "Atp5l"     "Bcl2l1"    "Bmpr1b"    "Btbd3"    
    ##  [13] "Bub3"      "Casq2"     "Cdh2"      "Cited2"    "Cldn5"     "Cltc"     
    ##  [19] "Cmklr1"    "Cntn1"     "Coch"      "Col12a1"   "Col1a2"    "Crat"     
    ##  [25] "Cstf3"     "Daam1"     "Ddah2"     "Ddx4"      "Def6"      "Dmd"      
    ##  [31] "Dmrt1"     "Dnajb4"    "Dusp4"     "Dyrk1a"    "Efna1"     "Eif5"     
    ##  [37] "Errfi1"    "Esrp2"     "Fgd4"      "Fgfr1"     "Foxo4"     "Fst"      
    ##  [43] "Gabra3"    "Gdnf"      "Gfi1"      "Hgf"       "Hsd3b7"    "Icam4"    
    ##  [49] "Id1"       "Igf1"      "Ip6k2"     "Itga10"    "Jpt1"      "Kcnip2"   
    ##  [55] "Kcnk5"     "Klf5"      "Lhx9"      "Ltbp1"     "Maml1"     "Map6"     
    ##  [61] "Meis1"     "Mgat1"     "Mid1"      "Mpped2"    "Myf6"      "Myh1"     
    ##  [67] "Myh10"     "Nav2"      "Ndufs4"    "Nek2"      "Nfyb"      "Nmi"      
    ##  [73] "Nmt1"      "Pacrg"     "Pdgfb"     "Pds5b"     "Pfn2"      "Phex"     
    ##  [79] "Plec"      "Plpp3"     "Prkcg"     "Prkn"      "Prok2"     "Pth1r"    
    ##  [85] "Ptpa"      "Pxn"       "Rplp0"     "Rtn1"      "Scn8a"     "Sdcbp"    
    ##  [91] "Sgcd"      "Sirpa"     "Slc4a1"    "Snx9"      "Sorbs2"    "Spock2"   
    ##  [97] "Stag1"     "Stc2"      "Tfdp2"     "Tnfaip8"   "Tns1"      "Trp53inp2"
    ## [103] "Trim28"    "Tspan7"    "Usp1"      "Vdr"       "Vip"       "Vwa5a"    
    ## [109] "Zfpm2"     "Zmat3"

``` r
# convert human gene set collection to a list of gene symbols
gs_name <- hs_gs$gs_name
hs_gs_list <- hs_gs %>%
    select(-c(collection, gs_name)) %>%
    .$gene_symbol %>%
    split(., gs_name)
hs_gs_list %>% head(.)
```

    ## $AAACCAC_MIR140
    ##   [1] "ABCC4"    "ABRAXAS2" "ACTN4"    "ACVR1"    "ADAM9"    "ADAMTS5" 
    ##   [7] "AMER2"    "ANK2"     "API5"     "ATOSA"    "BACH1"    "BAZ2B"   
    ##  [13] "BCL11A"   "BCL2L2"   "BCL9"     "C1orf21"  "CACNA1C"  "CEBPA"   
    ##  [19] "CHD4"     "CIT"      "COL23A1"  "CSK"      "CSNK1G3"  "CTCF"    
    ##  [25] "CUL3"     "DAZL"     "DBNDD2"   "DCUN1D4"  "DDX3X"    "DDX3Y"   
    ##  [31] "DHX57"    "DIPK2A"   "DPP4"     "DSCAM"    "DTNA"     "E2F3"    
    ##  [37] "EHD1"     "EPHB1"    "ERC2"     "ETV3"     "EYA2"     "GABARAP" 
    ##  [43] "GALNT16"  "GDF6"     "GIT1"     "GYS1"     "HDAC4"    "HNRNPH3" 
    ##  [49] "HSPA13"   "IGFBP5"   "KATNBL1"  "KCND2"    "LOXL3"    "LRRC4"   
    ##  [55] "LRRC8E"   "MAP3K8"   "MDGA2"    "MEX3C"    "MGAT1"    "MMD"     
    ##  [61] "MOK"      "NAV3"     "NKIRAS2"  "NR3C1"    "NSD3"     "NUTF2"   
    ##  [67] "OGT"      "OSTM1"    "PDGFRA"   "PFN1"     "PHF20L1"  "PHYHIP"  
    ##  [73] "PITX2"    "PPP1CC"   "PRIMA1"   "R3HDM1"   "REEP1"    "RNF19A"  
    ##  [79] "RTKN2"    "SAMTOR"   "SENP1"    "SIAH1"    "SLC25A13" "SLC38A2" 
    ##  [85] "SLC41A2"  "SLF2"     "SLMAP"    "SNX2"     "SOX4"     "SRR"     
    ##  [91] "STAG1"    "STRADB"   "SYT6"     "TAF9B"    "TBX3"     "TP53INP2"
    ##  [97] "TSHZ1"    "TSPAN2"   "TSSK2"    "TTYH2"    "UBASH3B"  "USP6"    
    ## [103] "VEGFA"    "WNT1"     "YES1"     "ZBED4"    "ZBTB10"   "ZNF182"  
    ## [109] "ZNF608"   "ZNF654"  
    ## 
    ## $AAAGACA_MIR511
    ##   [1] "ABCG8"     "ACE"       "ADAMTSL3"  "ADGRF5"    "ADSS2"     "AGBL3"    
    ##   [7] "AGO1"      "AGO2"      "AGO4"      "ALCAM"     "ANAPC15"   "ANKRD40CL"
    ##  [13] "ANKZF1"    "AQP6"      "ARHGEF17"  "ATL2"      "ATP2B2"    "ATRX"     
    ##  [19] "BCL11A"    "BTG1"      "BUB3"      "C1QL2"     "C1orf21"   "CALM1"    
    ##  [25] "CAMK2N1"   "CAMTA1"    "CAPRIN1"   "CCDC178"   "CCND1"     "CCNT2"    
    ##  [31] "CDH2"      "CDK14"     "CDK19"     "CELF1"     "CELF6"     "CEP350"   
    ##  [37] "CFAP298"   "CLK2"      "CLTC"      "CNOT4"     "CORIN"     "CREBRF"   
    ##  [43] "CREM"      "CRIM1"     "DCTN4"     "DDX3X"     "DDX3Y"     "DEDD"     
    ##  [49] "DNAJB12"   "DNAJC13"   "DSC1"      "DUSP6"     "DYRK1B"    "E2F3"     
    ##  [55] "EDEM3"     "EFR3A"     "ELAVL3"    "EMILIN2"   "EML4"      "ENPP1"    
    ##  [61] "ENPP4"     "EPHA4"     "ESRRG"     "EYA1"      "EYA4"      "FAM117A"  
    ##  [67] "FGF13"     "FIP1L1"    "FMR1"      "FN1"       "FNDC1"     "FNDC5"    
    ##  [73] "FOXK2"     "FOXN3"     "GAD2"      "GEMIN2"    "GFAP"      "GJA1"     
    ##  [79] "GLRA2"     "HAS2"      "HCN4"      "HLF"       "HLTF"      "HOXA13"   
    ##  [85] "IGF2BP1"   "IGF2BP3"   "ILRUN"     "KCNE1"     "KCNMA1"    "KHDRBS2"  
    ##  [91] "KLF9"      "KLHL18"    "KLHL24"    "LATS1"     "LMCD1"     "LPP"      
    ##  [97] "LRCH4"     "LUC7L3"    "MAP3K2"    "MAP4K4"    "MAPK1IP1L" "MBD2"     
    ## [103] "MBD6"      "MDGA2"     "METAP2"    "MIB1"      "MINK1"     "MRPL21"   
    ## [109] "MSTN"      "MTAP"      "MYCBP"     "MYO19"     "NAA50"     "NACC1"    
    ## [115] "NEUROD6"   "NHLH2"     "NLK"       "NR4A2"     "NRXN3"     "NTRK2"    
    ## [121] "NXPH1"     "ONECUT2"   "PAX8"      "PCARE"     "PCDH10"    "PCDH17"   
    ## [127] "PELI1"     "PHAF1"     "PHLPP1"    "PIK3R3"    "PMEPA1"    "POGK"     
    ## [133] "POU4F2"    "PPARGC1A"  "PRELP"     "PRP4K"     "PSMA1"     "PSMD10"   
    ## [139] "PTGR3"     "QKI"       "RAB22A"    "RAB2A"     "RBM15B"    "RBM26"    
    ## [145] "RECK"      "REV3L"     "RGL1"      "RHOJ"      "RHOT1"     "RNF19A"   
    ## [151] "ROBO2"     "RPS6KB1"   "RPS6KL1"   "SATB2"     "SCN4B"     "SELENOP"  
    ## [157] "SEMA3F"    "SEMA6D"    "SINHCAF"   "SLC22A17"  "SLC25A26"  "SLC6A6"   
    ## [163] "SLITRK1"   "SMARCE1"   "SOCS2"     "SORCS3"    "SOST"      "SOX12"    
    ## [169] "SPTBN4"    "SPTLC2"    "SRGAP3"    "SS18"      "ST18"      "SYT11"    
    ## [175] "TAF5"      "TBXT"      "THOC5"     "TIAL1"     "TMEM196"   "TMEM243"  
    ## [181] "TMEM248"   "TNRC6A"    "TNRC6B"    "TOB1"      "TRAPPC3"   "TRAPPC8"  
    ## [187] "TRIM2"     "TRIM24"    "TSPOAP1"   "TXNL1"     "UBE2H"     "VANGL2"   
    ## [193] "VAV3"      "VIRMA"     "VKORC1L1"  "VMP1"      "WNT16"     "YTHDF2"   
    ## [199] "YY1"       "ZCCHC24"   "ZDHHC21"   "ZNF319"    "ZNF654"    "ZNF706"   
    ## 
    ## $AAAGGAT_MIR501
    ##   [1] "ACACA"    "ACADSB"   "ADCYAP1"  "ADIPOR2"  "ALS2"     "ALX1"    
    ##   [7] "AMMECR1"  "APOLD1"   "ATP6V1H"  "BCL6"     "BCLAF1"   "C8orf82" 
    ##  [13] "CACHD1"   "CAMTA1"   "CBLIF"    "CCDC140"  "CD164"    "CDK2AP1" 
    ##  [19] "CELF2"    "CELSR2"   "CHODL"    "CLK1"     "CLK2"     "CTDSP1"  
    ##  [25] "CTDSPL2"  "CUL1"     "CUX2"     "DCX"      "DNAJB12"  "ELAVL4"  
    ##  [31] "ERRFI1"   "GRB10"    "H2AX"     "HAS2"     "HES5"     "HOXB8"   
    ##  [37] "JADE3"    "JUN"      "KCND2"    "KCNRG"    "KIF1C"    "KIF2A"   
    ##  [43] "KLHL14"   "KRR1"     "LARP1"    "LEPROTL1" "LPGAT1"   "LPIN1"   
    ##  [49] "LRRC1"    "MAP2K1"   "MAP3K8"   "MCU"      "MEF2C"    "MYB"     
    ##  [55] "MYCL"     "MYLK"     "NEXMIF"   "NFASC"    "NFIL3"    "NFIX"    
    ##  [61] "NPR3"     "NR2F2"    "NR4A3"    "PCDH19"   "PDK1"     "PHC1"    
    ##  [67] "PHF6"     "PIK3AP1"  "PITX2"    "PLP1"     "PLXNB1"   "PNN"     
    ##  [73] "PPP1CB"   "PPP2R5E"  "PPP4R3A"  "PPP6R3"   "PRKCE"    "PURA"    
    ##  [79] "QKI"      "RAB22A"   "RABGEF1"  "RASL10B"  "RCN1"     "RDX"     
    ##  [85] "RET"      "RGL1"     "RNF11"    "ROBO2"    "RPGRIP1L" "RSBN1"   
    ##  [91] "SATB2"    "SCN3A"    "SENP3"    "SEPHS1"   "SGPP1"    "SLC25A3" 
    ##  [97] "SLC35B3"  "SLITRK5"  "SMC1A"    "SNAP29"   "SOX11"    "SOX4"    
    ## [103] "SPOPL"    "SRR"      "SRSF2"    "SYNC"     "SYNJ1"    "SYT7"    
    ## [109] "TAF5L"    "TAPT1"    "TNNI2"    "TOGARAM1" "TOMM70"   "TRIM39"  
    ## [115] "TSC22D3"  "UBAP1"    "UBE2Q1"   "UBE4B"    "USP12"    "VDAC2"   
    ## [121] "WDFY3"    "WIPF2"    "WT1-AS"   "ZBTB18"   "ZC3H7A"   "ZIC4"    
    ## [127] "ZMYM5"   
    ## 
    ## $AAAGGGA_MIR204_MIR211
    ##   [1] "ABRAXAS2" "ADAMTS9"  "ADCY6"    "ADPRM"    "AGO4"     "AKAP1"   
    ##   [7] "ALPL"     "ANGPT1"   "ANKRD13A" "ANXA11"   "AP1S1"    "AP1S3"   
    ##  [13] "AP2A2"    "AP3M1"    "APH1A"    "ARAP2"    "ARCN1"    "ARGLU1"  
    ##  [19] "ARHGAP29" "ARL8B"    "ATF2"     "ATP2B1"   "AUP1"     "BAZ2A"   
    ##  [25] "BCL11B"   "BCL2"     "BCL9"     "BCL9L"    "BRD4"     "BRPF3"   
    ##  [31] "CAPRIN1"  "CCDC120"  "CCNT2"    "CCPG1"    "CDC25B"   "CDC42"   
    ##  [37] "CDH2"     "CELSR3"   "CHD5"     "CHN2"     "CHP1"     "CLIP1"   
    ##  [43] "CORO1C"   "COX5A"    "CPD"      "CPNE8"    "CREB5"    "CRKL"    
    ##  [49] "CTDNEP1"  "DAG1"     "DCAF5"    "DCUN1D3"  "DENND5A"  "DHH"     
    ##  [55] "DLG5"     "DMTF1"    "DNAJC13"  "DNM2"     "DTX1"     "DVL3"    
    ##  [61] "DYRK1A"   "EDEM1"    "EEF1E1"   "EFNB3"    "ELAVL3"   "ELF2"    
    ##  [67] "ELL2"     "ELMOD3"   "ELOVL6"   "EPHA7"    "EPHB6"    "ESR1"    
    ##  [73] "ESRRG"    "EVA1C"    "EZR"      "FAM117B"  "FAM120C"  "FARP1"   
    ##  [79] "FBN2"     "FBXW7"    "FHIP1B"   "FJX1"     "FNIP1"    "FRAS1"   
    ##  [85] "FREM1"    "FRY"      "GABRB3"   "GAPVD1"   "GGA2"     "GLIS3"   
    ##  [91] "GPM6A"    "GRM1"     "HAPSTR1"  "HIC2"     "HMGA2"    "HOOK3"   
    ##  [97] "HOXC8"    "HS2ST1"   "IGF2R"    "ING4"     "ITPR1"    "JPH3"    
    ## [103] "KCNA3"    "KCTD1"    "KDM2A"    "KHDRBS1"  "KHDRBS3"  "KITLG"   
    ## [109] "KLF12"    "KLHL13"   "KMT2A"    "KMT5A"    "LATS1"    "LPAR1"   
    ## [115] "LRRC8D"   "MALL"     "MAML3"    "MAP1LC3B" "MAP3K3"   "MBNL1"   
    ## [121] "MED13L"   "METAP1"   "MIR600HG" "MLLT3"    "MMGT1"    "MON2"    
    ## [127] "MRPL35"   "MRPL52"   "MYO10"    "NAA15"    "NBEA"     "NCOA7"   
    ## [133] "NEUROG1"  "NOVA1"    "NPTX1"    "NR3C1"    "NR4A2"    "NRBF2"   
    ## [139] "NTRK2"    "P4HB"     "PABIR2"   "PCDH9"    "PHF13"    "PID1"    
    ## [145] "PLAG1"    "POU3F2"   "PPARGC1A" "PPP3R1"   "PRDM2"    "PRPF38B" 
    ## [151] "PRRX1"    "RAB10"    "RAB14"    "RAB1A"    "RAP2C"    "RBSN"    
    ## [157] "REEP1"    "RERE"     "RHOBTB3"  "RHOT1"    "RICTOR"   "RPS6KA3" 
    ## [163] "RPS6KA5"  "RPS6KC1"  "RSPO3"    "RTKN2"    "RUNX2"    "SATB2"   
    ## [169] "SCRT2"    "SEC24D"   "SEC61A2"  "SERINC3"  "SF3B1"    "SGCZ"    
    ## [175] "SGIP1"    "SHC1"     "SIN3A"    "SIRT1"    "SLC17A7"  "SLC22A2" 
    ## [181] "SLC37A3"  "SLITRK4"  "SLTM"     "SMOC1"    "SOCS6"    "SOX11"   
    ## [187] "SOX4"     "SPOP"     "SPRED1"   "SPRYD7"   "SSRP1"    "ST7"     
    ## [193] "STXBP5"   "SUMO2"    "SUMO4"    "SZRD1"    "TAF5"     "TCF12"   
    ## [199] "TCF7L1"   "TGFBR2"   "TMEM30A"  "TMOD3"    "TNRC6B"   "TP53INP1"
    ## [205] "TRIAP1"   "TRIP12"   "TRPC5"    "TTYH1"    "UBE2R2"   "UHRF2"   
    ## [211] "USP6"     "WEE1"     "WNT4"     "WSB1"     "XRN1"     "YTHDF3"  
    ## [217] "YWHAG"    "ZCCHC14"  "ZCCHC24"  "ZDHHC17"  "ZFC3H1"   "ZFP91"   
    ## [223] "ZNF282"   "ZNF335"   "ZNF423"  
    ## 
    ## $AAANWWTGC_UNKNOWN
    ##   [1] "ACTB"      "ADHFE1"    "AFF4"      "ANK2"      "ANK3"      "APP"      
    ##   [7] "ASPA"      "ATOH7"     "ATP1B1"    "ATP2B4"    "ATXN7L1"   "BCL11A"   
    ##  [13] "BCL6"      "BNC2"      "C11orf87"  "CACNA1D"   "CACNG3"    "CALM1"    
    ##  [19] "CD14"      "CDC42EP3"  "CDC42EP5"  "CDH13"     "CDK2AP1"   "CEPT1"    
    ##  [25] "CHD2"      "CITED2"    "CNMD"      "CNTFR"     "DAB1"      "DCAF11"   
    ##  [31] "DCHS2"     "DDIT3"     "DIS3L"     "DLG2"      "DLGAP4"    "DMD"      
    ##  [37] "DNAJB5"    "DPYSL5"    "DRD3"      "DSCAM"     "DSEL"      "DSTN"     
    ##  [43] "DTX3L"     "DUSP1"     "DYNC1I2"   "EBF1"      "EFNA5"     "EGFLAM"   
    ##  [49] "EIF4EBP2"  "ELAVL4"    "ELF4"      "EPHA7"     "EPHB2"     "ESR1"     
    ##  [55] "FAM78A"    "FBXW7"     "FGF7"      "FGFR2"     "FN1"       "FOXN3"    
    ##  [61] "FOXP1"     "FOXP2"     "FTHL17"    "FZD7"      "GANAB"     "GATA3"    
    ##  [67] "GLRA2"     "GPC3"      "GPC6"      "GPR21"     "GPRIN3"    "GRHL3"    
    ##  [73] "GRIN2B"    "GTF2E2"    "HEPACAM"   "HHEX"      "HOXA2"     "HOXA3"    
    ##  [79] "HOXB2"     "HOXB6"     "HOXC4"     "IGF2BP1"   "INHBA"     "ITM2C"    
    ##  [85] "JADE2"     "KANK1"     "KCNJ13"    "KLF12"     "KLF14"     "KRTAP8-1" 
    ##  [91] "LEAP2"     "LINC00314" "LIPG"      "LOX"       "LOXL4"     "LRRC3B"   
    ##  [97] "LRRN1"     "LSAMP"     "LUC7L3"    "MAML3"     "MAN2A2"    "MAP3K4"   
    ## [103] "MAPK3"     "MBNL1"     "MEF2C"     "MEIS1"     "MGLL"      "MID1"     
    ## [109] "MLLT6"     "MMP3"      "MPZL3"     "MRAP"      "MRPL24"    "MRPS18B"  
    ## [115] "MYCL"      "MYH2"      "MYLK"      "NCBP3"     "NEK6"      "NEUROG1"  
    ## [121] "NFE2L2"    "NNAT"      "NR2F2"     "NRAS"      "NTN1"      "NTRK3"    
    ## [127] "OLFM1"     "OLIG2"     "OMG"       "OTX2"      "PARP9"     "PATZ1"    
    ## [133] "PAX1"      "PAX6"      "PCSK1"     "PCTP"      "PDGFRB"    "PHOX2B"   
    ## [139] "PHTF1"     "PIK3R3"    "POU2F1"    "POU4F1"    "PPARGC1A"  "PPFIA2"   
    ## [145] "PPP1R10"   "PPP2R2A"   "PPP3CC"    "PRDM16"    "PRIMA1"    "PRP4K"    
    ## [151] "RAB10"     "RBMX"      "RORA"      "RRS1"      "RSPO2"     "S100PBP"  
    ## [157] "SALL3"     "SAMD12"    "SATB2"     "SEMA6C"    "SESN2"     "SFRP2"    
    ## [163] "SGCD"      "SHC3"      "SIX5"      "SKIL"      "SKP2"      "SLMAP"    
    ## [169] "SNCAIP"    "SNX25"     "SORT1"     "SOX13"     "SOX4"      "SOX5"     
    ## [175] "SPAG9"     "SPARCL1"   "SSBP3"     "STEAP2"    "TBC1D8B"   "TFAP4"    
    ## [181] "TFDP2"     "TGIF1"     "THAP12"    "THBS2"     "TLE4"      "TLK1"     
    ## [187] "TLX3"      "TRAM1"     "TRPM3"     "TSC22D4"   "ZFPM1"     "ZHX3"     
    ## [193] "ZNF462"    "ZNF827"    "ZW10"     
    ## 
    ## $AAAYRNCTG_UNKNOWN
    ##   [1] "ABT1"     "ACVR1"    "ADAM12"   "ADD3"     "ADGRB3"   "AGGF1"   
    ##   [7] "ANKRD12"  "ANKRD28"  "AP4S1"    "APBB2"    "APOBR"    "AQP2"    
    ##  [13] "ARHGAP44" "ARID1A"   "ARID4A"   "ARPC2"    "ARSG"     "ARX"     
    ##  [19] "ASB4"     "ASPH"     "ATOH8"    "ATP1A2"   "ATP5IF1"  "ATP5MG"  
    ##  [25] "AXDND1"   "B4GALT6"  "BAMBI"    "BCL2L1"   "BCL9"     "BMPR1B"  
    ##  [31] "BMX"      "BRSK2"    "BTBD3"    "BUB3"     "C7orf33"  "CA3"     
    ##  [37] "CACNA2D3" "CACNB2"   "CAPN1"    "CAPZA1"   "CASQ2"    "CBX2"    
    ##  [43] "CCDC174"  "CCNJ"     "CCNY"     "CDC23"    "CDH2"     "CER1"    
    ##  [49] "CFAP161"  "CHRM1"    "CITED2"   "CLDN5"    "CLTC"     "CLTRN"   
    ##  [55] "CMKLR1"   "CNTLN"    "CNTN1"    "COCH"     "COL12A1"  "COL1A2"  
    ##  [61] "COL4A5"   "COL4A6"   "COLCA1"   "COLEC10"  "CRAT"     "CRH"     
    ##  [67] "CRKL"     "CRYGD"    "CRYGS"    "CSNK1A1"  "CSRNP3"   "CSTF3"   
    ##  [73] "CYBRD1"   "CYRIA"    "DAAM1"    "DBNDD2"   "DCAKD"    "DDAH2"   
    ##  [79] "DDX4"     "DEF6"     "DENND4A"  "DGKB"     "DHH"      "DHRS4"   
    ##  [85] "DHRS4L2"  "DIDO1"    "DMD"      "DMRT1"    "DNAJA2"   "DNAJB3"  
    ##  [91] "DNAJB3"   "DNAJB4"   "DSCAML1"  "DUSP4"    "DYNC1I1"  "DYRK1A"  
    ##  [97] "EDA"      "EFNA1"    "EGFLAM"   "EIF5"     "EMX2"     "EPC1"    
    ## [103] "EPHA7"    "ERBB4"    "ERG28"    "ERRFI1"   "ESRP2"    "ESRRB"   
    ## [109] "ESRRG"    "EYA1"     "FAM216B"  "FAM83F"   "FCER1A"   "FGD4"    
    ## [115] "FGF10"    "FGF12"    "FGFR1"    "FGFR1OP2" "FIZ1"     "FKRP"    
    ## [121] "FMNL3"    "FNDC9"    "FOXA1"    "FOXG1"    "FOXO4"    "FOXP2"   
    ## [127] "FSIP2"    "FST"      "GABRA3"   "GDNF"     "GFI1"     "GGNBP2"  
    ## [133] "GJB4"     "GLDN"     "GNAQ"     "GPR85"    "GPRC5D"   "GRIN2B"  
    ## [139] "GSE1"     "H3-3A"    "HDAC8"    "HESX1"    "HEXIM2"   "HGF"     
    ## [145] "HIC2"     "HID1"     "HIP1R"    "HOXA10"   "HOXA5"    "HOXB8"   
    ## [151] "HPSE2"    "HSD3B7"   "ICAM4"    "ID1"      "IGF1"     "IL1RAPL1"
    ## [157] "INHBC"    "IP6K2"    "ITGA10"   "ITGA8"    "JADE2"    "JPH1"    
    ## [163] "JPT1"     "KANK2"    "KCNIP2"   "KCNK5"    "KCNN3"    "KCNQ1DN" 
    ## [169] "KITLG"    "KLF5"     "KLHDC10"  "KLHL20"   "KLHL3"    "KMT2A"   
    ## [175] "LARS2"    "LBX2-AS1" "LENG9"    "LHFPL6"   "LHX9"     "LIX1"    
    ## [181] "LMO7"     "LRP5"     "LRRC4"    "LRRN4CL"  "LTBP1"    "MAML1"   
    ## [187] "MANF"     "MAP2"     "MAP3K5"   "MAP6"     "MEIS1"    "MGAT1"   
    ## [193] "MGAT4A"   "MID1"     "MIR646HG" "MOAP1"    "MPPED2"   "MRPL13"  
    ## [199] "MTA2"     "MTBP"     "MTRFR"    "MYF6"     "MYH1"     "MYH10"   
    ## [205] "MYO18A"   "NAGLU"    "NAPB"     "NAV2"     "NAV3"     "NCDN"    
    ## [211] "NDNF"     "NDST4"    "NDUFS4"   "NEK1"     "NEK2"     "NFATC4"  
    ## [217] "NFYB"     "NMI"      "NMT1"     "NR2F1"    "NRG1"     "NTRK2"   
    ## [223] "NUP54"    "NXPH4"    "OMA1"     "OMG"      "OR2L13"   "OTX2"    
    ## [229] "PACRG"    "PALS2"    "PART1"    "PCDH17"   "PCDH18"   "PCF11"   
    ## [235] "PCYT1B"   "PDGFB"    "PDGFRA"   "PDLIM2"   "PDS5B"    "PDZRN4"  
    ## [241] "PFN2"     "PHC2"     "PHEX"     "PHF1"     "PHF6"     "PHOX2B"  
    ## [247] "PLAGL2"   "PLEC"     "PLEKHM1"  "PLP2"     "PLPP3"    "PMCH"    
    ## [253] "PMCHL1"   "PODXL2"   "POFUT1"   "POU2AF1"  "POU4F1"   "PPP1R9B" 
    ## [259] "PPP2R3A"  "PPP2R5E"  "PPP3CA"   "PRELP"    "PRKCG"    "PRKCQ"   
    ## [265] "PRKN"     "PROK2"    "PTCHD4"   "PTGR3"    "PTH1R"    "PTPA"    
    ## [271] "PXN"      "R3HDM1"   "RAB30"    "RAB5B"    "RAB5C"    "RAPGEF4" 
    ## [277] "RBMS3"    "RGS17"    "RNF146"   "ROBO4"    "ROR1"     "RPLP0"   
    ## [283] "RTN1"     "S1PR2"    "SCN3B"    "SCN5A"    "SCN8A"    "SCOC"    
    ## [289] "SDCBP"    "SEMA6D"   "SEPTIN7"  "SESN3"    "SGCD"     "SH2D6"   
    ## [295] "SHC3"     "SHCBP1L"  "SIPA1"    "SIRPA"    "SLC26A6"  "SLC4A1"  
    ## [301] "SLC6A1"   "SMARCA2"  "SNX9"     "SORBS2"   "SOX12"    "SOX21"   
    ## [307] "SOX30"    "SOX5"     "SPINDOC"  "SPOCK2"   "SPTLC2"   "SRGAP2"  
    ## [313] "SRSF8"    "SSBP2"    "ST7L"     "STAC3"    "STAG1"    "STAG2"   
    ## [319] "STC2"     "STRN3"    "STRN4"    "TAS1R2"   "TEF"      "TENT4B"  
    ## [325] "TFAP4"    "TFDP2"    "TM2D3"    "TMEM182"  "TMEM69"   "TMSB4X"  
    ## [331] "TMSB4XP1" "TMSB4XP6" "TMSB4XP8" "TNFAIP8"  "TNS1"     "TNXB"    
    ## [337] "TP53INP2" "TRDN"     "TREML1"   "TRIM28"   "TRIM68"   "TRIM8"   
    ## [343] "TRIML1"   "TRPS1"    "TSC22D3"  "TSPAN7"   "TSSK3"    "TTC17"   
    ## [349] "TUSC2"    "UBE2W"    "UBXN10"   "USP1"     "VDR"      "VIP"     
    ## [355] "VKORC1L1" "VWA5A"    "WBP1"     "WNT2B"    "WT1"      "WT1-AS"  
    ## [361] "XRCC1"    "ZBTB11"   "ZBTB18"   "ZFP91"    "ZFPM2"    "ZIC1"    
    ## [367] "ZIC4"     "ZMAT3"    "ZNF296"   "ZNF503"   "ZNF521"   "ZNF524"  
    ## [373] "ZNF654"   "ZNF687"   "ZNF710"

``` r
# create GMT files for mouse and human gene set collections
write.gmt <- function(genelist, file) {
  file_conn <- file(file, open = "wt")
  for (set_name in names(genelist)) {
    line <- c(set_name, set_name, genelist[[set_name]])
    writeLines(paste(line, collapse = "\t"), file_conn)}
  close(file_conn)}
  
# Save GMT files to temporary directory
write.gmt(mm_gs_list, file.path(temp_dir, "20250319_mouse_msgidbr.gmt"))
write.gmt(hs_gs_list, file.path(temp_dir, "20250319_human_msgidbr.gmt"))
```
