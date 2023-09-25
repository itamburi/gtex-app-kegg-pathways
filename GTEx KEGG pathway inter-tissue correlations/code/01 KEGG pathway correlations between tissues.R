


here::i_am("code/01 KEGG pathway correlations between tissues.R")
source(here::here('code/Z-source.R'))



#### 0.0 - Load Data ################################
# For Natalie
# produce the initial tables to upload for SQL


### Old way to load and prepare working dataset
# load(here('data/raw files/GTEx NA included env.RData'))
# rm(GTEx_full)
# working_dataset=GTEx_subfiltered
# row.names(working_dataset) = working_dataset$gene_tissue
# working_dataset$gene_tissue=NULL
# working_dataset = as.data.frame(t(working_dataset))
# rm(GTEx_subfiltered)


#### 1.0 - Saved the processed working dataset here to reduce memory usage. Load into enviornment
load(here('data/GTEx working_dataset.RData'))
load(here('data/tissue_factors.RData')) # load tissue factors which were established at the end of this code after running hsa04062 individualy

print(tissue_factors)




#### 2.0 - Loop for KEGG pathways and compare against randomly selected genes ####


## 2.1 - KEGG pathways
library(KEGGREST)

# pathway IDs can be found on KEGG here: https://www.genome.jp/kegg/pathway.html#organismal

# metabolic pathways
kp1 = c( 
  #"hsa00120", # bile acids
  "hsa00010", # glycolysis gluconeogenesis
  "hsa00020", # TCA
  "hsa00190", # Oxidative phosphorylation
  "hsa04910", # Insulin signaling
  "hsa00071" # fatty acid degradation
  )

# Immune pathways
kp2 = c(
    "hsa04062", # Chemokine signaling
    "hsa04620", # TLRs
    "hsa04640"  # Hematopoietic cell lineage
    ) 

# Other
kp3 = c(
  "hsa04911", # insulin secretion 
  "hsa04710", # circadia rhythmn
  "hsa03010" # ribosome
  ) 



keggpaths = c(kp1, kp2, kp3)

# keggpaths = ("hsa04141")






## 2.2 - Plot function 

plot_cor_counts = function(pthres, table){
  
  
  sig_table = table
  
  
  
  ##set the color scheme up front
  
  #col_scheme = rev(met.brewer('Austria', length(unique(sig_table$tissue_2))))
  #names(col_scheme) = unique(sig_table$tissue_2)
  
  col_scheme = met.brewer('Austria', length(unique(tissue_factors)))
  names(col_scheme) = unique(tissue_factors)
  
  
  
  
  #Up front I think it will be helpful to have the user see the top-ranked genes correlating both within and excluding origin.  I set 30 as max_gene_length but this might be cool to let user enter.  We should max out at 40 or 50
  sig_table$logq = ifelse(sig_table$pvalue<pthres, 'S', 'NS')
  sig_table = sig_table[sig_table$logq=='S',]
  sig_table1 = sig_table %>%
    dplyr::group_by(tissue_2, Category, kegg_ENTRY, kegg_NAME) %>%
    dplyr::summarise(count = n())
  
  #top_genes = sig_table[1:max_gene_length,]
  
  sig_table1$color = col_scheme[match(sig_table1$tissue_2, names(col_scheme))]
  
  #sig_table = sig_table[sig_table$pvalue<0.05,]
  #pdf(file = here( paste0('plots/kegg pathway corrs/Number of Significant co-correlated ', NAME, " " , ENTRY , ' genes P less ', pthres, '.pdf')) ) 
  
  g2 = ggplot(sig_table1, aes(x=fct_reorder(tissue_2, count, .desc=T), y=count,  fill=tissue_2)) + geom_col() +
    #geom_col(position = position_dodge2(reverse = TRUE)) +  
    facet_grid(~Category, scales = "free") +
    theme_classic() + theme(legend.position = "none") +
    scale_fill_manual(values = sig_table1$color[order(sig_table1$count, decreasing = T)]) +
    labs(subtitle = paste0("genes from ", ENTRY, " - ", NAME)) +
    ggtitle(paste0('Number of Significant co-correlated genes P< ', pthres )) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    ) +
    scale_x_discrete(limits = tissue_factors)
  print(g2)
  
}




#### **** 2.3 - Loop over KEGG pathways and make tissue distribution plots ****

# establish specific tissue order
# show different distribution of these gene corrs over the tissues per pathway




for( path in keggpaths ){

  
  
  #### **** Make list of genes from KEGG pathway ****
  
  pathway = keggGet(path)
  print(paste0("Starting ", path))
  
  NAME = pathway[[1]]$NAME
  forbidden = c("/") #characters incompatible with filenames
  NAME = gsub(paste0("(", paste(forbidden, collapse = "|"), ")"), "", NAME)
  
  ENTRY = pathway[[1]]$ENTRY
  
  
  # extract the gene IDs from the pathway object
  gene_ids = pathway[[1]]$GENE
  gene_ids = ifelse(grepl("^\\d*$", gene_ids), NA, gene_ids)
  gene_ids = subset(gene_ids, is.na(gene_ids) == FALSE)
  gene_ids = sub(";.*", "", gene_ids)
  gene_set1 = paste(gene_ids, collapse = '|')
  
  
  # update tissue list of Gtex tissues in app
  tissue_list <- c('Adipose - Subcutaneous|Adipose - Visceral (Omentum)|Brain - Hypothalamus|Brain - Hippocampus|Small Intestine - Terminal Ileum|Stomach|Thyroid|Pancreas|Spleen|Muscle - Skeletal|Pituitary|Artery - Coronary|Liver|Kidney - Cortex|Heart - Left Ventricle|Colon - Transverse|Colon - Sigmoid|Adrenal Gland|Artery - Aorta')
  
  
  
  
  #### *** RAND CORS *** - Select a random number of genes equal to the number in the pathway
  
  gene_count = ( stringr::str_count(gene_set1, "\\|") + 1 )*18
    # 49 unique tissues determined by length(unique(cor_table2$tissue_1)) if we use everyhting in gtex
    # 18 unique tissues if we use the vector tissue_list from above
  #rand = working_dataset[,grepl(tissue_list, colnames(working_dataset))]
  rand = working_dataset
  rand = rand[,grepl(tissue_list, colnames(rand))]
  rand = rand[,sample(ncol(rand), gene_count) ]
  rand = rand[,!grepl('AS1', colnames(rand))]
  rand = rand[,!grepl('AS2', colnames(rand))]
  rand = rand[,!grepl('CATSPER', colnames(rand))]
  rand = rand[,!grepl('P1', colnames(rand))]
  colnames(rand)[1:100]
  rand_cors = bicorAndPvalue(rand, rand, use = 'p')
  rndcor_table = reshape2::melt(rand_cors$bicor)
  rand_new_p = reshape2::melt(rand_cors$p)
  
  colnames(rndcor_table) = c('gene_tissue_1', 'gene_tissue_2', 'bicor')
  #can drop here to clear CPU
  rand_cors=NULL
  
  rndcor_table$pvalue = signif(rand_new_p$value, 3)
  rndcor_table$bicor = round(rndcor_table$bicor, 3)
  rndcor_table$qvalue = signif(p.adjust(rndcor_table$pvalue, "BH"), 3)
  rndcor_table = rndcor_table[order(rndcor_table$qvalue, decreasing=F),]
  rndcor_table = na.omit(rndcor_table)
  rndcor_table$gene_symbol_1 = gsub("\\_.*","",rndcor_table$gene_tissue_1)
  rndcor_table$tissue_1 = gsub(".*_","",rndcor_table$gene_tissue_1)
  rndcor_table = rndcor_table[!is.na(rndcor_table$tissue_1),]
  
  rndcor_table$gene_symbol_2 = gsub("\\_.*","",rndcor_table$gene_tissue_2)
  rndcor_table$tissue_2 = gsub(".*_","",rndcor_table$gene_tissue_2)
  rndcor_table = rndcor_table[!is.na(rndcor_table$tissue_2),]
  
  rndcor_table$Category = "random.genes"
  rndcor_table$kegg_ENTRY = ENTRY
  rndcor_table$kegg_NAME = NAME
  # remove like-like corrs of the same gene_tissue object
  rndcor_table = rndcor_table[!rndcor_table$gene_tissue_1==rndcor_table$gene_tissue_2,]
  
  
  
  
  #### *** PATHWAY CORS **** - Target genes / pathways
  
  tissue1 = working_dataset[,grepl(gene_set1, colnames(working_dataset))]
  tissue1 = tissue1[,grepl(tissue_list, colnames(tissue1))]
  tissue1 = tissue1[,!grepl('AS1', colnames(tissue1))] # I'm not sure what the point of these is exactly, but im guessing they remve extraneous entries
  tissue1 = tissue1[,!grepl('AS2', colnames(tissue1))]
  tissue1 = tissue1[,!grepl('CATSPER', colnames(tissue1))]
  tissue1 = tissue1[,!grepl('P1', colnames(tissue1))]
  colnames(tissue1)[1:100]
  print(paste0("Computing bicors for ", path))
  full_cors = bicorAndPvalue(tissue1, tissue1, use = 'p')
  cor_table = reshape2::melt(full_cors$bicor)
  new_p = reshape2::melt(full_cors$p)
  
  colnames(cor_table) = c('gene_tissue_1', 'gene_tissue_2', 'bicor')
  #can drop here to clear CPU
  full_cors=NULL
  
  cor_table$pvalue = signif(new_p$value, 3)
  cor_table$bicor = round(cor_table$bicor, 3)
  cor_table$qvalue = signif(p.adjust(cor_table$pvalue, "BH"), 3)
  cor_table = cor_table[order(cor_table$qvalue, decreasing=F),]
  cor_table = na.omit(cor_table)
  cor_table$gene_symbol_1 = gsub("\\_.*","",cor_table$gene_tissue_1)
  cor_table$tissue_1 = gsub(".*_","",cor_table$gene_tissue_1)
  cor_table = cor_table[!is.na(cor_table$tissue_1),]
  
  cor_table$gene_symbol_2 = gsub("\\_.*","",cor_table$gene_tissue_2)
  cor_table$tissue_2 = gsub(".*_","",cor_table$gene_tissue_2)
  cor_table = cor_table[!is.na(cor_table$tissue_2),]
  
  cor_table$Category = "pathway.genes"
  cor_table$kegg_ENTRY = ENTRY
  cor_table$kegg_NAME = NAME
  cor_table = cor_table[!cor_table$gene_tissue_1==cor_table$gene_tissue_2,]
  
  # Bind pathway cors and rand cors
  plot_table = rbind(cor_table, rndcor_table)
 
  
  # free up some memory 
  rm(cor_table)
  rm(rndcor_table)
  
  
  
  ### ***** PLOT ****
  print(paste0("Creating plots for ", path))
  
  gg.object1 = plot_cor_counts(.01, plot_table)
  gg.object2 = plot_cor_counts(.0001, plot_table)
  
  #save(gg.object, file = here( paste0("plots/kegg pathway corrs vs random genes/ordered tissues/gg objects/p 0.01 ", ENTRY, ".Rdata")) )
  
  ggsave(plot = gg.object1, here( paste0("plots/kegg pathway corrs vs random genes/ordered tissues/pdfs/p 0.01/p 0.01 ", ENTRY, ".pdf")), width = 8, height = 6)
  ggsave(plot = gg.object2, here( paste0("plots/kegg pathway corrs vs random genes/ordered tissues/pdfs/p 0.0001/p 0.0001 ", ENTRY, ".pdf")), width = 8, height = 6)
  
  print(paste0("Completed ", ENTRY, "!"))

  # shouldnt rnd and cor have the same number of rows?
  
}





#### 3.0 - Code used to set tissue factor levels for the figure ###
# Saved it into memory as an .rdata file after determining the rank order for one of the KEGG pathways

# set factor levels - established this based on hsa04062. 
# saved the levels as an r.data and we will use these moving forard when we loop the script
tissue_factors = cor_table %>%
  filter(pvalue < .01, bicor!= 0) %>% 
  dplyr::group_by(tissue_1) %>%
  dplyr::summarise(count = n()) %>%
  arrange(-count)
  
tissue_factors = factor( tissue_factors$tissue_1, levels = tissue_factors$tissue_1)

print(tissue_factors)
# save(tissue_factors, file = "./data/tissue_factors.Rdata")









