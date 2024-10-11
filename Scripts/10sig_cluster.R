library(here)
library(gdscIC50) #for cleaning and normalizing the raw drug response data
library(pheatmap) #for plotting survival and sig scores
library(gtools) #for generating the permutations/scoring table
library(plyr) #for match_df to get acc score from rank order
library(patchwork) #for plotting
library(tidyverse)

conflicted::conflict_prefer("here", "here")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("mutate", "dplyr")
conflicted::conflict_prefer("select", "dplyr")

cell_lines = readRDS(here("Data", "cleaned_cell_line_meta.rds"))
subtype_key = readRDS(here("Data", "subtype_key.rds"))
# exprs_norm_all = readRDS(here("Data", "normalized_exprs.rds"))
exprs_norm = readRDS(here("Data", "normalized_exprs_to_gapdh.rds"))
drug_id_key = readRDS(here("Data", "drug_id_key.rds"))
norm_dose_resp = readRDS(here("Data", "prepped_drug_resp.rds"))


sig_files = list.files(here("Data", "Signatures"))
sigs = lapply(sig_files, function(x) readRDS(here("Data", "Signatures", x)))
names(sigs) = gsub("\\_.*", "", sig_files) #cut off everything after the first underscore

## we also need the full drug names to match with the correct drug response data. those need to be manually entered in the *same order* as the list shows
#make sure theyre spelled right! can check the `drug_id_key` table for spelling
names(sigs) #print this to make sure the order is correct
sig_drugs = c("Cisplatin", "Cytarabine", "5-Fluorouracil", "Gemcitabine", "Irinotecan",
              "Luminespib", "Paclitaxel", "Topotecan", "Vinblastine", "Vorinostat")

rm(sig_files)

sigs[["cis"]] = c("LRRC8C", "LY6K", "MMP10","SLFN11", "STOML2","USP31","WDR3", "ZNF750")
sigs[["gem"]] = c("CRLF3", "MLKL", "POLR3G", "PROCR", "RELB", "SLFN11")


fit_results = readRDS(here("Results",
                           paste0(length(sigs), "sig_nlme_fitted_data.rds")))

cl_order = subtype_key[order(subtype_key$subtype), "cell_line"]





print("##### PREDICTING #####")

#F(X) to subset genes within each drug sensitivity signature, then derive the sensitivity score per cell line
sig_score = function(sig_genes, exprs, drug_abbrev){
  
  exprs_norm2 = exprs %>%
    dplyr::select(COSMIC_ID, all_of(sig_genes))
  
  #Find signature score for each cell line
  sig_exprs = exprs_norm2[,-1] %>% 
    mutate(score = apply(.,1,median), cosmic_id = exprs$COSMIC_ID)
  
  #clean up result
  score = sig_exprs %>% dplyr::select(cosmic_id, score)
  rownames(score) = sig_exprs$cosmic_id
  colnames(score) = c("cell_line", "score")
  
  score$drug_abbrev = drug_abbrev
  
  return(score)
}


#rm "DATA." before each cosmic id
exprs_norm = exprs_norm %>%
  mutate(COSMIC_ID = substr(COSMIC_ID, 6, nchar(COSMIC_ID)))

# #Run the above function for each sig
scores = list()

for(i in 1:length(sigs)){
  scores[[i]] = sig_score(sig_genes=sigs[[i]],
                          exprs=exprs_norm,
                          drug_abbrev=names(sigs)[i])
  names(scores)[i] = names(sigs)[i]
}



#Label cancer type of each cell line
for(i in 1:length(scores)){
  #label subtypes
  scores[[i]] = full_join(scores[[i]], subtype_key, by="cell_line")
  
  #replace NAs with "unclassified"
  scores[[i]]$subtype = scores[[i]]$subtype %>%
    replace_na('UNCLASSIFIED') #replace any more cancer subtype NAs w/ unclassif.
}
scores_df = do.call(rbind, scores)

scores_df_wide = pivot_wider(scores_df,
                             id_cols = c(cell_line, subtype),
                             names_from = drug_abbrev,
                             values_from = score,
                             values_fn = mean #for duplicate scores on the same cl+drug, use the mean
)

scores = scores_df_wide
rm(scores_df, scores_df_wide)
scores = scores %>% select(-`NA`)

#Save
saveRDS(scores, here("Results",
                     paste0(length(sigs), "sig_scores.rds")))

scores = na.omit(scores)

scores_long = pivot_longer(scores,
                           cols = c(names(sigs)),
                           names_to = "drug_abbrev",
                           values_to = "score")


plot_score_ht = left_join(cl_order,scores,by="cell_line") %>% 
  select(-cell_line, -subtype)
pal = colorRampPalette(c("darkred", "#bb3e03", "#e9d8a6", "#94d2bd", "#0a9396", "#005f73"))

pheatmap::pheatmap(na.omit(plot_score_ht), #data to use
                   cluster_cols=FALSE, cluster_rows=FALSE, #specify if rows/cols should be clustered
                   angle_col=0, show_rownames = FALSE,
                   color = pal(100),
                   main="Signature Scores Across All Cell Lines", #title
                   legend_labels = "Sig Score", 
                   filename=here("Plots", "Data Prep", 
                                 paste0(length(sig_drugs), "sig_scores_CL.png"))) #Where to save the file

rm(plot_score_ht)


ggplot(scores_long, aes(x=score, col = drug_abbrev,)) +
  geom_density(size=1) +
  scale_color_manual(values = c("#9b2226","#f94149","#f3722c", "#f9c74f", "#90be6d",
                                "#43aa8b", "#277da1",  "#4361ee", "#7209b7", "#f72585")) +
  labs(x="Sensitivity Score",
       y="Frequency",
       color="Drug",
       title="Distribution of Signature-derived Sensitivity Scores",
       subtitle = paste0(length(unique(scores_long$cell_line)), " Cell Lines"))
ggsave(here("Plots" ,"Data Prep", 
            paste0(length(sigs), "drug_sigscore_dist.png")),
       width=8, height=6)

ggplot(scores_long, aes(x=drug_abbrev, y=score)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  labs(title="Distribution of Signature-derived Sensitivity Scores",
       subtitle = paste0(length(unique(scores_long$cell_line)), " Cell Lines"))
ggsave(here("Plots" ,"Data Prep", 
            paste0(length(sigs), "drug_sigscore_box.png")),
       width=8, height=6)




#Score

print("##### SCORING #####")

n = length(sigs)
acc_scores = readRDS(here("Results", paste0(length(sigs),"sig_score_table.rds")))


find_score = function(data, cl){
  
  subset = data %>% filter(cell_line==cl) #filter to just one cl at a time
  
  subset = subset[order(subset$surv, decreasing = FALSE),] #put drugs in order from most response to least response
  if(nrow(subset)==length(sigs)){
    subset$true_rank = (seq(1, length(sigs))) #label the observed sensitivity rank
    
    subset = subset[order(subset$score, decreasing = TRUE),] #put drugs in order from most response to least response
    
    score = match_df(acc_scores, as.data.frame(t(subset$true_rank))) #find which score goes with the prediction ranking
    score=as.numeric(score$score)
  } else { #there may be cell lines without data for all the drugs, so those need to be filtered out here
    print("Missing drug")
    score=NA
  }
  
  return(score)
}


all_cl = unique(fit_results$cell_line)
unscored_ranks = left_join(fit_results, scores_long,
                           by=c("cell_line", "subtype", "drug_abbrev"))

scored_ranks = data.frame(cell_line = all_cl, subtype = NA, accuracy = 0)

for(i in 1:length(all_cl)){
  print(paste(i, "out of", length(all_cl), "--", ((i/length(all_cl))*100), "%"))
  
  scored_ranks$cell_line[i] = all_cl[i]
  
  scored_ranks$subtype[i] = as.character(unique(unscored_ranks %>%
                                                  filter(cell_line==all_cl[i]) %>%
                                                  select(subtype)))
  
  value = find_score(unscored_ranks, all_cl[i])
  print(value)
  scored_ranks$accuracy[i] = value
}


saveRDS(scored_ranks, here("Results", paste0(length(sigs),"sig_accuracy_scores.rds")))