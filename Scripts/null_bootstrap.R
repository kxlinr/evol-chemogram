# Null Bootstrap

library(here)
library(plyr) #for match_df to get acc score from rank order
library(see) #for the half violins
library(tidyverse)

conflicted::conflict_prefer("here", "here")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("mutate", "dplyr")
conflicted::conflict_prefer("select", "dplyr")

#if running in the cluster, uncomment this
# n_sigs = args[1]

#if running locally, use this
n_sigs = 10

# Load Sigs
sig_files = list.files(here("Data", "Signatures"))
sigs = lapply(sig_files, function(x) readRDS(here("Data", "Signatures", x)))
names(sigs) = gsub("\\_.*", "", sig_files)
names(sigs)
sig_drugs = c("Cisplatin", "Cytarabine", "5-Fluorouracil", "Gemcitabine", "Irinotecan",
              "Luminespib", "Paclitaxel", "Topotecan", "Vinblastine", "Vorinostat")
rm(sig_files)

if(n_sigs == 3){
  #just the 3
  sigs = sigs[c(1,3,4)]
  sig_drugs=sig_drugs[c(1,3,4)]
  print(sig_drugs)
} else {
  print(sig_drugs)
}




# Load and Prep Score Table
#this can be generated from `chemogram_pipeline.Rmd` under the "Score Accuracy" section
acc_scores = readRDS(here("Results", paste0(length(sigs),"sig_score_table.rds")))


# Load Expression Data
exprs_norm = readRDS(here("Data", "normalized_exprs.rds"))


# Load Chemogram Results
fit_results = readRDS(here("Results",
                           paste0(length(sigs), "sig_nlme_fitted_data.rds")))

scored_ranks = readRDS(here("Results", 
                            paste0(length(sigs),"sig_accuracy_scores.rds"))) %>%
  mutate(s_score = accuracy) %>%
  select(-accuracy)

chemogram_results = na.omit(left_join(fit_results, scored_ranks, by="cell_line"))
unique(chemogram_results$subtype.x==chemogram_results$subtype.y) #double check this is true
chemogram_results = chemogram_results %>%
  mutate(subtype=subtype.x) %>%
  select(-subtype.x, -subtype.y)



# Bootstrap

## Helper Functions

### Generate Null Signatures
set_null_sigs = function(n, i, sigs, exprs){
  #set seed to n + i so its different each bootstrap iteration
  set.seed(n+i)
  
  #initialize list to store null sigs in
  null_sigs = list()
  
  #for each of the 10 sigs, take length and randomly select genes from colnames of exprs_norm[,-1]
  for(j in 1:length(sigs)){
    null_sigs[[j]] = sample(colnames(exprs)[2:ncol(exprs)], size = length(sigs[[j]]), 
                            replace = FALSE, prob = NULL)
    names(null_sigs)[j] = names(sigs)[j]
  }
  return(null_sigs)
}


### Generate Sig Scores
calc_sig_score = function(cl_surv, exprs, null_sigs){
  #loop through all cell lines in `survivals_full`
  cl = unique(cl_surv$cell_line)
  #initialize column to store null sig scores
  cl_surv = cl_surv %>% select(-s_score) #just using this df as a structural base, dont need this info
  cl_surv$n_score=NA #initialize a column to store null sig scores in
  cl_surv = do.call("rbind", replicate(length(null_sigs), #rep. all rows for the number of drugs we have
                                       cl_surv, simplify = FALSE)) 
  cl_surv = cl_surv[order(cl_surv$cell_line),]
  cl_surv$drug = rep(names(null_sigs), length(unique(cl_surv$cell_line)))
  
  
  
  #isolate row for the cell line
  for(j in 1:length(cl)){
    # print(j)
    exprs_subset = exprs %>%
      filter(COSMIC_ID==cl[j])
    
    if(nrow(exprs_subset) != 0){ #make sure there is expression data for the cell line
      #calc sig score for each of the 10 drugs
      if(table(cl_surv$cell_line)[j]==length(sigs)){ #make sure there are results for 10 drugs, no less
        for(k in 1:length(null_sigs)){
          #select just the columns with the sig genes
          exprs_sig_subset = exprs_subset[colnames(exprs_subset) %in% null_sigs[[k]]]
          
          #calc sig score (median of expression z-score)
          cl_surv$n_score[(cl_surv$cell_line==cl[j] & 
                             cl_surv$drug==names(null_sigs)[k])] = as.numeric(apply(exprs_sig_subset, 1, median))
        }#repeat for all 10 sigs
      } else {
        for(k in 1:length(null_sigs)){
          cl_surv$n_score[(cl_surv$cell_line==cl[j] &
                             cl_surv$drug==names(null_sigs)[k])] = NA
        }
        
      }
    } else {
      for(k in 1:length(null_sigs)){
        cl_surv$n_score[(cl_surv$cell_line==cl[j] &
                           cl_surv$drug==names(null_sigs)[k])] = NA
      }
    }
    
    
    
  }#repeat for all cell lines
  return(cl_surv)
}


### Calculate Accuracy Score
calc_acc_score = function(n_scores, acc_scores, acc_results){
  #initialize the column to store accuracy of null-sig chemogram
  colnames(acc_results)[which(colnames(acc_results)=="s_score")] = "s_accuracy"
  colnames(acc_results)[which(colnames(acc_results)=="drug_abbrev")] = "drug"
  acc_results$n_accuracy = NA
  
  n_scores = na.omit(n_scores)
  n_scores = left_join(n_scores, 
                       acc_results%>%select(cell_line, drug, surv),
                       by=c("cell_line", "drug"))
  
  
  #loop through 1 cell line at a time and subset data for each per loop iteration
  cl = unique(n_scores$cell_line)
  for(j in 1:length(cl)){
    n_scores_subset = n_scores %>% filter(cell_line==cl[j])
    
    #set DF in order of highest to lowest null score
    n_scores_subset = n_scores_subset[order(n_scores_subset$n_score, decreasing=TRUE),]
    
    #make column of n_rank and number from 1:10
    n_scores_subset$n_rank = seq(1:length(sigs))
    
    #set the DF in order of lowest to highest survival
    n_scores_subset = n_scores_subset[order(n_scores_subset$surv, decreasing=FALSE),]
    
    #store the new order of n_rank in a vector
    n_rank = as.data.frame(t(n_scores_subset$n_rank))
    
    #compare it to the acc_score matrix and retrieve the associated accuracy score
    n_score = suppressMessages(inner_join(acc_scores, n_rank)) %>%
      select(score)
    
    #store the acc score in a DF with cell line, cancer subtype, and n_acc
    acc_results$n_accuracy[acc_results$cell_line==cl[j]] = n_score
    
  }#reiterate for each cell line
  return(acc_results)
}



## Wrapper function
sig_bootstrap = function(n, i, sigs, exprs, cl_surv, acc_results, acc_scores){
  print(paste("Generating null sigs for run", i))
  null_sigs = set_null_sigs(n, i, sigs, exprs)
  
  print(paste("Calculating null sig scores for run", i))
  n_sig_scores = calc_sig_score(cl_surv, exprs, null_sigs)
  
  print(paste("Calculating null accuracy scores for run", i))
  scored_boot = calc_acc_score(n_scores=n_sig_scores, acc_scores, acc_results)
  
  return(scored_boot)
}





## Run bootstrap
n=1000

for(i in 1:n){
  print(i)
  if(i == 1){
    boot_results = sig_bootstrap(n, i,
                                 sigs=sigs,
                                 exprs=exprs_norm,
                                 cl_surv=scored_ranks,
                                 acc_results=chemogram_results,
                                 acc_scores=acc_scores)
    boot_results_full = boot_results[!grepl("NA", boot_results$n_accuracy),]
    boot_results_full$run = i
    
    boot_results_all = boot_results_full
  } else {
    boot_results = sig_bootstrap(n, i,
                                 sigs=sigs,
                                 exprs=exprs_norm,
                                 cl_surv=scored_ranks,
                                 acc_results=chemogram_results,
                                 acc_scores=acc_scores)
    boot_results_full = boot_results[!grepl("NA", boot_results$n_accuracy),]
    boot_results_full$run = i
    
    boot_results_all = rbind(boot_results_all, boot_results_full)
  }
  rm(boot_results, boot_results_full)
} #next iteration, i of n




# Calculate subtype average scores
subtype_vec = unique(boot_results_all$subtype)

boot_subtype_scores = data.frame(matrix(ncol=4, nrow=(n)))
colnames(boot_subtype_scores) = c("subtype", "n_accuracy", "s_accuracy", "run")

for(i in 1:length(subtype_vec)){ #go through each subtype
  result_subset = boot_results_all %>%
    filter(subtype==as.character(subtype_vec[i]))
  for(j in 1:n){ #go through each bootstrap run
    result_subset_run = result_subset %>%
      filter(run==j)
    boot_subtype_scores$subtype[j] = as.character(subtype_vec[i])
    boot_subtype_scores$run[j] = j
    boot_subtype_scores$n_accuracy[j] = mean(as.numeric(result_subset_run$n_accuracy))
    boot_subtype_scores$s_accuracy[j] = mean(result_subset_run$s_accuracy)
  }
  if(i==1){
    boot_subtype_results = boot_subtype_scores
  } else {
    boot_subtype_results = rbind(boot_subtype_results, boot_subtype_scores)
  }
}

rm(result_subset, result_subset_run, boot_subtype_scores)


saveRDS(boot_results_all, here("Results", paste0("null_", n, "_bootstrap_", 
                                                 length(sigs), "sig_indiv_scores.rds")))

saveRDS(boot_subtype_results, here("Results", paste0("null_", n, "_bootstrap_", 
                                                     length(sigs), "sig_subtype_scores.rds")))


#reload bootstrap data
boot_results_all = readRDS(here("Results", paste0("null_", n, "_bootstrap_", 
                                                  length(sigs), "sig_indiv_scores.rds")))
boot_subtype_results = readRDS(here("Results", paste0("null_", n, "_bootstrap_", 
                                                      length(sigs), "sig_subtype_scores.rds")))



################################################################################
# Reformat Chemogram Data
#create a new df with one row per subtype, and columns = subtype, proportion_correct
subtype = unique(scored_ranks$subtype)
subtype_scores = data.frame(subtype=subtype, accuracy=0, n = 0)

for (i in 1:length(subtype)){
  #Subset data per subtype
  subset = na.omit(scored_ranks[which(scored_ranks$subtype == subtype[i]),])
  
  #Store subtype name
  subtype_scores$subtype[i] = subtype[i]
  
  #calc and store % of cell line survivals correctly predicted
  subtype_scores$accuracy[i] = round((sum(subset$s_score)/nrow(subset)), digits=3)
  
  #Indicate number of cell lines in subtype (n)
  subtype_scores$n[i] = nrow(subset)
}
rm(subset)

#remove unclassified cell lines
chemogram_results = chemogram_results %>% filter(subtype != "UNCLASSIFIED") #-131
scored_ranks = na.omit(scored_ranks) %>% filter(subtype != "UNCLASSIFIED") #0
boot_results_all = boot_results_all %>% filter(subtype != "UNCLASSIFIED") #-130
boot_subtype_results = boot_subtype_results %>% filter(subtype != "UNCLASSIFIED")
subtype_scores = subtype_scores %>% filter(subtype != "UNCLASSIFIED")

#order the rows by highest to lowest correct
subtype_scores = subtype_scores[order(subtype_scores$accuracy),]

#To maintain this order when we plot later, factorize the column
subtype_scores$subtype <- factor(subtype_scores$subtype,
                                 ordered = TRUE,
                                 levels = subtype_scores$subtype)


#Set order to = subtype_scores
cancer_order = (as.character(subtype_scores$subtype))

scored_ranks$subtype <- factor(scored_ranks$subtype,
                               ordered = TRUE,
                               levels = c(cancer_order))
scored_ranks = scored_ranks[order(scored_ranks$subtype),]

#reorder cancer types
subtype_scores$subtype <- factor(subtype_scores$subtype,
                                 ordered = TRUE,
                                 levels = c(cancer_order))
subtype_scores = subtype_scores[order(subtype_scores$subtype),]

boot_results_all$subtype <- factor(boot_results_all$subtype,
                                   ordered = TRUE,
                                   levels = cancer_order)
boot_results_all = boot_results_all[order(boot_results_all$subtype),]

boot_subtype_results$subtype <- factor(boot_subtype_results$subtype,
                                       ordered = TRUE,
                                       levels = cancer_order)
boot_subtype_results = boot_subtype_results[order(boot_subtype_results$subtype),]




# Plots

## Reformat Data
half_box_data = boot_results_all %>% #label all the null results
  mutate(`Prediction Method`="Random Signatures") %>%
  mutate(accuracy = n_accuracy) %>%
  select(-run, -s_accuracy, -n_accuracy, -surv, -drug)
half_box_data_temp = scored_ranks %>%
  mutate(`Prediction Method`=paste0(length(sigs), "-sig Chemogram"),
         accuracy = s_score) %>%
  select(-s_score)
half_box_data = rbind(half_box_data, half_box_data_temp)
rm(half_box_data_temp)

half_box_data$subtype <- factor(half_box_data$subtype,
                                ordered = TRUE,
                                levels = cancer_order)

half_box_data$accuracy = as.numeric(half_box_data$accuracy)


## Plot Half-boxplot beeswarms
ggplot(data=half_box_data, aes(x=subtype, y=accuracy)) +
  
  geom_boxplot(width=0.6, aes(fill=`Prediction Method`, col=`Prediction Method`), alpha=0.8, fatten=NULL) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y.., group=`Prediction Method`, col=`Prediction Method`), #add in a line to indicate the mean
               width = 0.6, size = 1.5, linetype = "solid",
               position=position_dodge(preserve="total")) +
  
  scale_fill_manual(values=c("#e63946", "#a8dadc")) + #change boxplot color
  scale_color_manual(values=c("#6a040f", "#03045e")) + #change mean bar color
  
  geom_label(data=subtype_scores, aes(x = subtype, y = 0, label = n), #add n
             label.padding = unit(0.15, "lines"), label.size = 0.4, size = 8) +
  
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.25,.5,.75,1))+ #set y axis breaks
  scale_x_discrete(limits = cancer_order) +
  
  labs(title = paste0("Predictive Accuracy of Random Signatures vs ",length(sigs),"-sig Chemogram"), #titles and axis labels
       subtitle = paste0(length(unique(half_box_data$cell_line))," Cell Lines, ", n, " Bootstrap Iterations"),
       y = "Predictive Accuracy", x = "Disease Site") +
  theme_bw(base_size = 15) + #theme and sizing
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle=30, hjust=0.9),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 15))

ggsave(here("Plots", paste0("null", n, "_", length(sigs),"sig_halfbox.png")),
       width=15, height=9)


#Plot all results w/o subtype categorization
##Format data to plot
plot_null = as.data.frame(boot_results_all %>% select(n_accuracy))
plot_null$`Prediction Method` = "Random Signatures"
colnames(plot_null)[1] = "Predictive Accuracy"
plot_null$`Predictive Accuracy` = as.numeric(plot_null$`Predictive Accuracy`)

plot_real = as.data.frame(scored_ranks$s_score)
plot_real$`Prediction Method` = paste0(length(sigs), "-sig Chemogram")
colnames(plot_real)[1] = "Predictive Accuracy"

#t.test b/w null and real
p=wilcox.test(plot_null$`Predictive Accuracy`, plot_real$`Predictive Accuracy`)

##Plot
ggplot(plot_null, aes(x=`Prediction Method`, y=`Predictive Accuracy`)) +
  geom_violin(fill = "#a8dadc") +
  geom_boxplot(width = 0.15, fatten=NULL) +
  geom_violin(data=plot_real, aes(x=`Prediction Method`, y=`Predictive Accuracy`), fill="#e63946") +
  geom_boxplot(data=plot_real, aes(x=`Prediction Method`, y=`Predictive Accuracy`),
               width = 0.15, fatten=NULL) +
  
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y.., group=`Prediction Method`), #add in a line to indicate the mean
               width = 0.15, size = 1, linetype = "solid",
               position=position_dodge(preserve="total")) +
  stat_summary(data=plot_real, fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y.., group=`Prediction Method`), #add in a line to indicate the mean
               width = 0.15, size = 1, linetype = "solid",
               position=position_dodge(preserve="total")) +
  
  annotate("text", x = 1.5, y = 1.05, label = paste0("Wilcoxon test: p = ", round(p[["p.value"]], 6)), size = 5) +
  
  labs(title = paste0("Predictive Accuracy of Random Signatures vs ", length(sigs), "-sig Chemogram"), #titles and axis labels
       subtitle = paste0(length(unique(chemogram_results$cell_line))," Cell Lines, ", n, " Bootstrap Iterations"),
       y = "Predictive Accuracy", x = "Prediction Method") +
  
  theme_bw(base_size = 15) + #theme and sizing
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 10))

ggsave(here("Plots", paste0("null", n, "_", length(sigs),"sig_violins.png")),
       width = 15, height = 9) #save plot


#random prediction acc scores are uniformly distributed
boot_results_all$accuracy = as.numeric(boot_results_all$n_accuracy)
ggplot(boot_results_all, aes(x=accuracy)) +
  geom_histogram()

# plot violins with mean red points
ggplot(data=boot_subtype_results, aes(x=subtype, y=n_accuracy)) +
  geom_violin(fill="#71AEC1", trim=FALSE, alpha=1) +
  geom_boxplot(width=0.1) +
  geom_point(data=subtype_scores, aes(x=subtype, y=accuracy), color="#e63946", size=5, alpha = 0.9) +
  geom_label(data=subtype_scores, aes(x=subtype, y = 0, label = n),
             label.padding = unit(0.15, "lines"), label.size = 0.4, size = 8) +
  #coord_flip() +
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.25,.5,.75,1))+
  labs(title = paste0("Predictive Accuracy of Random Signatures vs ",length(sigs),"-sig Chemogram"), #titles and axis labels
       subtitle = paste0(length(unique(chemogram_results$cell_line))," Cell Lines, ", n, " Bootstrap Iterations"),
       y = "Average Predictive Accuracy", x = "Disease Site") +
  
  theme_bw(base_size = 15) + #theme and sizing
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle=30, hjust=0.9),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 10))

ggsave(here("Plots", paste0("null", n, "_", length(sigs),"sig_halfbox_avgs.png")),
       width = 15, height = 9) #save plot


## Reformat to Plot Epithelial Cancers Only

epi = c("HNSC", "ESCA", "BRCA",
        "COREAD", "LIHC", #"ACC", #leave this out bc it only has one cell line
        "STAD", "KIRC", "LUAD",
        "LUSC","MESO",
        "PAAD", "THCA", "BLCA",
        "CESC", "UCEC", "OV", "PRAD")

half_box_data_epi = half_box_data %>% filter(subtype %in% epi)
subtype_scores_epi = subtype_scores %>% filter(subtype %in% epi)
scored_ranks_full_epi = scored_ranks %>% filter(subtype %in% epi)
boot_subtype_scores_epi = boot_subtype_results %>% filter(subtype %in% epi)
boot_results_epi = boot_results_all %>% filter(subtype %in% epi)

cancer_order_epi = subtype_scores_epi$subtype

half_box_data_epi$subtype <- factor(half_box_data_epi$subtype,
                                    ordered = TRUE,
                                    levels = cancer_order_epi)
subtype_scores_epi$subtype <- factor(subtype_scores_epi$subtype,
                                     ordered = TRUE,
                                     levels = cancer_order_epi)
scored_ranks_full_epi$subtype <- factor(scored_ranks_full_epi$subtype,
                                        ordered = TRUE,
                                        levels = cancer_order_epi)
boot_subtype_scores_epi$subtype <- factor(boot_subtype_scores_epi$subtype,
                                          ordered = TRUE,
                                          levels = cancer_order_epi)
boot_results_epi$subtype <- factor(boot_results_epi$subtype,
                                   ordered = TRUE,
                                   levels = cancer_order_epi)

## Plot epi
ggplot(data=half_box_data_epi, aes(x=subtype, y=accuracy)) +
  
  geom_boxplot(width=0.6, aes(fill=`Prediction Method`, col=`Prediction Method`), alpha=0.8, fatten=NULL) +
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y.., group=`Prediction Method`, col=`Prediction Method`), #add in a line to indicate the mean
               width = 0.6, size = 1.5, linetype = "solid",
               position=position_dodge(preserve="total")) +
  
  scale_fill_manual(values=c("#e63946", "#a8dadc")) + #change boxplot color
  scale_color_manual(values=c("#6a040f", "#03045e")) + #change mean bar color
  
  geom_label(data=subtype_scores_epi, aes(x = subtype, y = 0, label = n), #add n
             label.padding = unit(0.15, "lines"), label.size = 0.4, size = 8) +
  
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.25,.5,.75,1))+ #set y axis breaks
  scale_x_discrete(limits = cancer_order_epi) +
  
  labs(title = paste0("Predictive Accuracy of Random Signatures vs ",length(sigs),"-sig Chemogram"), #titles and axis labels
       subtitle = paste0(length(unique(half_box_data_epi$cell_line))," Epithelial Cell Lines, ", n, " Bootstrap Iterations"),
       y = "Predictive Accuracy", x = "Disease Site") +
  theme_bw(base_size = 15) + #theme and sizing
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle=30, hjust=0.9),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 15))

ggsave(here("Plots", "Epithelial Cancers", 
            paste0("null", n, "_", length(sigs),"sig_halfbox_epi.png")),
       width=15, height=9)


#Plot all results w/o subtype categorization
##Format data to plot
plot_null_epi = boot_results_epi %>%
  select(n_accuracy)
plot_null_epi$`Prediction Method` = "Random Signatures"
colnames(plot_null_epi)[1] = "Predictive Accuracy"
plot_null_epi$`Predictive Accuracy` = as.numeric(plot_null_epi$`Predictive Accuracy`)

plot_real_epi = as.data.frame(scored_ranks_full_epi$s_score)
plot_real_epi$`Prediction Method` = "10-sig Chemogram"
colnames(plot_real_epi)[1] = "Predictive Accuracy"

#t.test b/w null and real
p_epi = wilcox.test(plot_null_epi$`Predictive Accuracy`, plot_real_epi$`Predictive Accuracy`)


##Plot
ggplot(plot_null_epi, aes(x=`Prediction Method`, y=`Predictive Accuracy`)) +
  geom_violin(fill = "#a8dadc") +
  geom_boxplot(width = 0.15, fatten=NULL) +
  geom_violin(data=plot_real_epi, aes(x=`Prediction Method`, y=`Predictive Accuracy`), fill="#e63946") +
  geom_boxplot(data=plot_real_epi, aes(x=`Prediction Method`, y=`Predictive Accuracy`),
               width = 0.15, fatten=NULL) +
  
  stat_summary(fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y.., group=`Prediction Method`), #add in a line to indicate the mean
               width = 0.15, size = 1, linetype = "solid",
               position=position_dodge(preserve="total")) +
  stat_summary(data=plot_real_epi, fun = mean, geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y.., group=`Prediction Method`), #add in a line to indicate the mean
               width = 0.15, size = 1, linetype = "solid",
               position=position_dodge(preserve="total")) +
  
  annotate("text", x = 1.5, y = 1, label = paste0("Wilcoxon test: p = ", round(p_epi[["p.value"]], 9)), size = 5) +
  
  labs(title = paste0("Predictive Accuracy of Random Signatures vs ",length(sigs),"-sig Chemogram"), #titles and axis labels
       subtitle = paste0(nrow(scored_ranks_full_epi)," Epithelial Cell Lines, ", n, " Bootstrap Iterations"),
       y = "Predictive Accuracy", x = "Prediction Method") +
  
  theme_bw(base_size = 15) + #theme and sizing
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 10))

ggsave(here("Plots", "Epithelial Cancers", 
            paste0("null", n, "_", length(sigs),"sig_violins_epi.png")),
       width = 15, height = 9) #save plot


# plot violins with mean red points
boot_subtype_scores_epi$accuracy = as.numeric(boot_subtype_scores_epi$n_accuracy)

ggplot(data=boot_subtype_scores_epi, aes(x=subtype, y=accuracy)) +
  geom_violin(fill="#71AEC1", trim=FALSE, alpha=1) +
  geom_boxplot(width=0.1) +
  geom_point(data=subtype_scores_epi, aes(x=subtype, y=accuracy), color="#e63946", size=5, alpha = 0.9) +
  geom_label(data=subtype_scores_epi, aes(x=subtype, y = 0, label = n),
             label.padding = unit(0.15, "lines"), label.size = 0.4, size = 8) +
  #coord_flip() +
  scale_y_continuous(limits=c(0,1), breaks = c(0,0.25,.5,.75,1))+
  labs(title = paste0("Predictive Accuracy of Random Signatures vs ",length(sigs),"-sig Chemogram"), #titles and axis labels
       subtitle = paste0(nrow(scored_ranks_full_epi)," Cell Lines, ", n, " Bootstrap Iterations"),
       y = "Average Predictive Accuracy", x = "Disease Site") +
  
  theme_bw(base_size = 15) + #theme and sizing
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle=30, hjust=0.9),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(size = 10))

ggsave(here("Plots", "Epithelial Cancers", 
            paste0("null", n, "_", length(sigs),"sig_halfbox_avgs_epi.png")),
       width = 15, height = 9) #save plot






# Significance tests
all_p = as.data.frame(matrix(nrow=length(cancer_order), ncol=2))
colnames(all_p) = c("subtype", "p_val")
class(boot_results_all$n_accuracy) = "numeric"

for(i in 1:length(cancer_order)){
  type = cancer_order[i]
  if(nrow(boot_results_all %>% filter(subtype == type)) > 1 & nrow(scored_ranks %>% filter(subtype == type)) > 1){
    result = wilcox.test(x=as.matrix(boot_results_all[which(boot_results_all$subtype==type), "n_accuracy"]),
                    y=as.matrix(scored_ranks[which(scored_ranks$subtype==type), "s_score"]),
                    paired=FALSE, alternative="two.sided")
    
    
    
    all_p$subtype[i] = type
    all_p$p_val[i] = result$p.value
  } else {
    all_p$subtype[i] = type
    all_p$p_val[i] = NA
  }
  
}

#label if subtype is epithelial origin or not
all_p = all_p %>%
  mutate(class = case_when(subtype %in% epi ~ "epi", 
                           TRUE ~ "non-epi"))
#add number of cell lines in each subtype
all_p = merge(all_p, subtype_scores[,c(1,3)], by="subtype")
#add if results are statistically significant +- bonferroni-adjusted alpha = 0.05
all_p = all_p %>%
  mutate(standard_signif = case_when(p_val <= 0.05 ~ TRUE,
                                     TRUE~FALSE),
         p_adj = p.adjust(p_val, method = "fdr"), #calculate adjusted pval using fdr
         fdr_signif = case_when(p_adj <= 0.05 ~ TRUE, #determine if padj is significant or not
                                       TRUE ~ FALSE),
         bonferroni_signif = case_when(p_val <= 0.05/n ~ TRUE, #determine if pval is significant or not using bonferroni correction
                                    TRUE ~ FALSE),)

#organize by decreasing significance
all_p = all_p[order(all_p$p_val, decreasing=FALSE),]
#reorder columns
all_p = all_p %>%
  select(subtype, class, n, 
         p_val, p_adj, 
         standard_signif, fdr_signif, bonferroni_signif)

#filter epithelial results
epi_p = all_p %>% filter(subtype %in% epi)

#save results
readr::write_csv(all_p, here("Results",
                            paste0("null", n, "_",
                              length(sigs), "_pvals.csv")))
