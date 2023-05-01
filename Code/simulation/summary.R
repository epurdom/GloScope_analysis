library(stringr)
library(vegan)
myfun_clus = function(dir, distmat){
  dist_list = list()
  for(i in list.files(dir)){
    load(paste0(dir, i))
    if(distmat == "dist_mat_GMM_pca"){
      dist_list[[i]] = dist_mat_GMM_pca
    }else if(distmat == "dist_mat_KNN_pca"){
      dist_list[[i]] = dist_mat_KNN_pca
    }else if(distmat == "dist_mat_GMM_scvi"){
      dist_list[[i]] = dist_mat_GMM_scvi
    }else if(distmat == "dist_mat_KNN_scvi"){
      dist_list[[i]] = dist_mat_KNN_scvi
    }else if(distmat == "clus_seurat_pca"){
      dist_list[[i]] == clus_seurat_pca
    }else if(distmat == "clus_seurat_scvi"){
      dist_list[[i]] == clus_seurat_scvi
    }else if(distmat == "clus_true"){
      dist_list[[i]] == clus_true
    }
  }

  
  
  set.seed(1)
  pvals = c()
  astat = c()
  ave_b = c()
  ave_w = c()
  for(i in names(dist_list)){
    meta = data.frame(sample_id = rownames(dist_list[[i]]))
    meta$group_id = substr(meta$sample_id, nchar(meta$sample_id), nchar(meta$sample_id)+1)
    test = anosim(dist_list[[i]], meta$group_id)
    pvals = c(pvals, test$signif)
    astat = c(astat, test$statistic)
    
    
    
    dist_mat = as.dist(dist_list[[i]])
    
    dist_test = c(dist_mat)
    mycompare = c()
    for(j in 1:length(attr(dist_mat, "Labels"))){
      mycompare = c(mycompare, paste0(attr(dist_mat, "Labels")[j], "-", attr(dist_mat, "Labels")[-(1:j)]))
      }
    
    mycompare = mycompare[-length(mycompare)]
    
    compare_group = str_replace_all(mycompare, paste0("sample", "[0-9]+", "[.](?!\\d+$)"), "")
    compare_group = gsub("[.](?!\\d+$)", "", compare_group, perl=TRUE)
    
    group = ifelse(substr(compare_group, 1,1) == substr(compare_group, 3,3), "same", "diff")
    ave_b = c(ave_b, mean(dist_test[group == "diff"]))
    ave_w = c(ave_w, mean(dist_test[group == "same"]))
  }

  return(list("p" = pvals,
              "t" = astat,
              "b" = ave_b,
              "w" = ave_w))
}



power_data = data.frame("n" = 10,"m" = 5000, "lfc" = 0.05, "sigma" = 0.13, "rou_DE" = 0.1,
             "omega_k" = "(1,...,1)", "lambda" = "~8.25", "tau" = 0, "delta" = 0, "alpha" = 100,
             "Pi_A" = "(0.4, 0.05, 0.2, 0.05, 0.05, 0.05, 0.1, 0.1)",
             "Pi_B" = "(0.15, 0.05, 0.1, 0.4, 0.05, 0.1, 0.05, 0.1)",
             "Power_GMM_PCA" = NA,
             "Power_KNN_PCA" = NA,
             "Power_GMM_scvi" = NA,
             "Power_KNN_scvi" = NA, 
             "ANOSIM_GMM_PCA" = NA,
             "ANOSIM_KNN_PCA" = NA,
             "ANOSIM_GMM_scvi" = NA,
             "ANOSIM_KNN_scvi" = NA, 
             "AveB_GMM_PCA" = NA,
             "AveB_KNN_PCA" = NA,
             "AveB_GMM_scvi" = NA,
             "AveB_KNN_scvi" = NA, 
             "AveW_GMM_PCA" = NA,
             "AveW_KNN_PCA" = NA,
             "AveW_GMM_scvi" = NA,
             "AveW_KNN_scvi" = NA)


################# baseline

baseline_list = list.files("../../results/simulation/distmat/baseline/")
power_baseline <- power_data[rep(1,length(baseline_list)),]
rownames(power_baseline) = baseline_list

#power_baseline = power_baseline[c(1,2,5),]

for(i in rownames(power_baseline)){
  power_GMM_PCA = try(myfun_clus(paste0("../../results/simulation/distmat/baseline/",i, "/"), "dist_mat_GMM_pca"))
  Power_KNN_PCA = try(myfun_clus(paste0("../../results/simulation/distmat/baseline/",i, "/"), "dist_mat_KNN_pca"))
  Power_GMM_scvi = try(myfun_clus(paste0("../../results/simulation/distmat/baseline/",i, "/"), "dist_mat_GMM_scvi"))
  Power_KNN_scvi = try(myfun_clus(paste0("../../results/simulation/distmat/baseline/",i, "/"), "dist_mat_KNN_scvi"))
#  Power_cluster = try(myfun_clus(paste0("baseline/update/",i, "_scvi/Step2/"), "clusprop_diss"))
  
  power_baseline[i,"Power_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$p<0.05))
  power_baseline[i,"Power_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$p<0.05))
  power_baseline[i,"Power_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$p<0.05))
  power_baseline[i,"Power_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$p<0.05))
  
  
  power_baseline[i,"ANOSIM_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$t))
  power_baseline[i,"ANOSIM_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$t))
  power_baseline[i,"ANOSIM_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$t))
  power_baseline[i,"ANOSIM_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$t))

  power_baseline[i,"AveB_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$b))
  power_baseline[i,"AveB_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$b))
  power_baseline[i,"AveB_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$b))
  power_baseline[i,"AveB_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$b))

  power_baseline[i,"AveW_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$w))
  power_baseline[i,"AveW_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$w))
  power_baseline[i,"AveW_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$w))
  power_baseline[i,"AveW_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$w))
  
  
  if(str_detect(i, "same")){
    power_baseline[i,"Pi_B"] = "(0.4, 0.05, 0.2, 0.05, 0.05, 0.05, 0.1, 0.1)"
  }
  power_baseline[i,"sigma"] = as.numeric(gsub('\\D', '', unlist(str_split(i, "_"))[1]))/100
  power_baseline[i,"lfc"] = as.numeric(gsub('\\D', '', unlist(str_split(i, "_"))[2]))/100
}


write.csv(power_baseline, file = "../../results/simulation/power/baseline.csv")
################# offset variation

offvar_list = list.files("../../results/simulation/distmat/offvar/")
power_offvar <- power_data[rep(1,length(offvar_list)),]

rownames(power_offvar) = offvar_list

power_offvar$lfc = 0.15

for(i in rownames(power_offvar)){
  power_GMM_PCA = try(myfun_clus(paste0("../../results/simulation/distmat/offvar/",i, "/"), "dist_mat_GMM_pca"))
  Power_KNN_PCA = try(myfun_clus(paste0("../../results/simulation/distmat/offvar/",i, "/"), "dist_mat_KNN_pca"))
  Power_GMM_scvi = try(myfun_clus(paste0("../../results/simulation/distmat/offvar/",i, "/"), "dist_mat_GMM_scvi"))
  Power_KNN_scvi = try(myfun_clus(paste0("../../results/simulation/distmat/offvar/",i, "/"), "dist_mat_KNN_scvi"))

  power_offvar[i,"Power_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$p<0.05))
  power_offvar[i,"Power_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$p<0.05))
  power_offvar[i,"Power_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$p<0.05))
  power_offvar[i,"Power_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$p<0.05))
  
  
  power_offvar[i,"ANOSIM_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$t))
  power_offvar[i,"ANOSIM_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$t))
  power_offvar[i,"ANOSIM_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$t))
  power_offvar[i,"ANOSIM_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$t))
  
  power_offvar[i,"AveB_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$b))
  power_offvar[i,"AveB_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$b))
  power_offvar[i,"AveB_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$b))
  power_offvar[i,"AveB_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$b))
  
  power_offvar[i,"AveW_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$w))
  power_offvar[i,"AveW_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$w))
  power_offvar[i,"AveW_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$w))
  power_offvar[i,"AveW_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$w))
  
  if(str_detect(i, "offvar")){
    power_offvar[i,"lambda"] =as.numeric(gsub('\\D', '', unlist(str_split(i, "_"))[1]))
    power_offvar[i,"delta"] =as.numeric(gsub('\\D', '', unlist(str_split(i, "_"))[2]))
    
  }else{
    power_offvar[i,"delta"] =as.numeric(gsub('\\D', '', unlist(str_split(i, "_"))[1]))
    power_offvar[i,"tau"] =as.numeric(gsub('\\D', '', unlist(str_split(i, "_"))[2]))
    
  }
}

write.csv(power_offvar, file = "../../results/simulation/power/offset.csv")


################################ parameters





parameters = list.files("../../results/simulation/distmat/parameter/")
parameter_list = list()
for(i in parameters){
  tmp_files = list.files(paste0("../../results/simulation/distmat/parameter/",i))
  power_tmp = power_data[rep(1,length(tmp_files)),]
  power_tmp$Pi_B = power_tmp$Pi_A
  rownames(power_tmp) = tmp_files
  
  
  for(j in rownames(power_tmp)){
    power_GMM_PCA = try(myfun_clus(paste0("../../results/simulation/distmat/parameter/",i,"/",j, "/"), "dist_mat_GMM_pca"))
    Power_KNN_PCA = try(myfun_clus(paste0("../../results/simulation/distmat/parameter/",i,"/",j, "/"), "dist_mat_KNN_pca"))
    Power_GMM_scvi = try(myfun_clus(paste0("../../results/simulation/distmat/parameter/",i,"/",j, "/"), "dist_mat_GMM_scvi"))
    Power_KNN_scvi = try(myfun_clus(paste0("../../results/simulation/distmat/parameter/",i,"/",j, "/"), "dist_mat_KNN_scvi"))
    #  Power_cluster = try(myfun_clus(paste0("baseline/update/",i, "_scvi/Step2/"), "clusprop_diss"))
    power_tmp[j,"Power_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$p<0.05))
    power_tmp[j,"Power_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$p<0.05))
    power_tmp[j,"Power_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$p<0.05))
    power_tmp[j,"Power_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$p<0.05))
    
    
    power_tmp[j,"ANOSIM_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$t))
    power_tmp[j,"ANOSIM_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$t))
    power_tmp[j,"ANOSIM_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$t))
    power_tmp[j,"ANOSIM_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$t))
    
    power_tmp[j,"AveB_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$b))
    power_tmp[j,"AveB_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$b))
    power_tmp[j,"AveB_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$b))
    power_tmp[j,"AveB_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$b))
    
    power_tmp[j,"AveW_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$w))
    power_tmp[j,"AveW_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$w))
    power_tmp[j,"AveW_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$w))
    power_tmp[j,"AveW_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$w))
    
    if(i == "more_lfc"){
      power_tmp[j,"lfc"] = as.numeric(gsub('\\D', '', j))/100
    }else if(i == "more_s"){
      power_tmp[j, "n"] = as.numeric(gsub('\\D', '', unlist(str_split(j, "_"))[1]))
      power_tmp[j,"lfc"] = as.numeric(gsub('\\D', '', unlist(str_split(j, "_"))[2]))/100
    }else if(i == "more_pde"){
      power_tmp[j,"rou_DE"] = as.numeric(gsub('\\D', '', unlist(str_split(j, "_"))[1]))/10
      power_tmp[j,"lfc"] = as.numeric(gsub('\\D', '', unlist(str_split(j, "_"))[2]))/100
      
    }else if(i == "more_m"){
      power_tmp[j,"m"] = as.numeric(gsub('\\D', '', unlist(str_split(j, "_"))[1]))
      power_tmp[j,"lfc"] = as.numeric(gsub('\\D', '', unlist(str_split(j, "_"))[2]))/100
      
    }else if(i == "L_L"){
      power_tmp[j,"omega_k"] =  paste0("(",as.numeric(gsub('\\D', '', j))/5, ",1...,1)")
    }else if(i == "S_L"){
      power_tmp[j, "omega_k"] =  paste0("(1,1,1,1,",paste(rep(as.numeric(gsub('\\D', '', j))/0.5,4), collapse = ","), ")")
    }
  }
  parameter_list[[i]] = power_tmp
}

lfcpde_power = rbind(parameter_list[["more_lfc"]], parameter_list[["more_pde"]])
omega_power = rbind(parameter_list[["L_L"]], parameter_list[["S_L"]])
write.csv(omega_power, file = "../../results/simulation/power/omega_k.csv")
write.csv(lfcpde_power, file = "../../results/simulation/power/lfc_pde.csv")
write.csv(parameter_list[["more_s"]], file = "../../results/simulation/power/S_power.csv")


################# cluster proportion test

myfun_clus = function(dir, distmat){
  dist_list = list()
  for(i in list.files(dir)){
    load(paste0(dir, i))
    if(distmat == "dist_mat_GMM_pca"){
      dist_list[[i]] = dist_mat_GMM_pca
    }else if(distmat == "dist_mat_KNN_pca"){
      dist_list[[i]] = dist_mat_KNN_pca
    }else if(distmat == "dist_mat_GMM_scvi"){
      dist_list[[i]] = dist_mat_GMM_scvi
    }else if(distmat == "dist_mat_KNN_scvi"){
      dist_list[[i]] = dist_mat_KNN_scvi
    }else if(distmat == "clus_seurat_pca"){
      dist_list[[i]] = clus_seurat_pca
    }else if(distmat == "clus_seurat_scvi"){
      dist_list[[i]] = clus_seurat_scvi
    }else if(distmat == "clus_true"){
      dist_list[[i]] = clus_true
    }
  }
  
  
  
  set.seed(1)
  pvals = c()
  astat = c()
  ave_b = c()
  ave_w = c()
  for(i in names(dist_list)){
    meta = data.frame(sample_id = rownames(dist_list[[i]]))
    meta$group_id = substr(meta$sample_id, nchar(meta$sample_id), nchar(meta$sample_id)+1)
    test = anosim(dist_list[[i]], meta$group)
    pvals = c(pvals, test$signif)
    astat = c(astat, test$statistic)
    
    
    
    dist_mat = as.dist(dist_list[[i]])
    
    dist_test = c(dist_mat)
    mycompare = c()
    for(j in 1:length(attr(dist_mat, "Labels"))){
      mycompare = c(mycompare, paste0(attr(dist_mat, "Labels")[j], "-", attr(dist_mat, "Labels")[-(1:j)]))
    }
    
    mycompare = mycompare[-length(mycompare)]
    
    compare_group = str_replace_all(mycompare, paste0("sample", "[0-9]+", "[.](?!\\d+$)"), "")
    compare_group = gsub("[.](?!\\d+$)", "", compare_group, perl=TRUE)
    
    group = ifelse(substr(compare_group, 1,1) == substr(compare_group, 3,3), "same", "diff")
    ave_b = c(ave_b, mean(dist_test[group == "diff"]))
    ave_w = c(ave_w, mean(dist_test[group == "same"]))
  }
  
  return(list("p" = pvals,
              "t" = astat,
              "b" = ave_b,
              "w" = ave_w))
}

power_data = data.frame("n" = 10,"m" = 5000, "lfc" = 0.05, "sigma" = 0.13, "rou_DE" = 0.1,
                        "omega_k" = "(1,...,1)", "lambda" = "~8.25", "tau" = 0, "delta" = 0, "alpha" = 100,
                        "Pi_A" = "(0.4, 0.05, 0.2, 0.05, 0.05, 0.05, 0.1, 0.1)",
                        "Pi_B" = "(0.15, 0.05, 0.1, 0.4, 0.05, 0.1, 0.05, 0.1)",
                        "Pi_A1" = NA,
                        "Pi_B1" = NA,
                        "Power_GMM_PCA" = NA, "Power_KNN_PCA" = NA,
                        "Power_GMM_scvi" = NA, "Power_KNN_scvi" = NA, 
                        "Power_cluster_true" = NA,"Power_cluster_scvi" = NA,
                        "Power_cluster_PCA" = NA,
                        "ANOSIM_GMM_PCA" = NA,"ANOSIM_KNN_PCA" = NA,
                        "ANOSIM_GMM_scvi" = NA,"ANOSIM_KNN_scvi" = NA, 
                        "ANOSIM_cluster_true" = NA, "ANOSIM_cluster_scvi" = NA,
                        "ANOSIM_cluster_PCA" = NA,
                        "AveB_GMM_PCA" = NA, "AveB_KNN_PCA" = NA,
                        "AveB_GMM_scvi" = NA, "AveB_KNN_scvi" = NA, 
                        "AveB_cluster_true" = NA, "AveB_cluster_scvi" = NA, 
                        "AveB_cluster_PCA" = NA, 
                        "AveW_GMM_PCA" = NA, "AveW_KNN_PCA" = NA,
                        "AveW_GMM_scvi" = NA,"AveW_KNN_scvi" = NA,
                        "AveW_cluster_true" = NA, "AveW_cluster_scvi" = NA,
                        "AveW_cluster_PCA" = NA)

clusprop = list.files("../../results/simulation/distmat/clusprop_V3/")
clus_list = list()
for(i in clusprop){
  tmp_files = list.files(paste0("../../results/simulation/distmat/clusprop_V3/",i))
  if(i == "Setting1"){
    tmp_files = tmp_files[-which(str_detect(tmp_files, "alpha|novar"))]
  }
  power_tmp = power_data[rep(1,length(tmp_files)),]
  power_tmp$lfc = NA
  power_tmp$rou_DE = 0
  
  rownames(power_tmp) = tmp_files
  if(i == "Setting1"){
    p_i = c(10,5,10,15,10,5,10)
    power_tmp$Pi_B = paste(round(c(25,p_i)/sum(c(25,p_i)),2), collapse = ",")
    power_tmp$Pi_B1 = round(25/sum(c(25,p_i)),2)
    
  }else{
    p_i = c(100,90,100,110,100,95,100)
    power_tmp$Pi_B = paste(round(c(20,p_i)/sum(c(20,p_i)),2), collapse = ",")
    power_tmp$Pi_B1 = round(20/sum(c(20,p_i)),2)
  }
  
  for(j in tmp_files){
    pi_A1 = as.numeric(unlist(str_split(j, "_"))[2])
    power_tmp[j, "Pi_A"] = paste(round(c(pi_A1,p_i)/sum(c(pi_A1,p_i)),2), collapse = ",")
    power_tmp[j, "Pi_A1"] = round(pi_A1/sum(c(pi_A1,p_i)),2)
    
    power_GMM_PCA = try(myfun_clus(paste0("../../results/simulation/distmat/clusprop_V3/",i, "/",j, "/Glos/"), "dist_mat_GMM_pca"))
    Power_KNN_PCA = try(myfun_clus(paste0("../../results/simulation/distmat/clusprop_V3/",i, "/",j, "/Glos/"), "dist_mat_KNN_pca"))
    Power_GMM_scvi = try(myfun_clus(paste0("../../results/simulation/distmat/clusprop_V3/",i, "/",j, "/Glos/"), "dist_mat_GMM_scvi"))
    Power_KNN_scvi = try(myfun_clus(paste0("../../results/simulation/distmat/clusprop_V3/",i, "/",j, "/Glos/"), "dist_mat_KNN_scvi"))
    Power_cluster_true = try(myfun_clus(paste0("../../results/simulation/distmat/clusprop_V3/",i, "/",j, "/cluster/"), "clus_true"))
    Power_cluster_scvi = try(myfun_clus(paste0("../../results/simulation/distmat/clusprop_V3/",i, "/",j, "/cluster/"), "clus_seurat_scvi"))
    Power_cluster_pca = try(myfun_clus(paste0("../../results/simulation/distmat/clusprop_V3/",i, "/",j, "/cluster/"), "clus_seurat_pca"))
   
    
    
    power_tmp[j,"Power_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$p<0.05))
    power_tmp[j,"Power_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$p<0.05))
    power_tmp[j,"Power_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$p<0.05))
    power_tmp[j,"Power_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$p<0.05))
    power_tmp[j,"Power_cluster_true"] =  ifelse(inherits(Power_cluster_true, "try-error"), NA, mean(Power_cluster_true$p<0.05))
    power_tmp[j,"Power_cluster_scvi"] =  ifelse(inherits(Power_cluster_scvi, "try-error"), NA, mean(Power_cluster_scvi$p<0.05))
    power_tmp[j,"Power_cluster_PCA"] =  ifelse(inherits(Power_cluster_pca, "try-error"), NA, mean(Power_cluster_pca$p<0.05))
    
    
    power_tmp[j,"ANOSIM_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$t))
    power_tmp[j,"ANOSIM_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$t))
    power_tmp[j,"ANOSIM_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$t))
    power_tmp[j,"ANOSIM_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$t))
    power_tmp[j,"ANOSIM_cluster_true"] =  ifelse(inherits(Power_cluster_true, "try-error"), NA, mean(Power_cluster_true$t))
    power_tmp[j,"ANOSIM_cluster_scvi"] =  ifelse(inherits(Power_cluster_scvi, "try-error"), NA, mean(Power_cluster_scvi$t))
    power_tmp[j,"ANOSIM_cluster_PCA"] =  ifelse(inherits(Power_cluster_pca, "try-error"), NA, mean(Power_cluster_pca$t))
    
    power_tmp[j,"AveB_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$b))
    power_tmp[j,"AveB_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$b))
    power_tmp[j,"AveB_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$b))
    power_tmp[j,"AveB_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$b))
    power_tmp[j,"AveB_cluster_true"] =  ifelse(inherits(Power_cluster_true, "try-error"), NA, mean(Power_cluster_true$b))
    power_tmp[j,"AveB_cluster_scvi"] =  ifelse(inherits(Power_cluster_scvi, "try-error"), NA, mean(Power_cluster_scvi$b))
    power_tmp[j,"AveB_cluster_PCA"] =  ifelse(inherits(Power_cluster_pca, "try-error"), NA, mean(Power_cluster_pca$b))
    
    power_tmp[j,"AveW_GMM_PCA"] =  ifelse(inherits(power_GMM_PCA, "try-error"), NA, mean(power_GMM_PCA$w))
    power_tmp[j,"AveW_KNN_PCA"] =  ifelse(inherits(Power_KNN_PCA, "try-error"), NA, mean(Power_KNN_PCA$w))
    power_tmp[j,"AveW_GMM_scvi"] =  ifelse(inherits(Power_GMM_scvi, "try-error"), NA, mean(Power_GMM_scvi$w))
    power_tmp[j,"AveW_KNN_scvi"] =  ifelse(inherits(Power_KNN_scvi, "try-error"), NA, mean(Power_KNN_scvi$w))
    power_tmp[j,"AveW_cluster_true"] =  ifelse(inherits(Power_cluster_true, "try-error"), NA, mean(Power_cluster_true$w))
    power_tmp[j,"AveW_cluster_scvi"] =  ifelse(inherits(Power_cluster_scvi, "try-error"), NA, mean(Power_cluster_scvi$w))
    power_tmp[j,"AveW_cluster_PCA"] =  ifelse(inherits(Power_cluster_pca, "try-error"), NA, mean(Power_cluster_pca$w))
  
    }
  clus_list[[i]] = power_tmp

  }

cluster_power = do.call(rbind, clus_list)
write.csv(cluster_power, file = "../../results/simulation/power/cluster.csv")

