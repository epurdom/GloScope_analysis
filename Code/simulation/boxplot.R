library(stringr)
between_fun = function(ddir, type){
  ave = c()
  for(i in list.files(ddir)){
    load(paste0(ddir,i))
    if(type == "gmmpca"){
      dist_mat = as.dist(dist_mat_GMM_pca)
    }else if(type == "knnpca"){
      dist_mat = as.dist(dist_mat_KNN_pca)
    }else if(type == "gmmscvi"){
      dist_mat = as.dist(dist_mat_GMM_scvi)
    }else if(type == "knnscvi"){    
      dist_mat = as.dist(dist_mat_KNN_scvi)
      }
    dist_test = c(dist_mat)
    mycompare = c()
    for(i in 1:length(attr(dist_mat, "Labels"))){
      mycompare = c(mycompare, paste0(attr(dist_mat, "Labels")[i], "-", attr(dist_mat, "Labels")[-(1:i)]))
    }
    
    mycompare = mycompare[-length(mycompare)]
    
    compare_group = str_replace_all(mycompare, paste0("sample", "[0-9]+", "[.](?!\\d+$)"), "")
    compare_group = gsub("[.](?!\\d+$)", "", compare_group, perl=TRUE)
    
    group = ifelse(substr(compare_group, 1,1) == substr(compare_group, 3,3), "same", "diff")
    ave = c(ave, mean(dist_test[group == "diff"]))
  }
  return(ave)
  
}

within_fun = function(ddir, type){
  ave_within = c()
  for(i in list.files(ddir)){
    load(paste0(ddir,i))
    if(type == "gmmpca"){
      dist_mat = as.dist(dist_mat_GMM_pca)
    }else if(type == "knnpca"){
      dist_mat = as.dist(dist_mat_KNN_pca)
    }else if(type == "gmmscvi"){
      dist_mat = as.dist(dist_mat_GMM_scvi)
    }else if(type == "knnscvi"){    
      dist_mat = as.dist(dist_mat_KNN_scvi)
    }
    
    dist_test = c(dist_mat)
    mycompare = c()
    for(i in 1:length(attr(dist_mat, "Labels"))){
      mycompare = c(mycompare, paste0(attr(dist_mat, "Labels")[i], "-", attr(dist_mat, "Labels")[-(1:i)]))
    }
    
    mycompare = mycompare[-length(mycompare)]
    
    compare_group = str_replace_all(mycompare, paste0("sample", "[0-9]+", "[.](?!\\d+$)"), "")
    compare_group = gsub("[.](?!\\d+$)", "", compare_group, perl=TRUE)
    
    group = ifelse(compare_group == "A-A", "A", "other")
    ave_within = c(ave_within, mean(dist_test[group == "A"]))
  }
  return(ave_within)
  
}
type = c("gmmpca", "knnpca", "gmmscvi", "knnscvi")



######################## sigma within
sd = c("sd013", "sd020", "sd030")
sd_list = list()

for(i in sd){
  for(j in type){
    sd_list[[paste0(i,"_",j)]] = 
      within_fun(ddir = paste0("../../results/simulation/distmat/baseline/",
                                i,"_lfc005/"),
                  type = j)
    
  }
}



########################## alpha within

alpha = c("alpha0", "alpha1", "alpha2","alpha10","alpha100")
alpha_list = list()

for(i in alpha){
  for(j in type){
    if(i == "alpha100"){
      alpha_list[[paste0(i,"_",j)]] = 
        within_fun(ddir = paste0("../../results/simulation/distmat/clusprop_V3/Setting1/25_100/Glos/"),
                   type = j)
    }else if(i == "alpha0"){
      alpha_list[[paste0(i,"_",j)]] = 
        within_fun(ddir = paste0("../../results/simulation/distmat/clusprop_V3/Setting1/25_100_novar/Glos/"),
                   type = j)
    }else{
      alpha_list[[paste0(i,"_",j)]] = 
        within_fun(ddir = paste0("../../results/simulation/distmat/clusprop_V3/Setting1/25_100_",i,"/Glos/"),
                   type = j)
    }
  }
}


################################### cluster between
setting = c("Setting1", "Setting3")
Set1_clus = c("25vs35", "25vs50", "25vs75", "25vs100")
Set3_clus = c("20vs40", "20vs60", "20vs80")
clus_list = list()

for(i in setting){
  if(i == "Setting1"){
    for(j in Set1_clus){
      for(k in type){
        clus_name = str_replace(j, "vs", "_")
        clus_list[[paste0("set1_",j,"_",k)]] = 
          between_fun(ddir = paste0("../../results/simulation/distmat/clusprop_V3/Setting1/",clus_name,"/Glos/"),
                      type = k)
      }
    }
  }else if(i == "Setting3"){
    for(j in Set3_clus){
      for(k in type){
        clus_name = str_replace(j, "vs", "_")
        clus_list[[paste0("set3_",j,"_",k)]] = 
          between_fun(ddir = paste0("../../results/simulation/distmat/clusprop_V3/Setting3/",clus_name,"/Glos/"),
                      type = k)
      }
    }    
  }
}

############################################# lfcpde

lfc = c("lfc005", "lfc010", "lfc015", "lfc020")
pde = c("pde01", "pde02", "pde03")
lfc_pde = list()
for(i in pde){
  for(j in lfc){
    for(k in type){
      if(i == "pde01"){
        if(j == "lfc005"){
          lfc_pde[[paste0(i,"_",j,"_",k)]] = 
            between_fun(ddir = "../../results/simulation/distmat/baseline/sd013_lfc005_same/",
                        type = k)
        }else{
          lfc_pde[[paste0(i,"_",j,"_",k)]] = 
            between_fun(ddir = paste0("../../results/simulation/distmat/parameter/more_lfc/",j,"/"),
                        type = k)
        }
      }else{
        if(j == "lfc020") next
        lfc_pde[[paste0(i,"_",j,"_",k)]] = 
          between_fun(ddir = paste0("../../results/simulation/distmat/parameter/more_pde/",i,"_",j,"/"),
                      type = k)
      }
    }

  }
}

