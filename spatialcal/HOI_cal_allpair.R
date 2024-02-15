# ---
# title: "Forest_spatial_neighbourhood_calculation"
# author: "Hanlun Liu" on the top of Bin Wang's codes for "Zofin_HOI" project
# date: "2024/02/14"
# output: csv_table
# ---


# "some functions"

#```{r}

create_directory_FUN <- function(){
  
  gc()
  
  dir_name <- paste("big_mat", 
                    format(Sys.time(), "%Y_%m%d_%H%M_%OS"), 
                    paste(sample(letters, 8), collapse = ""), 
                    sep = "_")
  
  backingpath <- sprintf(fmt = "/tmp/%s", dir_name)
  backingfile <- sprintf(fmt = "%s.bin", dir_name)
  descriptorfile <- sprintf(fmt = "%s.desc", dir_name)
  
  dir.create(path = backingpath)
  
  # return
  list(backingpath = backingpath, 
       backingfile = backingfile, 
       descriptorfile = descriptorfile)
}

#---

create_big_matrix_FUN <- function(nrow_0, colnames_0) {
  
  dir_description <- create_directory_FUN()
  
  # return
  bigmemory::filebacked.big.matrix(
    nrow = nrow_0, 
    ncol = length(colnames_0), 
    type = "double", 
    init = NA, 
    dimnames = list(NULL, colnames_0), 
    separated = TRUE, 
    backingfile = dir_description[["backingfile"]], 
    backingpath = dir_description[["backingpath"]], 
    descriptorfile = dir_description[["descriptorfile"]], 
    binarydescriptor = FALSE)
}

#---

deepcopy_big_matrix_FUN <- function(big_matrix, col_0, row_0 = NA) {
  
  
  if(all(is.na(row_0) == T)) {row_0 <- NULL} else {row_0 <- row_0}
  
  dir_description <- create_directory_FUN()
  
  big_matrix <- bigmemory::deepcopy(
    x = big_matrix, 
    cols = col_0, 
    rows = row_0, 
    separated = TRUE, 
    backingfile = dir_description[["backingfile"]], 
    backingpath = dir_description[["backingpath"]], 
    descriptorfile = dir_description[["descriptorfile"]], 
    binarydescriptor = FALSE, 
    shared = TRUE)
  
  options(bigmemory.allow.dimnames = TRUE)
  
  colnames(big_matrix) <- col_0
  
  # return
  big_matrix
  
}
#---

parallel_lapply_FUN <- function(
  x, fun, varlist = NA, text_to_parse = NA) {
  
  num_cpu <- 20
  
  cl <- parallel::makeCluster(num_cpu)
  
  if(all(is.na(varlist) == T)) {NULL} else {
    
    parallel::clusterExport(
      cl = cl, varlist = varlist, 
      envir = sys.frame(sys.nframe() - 1))
    
  }
  
  if(all(is.na(text_to_parse)) == T) {NULL} else {
    
    parallel::clusterExport(
      cl = cl, varlist = "text_to_parse", 
      envir = sys.frame(sys.nframe()))
    
    parallel::clusterEvalQ(
      cl = cl, 
      expr = eval(parse(text = text_to_parse)))
  }
  
  parallel_results <- parallel::parLapply(
    cl = cl, X = x, fun = compiler::cmpfun(fun))####put variable here
  
  names(parallel_results) <- x
  
  parallel::stopCluster(cl)
  
  # return
  parallel_results
}

#

big_tapply_FUN <- function(
  x = id_pair_matrix, sub_row = NA, 
  ccols = "focal_id", splitcol = "DI_mp", FUN = sum) {
  
  if(all(is.na(sub_row) == T)) {
    x <- x[, c(ccols, splitcol)]
  } else {
    x <- x[sub_row, c(ccols, splitcol)]
  }
  
  dat_split <- bigtabulate::bigsplit(
    x = x, ccols = ccols, splitcol = splitcol)
  
  gc()
  
  result <- sapply(X = dat_split, FUN = FUN)
  
  # return
  result
}


#```






# begin ~~~~~~~~~~~~~~~~~~~~




## df_FUN 
#```{r}
df_FUN <- function (df, r_max) {
  
  df <- as.data.frame(df[, c("tag", "Accepted_species_tnrs", "dbh1", "dbh2", "status1", "status2", "gx", "gy")])
  colnames(df) <- c("tree_tag", "species", "dbh", "dbh2", "STATUS1", "STATUS2", "gx", "gy")
  
  df[, "sp_name"] <- as.integer(as.factor(df[, "species"]))
  #print(df[2, "sp_name"])
  #print(as.factor(df[2, "species"]))
  colnames_add <- c("gx_square_num", "gy_square_num", 
                    "gxy_square_name", paste("square_name", 1:8, sep = "_"))
  
  df[, colnames_add] <- 0L
  
  df[, "id_name"] <- 1:nrow(df)
  
  rownames(df) <- NULL
  
  # return
  df
} 

#```




## df_add_square_FUN
#```{r}
df_add_square_FUN <- function(df, r_max){
  #
  df[, "gx_square_num"] <- as.integer(factor(ceiling(
    df[, "gx"] / r_max) * r_max)); gc()
  
  df[, "gy_square_num"] <- as.integer(factor(ceiling(
    df[, "gy"] / r_max) * r_max)); gc()
  
  df[, "gxy_square_name"] <- df[, "gx_square_num"] * 10000L + 
    df[, "gy_square_num"]; gc()
  
  df[, "square_name_1"] <- df[, "gx_square_num"] * 10000L + 
    (df[, "gy_square_num"] + 1L); gc()
  
  df[, "square_name_2"] <- df[, "gx_square_num"] * 10000L + 
    (df[, "gy_square_num"] - 1L); gc()
  
  df[, "square_name_3"] <- (df[, "gx_square_num"] + 1L) *
    10000L + df[, "gy_square_num"]; gc()
  
  df[, "square_name_4"] <- (df[, "gx_square_num"] + 1L) *
    10000L + (df[, "gy_square_num"] + 1L); gc()
  
  df[, "square_name_5"] <- (df[, "gx_square_num"] + 1L) *
    10000L + (df[, "gy_square_num"] - 1L); gc()
  
  df[, "square_name_6"] <- (df[, "gx_square_num"] - 1L) *
    10000L + df[, "gy_square_num"]; gc()
  
  df[, "square_name_7"] <- (df[, "gx_square_num"] - 1L) *
    10000L + (df[, "gy_square_num"] + 1L); gc()
  
  df[, "square_name_8"] <- (df[, "gx_square_num"] - 1L) *
    10000L + (df[, "gy_square_num"] - 1L); gc()
  
  #return
  df
}
#```



## annular_id_names_of_focal_ids_FUN
#```{r }
annular_id_names_of_focal_ids_FUN <- function (df, r_max) { 
  #
  annular_ids_in_gxy_squares <- split(
    x = df[, "id_name"], 
    f = df[, "gxy_square_name"])
  
  gxy_square_names <- as.character(names(annular_ids_in_gxy_squares))
  
  #
  gxy_squares_of_focal_ids <- split(
    x = df[, grep(pattern = "square_name", x = colnames(df))], 
    f = df[, "id_name"])
  
  index_focal_ids <- as.integer(as.character(names(gxy_squares_of_focal_ids)))
  
  #
  gx_0 <- df[, "gx"]
  gy_0 <- df[, "gy"]
  #
  
  fun <- function (index_focal_i) {
    
    index_gxy_squares <- gxy_square_names %in% gxy_squares_of_focal_ids[[index_focal_i]]
    
    annular_ids_of_focal_id <- unlist(x = annular_ids_in_gxy_squares[index_gxy_squares], 
                                      use.names = FALSE)
    
    gr <- ((gx_0[annular_ids_of_focal_id] - gx_0[index_focal_i]) ^ 2 + 
             (gy_0[annular_ids_of_focal_id] - gy_0[index_focal_i]) ^ 2) ^ 0.5
    # 
    as.integer(annular_ids_of_focal_id[gr <= r_max])
  }
  
  x <- index_focal_ids
  
  varlist <- c("annular_ids_in_gxy_squares", "gxy_square_names", "gxy_squares_of_focal_ids",  "gx_0", "gy_0", "r_max")
  
  annular_id_names_of_focal_ids <- 
    parallel_lapply_FUN(x = x, fun = fun, varlist = varlist)
  
  names(annular_id_names_of_focal_ids) <- index_focal_ids
  
  # return    
  annular_id_names_of_focal_ids
}

#```






## id_pair_matrix_FUN 
#```{r}
id_pair_matrix_FUN <- function (
  annular_id_names_of_focal_ids, df) 
{
  nrow_0 <- length(unlist(x = annular_id_names_of_focal_ids, recursive = TRUE))
  
  colnames_0 <- c("focal_id", "annular_id", "focal_sp", "annular_sp", "focal_dbh", "annular_dbh", "spatial_dist", "f_gx", "f_gy", "a_gx", "a_gy", "DI_mp", "DI_pm", "DI_pq", "DI_piq", "DI_pjq", "DI_pLq", "DI_pSq", "HOI_mpq", "HOI_mpiq", "HOI_mpjq", "HOI_mpLq", "HOI_mpSq")
  
  id_pair_matrix <- create_big_matrix_FUN(nrow_0, colnames_0)
  print("create_big_matrix done")
  # 
  #Nrow<-nrow(id_pair_matrix)
  #print(Nrow)
  #colnum<-which(colnames(id_pair_matrix)=="focal_id")
  #print(colnum)
  #I <- Nrow*(colnum-1)
  id_pair_matrix[, "focal_id"] <- as.integer(#id_pair_matrix[, "focal_id"]
    rep(x = as.integer(names(annular_id_names_of_focal_ids)), 
        times = sapply(X = annular_id_names_of_focal_ids, FUN = length))); gc()
  print("focal_id done")
  
  #colnum<-which(colnames(id_pair_matrix)=="annular_id")
  #I <- Nrow*(colnum-1)
  id_pair_matrix[, "annular_id"] <- as.integer(#id_pair_matrix[(I+1):(I+Nrow)]
    unlist(x = annular_id_names_of_focal_ids, use.names = FALSE)); gc()
  print("annual_id done")
  
  #colnum<-which(colnames(id_pair_matrix)=="focal_sp")
  #I <- Nrow*(colnum-1)
  id_pair_matrix[, "focal_sp"]<- # id_pair_matrix[(I+1):(I+Nrow)]
    df[id_pair_matrix[, "focal_id"], "sp_name"]; gc()
  print("focal_sp done")
  
  #colnum<-which(colnames(id_pair_matrix)=="annular_sp")
  #I <- Nrow*(colnum-1)
  id_pair_matrix[, "annular_sp"]<- #id_pair_matrix[(I+1):(I+Nrow)]
    df[id_pair_matrix[, "annular_id"], "sp_name"]; gc()
  print("annular_sp done")
  
  #colnum<-which(colnames(id_pair_matrix)=="focal_dbh")
  #I <- Nrow*(colnum-1)
  id_pair_matrix[, "focal_dbh"] <-#id_pair_matrix[(I+1):(I+Nrow)]
    df[id_pair_matrix[, "focal_id"], "dbh"]; gc()
  print("focal_dbh done")
  
  #colnum<-which(colnames(id_pair_matrix)=="annular_dbh")
  #I <- Nrow*(colnum-1)
  id_pair_matrix[,"annular_dbh"] <- #[(I+1):(I+Nrow)]
    df[id_pair_matrix[, "annular_id"], "dbh"]; gc()
  print("annular_dbh done")
  
  #colnum<-which(colnames(id_pair_matrix)=="f_gx")
  #I <- Nrow*(colnum-1)
  id_pair_matrix[, "f_gx"] <-#id_pair_matrix[(I+1):(I+Nrow)]
    df[id_pair_matrix[, "focal_id"], "gx"]; gc()
  print("f_gx done")
  
  #colnum<-which(colnames(id_pair_matrix)=="f_gy")
  #I <- Nrow*(colnum-1)
  id_pair_matrix[, "f_gy"]<- #id_pair_matrix[(I+1):(I+Nrow)]
    df[id_pair_matrix[, "focal_id"], "gy"]; gc()
  print("f_gy done")
  
  #colnum<-which(colnames(id_pair_matrix)=="a_gx")
  #I <- Nrow*(colnum-1)
  id_pair_matrix[, "a_gx"] <- #id_pair_matrix[(I+1):(I+Nrow)]
    df[id_pair_matrix[, "annular_id"], "gx"]; gc()
  print("a_gx done")
  
  #colnum<-which(colnames(id_pair_matrix)=="a_gy")
  #I <- Nrow*(colnum-1)
  id_pair_matrix[, "a_gy"] <- #id_pair_matrix[(I+1):(I+Nrow)]
    df[id_pair_matrix[, "annular_id"], "gy"]; gc()
  print("a_gy done")
  
  #colnum<-which(colnames(id_pair_matrix)=="spatial_dist")
  #I <- Nrow*(colnum-1)
  id_pair_matrix[, "spatial_dist"] <- #id_pair_matrix[(I+1):(I+Nrow)]
    ((id_pair_matrix[, "f_gx"] - id_pair_matrix[, "a_gx"]) ^ 2 +
       (id_pair_matrix[, "f_gy"] - id_pair_matrix[, "a_gy"]) ^ 2) ^ 0.5; gc()
  print("spatial_dist done")
  # return
  id_pair_matrix
}

#```




## single_big_matrix_FUN
#```{r}

single_big_matrix_FUN <- function(id_pair_matrix, R){
  #
  
  index_row <- bigmemory::mwhich(
    x = id_pair_matrix, cols = "spatial_dist", 
    vals = c(0, R), comps = c("gt", "le"), op = "AND")
  
  col_nam <- c("focal_id", "annular_id", "focal_sp", "annular_sp", "focal_dbh","annular_dbh", "spatial_dist", "DI_mp", "DI_pm", "DI_pq", "DI_piq", "DI_pjq", "DI_pLq", "DI_pSq", "HOI_mpq", "HOI_mpiq", "HOI_mpjq", "HOI_mpLq", "HOI_mpSq")
  
  single_big_matrix <- deepcopy_big_matrix_FUN(
    big_matrix = id_pair_matrix, 
    col_0 = col_nam, 
    row_0 = index_row)
  
  # return
  single_big_matrix
}




#```



## single_hoi_FUN
#```{r}

single_hoi_FUN <- function(single_big_matrix_0, df_0,splist) {
  
  #single_big_matrix_0 <- deepcopy_big_matrix_FUN(
    #big_matrix = single_big_matrix,
    #col_0 = colnames(single_big_matrix))
  
  #single_big_matrix_0[, "DI_mp"] <- (single_big_matrix_0[, "annular_dbh"] ^ U) / (single_big_matrix_0[, "spatial_dist"] ^ V)
  #single_big_matrix_0[, "DI_pm"] <- (single_big_matrix_0[, "focal_dbh"] ^ U) / (single_big_matrix_0[, "spatial_dist"] ^ V)
  
  
  #df_0 <- as.data.frame(df)
  
  #df_0 <- df_0[, c("tree_tag", "species", "dbh", "dbh2", "STATUS1", "STATUS2","gx","gy")]
  
  #names(df_0) <- c("TREE_TAG","SPECIES","DBH1","DBH2","STATUS1","STATUS2","X","Y")
  
  #df_0[, c("R", "u", "v", "DI", "DI_ii", "DI_ij", "DI_im")] <- 0L
  
  #df_0[,"R"] <- R; df_0[,"u"] <- U; df_0[,"v"] <- V 
  
  #----------------------
  print("part 1/4 ~~~")
  print(Sys.time())
  #----------------------

  DI <- big_tapply_FUN(
    x = single_big_matrix_0, sub_row = NA, 
    ccols = "focal_id", splitcol = "DI_mp", FUN = sum)
  df_0[as.numeric(names(DI)), "DI"] <- DI
  
  #
  index_row_sp <- single_big_matrix_0[, "annular_sp"] == single_big_matrix_0[, "focal_sp"]
  
#  DI_ii <- big_tapply_FUN(
#    x = single_big_matrix_0, sub_row = index_row_sp, 
#    ccols = "focal_id", splitcol = "DI_mp", FUN = sum)
#  df_0[as.numeric(names(DI_ii)), "DI_ii"] <- DI_ii
  
  

  funsp<- function(i){
    DI_new<-rep(0,nrow(df_0))
    names(DI_new)<- rownames(df_0)
    DI_ij <- big_tapply_FUN(
      x = single_big_matrix_0, sub_row = single_big_matrix_0[, "annular_sp"]==i, 
      ccols = "focal_id", splitcol = "DI_mp", FUN = sum)
    DI_new[as.numeric(names(DI_ij))] <- DI_ij
    DI_new<- data.frame(DI_new)
    colnames(DI_new)<- splist[i,"species"]
    return(DI_new)
  }
  
  #debug(funsp)
  #registerDoParallel(cores=10)
  print("loop starts")
  #DI_ijtab <- foreach(i=unique(as.integer(as.factor(df_0$SPECIES))),  .inorder = TRUE,.errorhandling="pass") %dopar%{funsp(i)}#.combine=cbind,
  #df_0<-cbind(df_0,DI_ijtab)
  DI_ijtab <-list()
  j<-1
  for (i in unique(as.integer(as.factor(df_0$SPECIES)))){
    print(j)
    DI_ijtab[[j]]<-funsp(i)
    j<-j+1
  }
  df_0<-list(df_0,DI_ijtab)
  
  print("DI part over")
  print(Sys.time())
  #hoi deleted
  ###
  path <- bigmemory::describe(single_big_matrix_0) @ description[["dirname"]]
  path <- gsub(pattern = "/$", replacement = "", x = path)
  rm(list = "single_big_matrix_0"); gc()
  unlink(x = path, recursive = T, force = T)
  print("delete part over")
  print(Sys.time())
  # return
  df_0
}

#```





#  parallel_all_hoi_FUN
#```{r}
parallel_all_hoi_FUN <- function(df, R_0, U_0, V_0, spdone=NULL)
{
  cat("\n")
  print("begin ~~~~")
  print(Sys.time())
  cat("\n")
  
  r_max <- max(R_0)
  
  df <- df_FUN(df, r_max)
  splist<-unique(df[,c("sp_name","species")])
  splist<-splist[order(splist[,"sp_name"]),]
  print("df_FUN done")
  df <- df_add_square_FUN(df, r_max)
  print("df_add_square done")
  annular_id_names_of_focal_ids <-
    annular_id_names_of_focal_ids_FUN(df, r_max)
  print("ids_FUN done")
  id_pair_matrix <-
    id_pair_matrix_FUN(annular_id_names_of_focal_ids, df)
  
  print("id_pair_matrix is finished.")
  print(Sys.time())
  cat("\n")
  
  dir()
  
  path_out <- file.path(getwd(), "out")
  if(dir.exists(paths = path_out) == T ) {NULL} else {
    dir.create(path = path_out)
  }
  
  
  
  for(R in R_0) {
    
    print(sprintf(fmt = "run single_hoi_FUN: R, %s", R))
    print(sprintf(fmt = "U, %s", paste(U_0, collapse = ",")))
    print(sprintf(fmt = "V, %s.", paste(V_0, collapse = ",")))
    print(Sys.time())
    cat("\n")
    
    single_big_matrix <- single_big_matrix_FUN(id_pair_matrix, R)
    
    describe_inf <- bigmemory::describe(single_big_matrix)
    #all loops use the same single big matrix
   single_big_matrix <- bigmemory::attach.big.matrix(describe_inf)    
      single_big_matrix_0 <- deepcopy_big_matrix_FUN(
      big_matrix = single_big_matrix,
      col_0 = colnames(single_big_matrix))
      
      single_big_matrix_0[, "DI_mp"] <- (single_big_matrix_0[, "annular_dbh"] ^ U_0) / (single_big_matrix_0[, "spatial_dist"] ^ V_0)
      single_big_matrix_0[, "DI_pm"] <- (single_big_matrix_0[, "focal_dbh"] ^ U_0) / (single_big_matrix_0[, "spatial_dist"] ^ V_0)
      
      df_00 <- as.data.frame(df)
      
      df_00 <- df_00[, c("tree_tag", "species", "dbh", "dbh2", "STATUS1", "STATUS2","gx","gy")]
      
      names(df_00) <- c("TREE_TAG","SPECIES","DBH1","DBH2","STATUS1","STATUS2","X","Y")
      
      df_00[, c("R", "u", "v", "DI", "DI_ii")] <- 0L
      
      df_00[,"R"] <- R; df_00[,"u"] <- U_0; df_00[,"v"] <- V_0 
      
    fun <- function(single0=single_big_matrix_0, df000=df_00) {
      
      #U = case[i, "U"]
      #V = case[i, "V"]
      #spj = case[i,"spj"]
      
      print("hoi start")
      hoi <- single_hoi_FUN(single0, df000,splist) 
      print("hoi done")
      
      write.csv(x = cbind(hoi[[1]],dplyr::bind_cols(hoi[[2]])),
                file = sprintf(fmt = "out/HOI_R%s_u%s_v%s.csv", R, U_0, V_0), row.names = F)
      #save(hoi,
       #             file = sprintf(fmt = "out/HOI_R%s_u%s_v%s.RData", R, U_0, V_0))
    }
    
    ###---
    fun()
  }
  print("finish ~~~~")
  print(Sys.time())
  cat("\n")
  
  #---------
  rm(list = c("id_pair_matrix", "single_big_matrix")); gc()
  big_mat_paths <- grep(pattern = "big_mat_[0-9]{4}_.*", x = dir(getwd()), value = T)
  lapply(X = big_mat_paths, FUN = unlink, recursive = TRUE, force = TRUE)
  #---------
  
  # return
  NULL
}
#```

#multiplot
di.cal <- function(n=NA,plotname,plotlist){#10 & 20 & 50 
  if(!is.na(n)){
  plotname<-names(plotlist)[n]}
  print(plotname)
  setwd(here(paste0(plotname)))
  plotdata<-plotlist[[plotname]]
  #10x10
  print("############## 10x10 begin")
  #### Cal (DBH and distance) 
  tryCatch(parallel_all_hoi_FUN(plotdata,R_0=10,U_0=1,V_0=1),error = function(e)e)#R--Radius; U--exponent on DBH; V--exponent on distance
  #20x20
  print("############## 20x20 begin")
  #### Cal (DBH and distance) 
  tryCatch(parallel_all_hoi_FUN(plotdata,R_0=20,U_0=1,V_0=1),error = function(e)e)
  
  #50x50
  print("############## 50x50 begin")
  tryCatch(parallel_all_hoi_FUN(plotdata,R_0=50,U_0=1,V_0=1),error = function(e)e)
}

di.cal2050 <- function(n=NA,plotname,plotlist){#10 & 20 & 50 
  if(!is.na(n)){
    plotname<-names(plotlist)[n]}
  print(plotname)
  #setwd(paste0("/data/gpfs/projects/punim1844/treeweb/multiplotcal/",plotname,"/split/",split))
  setwd(paste0("/data/gpfs/projects/punim1844/treeweb/multiplotcal/",plotname))
  plotdata<-plotlist[[plotname]]
  #20x20
  print("############## 20x20 begin")
  #### Cal (DBH and distance) 
  tryCatch(parallel_all_hoi_FUN(plotdata,R_0=20,U_0=1,V_0=1),error = function(e)e)
  
  #50x50
  print("############## 50x50 begin")
  tryCatch(parallel_all_hoi_FUN(plotdata,R_0=50,U_0=1,V_0=1),error = function(e)e)
}

di.cal50 <- function(n=NA,plotname,plotlist){#10 & 20 & 50 
  if(!is.na(n)){
    plotname<-names(plotlist)[n]}
  print(plotname)
  #setwd(paste0("/data/gpfs/projects/punim1844/treeweb/multiplotcal/",plotname,"/split/",split))
  setwd(paste0("/data/gpfs/projects/punim1844/treeweb/multiplotcal/",plotname))
  plotdata<-plotlist[[plotname]]
  
  #50x50
  print("############## 50x50 begin")
  tryCatch(parallel_all_hoi_FUN(plotdata,R_0=50,U_0=1,V_0=1),error = function(e)e)
}
##test
#debug(parallel_all_hoi_FUN)
#parallel_all_hoi_FUN(df=hsda1c500, R_0=10, U_0=1, V_0=1, spdone=NULL,spnum=1)
#di.cal(cbs)
#di.cal(hsda1)
