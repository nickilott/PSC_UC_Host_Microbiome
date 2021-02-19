


getOrphanModules <- function(features, method=c("absolute", "overlap.p")){
  
      all.modules <- features$module.module
      all.modules <- unlist(strsplit(as.character(all.modules), "_"))
      tissue1.modules <- as.character(all.modules[seq(1, length(all.modules), 2)])
      tissue2.modules <- as.character(all.modules[seq(2, length(all.modules), 2)])
      features$tissue1.modules <- tissue1.modules
      features$tissue2.modules <- tissue2.modules
      if (method == "absolute"){
          has.overlap.tissue1 <- features$tissue1.modules[features$overlap !=0]
          not.overlapping.tissue1 <- features$tissue1.modules[features$overlap == 0]
          orphan.tissue1 <- setdiff(not.overlapping.tissue1, has.overlap.tissue1)
          
          has.overlap.tissue2 <- features$tissue2.modules[features$overlap !=0]
          not.overlapping.tissue2 <- features$tissue2.modules[features$overlap == 0]
          orphan.tissue2 <- setdiff(not.overlapping.tissue2, has.overlap.tissue2)
      } else if (method == "overlap.p"){
        
        has.overlap.tissue1 <- features$tissue1.modules[features$overlap.p < 0.05]
        not.overlapping.tissue1 <- features$tissue1.modules[features$overlap.p > 0.05]
        orphan.tissue1 <- setdiff(not.overlapping.tissue1, has.overlap.tissue1)
        
        has.overlap.tissue2 <- features$tissue2.modules[features$overlap.p < 0.05]
        not.overlapping.tissue2 <- features$tissue2.modules[features$overlap.p > 0.05]
        orphan.tissue2 <- setdiff(not.overlapping.tissue2, has.overlap.tissue2)
      }
      orphans <- append(orphan.tissue1, orphan.tissue2)
      if (length(orphans) == 0){
        cat("No orphan modules found")}
      else{
        return(orphans)
      }
}
