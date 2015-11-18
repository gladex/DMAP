
fun_Pheatmap<-function(mymatrix){
  library(pheatmap)
  p=pheatmap(mymatrix, border_color='white',
           cluster_rows=F, cluster_cols=F, clustering_method="complete",
           #kmeans_k=3,
           color=colorRampPalette(c("skyblue","navy","firebrick3"))(50),#firebrick3
           #color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100),
           treeheight_row=5, treeheight_col=5,
           legend=F,
           number_color='white',fontsize_number=11,fontsize=11,
           #cellwidth=40, #cellheight=40,
           display_numbers=T, number_format="%.3f")   # "%.1f"
  p
}