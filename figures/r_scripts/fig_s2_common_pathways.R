library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(devtools)
library(pheatmap)

data_base_path = '../figure_data/'
pathways_dt <- fread(paste0(data_base_path,'fig2/common_pathways_mouse_human_within_tissues.csv'))


rownames(pathways_dt) <- pathways_dt[, Term]
pathways_dt[, Term:=NULL]

colorscale <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(100)

pheatmap(pathways_dt, na_col = "grey", cluster_rows=FALSE, cluster_cols=FALSE, fontsize=15,  color=colorscale)

# considering major pathways, discarding overlapping subsets

pathways_dt <- fread(paste0(data_base_path,'fig2/common_pathways_mouse_human_within_tissues.csv'))
major_pathways <- c('mitochondrial ATP synthesis coupled electron transport',
                    'respiratory electron transport chain',
                    'NADH dehydrogenase complex assembly', 
                    'mitochondrial translational elongation', 
                    'mitochondrial translational termination',
                    'translational termination')

pathways_dt <- pathways_dt[Term %in% major_pathways]
row_indices <- match(major_pathways, pathways_dt[, Term])
pathways_dt <- pathways_dt[row_indices,]
rownames(pathways_dt) <- pathways_dt[, Term]

pathways_dt[, Term:=NULL]



cairo_pdf('../plots/fig2/heatmap_common_pathways.pdf',width = 8, height = 3)
pheatmap(pathways_dt, legend=TRUE,na_col = "grey", 
         cluster_rows=FALSE, cluster_cols=FALSE, 
         fontsize=10,  color=colorscale, cexCol = 0.25)
dev.off() 
