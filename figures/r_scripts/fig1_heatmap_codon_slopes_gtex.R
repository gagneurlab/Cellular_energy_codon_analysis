library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(stringr)

data_path = '../figure_data/fig1/codon_subtissue_effects_rplot.csv'

data_heatmap_df <- fread(data_path) #estimated codon effects were saturated to 10 or -10
data_heatmap_df[,codon:=str_replace_all(codon, 'T','U')]
sat_value_decoding_rate<-9 #(decoding rate has some outliers higher than 9codon/s, to better visualize its distribution we will saturate them to 9codon/s )
setnames(data_heatmap_df, 'decoding rate (HEK293)', 'codon_decoding_rate')
col_fun_ann = colorRamp2(c(min(data_heatmap_df$codon_decoding_rate), sat_value_decoding_rate), c("white", "seagreen"))
col_ha = HeatmapAnnotation(`Decoding rate\n[codon/s]`=data_heatmap_df$codon_decoding_rate, 
                       col=list(`Decoding rate\n[codon/s]`=col_fun_ann), annotation_name_gp= gpar(fontsize = 0),
                       annotation_legend_param = list(col_fun = col_fun_ann, at = c(5, 6, 7, 8, 9), labels = c("5", "6", "7", "8", "> 9"), title_gp = gpar(fontsize = 8)))
codons_list <- data_heatmap_df[, codon]
data_heatmap_df[, codon:=NULL]
data_heatmap_df[, codon_decoding_rate:=NULL]
hmat = t(as.matrix(data_heatmap_df))
colnames(hmat) = codons_list



pheatmap_default_colorscale <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
filename <- "../plots/fig1/heatmap_codon_slope_subtissues"
#svg(paste0(filename,'.svg'))
cairo_pdf(paste0(filename,'.pdf'))

ht <- Heatmap(hmat,  bottom_annotation = col_ha, name='Estimated \ncodon effect', row_names_gp = gpar(fontsize = 6),row_names_side = "right",
        heatmap_legend_param = list(col_fun = pheatmap_default_colorscale, at = c(-10, -5, 0, 5, 10), labels = c("< -10", "-5", "0", "5", "> 10"),title_gp = gpar(fontsize = 8)),
        rect_gp = gpar(col = "gray", lwd = 1),col=pheatmap_default_colorscale,
        cluster_rows=FALSE, cluster_columns=FALSE, column_names_rot = 60, column_names_gp = gpar(fontsize = 5),
        width = unit(10, "cm"), height = unit(11, "cm"))

draw(ht, heatmap_legend_side = "left", merge_legend = TRUE,annotation_legend_side="bottom", legend_grouping="original")

dev.off() 









