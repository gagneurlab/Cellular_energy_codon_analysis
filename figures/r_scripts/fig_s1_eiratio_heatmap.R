library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

chosen.colors <- c("#2f4f4f", "#7f0000", "#006400", "#808000", "#483d8b",
                    "#3cb371", "#cd853f", "#4682b4", "#00008b", "#8b008b",
                    "#b03060", "#ff4500", "#ffa500", "#00ff00", "#8a2be2",
                    "#00ff7f", "#00ffff", "#0000ff", "#adff2f", "#ff00ff",
                    "#1e90ff", "#ffff54", "#add8e6", "#ff1493", "#ee82ee",
                    "#ffdead", "#ffb6c1")

data_base_path = '../figure_data/'

## columns as samples

ei_ratios_df <- fread(paste0(data_base_path,'fig1/samples_ei_ratios_major_iso_centered_2_3_non_na_comd_coef.csv'))
names(chosen.colors) <- ei_ratios_df[, sort(unique(tissue))]


col_ha = HeatmapAnnotation(Tissue=ei_ratios_df[, tissue], col=list(tissue=chosen.colors),annotation_name_side = "left")
ei_ratios_df[, tissue:=NULL]
ei_ratios_df[, comd_coef:=NULL]
ei_ratios_df[, sample.id:=NULL]
hmat = t(as.matrix(ei_ratios_df))

filename <- "../plots/fig1/ei_ratios_heatmap"
cairo_pdf(paste0(filename,'.pdf'))

#rows_to_samples <- sample(dim(hmat)[2], size=2000)
colors <- colorRamp2(c(-4,0,4), c("#00008B", "#FFFFFF", "#800000"))
Heatmap(hmat,  top_annotation = col_ha, name='log2 mRNA \nrelative half-life', col = colors,
        column_title='Samples', row_title='Transcripts', show_row_names=FALSE, 
        #row_names_gp = gpar(fontsize = 6), 
        cluster_rows=TRUE, cluster_columns=TRUE,
        width = unit(10, "cm"), height = unit(10, "cm"))
dev.off()







