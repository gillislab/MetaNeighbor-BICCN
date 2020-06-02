# Utility functions used in MetaNeighbor vignette for BICCN data


mn_heatmap = function(aurocs) {
    subclass_to_col = c(
        Lamp5='#DA808C', Sncg='#D633FF', Vip='#B864CC', Sst='#FF9900',
        Pvalb='#D93137', "L2/3 IT"='#C4EC04', "L5 IT"='#50B2AD',
        "L6 IT"='#A19922', "L5 ET"='#0D5B78', "L6 IT Car3"='#5100FF',
        "L5/6 NP"='#3E9E64', "L6 CT"='#2D8CB8', L6b='#53377D', CR='#00FF66',
        Astro='#665C47', OPC='#2E3E39', Oligo='#2E3E39', Endo='#8D6C62',
        SMC='#807059', Peri='#665547', VLMC='#697255', Macrophage='#94AF97'
    )
    dataset_to_col = c(scSS='#EC2028', scCv3='#F15A25',
                       scCv2='#F68921',
                       snSS='#7E9C3D', snCv3M='#B5D554',
                       snCv3Z='#D3E3A6', snCv2='#F2F5C3')
    cell_type = MetaNeighbor::getCellType(colnames(aurocs))
    subclass_label = cell_type
    for (s in names(subclass_to_col)) {subclass_label[startsWith(cell_type, s)] = s}
    row_order = order(factor(subclass_label, levels = names(subclass_to_col)),
                      na.last = TRUE)
    aurocs = aurocs[row_order, row_order]
    subclass_label = subclass_label[row_order]

    gplots::heatmap.2(
        x = aurocs, margins = c(11,11), key = TRUE, keysize = 1, key.xlab="AUROC",
        key.title="NULL", offsetRow=0.1, offsetCol=0.1, trace = "none",
        Rowv = FALSE, Colv = FALSE, labRow = NA, labCol = NA, na.color = gray(0.95), 
        dendrogram = "none",  density.info = "none",
        col = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(20)),
        breaks = seq(0, 1, length = 21),
        RowSideColors = subclass_to_col[subclass_label],
        ColSideColors = dataset_to_col[MetaNeighbor::getStudyId(colnames(aurocs))]
    )

    par(lend = 1)
    legend("topleft", inset = c(0.82, 0.28), legend = names(subclass_to_col),
           col = subclass_to_col, cex = 0.7, lwd = 10, bty="n")
    legend("topleft", inset = c(0.82, 0.03), legend = names(dataset_to_col),
           col = dataset_to_col, pt.cex = 1, cex = 0.7, lwd = 10, bty="n")
}

mn_heatmap_pretrained = function(aurocs) {
    subclass_to_col = c(
        Lamp5='#DA808C', Sncg='#D633FF', Vip='#B864CC', Sst='#FF9900',
        Pvalb='#D93137', "L2/3 IT"='#C4EC04', "L5 IT"='#50B2AD',
        "L6 IT"='#A19922', "L5 ET"='#0D5B78', "L6 IT Car3"='#5100FF',
        "L5/6 NP"='#3E9E64', "L6 CT"='#2D8CB8', L6b='#53377D', CR='#00FF66',
        Astro='#665C47', OPC='#2E3E39', Oligo='#2E3E39', Endo='#8D6C62',
        SMC='#807059', Peri='#665547', VLMC='#697255', Macrophage='#94AF97'
    )
    cell_type = MetaNeighbor::getCellType(colnames(aurocs))
    subclass_label = cell_type
    for (s in names(subclass_to_col)) {subclass_label[startsWith(cell_type, s)] = s}
    col_order = order(factor(subclass_label, levels = names(subclass_to_col)),
                      na.last = TRUE)
    aurocs = aurocs[, col_order]
    subclass_label = subclass_label[col_order]
    positive_aurocs = aurocs
    positive_aurocs[aurocs < 0.5] = NA
    row_order = order(colSums(t(positive_aurocs)*seq_len(ncol(positive_aurocs)), na.rm=TRUE)/rowSums(positive_aurocs, na.rm=TRUE))
    aurocs = aurocs[row_order,]

    gplots::heatmap.2(
        x = aurocs, margins = c(12,12), key = TRUE, keysize = 1, key.xlab="AUROC",
        key.title="NULL", offsetRow=0.1, offsetCol=0.1, trace = "none",
        density.info="none", labCol=NA, na.color=gray(0.95),
        col = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(20)),
        breaks = seq(0, 1, length = 21),
        ColSideColors = subclass_to_col[subclass_label],
        dendrogram = "none", Rowv = FALSE, Colv = FALSE
    )

    par(lend = 1)
    legend("bottom", inset = c(0, 0), legend = names(subclass_to_col),
           col = subclass_to_col, cex = 0.5, lwd = 10, bty="n", ncol = 6)
}
