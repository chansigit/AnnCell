library(SingleR)
BlacklistGenes<- function(seu){
    results <- c()
    # mito
    g <- grep("^mt-", rownames(seu), value = T); results <- c(results, g);
    # ribosomal
    g <- grep("^Rpl", rownames(seu), value = T); results <- c(results, g);
    g <- grep("^Rps", rownames(seu), value = T); results <- c(results, g);
    # mito rpl rps
    g <- grep("^Mrpl", rownames(seu), value = T); results <- c(results, g);
    g <- grep("^Mrps", rownames(seu), value = T); results <- c(results, g);
    
    # 	RNA Polymerase
    g <- grep("^Polr", rownames(seu), value = T); results <- c(results, g);
    
    # Translation factors
    g <- grep("^Eif", rownames(seu), value = T); results <- c(results, g);
    # Interferon stimulated genes
    g <- grep("^Ifi", rownames(seu), value = T); results <- c(results, g);
    g <- grep("^Isg", rownames(seu), value = T); results <- c(results, g);
    # Heat shock protein
    g <- grep("Hsp\\d+.*", rownames(seu), value = T); results <- c(results, g);
    g <- grep("^Hsbp", rownames(seu), value = T); results <- c(results, g);
    # Structural Maintenance Of Chromosomes
    g <- grep('Smc\\d+.*',rownames(seurat), value = T); results <- c(results, g);
    # Histone
    g <- grep('^Hist',rownames(seurat), value = T); results <- c(results, g);
    # Minichromosome Maintenance Complex Component
    g <- grep("^Mcm", rownames(seu), value = T); results <- c(results, g);
    # Cyclin
    g <- grep('^Ccn',rownames(seurat), value = T); results <- c(results, g);
    # Tryptophanyl-TRNA Synthetase
    g <- grep('^Wars',rownames(seurat), value = T); results <- c(results, g);
    # Atp
    g <- grep('^Atp',rownames(seurat), value = T); results <- c(results, g);
    # NADH dehydrogenase
    g <- grep('^Nduf',rownames(seurat), value = T); results <- c(results, g);
    
    # cc genes
    results<-c(results, "Anln","Anp32e","Aurka","Aurkb","Birc5","Blm","Bub1","Casp8ap2","Cbx5",
            "Ccnb2","Ccne2","Cdc20","Cdc25c","Cdc45","Cdc6","Cdca2","Cdca3","Cdca7",
            "Cdca8","Cdk1","Cenpa","Cenpb","Cenpe","Cenpf","Cenpu","Chaf1b","Ckap2","Ckap2l",
            "Ckap5","Cks1b","Cks2","Clspn","Ctcf","Dlgap5","Dscc1","Dtl",
            "E2f8","Ect2","Exo1","Fen1","G2e3","Gas2l3","Gins2","Gmnn","Gtse1","Hells",
            "Hjurp","Hmgb2","Hmmr","Jpt1","Kif11","Kif20b","Kif23","Kif2c","Lbr",
            "Mcm4","Mcm5","Mcm6","Mcm7","Mki67","Mrpl36","Msh2",
            "Nasp","Ncapd2","Ndc80","Nek2","Nuf2","Nusap1","Pcna","Pimreg","Pola1",
            "Pold3","Polr1b","Prim1","Psrc1","Rad51","Rad51ap1","Rangap1","Rfc2","Rrm1","Rrm2",
            "Slbp","Smc4","Tacc3","Tipin","Top2a","Tpx2","Ttk","Tubb4b","Tyms",
            "Ubr7","Uhrf1","Ung","Usp1","Wdr76")
    results <- c(results,"Jun","Junb","Jund")
    
    
    results <- unique(results)
    results <- intersect(rownames(seurat), results)
    results
}

genes.blacklist <- BlacklistGenes(seurat)




ref <- ImmGenData()

seu = seurat %>% filter(batch=='ABMaLP2')

ann_key = 'ImmGen_fine'
genes.restrict <- intersect( setdiff( rownames(ref), genes.blacklist ), rownames(seu))
mat            <- GetAssayData(seu, slot="data")[genes.restrict, ]
pred.fine      <- SingleR(test = mat, ref = ref, labels = ref$label.fine)
seu[[ann_key]] <- pred.fine$pruned.labels


library(dplyr)
newmeta <- seu@meta.data; rownames(newmeta) -> newmeta$bc
newmeta = merge(x=newmeta, y=read.csv('./immgen_label_v221009.tsv', sep ='\t', 
                                     col.names = c('immgen','ImmGen_fine_short')),
                by.x=ann_key, by.y='immgen', all.x=T) %>% `rownames<-`(.$bc)  
newmeta = newmeta[rownames(seu@meta.data), ]
newmeta -> seu@meta.data

seu-> seu.abm2
