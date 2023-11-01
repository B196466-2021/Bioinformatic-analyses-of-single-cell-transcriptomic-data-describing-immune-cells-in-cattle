#! /usr/bin/Rscript
# Go enrichment analysis of dendritic cells group 1
DC_g1.response <- data.frame(gene = rownames(DC_g1.response), DC_g1.response)
head(DC_g1.response)

logFC_t=0.5
P.Value_t = 0.05
k1 = (DC_g1.response$p_val_adj < P.Value_t)&(DC_g1.response$avg_log2FC < -logFC_t)
k2 = (DC_g1.response$p_val_adj < P.Value_t)&(DC_g1.response$avg_log2FC > logFC_t)
table(k1)
table(k2)

change = ifelse(k1,"down",ifelse(k2,"up","stable"))
DC_g1.response$change <- change
head(DC_g1.response)
s2e <- bitr(DC_g1.response$gene, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

DC_g1.response <- inner_join(DC_g1.response,s2e,by=c("gene"="SYMBOL"))
head(DC_g1.response)

#GO[Symbol]
gene_up = DC_g1.response[DC_g1.response$change == 'up','gene'] 
gene_down = DC_g1.response[DC_g1.response$change == 'down','gene'] 
gene_diff = c(gene_up,gene_down)

# up regulated genes
ego_CC <- enrichGO(gene = gene_up,keyType = 'SYMBOL', OrgDb = org.Hs.eg.db,  ont= "CC",pAdjustMethod = "BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_CC))

ego_BP <- enrichGO(gene= gene_up,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_BP))

ego_MF <- enrichGO(gene= gene_up,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_MF))

go <- enrichGO(gene = gene_up, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL',ont="all")
head(go)

p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p
ggsave("p.png",p,height=20)

go<- cnetplot(go, showCategory = 5) 
go
ggsave("go.png",go)

# down regulated genes
ego_CC <- enrichGO(gene = gene_down,keyType = 'SYMBOL', OrgDb = org.Hs.eg.db,  ont= "CC",pAdjustMethod = "BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_CC))

ego_BP <- enrichGO(gene= gene_down,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_BP))

ego_MF <- enrichGO(gene= gene_down,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_MF))

go <- enrichGO(gene = gene_down, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL',ont="all")
head(go)

p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p
ggsave("p.png",p,height=20)

go<- cnetplot(go, showCategory = 5) 
go
ggsave("go.png",go)




#Go enrichment analysis of dendritic cells group 3
DC_g3.response <- data.frame(gene = rownames(DC_g3.response), DC_g3.response)
head(DC_g3.response)

logFC_t=0.5
P.Value_t = 0.05
k1 = (DC_g3.response$p_val_adj < P.Value_t)&(DC_g3.response$avg_log2FC < -logFC_t)
k2 = (DC_g3.response$p_val_adj < P.Value_t)&(DC_g3.response$avg_log2FC > logFC_t)
table(k1)

table(k2)

change = ifelse(k1,"down",ifelse(k2,"up","stable"))
DC_g3.response$change <- change
head(DC_g3.response)

s2e <- bitr(DC_g3.response$gene, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

DC_g3.response <- inner_join(DC_g3.response,s2e,by=c("gene"="SYMBOL"))
head(DC_g3.response)

#GO[Symbol]
gene_up = DC_g3.response[DC_g3.response$change == 'up','gene'] 
gene_down = DC_g3.response[DC_g3.response$change == 'down','gene'] 
gene_diff = c(gene_up,gene_down)

# up regulated genes
ego_CC <- enrichGO(gene = gene_up,keyType = 'SYMBOL', OrgDb = org.Hs.eg.db,  ont= "CC",pAdjustMethod = "BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_CC))

ego_BP <- enrichGO(gene= gene_up,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_BP))

ego_MF <- enrichGO(gene= gene_up,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_MF))

go <- enrichGO(gene = gene_up, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL',ont="all")
head(go)

p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p
ggsave("p.png",p,height=20)

go<- cnetplot(go, showCategory = 5) 
go
ggsave("go.png",go)

# down regulated genes
ego_CC <- enrichGO(gene = gene_down,keyType = 'SYMBOL', OrgDb = org.Hs.eg.db,  ont= "CC",pAdjustMethod = "BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_CC))

ego_BP <- enrichGO(gene= gene_down,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_BP))

ego_MF <- enrichGO(gene= gene_down,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_MF))

go <- enrichGO(gene = gene_down, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL',ont="all")
head(go)

p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p

ggsave("p.png",p,height=20)

go<- cnetplot(go, showCategory = 5) 
go
ggsave("go.png",go)








#Go enrichment analysis of dendritic cells group 4
DC_g4.response <- data.frame(gene = rownames(DC_g4.response), DC_g4.response)
head(DC_g4.response)

logFC_t=0.5
P.Value_t = 0.05
k1 = (DC_g4.response$p_val_adj < P.Value_t)&(DC_g4.response$avg_log2FC < -logFC_t)
k2 = (DC_g4.response$p_val_adj < P.Value_t)&(DC_g4.response$avg_log2FC > logFC_t)
table(k1)
table(k2)

change = ifelse(k1,"down",ifelse(k2,"up","stable"))
DC_g4.response$change <- change
head(DC_g4.response)

s2e <- bitr(DC_g4.response$gene, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

DC_g4.response <- inner_join(DC_g4.response,s2e,by=c("gene"="SYMBOL"))
head(DC_g4.response)

#GO[Symbol]
gene_up = DC_g4.response[DC_g4.response$change == 'up','gene'] 
gene_down = DC_g4.response[DC_g4.response$change == 'down','gene'] 
gene_diff = c(gene_up,gene_down)

# up regulated gene
ego_CC <- enrichGO(gene = gene_up,keyType = 'SYMBOL', OrgDb = org.Hs.eg.db,  ont= "CC",pAdjustMethod = "BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_CC))

ego_BP <- enrichGO(gene= gene_up,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_BP))

ego_MF <- enrichGO(gene= gene_up,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_MF))

go <- enrichGO(gene = gene_up, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL',ont="all")
head(go)

p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p
ggsave("p.png",p,height=20)

go<- cnetplot(go, showCategory = 5) 
go
ggsave("go.png",go)

# down regulated gene
ego_CC <- enrichGO(gene = gene_down,keyType = 'SYMBOL', OrgDb = org.Hs.eg.db,  ont= "CC",pAdjustMethod = "BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_CC))

ego_BP <- enrichGO(gene= gene_down,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_BP))

ego_MF <- enrichGO(gene= gene_down,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_MF))

go <- enrichGO(gene = gene_down, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL',ont="all")
head(go)

p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p

ggsave("p.png",p,height=20)

go<- cnetplot(go, showCategory = 5) 
go
ggsave("go.png",go)




# Go enrichment analysis of macrophage
Macrophage.compare.response <- data.frame(gene = rownames(Macrophage.compare.response), Macrophage.compare.response)
head(Macrophage.compare.response)

logFC_t=0.5
P.Value_t = 0.05
k1 = (Macrophage.compare.response$p_val_adj < P.Value_t)&(Macrophage.compare.response$avg_log2FC < -logFC_t)
k2 = (Macrophage.compare.response$p_val_adj < P.Value_t)&(Macrophage.compare.response$avg_log2FC > logFC_t)
table(k1)

table(k2)

change = ifelse(k1,"down",ifelse(k2,"up","stable"))
Macrophage.compare.response$change <- change
head(Macrophage.compare.response)

s2e <- bitr(Macrophage.compare.response$gene, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

Macrophage.compare.response <- inner_join(Macrophage.compare.response,s2e,by=c("gene"="SYMBOL"))
head(Macrophage.compare.response)

#GO[Symbol]
gene_up = Macrophage.compare.response[Macrophage.compare.response$change == 'up','gene'] 
gene_down = Macrophage.compare.response[Macrophage.compare.response$change == 'down','gene'] 
gene_diff = c(gene_up,gene_down)

# up regulated gene
ego_CC <- enrichGO(gene = gene_up,keyType = 'SYMBOL', OrgDb = org.Hs.eg.db,  ont= "CC",pAdjustMethod = "BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_CC))

ego_BP <- enrichGO(gene= gene_up,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_BP))

ego_MF <- enrichGO(gene= gene_up,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_MF))

go <- enrichGO(gene = gene_up, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL',ont="all")
head(go)

p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p

ggsave("p.png",p,height=20)

go<- cnetplot(go, showCategory = 5) 
go
ggsave("go.png",go)

# down regulated gene
ego_CC <- enrichGO(gene = gene_down,keyType = 'SYMBOL', OrgDb = org.Hs.eg.db,  ont= "CC",pAdjustMethod = "BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_CC))

ego_BP <- enrichGO(gene= gene_down,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_BP))

ego_MF <- enrichGO(gene= gene_down,OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
head(as.data.frame(ego_MF))

go <- enrichGO(gene = gene_down, OrgDb = "org.Hs.eg.db", keyType= 'SYMBOL',ont="all")
head(go)

p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p

ggsave("p.png",p,height=20)

go<- cnetplot(go, showCategory = 5) 
go
ggsave("go.png",go)
