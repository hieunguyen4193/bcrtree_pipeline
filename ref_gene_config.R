# ref.genes <- list(
#   `10x` = file.path(path.to.main.src, "10x/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0/fasta/regions.fa"),
#   IMGT = list(
#     # heavy chain V J gene
#     heavy.J.gene = file.path(path.to.main.src, "FASTA", "IMGT/V_J_genes/IGHJ.fasta"),
#     heavy.V.gene = file.path(path.to.main.src, "FASTA", "IMGT/V_J_genes/IGHV.fasta"),
#     # light chain V J gene
#     light.J.geneK = file.path(path.to.main.src, "FASTA", "IMGT/V_J_genes/IGKJ.fasta"),
#     light.V.geneK = file.path(path.to.main.src, "FASTA", "IMGT/V_J_genes/IGKV.fasta"),
#     light.J.geneL = file.path(path.to.main.src, "FASTA", "IMGT/V_J_genes/IGLJ.fasta"),
#     light.V.geneL = file.path(path.to.main.src, "FASTA", "IMGT/V_J_genes/IGLV.fasta)"
#     )
#   )

ref.genes <- list(
  `10x` = file.path(path.to.main.src, "FASTA", "10x/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0/fasta/regions.fa"),
  IMGT = list(
    # heavy chain V J gene
    J.gene = file.path(path.to.main.src, "FASTA", "IMGT/V_J_genes/IGHJ.fasta"),
    V.gene = file.path(path.to.main.src, "FASTA", "IMGT/V_J_genes/IGHV.fasta")
    )
  )
