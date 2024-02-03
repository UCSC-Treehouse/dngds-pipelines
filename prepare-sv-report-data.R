## Load packages
library(sveval)
library(GenomicRanges)
library(dplyr)

## for debug:
cfg = c(pli_gene='gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz',
        simprep='simpleRepeat.txt.gz',
        hsvlr='hsvlr.vcf.gz',
        dgv='GRCh38_hg38_variants_2016-08-31.txt',
        clingen='iscaPathogenic.txt.gz',
        ctcf='ENCFF010WHH.bed.gz',
        gnomadsv='gnomad_v2_sv.sites.pass.lifted.vcf.gz',
        cre='ENCFF166QIT.lifted.bed.gz',
        cons='phastConsElements100way.txt.gz',
        gencode='gencode.v35.annotation.gtf.gz',
        cytoband='cytoBandIdeo.txt.gz')

message('Gene annotation...')
types.ranked = c('CDS', 'UTR', 'promoter', 'gene')
types.labels = c('coding', 'UTR', 'promoter', 'intronic')
genc = rtracklayer::import(cfg['gencode'])
genc = subset(genc, type %in% types.ranked)
prom = promoters(subset(genc, type=='gene'))
prom$type = 'promoter'
genc = c(genc, prom)
mcols(genc) = mcols(genc)[,c('type', 'gene_name', 'gene_type')]

message('Import pli scores...')
pli.df = read.table(cfg['pli_gene'], as.is=TRUE, header=TRUE, sep='\t')
pli.df = pli.df %>% select(gene, pLI) %>% unique

message('Overlap with simple repeats...')
sr = read.table(cfg['simprep'])
sr = sr[,c(2,3,4)]
colnames(sr) = c('chrom', 'start', 'end')
sr = makeGRangesFromDataFrame(sr, keep.extra.columns=TRUE)
sr = reduce(sr)

message('Overlap with DGV catalog...')
dgv = read.table(cfg['dgv'], as.is=TRUE, header=TRUE, sep='\t')
dgv = dgv[,c('chr','start','end', 'variantsubtype')]
dgv$chr = paste0('chr', dgv$chr)
dgv = makeGRangesFromDataFrame(dgv, keep.extra.columns=TRUE)
dgv.loss = subset(dgv, variantsubtype %in% c('loss', 'deletion', 'gain+loss'))
dgv.gain = subset(dgv, variantsubtype %in% c('duplication', 'gain', 'insertion',
                                             'mobile element insertion',
                                             'novel sequence insertion',
                                             'tandem duplication', 'gain+loss'))

message('Overlap with ClinGenn pathogenic CNVs...')
clingen = read.table(cfg['clingen'], as.is=TRUE, sep='\t')
clingen = GRanges(clingen$V2, IRanges(clingen$V3, clingen$V4))

message('Overlap with CTCF peaks...')
ctcf = read.table(cfg['ctcf'], as.is=TRUE)
ctcf = with(ctcf, GRanges(V1, IRanges(V2, V3), score=V5))
ctcf = reduce(ctcf)

message('Overlap with regulatory regions...')
cres = read.table(cfg['cre'], as.is=TRUE)
cres = with(cres, GRanges(V1, IRanges(V2, V3), score=V5))
cres = reduce(cres)

message('Overlap with conserved regions...')
cons.gr = read.table(cfg['cons'], as.is=TRUE)
cons.gr = with(cons.gr, GRanges(V2, IRanges(V3, V4)))
cons.gr = reduce(cons.gr)

message('Annotate frequency using gnomAD...')
gnomad = readSVvcf(cfg['gnomadsv'], other.field='AF')
gnomad$ac = 1

message('Overlap with SVs from public long-read studies...')
lr.gr = readSVvcf(cfg['hsvlr'])
lr.gr$ac = 1

rm.after.annot = c('genc', 'sr', 'dgv.loss', 'dgv.gain', 'clingen', 'cres',
                   'cons.gr', 'gnomad', 'lr.gr')

save(rm.after.annot, genc, types.ranked, types.labels,
     pli.df, sr, dgv.loss, dgv.gain, clingen, cres, ctcf,
     cons.gr, gnomad, lr.gr, file='sv-report-db.RData')
