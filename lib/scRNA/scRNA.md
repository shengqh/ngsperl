# pool_sample

Sometime, we will want to pool some sample together.

```
  pool_sample        => 1,
  pool_sample_groups => {
    "iSGS" => [ "iSGS1",  "iSGS2",  "iSGS3" ],
    "LTS"  => [ "LTS1_1", "LTS3_1", "LTS4_3" ],
    "GPA" => ["GPA5_3"],
  },
```

# batch_for_integration

We can perform integration based on batch. The samples in same batch will just be merged.

```
  batch_for_integration => 1,
  batch_for_integration_groups => {
    "3364" => [ "iSGS" ],
    "3855" => [ "LTS", "GPA" ]
  },
```

# Differential expression analysis

## perform_edgeR

We can perform differential expression analysis on different level, such as between samples in same celltype, or in same cluster. We can perform analysis on cell level using sample as covariance, or on sample level which merges all cells in same sample.

```
  #differential expression analysis
  perform_edgeR => 1,

  DE_by_celltype => 0,
  DE_by_cluster  => 1,
  DE_by_sample   => 0,
  DE_by_cell     => 1,

  DE_pvalue         => 0.05,
  DE_use_raw_pvalue => 0,
  DE_fold_change    => 1.5,

  groups => {
    "iSGS" => ["iSGS"],
    "LTS"  => ["LTS"],
    "GPA"  => ["GPA"],
  },

  pairs => {
    "iSGS_vs_LTS" => [ "LTS", "iSGS" ],
    "iSGS_vs_GPA" => [ "GPA", "iSGS" ],
    "GPA_vs_LTS"  => [ "LTS",  "GPA" ],
  },

```

## DE_cluster_pairs

By define DE_cluster_pairs, we can also perform differential expression analysis between different clusters.

```
  DE_cluster_groups => {
     "Microphage_01" => [ "1" ],
     "Microphage_06" => [ "6" ],
     "Microphage_20" => [ "20" ],
     "Microphage_Not01" => [ "6", "20" ],
     "Microphage_Not06" => [ "1", "20" ],
     "Microphage_Not20" => [ "1", "6" ],
  },
  
  DE_cluster_pairs => {
     "Microphage_01" => ["Microphage_Not01", "Microphage_01"],
     "Microphage_06" => ["Microphage_Not06", "Microphage_06"],
     "Microphage_20" => ["Microphage_Not20", "Microphage_20"],
  },
```

## perform_webgestalt

```
  #DE gene annotation
  perform_webgestalt  => 1,
  webgestalt_organism => "hsapiens",
```

## perform_gsea

```
  perform_gsea    => 1,
  gsea_jar        => "gsea-cli.bat",
  gsea_db         => "C:/projects/database/gsea/v7.1/",
  gsea_categories => "'h.all.v7.1.symbols.gmt', 'c2.all.v7.1.symbols.gmt', 'c5.all.v7.1.symbols.gmt', 'c6.all.v7.1.symbols.gmt', 'c7.all.v7.1.symbols.gmt'",
  gsea_makeReport => 0,
```

# perform_rename_cluster

Based on manual interpretation, we may need to rename cell type of some clusters.

```
  perform_rename_cluster => 1,
  rename_cluster     => {
    "0"  => "CD4_Naive",
    "2"  => "CD8_nonTRM",
    "3"  => "NK",
    "5"  => "CD4_Activated",
    "7"  => "CD4_Treg",
    "8"  => "gammaDelta",
    "9"  => "mDC",
    "10"  => "DN_Tcell",
    "13"  => "CD8_TRM",
    "14"  => "CD4_TFH",
    "16"  => "Erythrocyte",
    "17" => "CD8_EM",
    "19" => "pDC",
    "22" => "Ciliated"
  },
  
```

# perform_marker_dotplot

For specific clusters, we perform Seurat::FindAllMarkers for those cluster only, trying to find marker genes which can represent each cluster in those clusters.

```
  perform_marker_dotplot => 1,
  marker_dotplot_clusters => {
    "Microphage" => ["1","6","20"],
  },
  marker_dotplot_gene_number => 20,
```

# perform_curated_gene_dotplot

We can draw dotplot for specific geneset on specific clusters. You can have multiple definitions, each contains key clusters and genes.

```
  perform_curated_gene_dotplot => 1,
  curated_gene_dotplot => {
    "Microphage" => {
      clusters => ["1","6","20"],
      genes => {
        "IL-17A_pathway_1" => [qw(IL17A
    IL17RC
    IL17F
    TLR4)],
      "mTOR" => [qw(PIH1D1
    SPAAR
    DYRK3)],
      }
    }
  },
```