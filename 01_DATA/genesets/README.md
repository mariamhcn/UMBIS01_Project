# Gene set files

One CSV per signature (or collection), with two columns: `gs_name` (name of the gene set) and `gene_symbol` (human gene
symbol, one row per gene). The notebooks read these files directly: nothing needs to be downloaded.

## Q3 – pathway enrichment

| File | Content | Source |
|------|---------|--------|
| `Q3_hallmark.csv` | MSigDB Hallmark collection (50 sets) | MSigDB v7.5.1 [Liberzon et al. 2015, *Cell Systems*] |
| `Q3_reactome.csv` | Reactome pathways (MSigDB C2:CP:REACTOME, 1,615 sets) | MSigDB v7.5.1 |
| `Q3_syngo.csv` | SynGO: the two top-level terms ("synapse", "process in the synapse") and their direct children (13 sets); each term includes the genes of its sub-terms | SynGO release 1.3, table `ontologies.xlsx` [Koopmans et al. 2019, *Neuron*] |

GO Biological Process is read from the `org.Hs.eg.db` package, not from a file.

## Q5 – ageing signatures

| File | Signature | Source |
|------|-----------|--------|
| `Q5_senmayo.csv` | SenMayo, senescence (124 genes) | Saul et al. 2022, *Nature Communications*, Supplementary Data 1 (MSigDB gene set `SAUL_SEN_MAYO`) |
| `Q5_inflammatory_response.csv` | Hallmark Inflammatory Response (200 genes) | MSigDB v7.5.1 |
| `Q5_tnfa_nfkb.csv` | Hallmark TNFα signalling via NF-κB (200 genes) | MSigDB v7.5.1 |
| `Q5_oxphos.csv` | Hallmark Oxidative Phosphorylation (200 genes) | MSigDB v7.5.1 |
| `Q5_myelination.csv` | GO "ensheathment of neurons" (GO:0007272, contains myelination GO:0042552; 140 genes) | MSigDB v7.5.1, C5:GO:BP |
| `Q5_microglia_ageing_galatro.csv` | genes increased (`GALATRO_AGING_UP`, 212) or decreased (`GALATRO_AGING_DOWN`, 360) with age in purified human cortical microglia | Galatro et al. 2017, *Nature Neuroscience*, Supplementary Table 7 |
| `Q5_lifespan_ageing.csv` | genes rising (`AGING_UP_<subclass>`) or falling (`AGING_DOWN_<subclass>`) with age in adult DLPFC, per subclass: the authors' age-associated genes (adjusted p < 0.05) of the young, middle and late adulthood groups (20–39, 40–59, ≥ 60 years; development excluded); sets with fewer than 10 genes removed (19 sets) | Yang et al. 2026, *Nature*, Supplementary Data 4 |
