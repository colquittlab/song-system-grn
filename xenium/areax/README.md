# Locating Area X in the Xenium sections

Goal is anatomical: mark where Area X sits in these ten **sagittal** sections so
it can be drawn on a schematic. Not a cell-type claim.

## Gene list provenance

`zebra_areax_markers.csv` is a copy of a scrape of the **ZEBrA portal**
(Mello lab, `zebrafinchatlas.org`, *Markers of the Song System -> Area X*),
originally made for a separate project. It is public-database content, not
anyone's unpublished data, and it is copied here so this analysis does not
reach into another project's tree at run time.

**No snRNA-seq data from that project is used anywhere in this analysis** --
only the gene list. Everything else comes from this project's Xenium data and
its own hybrid-label RCTD calls.

ZEBrA's contrast is Area X **versus adjacent brain tissue** in non-singing
adult males. That is the right contrast for this purpose: the markers separate
Area X from the striatum around it, which is exactly the boundary a schematic
needs. It also means they are weak at separating striatum from non-striatum,
so the index is only computed inside striatal (LGE) calls.

## Panel overlap

17 of the 73 ZEBrA marker genes are on the 425-gene panel (13 general/Up,
3 general/Down, 3 sparse/Up; CHRNA4 and NEFL appear in two sets).

`NEFL`, `SV2B` and `UCHL1` are excluded from the primary index: they are
general neuronal-maturity genes rather than regional ones, so they would pull
the index toward whichever cells are simply more mature. Both versions are
computed; see `areax_index_gene_sensitivity.csv`.
