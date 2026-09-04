"""Anatomical-position row colours for the glut-only hybrid-label heatmaps, replacing the
earlier song(red)/non-song(blue) binary scheme with the lab's standard per-region palette
(group/colquittlab-utils/R/common_aesthetics.R's `position_colors`), keyed to each hybrid
cluster's actual dominant 'position' value in bf_adult_hybrid.h5ad (verified empirically --
e.g. Glut-DACH2-1..4 are ~80-90% 'nc', Glut-DACH2-5..8 are ~65-90% 'nr', Glut-CACNA1H-1..4
and Glut-CACNA1H-RA are 100% 'arco'/'ra').

common_aesthetics.R has no 'nc'/'nr' entries (those are this project's actual data labels,
not the R file's palette keys) -- 'nc' is mapped to that file's 'ncl' colour (nidopallium
caudolaterale, the closest named match to 'nc' = nidopallium caudale) and 'nr' to its
generic 'nido' colour (nidopallium, standing in for the rostral portion since no dedicated
'nr' swatch exists there). Legend labels use this project's own 'NC'/'NR' terms, not the
R file's 'NCL'/'Nido' pretty names, to stay consistent with how this analysis already
talks about these regions.
"""

# Subset of group/colquittlab-utils/R/common_aesthetics.R's `position_colors` (the
# non-reversed variant) needed here, copied verbatim.
POSITION_COLORS = {
    "hvc": "#ff7f00",
    "lman": "#1f78b4",
    "ra": "#e31a1c",
    "nc": "#fdbf6f",    # common_aesthetics.R's 'ncl'
    "nr": "#a6cee3",    # common_aesthetics.R's 'nido'
    "arco": "#fb9a99",
}
POSITION_LABELS = {"hvc": "HVC", "lman": "LMAN", "ra": "RA", "nc": "NC", "nr": "NR", "arco": "Arco"}
# Display order for the legend -- song nuclei first (matching ROW_ORDER's grouping), then
# the three non-song regions.
POSITION_ORDER = ["hvc", "lman", "ra", "nc", "nr", "arco"]

CLUSTER_POSITION = {
    "Glut-DACH2-HVCra": "hvc",
    "Glut-DACH2-HVCra-Pre": "hvc",
    "Glut-DACH2-HVCx": "hvc",
    "Glut-DACH2-LMANco": "lman",
    "Glut-DACH2-LMANsh": "lman",
    "Glut-CACNA1H-RA": "ra",
    "Glut-DACH2-1": "nc", "Glut-DACH2-2": "nc", "Glut-DACH2-3": "nc", "Glut-DACH2-4": "nc",
    "Glut-DACH2-5": "nr", "Glut-DACH2-6": "nr", "Glut-DACH2-7": "nr", "Glut-DACH2-8": "nr",
    "Glut-CACNA1H-1": "arco", "Glut-CACNA1H-2": "arco", "Glut-CACNA1H-3": "arco", "Glut-CACNA1H-4": "arco",
}

# Song-nucleus cluster membership (hvc/lman/ra positions), matching the set used
# throughout this analysis's heatmaps/stats.
SONG_CLUSTERS = {c for c, p in CLUSTER_POSITION.items() if p in ("hvc", "lman", "ra")}

# Glut-DACH2-HVCra-Pre: an edge-case putative HVCra precursor population, excluded
# everywhere in this analysis (heatmaps, stats, and now UMAPs) per explicit user request.
EXCLUDED_CLUSTERS = {"Glut-DACH2-HVCra-Pre"}


def cluster_color(c):
    return POSITION_COLORS[CLUSTER_POSITION[c]]


def song_group(c):
    return "song" if c in SONG_CLUSTERS else "non-song"
