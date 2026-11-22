#!/usr/bin/env bash
# Figure 6 -- every scripted panel (C-G here; B needs the SCENIC+ tables on lark). Inputs must exist (make_inputs.sh).
#   bash figures/figure6/make_all.sh [B] [C] [D] [E] [F] [G]      # default: C D E F G
here=$(dirname "$0")
panels=${*:-C D E F G}
for p in $panels; do
  f=$(ls "$here"/panel_${p}_*.sh)
  echo "== panel $p: $(basename "$f")"; bash "$f"
done
