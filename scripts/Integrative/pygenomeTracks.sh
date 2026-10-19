# nolint start
####################################################################################
#' Project: Repeated malaria exposure reprograms Vδ2+ γδ T cell chromatin accessibility and immune function toward a reversible hyporesponsive state
#' Script: pyegnomeTracks.R
#' Author: MIJM
#' Last Update: 27.08.26
#' Description: Plotting genome tracks

####################################################################################

#!/bin/bash

# Convert into bedGraph
for f in *.bed; do
  if [ -f "$f" ]; then
    awk '{printf "%s\t%d\t%d\t%2.3f\n", $1,$2,$3,$5}' "$f" > "${f%.bed}.bedgraph"
  fi
done


for p2g_file in *_p2g.links; do
  
  filename=$(basename "$p2g_file")
  gene="${filename%_p2g.links}"
  view_file="${gene}_view.links"

  echo "========================================"
  echo "Generating PyGenomeTracks for: $gene"

  if [ ! -s "$p2g_file" ]; then
    echo "No linkages found for $gene. Skipping."
    continue
  fi

  chr=$(head -n 1 "$view_file" | awk '{print $1}')


  # start=$(awk -F',' -v g="$gene" '$1 == g || $2 == g { s=$3; print (s - 250000 < 0 ? 0 : s - 250000) }' tss.csv)
  # end=$(awk -F',' -v g="$gene" '$1 == g || $2 == g { print $3 + 250000 }' tss.csv)

  # start=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$p2g_file" | sort -n | head -n 1 | awk '{print $1-10000}')
  # end=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$p2g_file" | sort -n | tail -n 1 | awk '{print $1+10000}')

  start=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | head -n 1 | awk '{print $1-10000}')
  end=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | tail -n 1 | awk '{print $1+10000}')
  
  # : "${start:=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | head -n 1 | awk '{print $1-10000}')}"
  # : "${end:=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | tail -n 1 | awk '{print $1+10000}')}"

  echo "Plotting region: ${chr}:${start}-${end}"

  max_y=$(awk -v c="$chr" -v s="$start" -v e="$end" '
    $1 == c && $3 >= s && $2 <= e {
      if ($4 > max) max = $4
    }
    END {
      if (max == "") max = 1
      print int(max * 1.1) + 1
    }
  ' *.bedgraph)

  echo "Calculated shared max y-axis: $max_y"

  sed -e "s/MAX_Y/$max_y/g" \
      -e "s|GENE_P2G_LINKS|$p2g_file|g" \
      -e "s|GENE_VIEW_LINKS|$view_file|g" \
      -e "s|GENE_VLINES|${gene}_vlines.bed|g" ./MAIN.ini > "MAIN_${gene}.ini"

  pyGenomeTracks --tracks "MAIN_${gene}.ini" --region "${chr}:${start}-${end}" -o "./plots/mean_full_${gene}.pdf"

unset start end

done

# zoomed in

ls *_p2g.links


p2g_file="PLCG2_p2g.links"
start=81770000
end=81800000

filename=$(basename "$p2g_file")
gene="${filename%_p2g.links}"
view_file="${gene}_view.links"

echo "========================================"
echo "Generating PyGenomeTracks for: $gene"

if [ ! -s "$p2g_file" ]; then
  echo "No linkages found for $gene. Skipping."
  continue
fi

chr=$(head -n 1 "$view_file" | awk '{print $1}')


: "${start:=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | head -n 1 | awk '{print $1-10000}')}"
: "${end:=$(awk '{print $2"\n"$3"\n"$5"\n"$6}' "$view_file" | sort -n | tail -n 1 | awk '{print $1+10000}')}"

echo "Plotting region: ${chr}:${start}-${end}"

max_y=$(awk -v c="$chr" -v s="$start" -v e="$end" '
  $1 == c && $3 >= s && $2 <= e {
    if ($4 > max) max = $4
  }
  END {
    if (max == "") max = 1
    print int(max * 1.1) + 1
  }
' *.bedgraph)

echo "Calculated shared max y-axis: $max_y"

sed -e "s/MAX_Y/$max_y/g" \
    -e "s|GENE_P2G_LINKS|$p2g_file|g" \
    -e "s|GENE_VIEW_LINKS|$view_file|g" \
    -e "s|GENE_VLINES|${gene}_vlines.bed|g" ./MAIN.ini > "MAIN_${gene}.ini"

pyGenomeTracks --tracks "MAIN_${gene}.ini" --region "${chr}:${start}-${end}" -o "./plots/mean_view_reduced_${gene}.pdf"

unset start end

