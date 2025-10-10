for gfa in ../mw_nanopore_out_all/flye/*/*_graph.gfa
do
    Bandage image "$gfa" "${gfa%.gfa}.png"
done