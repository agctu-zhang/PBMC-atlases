#default value
input_loom=out.loom
n_workers=20
#help function
function usage() {
echo -e "OPTIONS:\n-i|--input_loom:\t input loom file"
echo -e "-n|--n_workers:\t working core number"
echo -e "-h|--help:\t Usage information"
exit 1
}
#get value
while getopts :i:n:h opt
do
    case "$opt" in
        i) input_loom="$OPTARG" ;;
        n) n_workers="$OPTARG" ;;
        h) usage ;;
        :) echo "This option -$OPTARG requires an argument."
           exit 1 ;;
        ?) echo "-$OPTARG is not an option"
           exit 2 ;;
    esac
done
#需要更改路径
tfs=~/download/database/SCENIC/allTFs_hg38.txt
feather=~/download/database/SCENIC/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
tbl=~/download/database/SCENIC/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
pyscenic=~/miniconda3/envs/zsy_sc/bin/pyscenic

# grn
$pyscenic grn \
--num_workers $n_workers \
--output ${input_loom%.*}.grn.tsv \
--method grnboost2 \
$input_loom  $tfs

# cistarget
$pyscenic ctx \
${input_loom%.*}.grn.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output ${input_loom%.*}.ctx.csv \
--num_workers $n_workers   \
--mask_dropouts

# AUCell
$pyscenic aucell \
$input_loom \
${input_loom%.*}.ctx.csv \
--output ${input_loom%.*}.AUCell.loom \
--num_workers $n_workers
