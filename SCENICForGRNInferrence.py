# Download pyscenic prelim dataset:

cd ../
# download ranking database
mkdir -p data/ranking_db
cd data/ranking_db
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
# download motif database
cd ../
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl

# Install pyscenic using python virtual env
mkdir code # make sure the code dir is at the same level as the data dir
cd code
python3 -m venv --system-site-packages pyscenic_env #make sure you have python > 3.7
source pyscenic_env/bin/activate # activate the virtual env everytime
pip install pyscenic
pip install igraph

# To test the installation, try to run the following command in python
import os
import glob
import pickle
import pandas as pd
import numpy as np
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from scipy.io import mmread
import igraph as ig

# run SCENIC to infer GRN
import os
import glob
import pickle
import pandas as pd
import numpy as np
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from scipy.io import mmread
import igraph as ig

BASEFOLDER_NAME = '../data'
K562_SCRNASEQ_NAME = '../data/K562_scRNA_seq'
MOTIFS_HGNC_FNAME = os.path.join(BASEFOLDER_NAME, 'motifs-v9-nr.hgnc-m0.001-o0.0.tbl')
CURATED_TFS_HGNC_FNAME = os.path.join(BASEFOLDER_NAME, 'lambert2018.txt')
OUT_TFS_HGNC_FNAME = os.path.join(BASEFOLDER_NAME, 'hs_hgnc_curated_tfs.txt')
"https://github.com/aertslab/pySCENIC/tree/master/resources"
EXP_MAT_FNAME = os.path.join(K562_SCRNASEQ_NAME, 'GSM5687481_k562_rep1_matrix.mtx')
BARCODE_FNAME = os.path.join(K562_SCRNASEQ_NAME, 'GSM5687481_k562_rep1_barcodes.tsv.gz')
GENE_NAME_FNAME = os.path.join(K562_SCRNASEQ_NAME, 'GSM5687481_k562_rep1_features.tsv.gz')
DATABASE_FOLDER = '../data/ranking_db'
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "hg38__refseq-r80__*_tss.mc9nr.genes_vs_motifs.rankings.feather")
REGULONS_FNAME = os.path.join(BASEFOLDER_NAME, "regulons_two_dir.pk")
MOTIFS_FNAME = os.path.join(BASEFOLDER_NAME, "motifs_two_dir.csv")

# prepare the valid motif and gene list
with open(CURATED_TFS_HGNC_FNAME, 'rt') as f:
    hs_curated_tfs = list(map(lambda s: s.strip(), f.readlines()))
len(hs_curated_tfs)
df_motifs_hgnc = pd.read_csv(MOTIFS_HGNC_FNAME, sep='\t')
hs_tfs = df_motifs_hgnc.gene_name.unique()
hs_curated_tfs_with_motif = list(set(hs_tfs).intersection(hs_curated_tfs))
len(hs_curated_tfs_with_motif)
with open(OUT_TFS_HGNC_FNAME, 'wt') as f:
    f.write('\n'.join(hs_curated_tfs_with_motif) + '\n')

# format K562 scRNA-seq data
exp_mat = mmread(EXP_MAT_FNAME)
exp_mat = pd.DataFrame.sparse.from_spmatrix(exp_mat).T
barcode = pd.read_csv(BARCODE_FNAME,sep='\t',header=None)
gene_name = pd.read_csv(GENE_NAME_FNAME,sep='\t',header=None)
exp_mat.columns = gene_name[1].values.tolist()
exp_mat.insert(0,'barcode',barcode[0].values.tolist())


# pyscenic workflow
tf_names = load_tf_names(OUT_TFS_HGNC_FNAME)
adjacencies = grnboost2(exp_mat, tf_names=tf_names, verbose=True)
adjacencies = pd.DataFrame(adjacencies)
adjacencies.to_csv('../data/grnboost2_out.csv')

# load ranking database
db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

# add TF to the module
exp_mat = exp_mat.drop('barcode',axis=1)
modules = list(modules_from_adjacencies(adjacencies, exp_mat,keep_only_activating=False))
if __name__ == '__main__':
    with ProgressBar():
        df = prune2df(dbs, modules, MOTIFS_HGNC_FNAME)

    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)

    # Save the enriched motifs and the discovered regulons to disk.
    df.to_csv(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "wb") as f:
        pickle.dump(regulons, f)
        
# generate a heatmap
tf_names = []
target_gene = []
score = dict()
for i in regulons:
    tmp_tf_names = i.transcription_factor
    tmp_score = dict(i.gene2weight)
    if '-' in i.name:
        for k in tmp_score:
            tmp_score[k] *= -1
    score[tmp_tf_names] = tmp_score
    target_gene.extend(list(tmp_score.keys()))
    tf_names.append(tmp_tf_names)

target_gene = [*set(target_gene)]

scenic_score_mat = np.zeros((len(target_gene),len(tf_names)))
for idx,i in enumerate(tf_names):
    for j in score[i]:
        tmp_idx = target_gene.index(j)
        tmp_score = score[i][j]
        scenic_score_mat[tmp_idx,idx] = tmp_score

# visualize the graph
edge_list = []
node_list = set()
for i in regulons:
    tmp_tf = i.transcription_factor
    direction = -1 if '-' in i.name else 1
    node_list.add(tmp_tf)
    denorm = sum(list(i.gene2weight.values()))
    sign = 1 if denorm > 0 else -1
    for k in i.gene2weight:
        tmp_tg = k
        tmp_val = sign*i.gene2weight[k]/denorm
        edge_list.append([tmp_tf,tmp_tg,direction*tmp_val])
        node_list.add(tmp_tg)

node_id_mapping = dict()
for idx,i in enumerate(node_list):
    node_id_mapping[i] = idx



edge_frame = pd.DataFrame(edge_list)
edge_frame[0] = edge_frame[0].map(node_id_mapping)
edge_frame[1] = edge_frame[1].map(node_id_mapping)
filtered_edge_frame_pos = edge_frame[edge_frame[2]>0.05]
filtered_edge_frame_neg = edge_frame[edge_frame[2]<0]
filtered_edge_frame = pd.concat([filtered_edge_frame_pos,filtered_edge_frame_neg])

scenic_regulon = ig.Graph.TupleList(filtered_edge_frame.itertuples(index=False), directed=True, weights=False, edge_attrs="weight")
components = scenic_regulon.connected_components()

visual_style = {}
visual_style['vertex_size'] = 3
#visual_style['vertex_color'] = list(map(int, ig.rescale(components.membership, (0, 200), clamp=True)))
visual_style['edge_width'] = 0.01
visual_style['edge_arrow_size'] = 0.07
visual_style['edge_arrow_width'] = 0.07
edge_color = []
for i in scenic_regulon.es()['weight']:
    if i < 0:
        edge_color.append('black')
    else:
        edge_color.append('red')
visual_style['edge_color'] = edge_color


ig.plot(
    components,
    '../data/ig_plot.pdf',
    palette=ig.RainbowPalette(),
    **visual_style
)

# write out formated data
filtered_edge_frame.to_csv("../data/filtered_edge_list.csv", index=None)
tf_names = pd.DataFrame(tf_names)
tf_names.to_csv("../data/tf_list.csv",index=None)
node_id_mapping = pd.DataFrame(node_id_mapping.items())
node_id_mapping = node_id_mapping.to_csv("../data/gene_name_id_mapping.csv",index=None)