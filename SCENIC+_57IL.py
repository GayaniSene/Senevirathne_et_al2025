
#Scanpy, PycistTarget, Cis-topicTarget and pycistopic codes; SCENIC+ pipeline at the end.

#Scenic+ input files 
#Can get it from scanpy (refer to SCENIC+_E53_IL code for scanpy RNA-processing)

or from Seurat

#cell_data table with annotations and barcodes for the input file for SCENIC+ pre-ATAC analyses 
#to export barcodes and celltypes use R seurat dataset

library(tibble)
export_df <- pbmc@meta.data %>% 
  rownames_to_column("barcodes")

head(export_df)
write.csv(export_df, "cell_data.csv")

#the csv file neesds to be saved as .tsv after removing the first row
#After you export the cell_data.tsv make sure to remove the first column. 



#scATAC-seq pre-procssing for SCENIC+
#created a scenicplus enviornment with python 3.11

source activate scenicplus2


#copy the fragment.tz and fragment tx.tbi and the seurat output with cell annotation files into a data folder
#all the scenic plus commands are written in python. 
#To run the commands in the python environment, type
which python
#copy paste what you get from here and run again. 


#!/usr/bin/env python
import pycisTopic
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os 
import pandas as pd
import pickle
import pandas as pd
import scanpy as sc


#Let’s make some directories to store the output of pycisTopic
import os
out_dir = "outs"
os.makedirs(out_dir, exist_ok = True)


#Here "SeuratProject" is the name of the "orig.ident" name from the cell_data.tsv file. 
#pycisTopic will automatically append sample ids to barcodes to avoid barcode collisions between samples.
#make sure to copy paste the "atac_fragments.tsv.gz" file to the data foler along with the tsv.tbi.gz file and the cell_data.tsv file. 
fragments_dict = {"SeuratProject": "data/atac_fragments.tsv.gz"}


#First we read the barcode-to-cell type annotation as a pd.DataFrame
import pandas as pd
cell_data = pd.read_table("data/cell_data.tsv", index_col = 0)
cell_data.head()



#size of each chromosome, this can be downloaded from the UCSC databases.
chromsizes = pd.read_table(
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
chromsizes.head()

from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
os.makedirs(os.path.join(out_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)



bw_paths, bed_paths = export_pseudobulk(
    input_data = cell_data,
    variable = "celltype",
    sample_id_col = "orig.ident",
    chromsizes = chromsizes,
    bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    bigwig_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    path_to_fragments = fragments_dict,
    n_cpu = 1,
    normalize_bigwig = True,
    temp_dir = "/tmp",
    split_pattern = "-"
)


#We will need the paths to the bed files later on, so let’s save them to disk.
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")
        
        
with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv"), "wt") as f:
    for v in bw_paths:
        _ = f.write(f"{v}\t{bw_paths[v]}\n")
        


#Next we will use MACS to call peaks for each pseudobulk fragments.tsv.gz file.
bw_paths = {}
with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bw_paths.update({v: p})

bed_paths = {}
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})       


from pycisTopic.pseudobulk_peak_calling import peak_calling
macs_path = "macs2"

os.makedirs(os.path.join(out_dir, "consensus_peak_calling/MACS"), exist_ok = True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.path.join(os.path.join(out_dir, "consensus_peak_calling/MACS")),
    genome_size = 'hs',
    n_cpu = 10,
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05,
    _temp_dir = '/mscstemp'
)


#before removing the blacklisted genomic regions, make sure to copy paste your pycisTopic downloaded folder to your working directory. 

from pycisTopic.iterative_peak_calling import get_consensus_peaks

# Other param
peak_half_width=250
path_to_blacklist="pycisTopic/blacklist/hg38-blacklist.v2.bed"

# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = path_to_blacklist)


consensus_peaks.to_bed(
    path = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep =True,
    compression = 'infer',
    chain = False)



#NEW code 19th May 2024 for generating the CisTopic object: 

fragments_dict = {
    "SeuratProject": "data/atac_fragments.tsv.gz"
}


wget https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz
!tar -xf Mallet-202108-bin.tar.gz



#!/usr/bin/env python

out_dir="outs"
import os

fragments_dict = {"SeuratProject": "data/atac_fragments.tsv.gz"}

from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
import matplotlib.pyplot as plt

for sample_id in fragments_dict:
    fig = plot_sample_stats(
        sample_id = sample_id,
        pycistopic_qc_output_dir = "outs/qc"
    ) 

from pycisTopic.qc import get_barcodes_passing_qc_for_sample
sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in fragments_dict:
    (
        sample_id_to_barcodes_passing_filters[sample_id],
        sample_id_to_thresholds[sample_id]
    ) = get_barcodes_passing_qc_for_sample(
            sample_id = sample_id,
            pycistopic_qc_output_dir = "outs/qc",
            unique_fragments_threshold = None, # use automatic thresholding
            tss_enrichment_threshold = None, # use automatic thresholding
            frip_threshold = 0,
            use_automatic_thresholds = True,
    )    


for sample_id in fragments_dict:
    fig = plot_barcode_stats(
        sample_id = sample_id,
        pycistopic_qc_output_dir = "outs/qc",
        bc_passing_filters = sample_id_to_barcodes_passing_filters[sample_id],
        detailed_title = False,
        **sample_id_to_thresholds[sample_id]
    )


path_to_regions = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
path_to_blacklist = "pycisTopic/blacklist/hg38-blacklist.v2.bed"
pycistopic_qc_output_dir = "outs/qc"

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl

cistopic_obj_list = []
for sample_id in fragments_dict:
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[ sample_id_to_barcodes_passing_filters[sample_id] ]
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
        n_cpu = 1,
        project = sample_id,
        split_pattern = '-'
    )
    cistopic_obj_list.append(cistopic_obj)


cistopic_obj = cistopic_obj_list[0]
print(cistopic_obj)

import pickle
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)

import pandas as pd
cell_data = pd.read_table("data/cell_data.tsv", index_col = 0)
cell_data.head()
cistopic_obj.add_cell_data(cell_data, split_pattern='-')
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)

cistopic_obj.cell_data

import os

#running models
os.environ['MALLET_MEMORY'] = '200G'
from pycisTopic.lda_models import run_cgs_models_mallet
# Configure path Mallet
mallet_path="Mallet-202108/bin/mallet"
# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    n_cpu=12,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path="mallet/tutorial",
    save_path="mallet/tutorial",
    mallet_path=mallet_path,
)



pickle.dump(
    models,
    open(os.path.join(out_dir, "models.pkl"), "wb")
)


################# can open a .py file
#!/usr/bin/env python
import pickle
import os
out_dir = "outs"
#loading the model file should be in the current working directory
#change the working directory to the output folder

models = pickle.load (open("models.pkl", 'rb'))

from pycisTopic.lda_models import evaluate_models
model = evaluate_models(
    models,
    select_model = 40, save="model.pdf",
    return_model = True
)

cistopic_obj = pickle.load (open ("cistopic_obj.pkl", 'rb'))

cistopic_obj.add_LDA_model(model)

pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)

###########################################################

#read the cistopicObj once the model data is added. can run it on the laptop
#!/usr/bin/env python
import pickle
import os
cistopic_obj = pickle.load (open ("cistopic_obj.pkl", 'rb'))

from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)

find_clusters(
    cistopic_obj,
    target  = 'cell',
    k = 10,
    res = [0.6, 1.2, 3],
    prefix = 'pycisTopic_',
    scale = True,
    split_pattern = '-'
)

#run UMAP
run_umap(
    cistopic_obj,
    target  = 'cell', scale=True)

#run tsne
run_tsne(
    cistopic_obj,
    target  = 'cell', scale=True)

#mapping the UMAP

plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['celltype', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],
    target='cell', num_columns=4,
    text_size=10, save="UMAP_pycisTopic.pdf",
    dot_size=5)


#annotate the dataset
annot_dict = {}
for resolution in [0.6, 1.2, 3]:
    annot_dict[f"pycisTopic_leiden_10_{resolution}"] = {}
    for cluster in set(cistopic_obj.cell_data[f"pycisTopic_leiden_10_{resolution}"]):
        counts = cistopic_obj.cell_data.loc[
            cistopic_obj.cell_data.loc[cistopic_obj.cell_data[f"pycisTopic_leiden_10_{resolution}"] == cluster].index,
            "celltype"].value_counts()
        annot_dict[f"pycisTopic_leiden_10_{resolution}"][cluster] = f"{counts.index[counts.argmax()]}({cluster})"

annot_dict

for resolution in [0.6, 1.2, 3]:
    cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'] = [
        annot_dict[f'pycisTopic_leiden_10_{resolution}'][x] for x in cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'].tolist()
    ]

 #plotting UMAP with annotated cell clusters  

plot_metadata(cistopic_obj, reduction_name='UMAP', variables=['celltype', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'], target='cell', num_columns=4, text_size=10, save="UMAP_annoted_cell_Clusters.pdf", dot_size=5)


#plotting contnuos data
plot_metadata(cistopic_obj,reduction_name='UMAP',variables=['log10_unique_fragments_count', 'tss_enrichment', 'fraction_of_fragments_in_peaks'], target='cell', num_columns=4, save="UMAP_continuous_data.pdf", text_size=10, dot_size=5)

#visualizing cell topic contributions

plot_topic(
    cistopic_obj,
    reduction_name = 'UMAP',
    target = 'cell', save="UMAP_celltype_contribution.pdf",
    num_columns=5
)

#heatmap with cell topic and annotations

cell_topic_heatmap(
    cistopic_obj,
    variables = ['celltype'],
    scale = False,
    legend_loc_x = 1.0,
    legend_loc_y = -1.2,
    legend_dist_y = -1,
    figsize = (10, 10), save="celltopic_Heatmap.pdf",
)


#Topic binarization & QC

#!/usr/bin/env python
import pycisTopic
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import pickle
import pandas as pd
import scanpy as sc


import pickle
import os
cistopic_obj = pickle.load (open ("cistopic_obj.pkl", 'rb'))
out_dir = "outs"


from pycisTopic.topic_binarization import binarize_topics

region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=True, save="topics3k.pdf", num_columns=5
)

region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, save="topics_otsu.pdf", num_columns=5
)

binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True, save="bonarized_topic.pdf",
    num_columns=5, nbins=100)


from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img


topic_qc_metrics = compute_topic_metrics(cistopic_obj)

fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)


topic_annot = topic_annotation(
    cistopic_obj,
    annot_var='celltype',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)






#Finding DARs/differentially accessible regions


from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
import numpy as np

imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6)

normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)

markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='celltype',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=1,
    _temp_dir='/tmp',
    split_pattern = '-'
)





##### saving the gene sets

os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok = True)


from pycisTopic.utils import region_names_to_coordinates

for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False)


for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )



  input_data:
  cisTopic_obj_fname: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/cistopic_obj.pkl"
  GEX_anndata_fname: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/adata.h5ad"
  region_set_folder: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/region_sets"
  ctx_db_fname: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
  dem_db_fname: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
  path_to_motif_annotations: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

output_data:
  # output for prepare_GEX_ACC .h5mu
  combined_GEX_ACC_mudata: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/ACC_GEX.h5mu"
  # output for motif enrichment results .hdf5
  dem_result_fname: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/dem_results.hdf5"
  ctx_result_fname: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/ctx_results.hdf5"
  # output html for motif enrichment results .html
  output_fname_dem_html: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/dem_results.html"
  output_fname_ctx_html: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/ctx_results.html"
  # output for prepare_menr .h5ad
  cistromes_direct: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/cistromes_direct.h5ad"
  cistromes_extended: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/cistromes_extended.h5ad"
  # output tf names .txt
  tf_names: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/tf_names.txt"
  # output for download_genome_annotations .tsv
  genome_annotation: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/genome_annotation.tsv"
  chromsizes: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/chromsizes.tsv"
  # output for search_space .tsb
  search_space: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/search_space.tsv"
  # output tf_to_gene .tsv
  tf_to_gene_adjacencies: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/tf_to_gene_adj.tsv"
  # output region_to_gene .tsv
  region_to_gene_adjacencies: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/region_to_gene_adj.tsv"
  # output eGRN .tsv
  eRegulons_direct: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/eRegulon_direct.tsv"
  eRegulons_extended: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/eRegulons_extended.tsv"
  # output AUCell .h5mu
  AUCell_direct: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/AUCell_direct.h5mu"
  AUCell_extended: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/AUCell_extended.h5mu"
  # output scplus mudata .h5mu
  scplus_mdata: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E57_multiomics_September2023_Hg38/count/E57_IL_multi/outs/filtered_feature_bc_matrix/outs/outs/scplusmdata.h5mu"




#With the Scenic+ output data
#!/usr/bin/env python
import os
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os 
import pandas as pd
import pickle
import pandas as pd
import scanpy as sc


import mudata
scplus_mdata = mudata.read("scplusmdata.h5mu")

scplus_mdata.uns["direct_e_regulon_metadata"]

scplus_mdata.uns["extended_e_regulon_metadata"]
from scenicplus.RSS import(regulon_specificity_scores, plot_rss)
rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scATAC_counts:celltype",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC"])





eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)

eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]
sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")

sc.tl.umap(eRegulon_gene_AUC)

sc.pl.umap(eRegulon_gene_AUC, color = "scATAC_counts:celltype", save="UMAP_ATAC.pdf")


from scenicplus.RSS import(regulon_specificity_scores, plot_rss)
rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scATAC_counts:celltype",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC"]
)

plot_rss(
    data_matrix = rss,
    top_n = 10,
    num_columns = 5, save="plotrss.pdf"
)


#plotting the eregulons and heatmap
from scenicplus.plotting.dotplot import heatmap_dotplot

heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scATAC_counts:celltype",
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "horizontal",
    figsize = (16, 5), save="heatmap_eregulon.pdf"
)

#Cytoscape analysis
eRegulons_direct=pd.read_table("eRegulon_direct.tsv")
eRegulons_direct.sort_values('triplet_rank').head(200).loc[:,['TF','Gene','Region']].to_csv('eRegulon_direct_topTriplet.tsv', sep='\t')




                                                                                                                                                                                                                        
