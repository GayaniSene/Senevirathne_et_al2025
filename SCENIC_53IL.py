#!/bin/bash
#SBATCH -N 1
#SBATCH -p shared # partition (queue)
#SBATCH -c 10 # number of cores
#SBATCH --mem 60000 # memory pool for all cores
#SBATCH -t 2-8:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=msenevirathne@fas.harvard.edu
python s1.py



#running snakemake
#!/bin/bash
#SBATCH -N 1
#SBATCH -p sapphire # partition (queue)
#SBATCH -c 10 # number of cores
#SBATCH --mem 250000 # memory pool for all cores
#SBATCH -t 2-8:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=msenevirathne@fas.harvard.edu
snakemake --cores 10



python s1.py


Scenic+ input files 
#Can get it from scanpy 

or from Seurat

#cell_data table with annotations and barcodes for the input file for SCENIC+ 
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

#shell commands
ssh msenevirathne@login.rc.fas.harvard.edu
GandRo3147*

source activate scenicplus2 #this can be used for all the pycistopic 

#set the working directory
cd /n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/Multiomics_E53_29271_August_Replicate01_hg38/count/E53_IL/outs/outs/scplus_pipeline/Snakemake


cd /n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/Multiomics_E53_August_Replicate01_hg38/count/E53_IL/outs/outs/scplus_pipeline/Snakemake


cd /n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/E67_E72_2024_April_Hg38/RNA-seq/count/E67_IL/outs

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


#generating pseudobulk profiles for ATAC-seq
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
with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv"), "wt") as f:
    for v in bw_paths:
        _ = f.write(f"{v}\t{bw_paths[v]}\n")
        
        
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")
        



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
    n_cpu = 1,
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05,
    _temp_dir = '/mscstemp'
)


#before removing the blacklisted genomic regions, make sure to copy paste your pycisTopic software downloaded folder to your working directory. 


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




----------------------------------------------------------------
#no need to run .py. These are shell commands
#Next step is to perform the QC for ATAC-seq 
#Perform these steps in a shell environment

pycistopic tss gene_annotation_list | grep Human

#Check if the "qc" folder is made under the "outs" folder
mkdir -p outs/qc

pycistopic tss get_tss \
    --output outs/qc/tss.bed \
    --name "hsapiens_gene_ensembl" \
    --to-chrom-source ucsc \
    --ucsc hg38

head outs/qc/tss.bed | column -t

pycistopic qc \
    --fragments data/atac_fragments.tsv.gz \
    --regions outs/consensus_peak_calling/consensus_regions.bed \
    --tss outs/qc/tss.bed \
    --output outs/qc/SeuratProject


#After the SHELL commands, run these in a python environment (no need for a .py command line)

from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
import matplotlib.pyplot as plt
fragments_dict = {
    "SeuratProject": "data/atac_fragments.tsv.gz"
}

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

---------------------------------



#NEW code 19th May 2024 for generating the CisTopic object: 


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

#Adding meta data for the object
import pandas as pd
cell_data = pd.read_table("data/cell_data.tsv", index_col = 0)
cell_data.head()
cistopic_obj.add_cell_data(cell_data, split_pattern='-')
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)

cistopic_obj.cell_data


------------------------------------------------------------------------------------
------------------------------------------------------------------

#Running models using Mallet
#Shell command to downlad Mallet

wget https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz
!tar -xf Mallet-202108-bin.tar.gz

mkdir -p /mallet/tutorial/


#run these on python

#!/usr/bin/env python
from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
import matplotlib.pyplot as plt

out_dir="outs"
import os
fragments_dict = {
    "SeuratProject": "data/atac_fragments.tsv.gz"
}   



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

cistopic_obj.add_LDA_model(model)

pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir,"cistopic_obj.pkl"), "wb")
)




#### 
###########################################################

#read the cistopicObj once the model data is added. can run it on the laptop
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

plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['celltype', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],
    target='cell', num_columns=4,
    text_size=10, save="UMAP_annoted_cell_Clusters.pdf",
    dot_size=5)


#plotting contnuos data
plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['log10_unique_fragments_count', 'tss_enrichment', 'fraction_of_fragments_in_peaks'],
    target='cell', num_columns=4, save="UMAP_continuous_data.pdf",
    text_size=10,
    dot_size=5)

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

import pickle
import os
cistopic_obj = pickle.load (open ("cistopic_obj.pkl", 'rb'))

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


# Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1
plt.subplots_adjust(wspace=0, hspace=-0.70)
plt.show()



topic_annot = topic_annotation(
    cistopic_obj,
    annot_var='celltype',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)




#Finding DARs/differentially accessible regions

#!/usr/bin/env python

import pickle
import os
cistopic_obj = pickle.load (open ("cistopic_obj.pkl", 'rb'))

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
    scale_factor=10**6
)

normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True, save="normalized_dispersion_of_features.pdf",
)

len(variable_regions)

markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='celltype',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=5,
    _temp_dir="/tmp",
    split_pattern = '-'
)

from pycisTopic.clust_vis import plot_imputed_features


print("Number of DARs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(markers_dict[x])}")




##### saving the gene sets

os.makedirs(os.path.join("region_sets"), exist_ok = True)
os.makedirs(os.path.join("region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join("region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join("region_sets", "DARs_cell_type"), exist_ok = True)


from pycisTopic.utils import region_names_to_coordinates


for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join("region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

  for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join("region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join("region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False



#edit the .yaml

input_data:
  cisTopic_obj_fname: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/Multiomics_E53_August_Replicate01_hg38/count/E53_IL/outs/outs/cistopic_obj.pkl"
  GEX_anndata_fname: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/Multiomics_E53_August_Replicate01_hg38/count/E53_IL/outs/outs/adata.h5ad"
  region_set_folder: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/Multiomics_E53_August_Replicate01_hg38/count/E53_IL/outs/outs/region_sets"
  ctx_db_fname: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/Multiomics_E53_August_Replicate01_hg38/count/E53_IL/outs/outs/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
  dem_db_fname: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/Multiomics_E53_August_Replicate01_hg38/count/E53_IL/outs/outs/hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
  path_to_motif_annotations: "/n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/Multiomics_E53_August_Replicate01_hg38/count/E53_IL/outs/outs/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"




-------------------------------------

source activate scenicplus2


#shell commands
ssh msenevirathne@login.rc.fas.harvard.edu
GandRo3147*

source activate scenicplus2 #this can be used for all the pycistopic 

#set the working directory

cd /n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/Multiomics_E53_29271_August_Replicate01_hg38/count/E53_IL/outs/outs/scenic/scplus_pipeline/Snakemake



----------python file s2.py -----------------

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



eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)


eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]
sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")

sc.tl.umap(eRegulon_gene_AUC)

sc.pl.umap(eRegulon_gene_AUC, color = "scATAC_counts:celltype", save="UMAP_ATAC_new.pdf")


from scenicplus.RSS import (regulon_specificity_scores, plot_rss)
rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scATAC_counts:celltype",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC"]
)


plot_rss(
    data_matrix = rss,
    top_n = 10,
    num_columns = 5, figsize = (5, 10), save="plot_eregulon_top10_new.pdf"
)


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
    scale_size_matrix = True,
        scale_color_matrix = True,figsize = (20,20),
        orientation = 'vertical',limitsize=False,
    save="heatmap_eregulon_new.pdf"
)











------------end of s2.py------------------------------------------



#input for the Cytoscape



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



eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)


eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]
sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")

sc.tl.umap(eRegulon_gene_AUC)

sc.pl.umap(eRegulon_gene_AUC, color = "scATAC_counts:celltype", save="UMAP_ATAC_new.pdf")


from scenicplus.RSS import (regulon_specificity_scores, plot_rss)







------------- beginning of s4.py------------------



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
import os
import scanpy as sc
import anndata
import matplotlib
import matplotlib.pyplot as plt
import adjustText
import numpy as np
import pandas as pd




from scenicplus.simulation import (
    train_gene_expression_models,
    simulate_perturbation,
    plot_perturbation_effect_in_embedding)

%matplotlib inline

scplus_mdata = mudata.read("scplusmdata.h5mu")

eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)
eRegulon_gene_AUC.obs = scplus_mdata.obs


sc.pp.pca(eRegulon_gene_AUC)


color_dict_line = {
    'ChondroProg': '#9A031E',
    'EndothelialCells': '#C75146',
    'Fibro1': '#FFA987',
    'Fibro2': '#222E50',
    'Fibro3': '#8BB174',
    'Mes1': '#2A4849',
    'Mes2': '#3E5641',
    'MyoProg': '#59A96A',
    'Perichondro+Osteoblasts': '#56E39F'}



def plot_mm_line_pca(ax):
    texts = []
    # Plot PCA
    ax.scatter(
        eRegulon_gene_AUC.obsm["X_pca"][:, 0],
        eRegulon_gene_AUC.obsm["X_pca"][:, 1],
        color = [color_dict_line[line] for line in eRegulon_gene_AUC.obs["scATAC_counts:MMline"]]
    )

    # Plot labels
    for line in set(eRegulon_gene_AUC.obs["scATAC_counts:MMline"]):
        line_bc_idc = np.arange(len(eRegulon_gene_AUC.obs_names))[eRegulon_gene_AUC.obs["scATAC_counts:MMline"] == line]
        avg_x, avg_y = eRegulon_gene_AUC.obsm["X_pca"][line_bc_idc, 0:2].mean(0)
        texts.append(
            ax.text(
                avg_x,
                avg_y,
                line,
                fontweight = "bold"
            )
        )
    adjustText.adjust_text(texts)

fig, ax = plt.subplots()
plot_mm_line_pca(ax)



--------------------
-------------------
-----------------------
-----------------------

######### How to compare with the Hi-c file (Original file: GSE200345_C28_HiC_megaMap_SIP_20kb_loops.bedpe;Diekman’s paper 
### file name "chondrocyte_hic_interaction_data_anchors_HITTING_PROMOTERS_RES_WITH_SCORES.txt" V1:V3 being the left anchor, ‘chr’, ‘start’ and ‘stop’ referring to the right anchor. V6-V10 refers to the gene promoter I intersected with the anchor)

#shell commands
ssh msenevirathne@login.rc.fas.harvard.edu
GandRo3147*

source activate scenicplus2 #this can be used for all the pycistopic 

#set the working directory
cd /n/boslfs02/LABS/capellini_lab/users/msenevirathne/Multiomics/Multiomics_E53_29271_August_Replicate01_hg38/count/E53_IL/outs/outs/scenic/scplus_pipeline/Snakemake



#HiC-links


import pandas as pd
import pandas as pd

#Extract the Promoter Anchor: I also copy pasted the promoter columns 6-9 and saved as a .bed file
cut -f6-9 chondrocyte_hic_interaction.txt > hic_promoters.bed

#Extract the enhancer 
#Extract the Distal Anchor (Candidate Enhancer):
#Choose the non-promoter anchor—in many cases this is the right anchor (columns 11–13). Create a BED file for these regions:
#Identifying Candidate Enhancers
    #•   Distal Anchors: In the context of Hi-C data, a distal anchor refers to a region that is not immediately adjacent to a promoter 
      #but is linked via Hi-C interaction. These regions can potentially serve as enhancers.
    #•   ATAC-seq Overlap: By overlapping distal Hi-C anchors with ATAC-seq data, which highlights regions of open chromatin, 
    #researchers can identify which of these distal anchors are likely to be functional enhancers. 
    #Open chromatin in distal anchors suggests that these regions are accessible to transcription factors and other regulatory proteins, which is a hallmark of enhancer activity.

cut -f11-13 chondrocyte_hic_interaction_data_anchors_HITTING_PROMOTERS_RES_WITH_SCORES.txt > hic_distal.bed


#ATAC‑seq .bed files mark regions of open chromatin, which are candidate enhancers. 
#To test if the distal Hi‑C anchor is in an accessible region, use BEDTools to find overlaps between your ATAC‑seq peaks and the distal anchor BED file.

bedtools intersect -a hic_distal.bed -b RestingChondro2_HARs_hg19.bed -wa -wb > HIC_RestingChondro2_HARs_hg19.txt

bedtools intersect -a hic_distal.bed -b RestingChondro_HARs_hg19.bed -wa -wb > HIC_RestingChondro_HARs_hg19.txt

bedtools intersect -a hic_distal.bed -b Perichondro+Osteoblasts_HARs_hg19.bed -wa -wb > HIC_Perichondro+Osteoblasts_HARs_hg19.txt

bedtools intersect -a hic_distal.bed -b Perichondro+Osteoblasts_unipeaks_hg19.bed -wa -wb > HIC_Perichondro+Osteoblasts_unipeaks_hg19.txt







bedtools intersect -a hic_distal.bed -b E53RestingChondro2_peaks.narrowPeak -wa -wb > E53HIC_RestingChondro2_distal_ATAC_overlap.txt

bedtools intersect -a hic_distal.bed -b E53ChondroProg_commonpeaks.bed -wa -wb > E53HIC_ChondroProg_distal_HARs.txt

bedtools intersect -a hic_distal.bed -b E53RestingChondro2_commonpeaks.bed -wa -wb > E53HIC_RestingChondro2_distal_HARs.txt

bedtools intersect -a hic_distal.bed -b E53RestingChondro_commonpeaks.bed -wa -wb > E53HIC_RestingChondro_distal_HARs.txt


#Next, I used the HiC overlapped files, and imported them to Juicebox, to visualize the loops. And mapped the cell-type peak files and HAR-overllaped ones on to the HiC file that I downloaded. 


import pandas as pd


hic_cols = [
    "L_chr", "L_start", "L_end", "L_side", "L_ID",       # Left anchor
    "P_chr", "P_start", "P_end", "Gene", "Window",         # Promoter anchor (with gene)
    "R_chr", "R_start", "R_end", "R_side", "R_ID", "Score", "Orig"  # Right anchor (candidate enhancer)
]


#load the Hi-C file
hic_df = pd.read_csv("chondrocyte_hic_interaction_data_anchors_HITTING_PROMOTERS_RES_WITH_SCORES.txt",
                     sep="\t", names=hic_cols)


# Define column names for the distal ATAC-seq overlap file.
# Here we assume that the overlap file was created using:
# bedtools intersect -a hic_distal.bed -b ATAC_peaks.bed -wa -wb > distal_ATAC_overlap.txt
# And that hic_distal.bed has three columns: R_chr, R_start, R_end.
# Followed by, say, four columns from the ATAC file: ATAC_chr, ATAC_start, ATAC_end, ATAC_info.
overlap_cols = ["R_chr", "R_start", "R_end", "ATAC_chr", "ATAC_start", "ATAC_end", "ATAC_info"]


# Load the overlap file
overlap_df = pd.read_csv("distal_ATAC_overlap_PerichondroHARs.txt", sep="\t", names=overlap_cols)



# Merge the Hi-C data with the ATAC overlap file on the distal anchor coordinates
merged_df = pd.merge(hic_df, overlap_df, on=["R_chr", "R_start", "R_end"])













# Load the SCENIC+ output file (adjust the file path as necessary)
scenic_file = "eRegulon_direct.tsv"
df = pd.read_csv(scenic_file, sep="\t")

# Inspect column names to identify those for TFs and targets
print(df.columns)  # Look for columns like 'TF' and 'Target'

# Extract the TF-target interactions and remove duplicates if needed
tf_target_df = df[['TF', 'Gene']].drop_duplicates()

# Save as an edge list suitable for network analysis or integration with other data (e.g., Hi-C)
tf_target_df.to_csv("scenic_gene_links.txt", sep="\t", index=False, header=False)



import pandas as pd

# Load SCENIC+ TF-target file (adjust the path and delimiter as needed)
scenic_df = pd.read_csv("scenic_gene_links.txt", sep="\t", names=["TF", "Target"])

# Define column names for the Hi-C file based on its structure
hic_cols = ["L_chr", "L_start", "L_end", "L_side", "L_ID",
            "P_chr", "P_start", "P_end", "Gene", "Window",
            "R_chr", "R_start", "R_end", "R_side", "R_ID", "Score", "Orig"]

# Load the Hi-C data file (adjust path and delimiter as needed)
hic_df = pd.read_csv("chondrocyte_hic_interaction.txt",
                     sep="\t", names=hic_cols)

# Merge SCENIC+ links with Hi-C data based on the target gene matching the promoter hit
integrated_df = pd.merge(scenic_df, hic_df, left_on="Target", right_on="Gene", how="inner")

# Optionally, you can filter by the Hi-C interaction score to retain high-confidence links:
score_threshold = 0.02  # adjust threshold as appropriate
integrated_df = integrated_df[integrated_df["Score"] >= score_threshold]

# Save the integrated output for downstream analysis or visualization
integrated_df.to_csv("integrated_TF_target_HiC.txt", sep="\t", index=False)



#REMAIN THE UNIQUE INTERSECTIONS

# Load the integrated file (adjust delimiter if needed)
df = pd.read_csv("integrated_TF_target_HiC.txt", sep="\t")

# Remove duplicates (if you want to consider all columns)
df_unique = df.drop_duplicates()

# Alternatively, to remove duplicates based on just the TF and Target columns:
# df_unique = df.drop_duplicates(subset=["TF", "Target"])

# Save the result to a new file
df_unique.to_csv("integrated_TF_target_HiC_unique.txt", sep="\t", index=False)






