"""GRN analysis script.

Consolidated imports placed at top; duplicate imports removed.
"""

# %% Consolidated imports
# Standard library
import sys
import os
import re
import json
import math
import glob
import itertools
import collections
import pickle
import warnings

# Third-party
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.patches as mpatches
import adjustText as adjust_text
import celloracle as co
import scipy
import scipy.stats  # noqa: F401
from scipy.sparse import csc_matrix
from scipy.stats import zscore
# import hotspot as hs
import plotnine as pn
import gseapy as gp
from gseapy.plot import barplot, dotplot
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from statannotations.Annotator import Annotator
from collections import Counter
from matplotlib_venn import venn2
from statannotations.stats.StatTest import StatTest
from collections import Counter
import networkx as nx
import math

# Optional dependency
try:
    import mygene  # type: ignore
    from mygene import MyGeneInfo  # type: ignore
    _MG_AVAILABLE = True
except ImportError:  # pragma: no cover
    mygene = None  # type: ignore
    MyGeneInfo = None  # type: ignore
    _MG_AVAILABLE = False
    print("mygene not installed: pip install mygene (falling back to empty annotations)")

warnings.filterwarnings('ignore')

# Custom imports
sys.path.append('/storage/sahuarea/samuele_cancellieri/code_for_projects/import')
import celloracle_utils as co_utils
# %%
co_utils.myCellOracle.test_import()

# %%
co.__version__,sc.__version__,sns.__version__

# %%
sns.reset_defaults()
# plt.rcdefaults()
# sns.set_theme(style="whitegrid", context="paper", font_scale=1.5, rc={"figure.figsize":(7,7)}, palette="deep")
# sc.set_figure_params(dpi=100, dpi_save=300, format='pdf', frameon=False, fontsize=12,color_map="inferno", vector_friendly=True)
mpl.rcParams["pdf.fonttype"] = 42

%matplotlib inline
# %% [markdown]
# ## SAVE FOLDERS

# %%
# save folders
filename="regev"

save_path_matrices=f"/storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/{filename}/matrices/"
save_path_metadata=f"/storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/{filename}/metadata/"
save_path_figures=f"/storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/{filename}/figures/"
save_path_networks=f"/storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/{filename}/networks/"

TG_to_TF_dictionary="/storage/sahuarea/samuele_cancellieri/reference_data/TFs/TG_to_TF_dictionary.pkl"

save_folders={'save_path_matrices':save_path_matrices, 'save_path_metadata':save_path_metadata, 'save_path_figures':save_path_figures, 'save_path_networks':save_path_networks}

# import os
# for folder in save_folders.values():
#     os.makedirs(folder, exist_ok=True)

# %% [markdown]
# # OPEN DATA
tf_list=open("/storage/sahuarea/samuele_cancellieri/reference_data/TFs/TF_names_v_1.01.txt", 'r').readlines()
tf_list=[i.strip() for i in tf_list]
tf_list[:10]
# %%
# adata = sc.read_10x_mtx(save_path_matrices+"filtered_gene_bc_matrices_h5.h5ad")

# %%
adata=sc.read_h5ad(f"{save_path_matrices}/analyzed/regev_caf-epithelial(malignant-non_malignant).h5ad")
adata
# %%
adata.obs.groupby('pid')['response'].unique().value_counts()
#%%
order=['Untreated','Poor response','Minimal response','Moderate response']
adata.obs['response']=adata.obs['response'].cat.reorder_categories(order,ordered=True)
# %%
adata.obs['level_1_annotation']=adata.obs['Level 1 Annotation'].astype(str)
adata.obs['level_2_annotation']=adata.obs['Level 2 Annotation'].astype(str)
adata.obs['level_3_annotation']=adata.obs['Level 3 Annotation'].astype(str)
adata.obs['new_celltypes']=adata.obs['new_celltypes'].astype(str)
adata.obs_keys()
# %%
cell_type_annotation='Epithelial (malignant)'
adata_malignant_original=adata[adata.obs['level_1_annotation']==cell_type_annotation].copy()
adata_malignant_original
# %%
sc.pl.umap(adata, color=['level_1_annotation'], title="Cell types",save="regev_level_1_annotation.pdf")
sc.pl.umap(adata, color=['level_2_annotation'], title="Cell types",save="regev_level_2_annotation.pdf")
sc.pl.umap(adata, color=['level_3_annotation'], title="Cell types",save="regev_level_3_annotation.pdf")

# %%
sc.pl.umap(adata, color=['NFATC2'],cmap='inferno')

# %%
sc.tl.rank_genes_groups(adata, 'new_celltypes', method='wilcoxon')

# %% [markdown]
# # MERGE TREATMENT AND RESPONSE INTO BINARY GROUPS
# %%
adata.obs['merged_treatment']=adata.obs['new_treatment'].apply(lambda x: "Untreated" if x=="Untreated" else "Treated")
adata.obs['merged_treatment'].unique()
# %%
adata.obs['merged_response']=adata.obs['response'].apply(lambda x: "U/PR" if ['Untreated','Poor response'].count(x) else "MM/MR")
adata.obs['merged_response'].unique()
# %%
print(adata.obs['merged_response'].value_counts())
print(adata.obs['response'].value_counts())

print(adata.obs['merged_treatment'].value_counts())
print(adata.obs['new_treatment'].value_counts())

# %%
adata.write_h5ad(f"{save_path_matrices}/regev_original_with_merged_metadata.h5ad")

# %%
for annot in adata.obs['level_1_annotation'].unique():
    print(annot)
    adata_subset=adata[adata.obs['level_1_annotation']==annot].copy()
    adata_subset.write_h5ad(f"{save_path_matrices}/regev_{annot}_merged.h5ad")
    

# %% [markdown]
# ### SPLIT DATASET & SUBSET CELLTYPES

# %%
cell_types_extract=["Epithelial (malignant)", "Cancer-associated fibroblast", "Epithelial (non-malignant)"]

adata_filtered=adata[adata.obs['level_1_annotation'].isin(cell_types_extract)].copy()
sc.pl.umap(adata_filtered,color=['level_2_annotation','merged_treatment','merged_response'], wspace=0.3,save="regev_filtered.pdf")

# %%
adata_response_treat=adata_filtered[adata_filtered.obs['merged_response']=="MM/MR"].copy()
adata_response_treat

# %%
adata_response_untreat=adata_filtered[adata_filtered.obs['merged_response']=="U/PR"].copy()
adata_response_untreat

# %% [markdown]
# # CLEAN DATA AND PROCESS TREATMENT BASE NETWORK

# %%
# adata_filtered=adata[adata.obs['level_1_annotation']=="Myeloid"].copy()
cell_type_annotation='Epithelial (malignant)'
adata_filtered=adata[adata.obs['level_1_annotation']==cell_type_annotation].copy()
adata_filtered
# %%
adata_sampled = adata_filtered.obs.groupby("response", group_keys=False).apply(
    lambda x: x.sample(n=1000, random_state=42) if len(x) >= 1000 else x
)
adata_sampled = adata_filtered[adata_sampled.index].copy()
# %%
sc.tl.rank_genes_groups(adata_sampled, groupby='merged_response', method='wilcoxon', use_raw=False)
# %%
df_rank = sc.get.rank_genes_groups_df(adata_sampled, group=None, pval_cutoff=0.05)
# %%
df_rank_tf_filtered = df_rank[df_rank['names'].isin(tf_list)]
df_rank_tf_filtered
# %%
sc.set_figure_params(dpi=120, dpi_save=200, format='pdf', frameon=True, fontsize=12, vector_friendly=True)
sc.pl.umap(adata_filtered, color=['level_2_annotation','merged_treatment','merged_response'],wspace=0.4)

# %%
uti.myScanpy.quality_control(adata=adata)

# %%
celldownsample=30000 if adata.n_obs>=30000 else adata.n_obs
co_utils.myCellOracle.necessary_preprocessing(adata, 
                                              cluster_column_name="merged_response", 
                                              cell_downsample=celldownsample,
                                              top_genes=2000)

# %%
adata_filtered

# %%
sc.pl.umap(adata,color=['level_2_annotation'], ncols=1)

# %%
sc.pl.draw_graph(adata_filtered, color=['level_1_annotation','level_2_annotation'], ncols=1)

# %% [markdown]
# ## CREATE NETWORKS

# %%
adata_filtered.obs.value_counts('level_2_annotation')

# %%
adata_filtered.layers['counts'].toarray()[1:25,1:25]

# %%
sc.pp.filter_cells(adata_filtered, min_genes=200)
sc.pp.filter_genes(adata_filtered, min_cells=3)

# %%
adata_filtered.X.toarray()[1:25,1:25]

# %%
co_utils.myCellOracle.necessary_preprocessing(adata_filtered, cluster_column_name="response", cell_downsample=20000, top_genes=2000)

# %%
oracle=co_utils.myCellOracle.create_oracle_object(adata=adata, 
                                                  TG_to_TF_dictionary=TG_to_TF_dictionary, 
                                                  cluster_column_name="merged_response", 
                                                  embedding_name="X_draw_graph_fa",
                                                  raw_count_layer="counts")

# %%
co_utils.myCellOracle.run_PCA(oracle)

# %%
co_utils.myCellOracle.run_KNN(oracle)

# %%
# Save oracle object.
# oracle.to_hdf5(f"{save_path_networks}/regev_caf-epithelial(malignant-non_malignant)_reclustered.celloracle.oracle")
oracle.to_hdf5(f"{save_path_networks}/regev_malignant_merged-response.celloracle.oracle")

# %% [markdown]
# ## PROCESS NETWORKS

# %%
# oracle=co.load_hdf5(f"{save_path_networks}/regev_endothelial.celloracle.oracle")
# oracle

# %%
links=co_utils.myCellOracle.run_links(oracle, cluster_column_name="merged_response")

# %%
# links.filter_links(p=0.001, weight="coef_abs")

# %%
# links.to_hdf5(file_path=f"{save_path_networks}/regev_endothelial.celloracle.links")

# %%
links.to_hdf5(file_path=f"{save_path_networks}/regev_malignant_merged-response.celloracle.links")

# %%
# links.plot_degree_distributions(plot_model=True)

# %%
# links.get_network_score()
# %%
# filename_treat="regev_epithelial-CAF_treated"
# filename_untreat="regev_epithelial-CAF_untreated"

filename_poor="regev_caf-epithelial_bad_response"
filename_good="regev_caf-epithelial_good_response"

oracle_good = co.load_hdf5(f"{save_path_networks}/{filename_good}.celloracle.oracle")
links_good = co.load_hdf5(f"{save_path_networks}/{filename_good}.celloracle.links")
oracle_poor = co.load_hdf5(f"{save_path_networks}/{filename_poor}.celloracle.oracle")
links_poor = co.load_hdf5(f"{save_path_networks}/{filename_poor}.celloracle.links")
#%%
links_treat.cluster
# %%
links_untreat.cluster
# %%
df_treat=links_treat.merged_score.sort_values("eigenvector_centrality", ascending=False).reset_index()
df_untreat=links_untreat.merged_score.sort_values("eigenvector_centrality", ascending=False).reset_index()
#%%
# Filter df to only include genes in tf_list
df_filtered = df_treat[df_treat['index'].isin(tf_list)]
df_filtered=df_filtered[df_filtered['index'].isin(df_rank_tf_filtered['names'].unique().tolist())]

# Group by cluster and get top 5 by eigenvector_centrality for each cluster
top5_by_cluster = (
    df_filtered
    .sort_values(['cluster', 'eigenvector_centrality'], ascending=[True, False])
    .groupby('cluster')
    .head(25)
)

top5_by_cluster
#%%
# Pivot the top5_by_cluster DataFrame to create a matrix for heatmap
heatmap_data = top5_by_cluster.pivot(index='index', columns='cluster', values='eigenvector_centrality')

# Optionally, fill NaN with 0 or another value
heatmap_data = heatmap_data.fillna(min(heatmap_data.min().min(), 0))

# plt.close('all')
# fig,ax=plt.subplots(figsize=(10, 10))
sns.clustermap(heatmap_data, figsize=(15,20),square=False,cmap='inferno', annot=False, fmt=".2f",row_cluster=True,col_cluster=True, linewidths=0.1)
plt.suptitle("Top TFs by Eigenvector Centrality per Cluster - Good response", fontsize=16)
plt.ylabel("Gene (TF)")
plt.xlabel("Cluster")
# plt.xticks(rotation=90)
plt.tight_layout()
# plt.savefig(f"{save_path_figures}/top25_TFs_eigenvector_centrality_per_cluster_good.pdf", bbox_inches='tight', dpi=300)
plt.show()
#%%
# Filter df to only include genes in tf_list
df_filtered = df_untreat[df_untreat['index'].isin(tf_list)]
df_filtered=df_filtered[df_filtered['index'].isin(df_rank_tf_filtered['names'].unique().tolist())]

# Group by cluster and get top 5 by eigenvector_centrality for each cluster
top5_by_cluster = (
    df_filtered
    .sort_values(['cluster', 'eigenvector_centrality'], ascending=[True, False])
    .groupby('cluster')
    .head(25)
)

top5_by_cluster
#%%
# Pivot the top5_by_cluster DataFrame to create a matrix for heatmap
heatmap_data = top5_by_cluster.pivot(index='index', columns='cluster', values='eigenvector_centrality')

# Optionally, fill NaN with 0 or another value
heatmap_data = heatmap_data.fillna(min(heatmap_data.min().min(), 0))

# plt.close('all')
# fig,ax=plt.subplots(figsize=(10, 10))
sns.clustermap(heatmap_data, figsize=(15,20),square=False,cmap='inferno', annot=False, fmt=".2f",row_cluster=True,col_cluster=True, linewidths=0.1)
plt.suptitle("Top TFs by Eigenvector Centrality per Cluster - Bad response", fontsize=16)
plt.ylabel("Gene (TF)")
plt.xlabel("Cluster")
# plt.xticks(rotation=90)
plt.tight_layout()
# plt.savefig(f"{save_path_figures}/top25_TFs_eigenvector_centrality_per_cluster_bad.pdf", bbox_inches='tight', dpi=300)
plt.show()
#%%
# Compare top TFs between treated and untreated: plot overlap, unique, and differences
df_treat=links_treat.merged_score.sort_values("eigenvector_centrality", ascending=False).reset_index()
df_untreat=links_untreat.merged_score.sort_values("eigenvector_centrality", ascending=False).reset_index()

# Create TF comparison folder and ensure it exists
tf_comparison_folder = f"{save_path_figures}/tf_comparison_treatment/"
os.makedirs(tf_comparison_folder, exist_ok=True)

N = 10000 ##take all the TF
clusters = sorted(set(df_treat['cluster']).intersection(df_untreat['cluster']))

for cluster in clusters:
    top_treat = (
        df_treat[df_treat['cluster'] == cluster]
        .sort_values('eigenvector_centrality', ascending=False)
        .loc[lambda df: df['index'].isin(tf_list)]
        .head(N)
        .set_index('index')
    )
    top_untreat = (
        df_untreat[df_untreat['cluster'] == cluster]
        .sort_values('eigenvector_centrality', ascending=False)
        .loc[lambda df: df['index'].isin(tf_list)]
        .head(N)
        .set_index('index')
    )

    # Find overlap and unique TFs
    tf_treat = set(top_treat.index)
    tf_untreat = set(top_untreat.index)
    tf_overlap = tf_treat & tf_untreat
    tf_unique_treat = tf_treat - tf_untreat
    tf_unique_untreat = tf_untreat - tf_treat

    # Prepare DataFrame for plotting differences
    common_tfs = list(tf_overlap)
    if common_tfs:
        diff_df = pd.DataFrame({
            'Treated': top_treat.loc[common_tfs, 'eigenvector_centrality'],
            'Untreated': top_untreat.loc[common_tfs, 'eigenvector_centrality'],
        })
        diff_df['Difference'] = diff_df['Treated'] - diff_df['Untreated']
        diff_df = diff_df.sort_values('Difference', ascending=False)
        diff_df_tmp_head=diff_df.head(10)
        diff_df_tmp_tail=diff_df.tail(10)
        diff_df=pd.concat([diff_df_tmp_head,diff_df_tmp_tail])

        # plt.figure(figsize=(7, 7))
        fig,ax=plt.subplots(figsize=(7, 5))
        diff_df['Difference'].plot(kind='bar', color='purple',ax=ax)
        plt.title(f"Difference in Eigenvector Centrality (Treated - Untreated)\nTop TFs in {cluster}")
        plt.ylabel("Difference")
        plt.xlabel("TF")
        sns.despine()
        ax.grid(False)
        plt.tight_layout()
        plt.savefig(f"{tf_comparison_folder}/top_TFs_difference_treated-untreated_{cluster}.pdf", bbox_inches='tight', dpi=300)
        plt.show()

    # Plot unique TFs for treated
    if tf_unique_treat:
        # plt.figure(figsize=(5, 5))
        fig,ax=plt.subplots(figsize=(5, 5))
        vals = top_treat.loc[list(tf_unique_treat), 'eigenvector_centrality'].sort_values(ascending=True)
        vals.plot(kind='barh', color="#0096FF",ax=ax)
        plt.title(f"Unique Top TFs in Treated ({cluster})")
        plt.xlabel("Eigenvector Centrality")
        plt.ylabel("TF")
        sns.despine()
        ax.grid(False)
        plt.tight_layout()
        plt.savefig(f"{tf_comparison_folder}/unique_TFs_treated_{cluster}.pdf", bbox_inches='tight', dpi=300)
        plt.show()

    # Plot unique TFs for untreated
    if tf_unique_untreat:
        # plt.figure(figsize=(5, 5))
        fig,ax=plt.subplots(figsize=(5, 5))
        vals = top_untreat.loc[list(tf_unique_untreat), 'eigenvector_centrality'].sort_values(ascending=True)
        vals.plot(kind='barh', color="#fb3310",ax=ax)
        plt.title(f"Unique Top TFs in Untreated ({cluster})")
        plt.xlabel("Eigenvector Centrality")
        plt.ylabel("TF")
        sns.despine()
        ax.grid(False)
        plt.tight_layout()
        plt.savefig(f"{tf_comparison_folder}/unique_TFs_untreated_{cluster}.pdf", bbox_inches='tight', dpi=300)
        plt.show()
# %%
df_plot = links_treat.merged_score.copy()
df_plot.query("degree_out > 0", inplace=True)
df_plot = df_plot.sort_values(['cluster', 'eigenvector_centrality'], ascending=[True, False])
df_plot['rank'] = df_plot.groupby('cluster')['eigenvector_centrality'].rank(ascending=False, method='first')
top_n = 20
df_plot = df_plot[df_plot['rank'] <= top_n]
# df_plot.index.name = 'gene'
df_plot["gene"] = df_plot.index
df_plot.reset_index(drop=True, inplace=True)
# df_plot = df_plot[df_plot['cluster'].isin(['ADM', 'CAF', 'Malignant', 'Ductal'])]

for cluster in df_plot['cluster'].unique():
    fig,ax=plt.subplots(figsize=(6, 6))
    
    df_plot_cluster = df_plot[df_plot['cluster'] == cluster]
    df_plot_cluster = df_plot_cluster.sort_values('eigenvector_centrality', ascending=False)
    sns.scatterplot(data=df_plot_cluster, x='eigenvector_centrality', y='gene',edgecolor="black", color="#0096FF", s=50, ax=ax)
    plt.title(f"Top {top_n} Eigenvector Centrality Treated\n{cluster}")
    plt.xlim(min(df_plot_cluster['eigenvector_centrality']) - 0.05, max(df_plot_cluster['eigenvector_centrality']) + 0.05)
    sns.despine()
    ax.grid(False)
    plt.tight_layout()
    fig.savefig(f"{tf_comparison_folder}{filename_treat}_{cluster}_top{top_n}_eigenvector_centrality.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{tf_comparison_folder}{filename_treat}_{cluster}_top{top_n}_eigenvector_centrality.svg", dpi=300, bbox_inches='tight')
    plt.show()
# %%
df_plot = links_untreat.merged_score.copy()
df_plot.query("degree_out > 0", inplace=True)
df_plot = df_plot.sort_values(['cluster', 'eigenvector_centrality'], ascending=[True, False])
df_plot['rank'] = df_plot.groupby('cluster')['eigenvector_centrality'].rank(ascending=False, method='first')
top_n = 20
df_plot = df_plot[df_plot['rank'] <= top_n]
# df_plot.index.name = 'gene'
df_plot["gene"] = df_plot.index
df_plot.reset_index(drop=True, inplace=True)
# df_plot = df_plot[df_plot['cluster'].isin(['ADM', 'CAF', 'Malignant', 'Ductal'])]

for cluster in df_plot['cluster'].unique():
    fig,ax=plt.subplots(figsize=(6, 6))
    
    df_plot_cluster = df_plot[df_plot['cluster'] == cluster]
    df_plot_cluster = df_plot_cluster.sort_values('eigenvector_centrality', ascending=False)
    sns.scatterplot(data=df_plot_cluster, x='eigenvector_centrality', y='gene',edgecolor="black", color="#0096FF", s=50, ax=ax)
    plt.title(f"Top {top_n} Eigenvector Centrality Treated\n{cluster}")
    plt.xlim(min(df_plot_cluster['eigenvector_centrality']) - 0.05, max(df_plot_cluster['eigenvector_centrality']) + 0.05)
    sns.despine()
    ax.grid(False)
    plt.tight_layout()
    fig.savefig(f"{tf_comparison_folder}{filename_untreat}_{cluster}_top{top_n}_eigenvector_centrality.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{tf_comparison_folder}{filename_untreat}_{cluster}_top{top_n}_eigenvector_centrality.svg", dpi=300, bbox_inches='tight')
    plt.show()
# %%
# Compare top TFs between treated and untreated: plot overlap, unique, and differences
df_treat=links_good.merged_score.sort_values("eigenvector_centrality", ascending=False).reset_index()
df_untreat=links_poor.merged_score.sort_values("eigenvector_centrality", ascending=False).reset_index()

N = 10000 ##take all the TF
clusters = sorted(set(df_treat['cluster']).intersection(df_untreat['cluster']))

# Create TF comparison folder and ensure it exists
tf_comparison_folder = f"{save_path_figures}/tf_comparison_response/"
os.makedirs(tf_comparison_folder, exist_ok=True)

for cluster in clusters:
    top_treat = (
        df_treat[df_treat['cluster'] == cluster]
        .sort_values('eigenvector_centrality', ascending=False)
        .loc[lambda df: df['index'].isin(tf_list)]
        .head(N)
        .set_index('index')
    )
    top_untreat = (
        df_untreat[df_untreat['cluster'] == cluster]
        .sort_values('eigenvector_centrality', ascending=False)
        .loc[lambda df: df['index'].isin(tf_list)]
        .head(N)
        .set_index('index')
    )

    # Find overlap and unique TFs
    tf_treat = set(top_treat.index)
    tf_untreat = set(top_untreat.index)
    tf_overlap = tf_treat & tf_untreat
    tf_unique_treat = tf_treat - tf_untreat
    tf_unique_untreat = tf_untreat - tf_treat

    # Prepare DataFrame for plotting differences
    common_tfs = list(tf_overlap)
    if common_tfs:
        diff_df = pd.DataFrame({
            'Treated': top_treat.loc[common_tfs, 'eigenvector_centrality'],
            'Untreated': top_untreat.loc[common_tfs, 'eigenvector_centrality'],
        })
        diff_df['Difference'] = diff_df['Treated'] - diff_df['Untreated']
        diff_df = diff_df.sort_values('Difference', ascending=False)
        diff_df_tmp_head=diff_df.head(10)
        diff_df_tmp_tail=diff_df.tail(10)
        diff_df=pd.concat([diff_df_tmp_head,diff_df_tmp_tail])

        # plt.figure(figsize=(7, 7))
        fig,ax=plt.subplots(figsize=(7, 5))
        diff_df['Difference'].plot(kind='bar', color='purple',ax=ax)
        plt.title(f"Difference in Eigenvector Centrality (Good - Poor)\nTop TFs in {cluster}")
        plt.ylabel("Difference")
        plt.xlabel("TF")
        sns.despine()
        ax.grid(False)
        plt.tight_layout()
        plt.savefig(f"{tf_comparison_folder}/top_TFs_difference_good-poor_{cluster}.pdf", bbox_inches='tight', dpi=300)
        plt.show()

    # Plot unique TFs for treated
    if tf_unique_treat:
        # plt.figure(figsize=(5, 5))
        fig,ax=plt.subplots(figsize=(5, 5))
        vals = top_treat.loc[list(tf_unique_treat), 'eigenvector_centrality'].sort_values(ascending=True).head(10)
        vals.plot(kind='barh', color="#0096FF",ax=ax)
        plt.title(f"Unique Top TFs in Good response ({cluster})")
        plt.xlabel("Eigenvector Centrality")
        plt.ylabel("TF")
        sns.despine()
        ax.grid(False)
        plt.tight_layout()
        plt.savefig(f"{tf_comparison_folder}/unique_TFs_good_{cluster}.pdf", bbox_inches='tight', dpi=300)
        plt.show()

    # Plot unique TFs for untreated
    if tf_unique_untreat:
        # plt.figure(figsize=(5, 5))
        fig,ax=plt.subplots(figsize=(5, 5))
        vals = top_untreat.loc[list(tf_unique_untreat), 'eigenvector_centrality'].sort_values(ascending=True).head(10)
        vals.plot(kind='barh', color="#fb3310",ax=ax)
        plt.title(f"Unique Top TFs in Poor response ({cluster})")
        plt.xlabel("Eigenvector Centrality")
        plt.ylabel("TF")
        sns.despine()
        ax.grid(False)
        plt.tight_layout()
        plt.savefig(f"{tf_comparison_folder}/unique_TFs_poor_{cluster}.pdf", bbox_inches='tight', dpi=300)
        plt.show()
# %%
df_plot = links_good.merged_score.copy()
df_plot.query("degree_out > 0", inplace=True)
df_plot = df_plot.sort_values(['cluster', 'eigenvector_centrality'], ascending=[True, False])
df_plot['rank'] = df_plot.groupby('cluster')['eigenvector_centrality'].rank(ascending=False, method='first')
top_n = 20
df_plot = df_plot[df_plot['rank'] <= top_n]
# df_plot.index.name = 'gene'
df_plot["gene"] = df_plot.index
df_plot.reset_index(drop=True, inplace=True)
# df_plot = df_plot[df_plot['cluster'].isin(['ADM', 'CAF', 'Malignant', 'Ductal'])]

for cluster in df_plot['cluster'].unique():
    fig,ax=plt.subplots(figsize=(6, 6))
    
    df_plot_cluster = df_plot[df_plot['cluster'] == cluster]
    df_plot_cluster = df_plot_cluster.sort_values('eigenvector_centrality', ascending=False)
    sns.scatterplot(data=df_plot_cluster, x='eigenvector_centrality', y='gene',edgecolor="black", color="#0096FF", s=50, ax=ax)
    plt.title(f"Top {top_n} Eigenvector Centrality Good Response\n{cluster}")
    plt.xlim(min(df_plot_cluster['eigenvector_centrality']) - 0.05, max(df_plot_cluster['eigenvector_centrality']) + 0.05)
    sns.despine()
    ax.grid(False)
    plt.tight_layout()
    fig.savefig(f"{tf_comparison_folder}{filename_good}_{cluster}_top{top_n}_eigenvector_centrality.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{tf_comparison_folder}{filename_good}_{cluster}_top{top_n}_eigenvector_centrality.svg", dpi=300, bbox_inches='tight')
    plt.show()
# %%
df_plot = links_poor.merged_score.copy()
df_plot.query("degree_out > 0", inplace=True)
df_plot = df_plot.sort_values(['cluster', 'eigenvector_centrality'], ascending=[True, False])
df_plot['rank'] = df_plot.groupby('cluster')['eigenvector_centrality'].rank(ascending=False, method='first')
top_n = 20
df_plot = df_plot[df_plot['rank'] <= top_n]
# df_plot.index.name = 'gene'
df_plot["gene"] = df_plot.index
df_plot.reset_index(drop=True, inplace=True)
# df_plot = df_plot[df_plot['cluster'].isin(['ADM', 'CAF', 'Malignant', 'Ductal'])]

for cluster in df_plot['cluster'].unique():
    fig,ax=plt.subplots(figsize=(6, 6))
    
    df_plot_cluster = df_plot[df_plot['cluster'] == cluster]
    df_plot_cluster = df_plot_cluster.sort_values('eigenvector_centrality', ascending=False)
    sns.scatterplot(data=df_plot_cluster, x='eigenvector_centrality', y='gene',edgecolor="black", color="#0096FF", s=50, ax=ax)
    plt.title(f"Top {top_n} Eigenvector Centrality Poor Response\n{cluster}")
    plt.xlim(min(df_plot_cluster['eigenvector_centrality']) - 0.05, max(df_plot_cluster['eigenvector_centrality']) + 0.05)
    sns.despine()
    ax.grid(False)
    plt.tight_layout()
    fig.savefig(f"{tf_comparison_folder}{filename_poor}_{cluster}_top{top_n}_eigenvector_centrality.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{tf_comparison_folder}{filename_poor}_{cluster}_top{top_n}_eigenvector_centrality.svg", dpi=300, bbox_inches='tight')
    plt.show()

# %% [markdown]
# # NETWORK ANALYSIS

# %% [markdown]
# ## LOAD ORACLE NETWORK AND LINKS

# %%
gene_list=open('/storage/sahuarea/samuele_cancellieri/regev/reference_n_results/personalized_gene_sets/gene_merged.txt','r').readlines()
# gene_list=open('/storage/sahuarea/samuele_cancellieri/regev/filtered_tfs.txt','r').readlines()
gene_list=[x.strip() for x in gene_list]
gene_list.append("GLIS3")
gene_list.sort()
gene_list[:5]

# %%
filtered_gene_list=open('/storage/sahuarea/samuele_cancellieri/regev/reference_n_results/personalized_gene_sets/filtered_TF_rummagenes.txt','r').readlines()
filtered_gene_list=[x.strip() for x in filtered_gene_list]
filtered_gene_list.sort()
filtered_gene_list[:5]

# %%
save_folders['save_path_figures']

# %%
save_path_networks

# %%
filename_treat="regev_epithelial-Malignant_treated"
filename_untreat="regev_epithelial-Malignant_untreated"

oracle_treat = co.load_hdf5(f"{save_path_networks}/{filename_treat}.celloracle.oracle")
links_treat = co.load_hdf5(f"{save_path_networks}/{filename_treat}.celloracle.links")
oracle_untreat = co.load_hdf5(f"{save_path_networks}/{filename_untreat}.celloracle.oracle")
links_untreat = co.load_hdf5(f"{save_path_networks}/{filename_untreat}.celloracle.links")
# %%
oracle_treat.cluster_column_name,oracle_untreat.cluster_column_name
#%%
links_treat.cluster,links_untreat.cluster
# %%
oracle_untreat.adata.obs['level_2_annotation'].value_counts()

# %%
# filename="caf-epithelial(malignant-non_malignant)"
filename="malignant"
# filename="endothelial"
oracle=co.load_hdf5(f"{save_path_networks}/regev_{filename}.celloracle.oracle")
links=co.load_hdf5(f"{save_path_networks}/regev_{filename}.celloracle.links")

# %%
links.filtered_links.keys()

# %%
links.filter_links(p=0.001, weight="coef_abs")

# %%
links.get_network_score()

# %%
links.merged_score

# %%
oracle.adata
#%%
# links.plot_scores_as_rank(cluster='Untreated',n_gene=10)
# links.plot_scores_as_rank(cluster='Poor response',n_gene=10)
links.plot_scores_as_rank(cluster='Moderate response',n_gene=10)
#%%
df_base_promoter=co.data.load_human_promoter_base_GRN()
#%%
oracle.TFdict['DDIT3']
# %%
links.merged_score.to_csv(f"{save_path_metadata}/regev_{filename}.merged_score.csv",index=True,index_label="gene")

# %%
# sc.pl.umap(oracle.adata, color=['response','NFATC2','CEBPD','DDIT3','ATF3','RASGEF1B','SLC26A3'], cmap='coolwarm',vcenter=0,show=True, use_raw=False)
# sc.pl.draw_graph(oracle.adata,color=['level_2_annotation'], cmap='coolwarm',vcenter=0)
sc.pl.paga(oracle.adata,color=['level_2_annotation'], title='Paga GRN',save='_grn.pdf')
# sc.pl.draw_graph(oracle_treat.adata,color=['level_2_annotation','NFATC2','MECOM'], cmap='Reds')
# sc.pl.draw_graph(oracle_untreat.adata,color=['level_2_annotation','NFATC2','MECOM'], cmap='Reds')

# %% [markdown]
# ## NETWORK STATISTICAL ANALYSIS

# %% [markdown]
# ### TREAT vs UNTREAT

# %%
gene_dict=dict()
for cluster in links_treat.cluster:
    print(cluster)
    df_treat=links_treat.filtered_links[cluster]
    df_untreat=links_untreat.filtered_links[cluster]
    # for gene in gene_list:
    for gene in filtered_gene_list:
        gene_dict[f'{cluster}_{gene}']=list()
        if gene in df_treat['source'].to_list() and gene in df_untreat['source'].to_list():
            df_treat_gene=df_treat[df_treat['source']==gene]
            df_untreat_gene=df_untreat[df_untreat['source']==gene]
            gene_dict[f'{cluster}_{gene}'].append(cluster)
            gene_dict[f'{cluster}_{gene}'].append(gene)
            gene_dict[f'{cluster}_{gene}'].append(df_treat_gene['coef_abs'].to_list())
            gene_dict[f'{cluster}_{gene}'].append(df_untreat_gene['coef_abs'].to_list())
            if df_treat_gene.shape[0] and df_untreat.shape[0]:
                # stats,pvalue=scipy.stats.mannwhitneyu(df_treat_gene['coef_abs'],df_untreat_gene['coef_abs'], alternative='two-sided')
                stats,pvalue=scipy.stats.ks_2samp(sorted(df_treat_gene['coef_abs']),sorted(df_untreat_gene['coef_abs']), alternative='two-sided')
                gene_dict[f'{cluster}_{gene}'].append(stats)
                gene_dict[f'{cluster}_{gene}'].append(pvalue)
            else:
                gene_dict[f'{cluster}_{gene}'].append(None)
                gene_dict[f'{cluster}_{gene}'].append(None)
        else:
            gene_dict[f'{cluster}_{gene}'].append(None)
            gene_dict[f'{cluster}_{gene}'].append(None)
            gene_dict[f'{cluster}_{gene}'].append(None)
            gene_dict[f'{cluster}_{gene}'].append(None)
            gene_dict[f'{cluster}_{gene}'].append(None)
            gene_dict[f'{cluster}_{gene}'].append(None)

# %% [markdown]
# ### RESPONDERS

# %%
gene_dict=dict()
for cell_t in links.merged_score.cell_type.unique():
    print(cell_t)
    df_treat=links.filtered_links[f"{cell_t}-MM/MR"]
    df_untreat=links.filtered_links[f"{cell_t}-U/PR"]
    # for gene in gene_list:
    for gene in filtered_gene_list:
        gene_dict_key=f'{cell_t}_{gene}'
        gene_dict[gene_dict_key]=list()
        if gene in df_treat['source'].to_list() and gene in df_untreat['source'].to_list():
            df_treat_gene=df_treat[df_treat['source']==gene]
            df_untreat_gene=df_untreat[df_untreat['source']==gene]
            gene_dict[gene_dict_key].append(cell_t)
            gene_dict[gene_dict_key].append(gene)
            gene_dict[gene_dict_key].append(df_treat_gene['coef_abs'].to_list())
            gene_dict[gene_dict_key].append(df_untreat_gene['coef_abs'].to_list())
            if df_treat_gene.shape[0]>0 and df_untreat.shape[0]>0:
                # stats,pvalue=scipy.stats.mannwhitneyu(df_treat_gene['coef_abs'],df_untreat_gene['coef_abs'], alternative='two-sided')
                stats,pvalue=scipy.stats.ks_2samp(sorted(df_treat_gene['coef_abs']),sorted(df_untreat_gene['coef_abs']), alternative='two-sided')
                gene_dict[gene_dict_key].append(stats)
                gene_dict[gene_dict_key].append(pvalue)
            else:
                gene_dict[gene_dict_key].append(None)
                gene_dict[gene_dict_key].append(None)
        else:
            gene_dict[gene_dict_key].append(None)
            gene_dict[gene_dict_key].append(None)
            gene_dict[gene_dict_key].append(None)
            gene_dict[gene_dict_key].append(None)
            gene_dict[gene_dict_key].append(None)
            gene_dict[gene_dict_key].append(None)

# %% [markdown]
# ### GENERAL ANALYSIS

# %%
df_gene_total=pd.DataFrame(gene_dict, index=["cluster","gene","treat","untreat","stats","pvalue"]).T
df_gene_total.dropna(inplace=True)
df_gene_total['pvalue']=df_gene_total['pvalue'].astype(float)
df_gene_total.sort_values(by=["stats","pvalue"],ascending=[False,True], inplace=True)
df_gene_total_no_filter=df_gene_total.copy()
df_gene_total.query("pvalue<0.05", inplace=True)
df_gene_total['pvalue']=df_gene_total['pvalue'].astype(float)
df_gene_total.groupby("cluster").head(20).to_csv(f"{save_path_metadata}{filename_treat}_{filename_untreat}_stats.csv",sep=",")
df_top20=df_gene_total.groupby("cluster").head(20)

# %%
df_top20.index.name="cluster-gene"
df_top20.to_csv(f"{save_path_metadata}{filename_treat}_{filename_untreat}_stats.csv",sep=",")
# df_top20.to_csv(f"{save_path_metadata}oracle_filtered_stats.csv",sep=",")

# %%
with open(f"{save_path_metadata}{filename_treat}_{filename_untreat}_stats.csv","w") as dd:
    for elem in df_top20.gene.unique().tolist():
        dd.write(f"{elem}\n")

# %%
for gene in df_gene_total.index.tolist():
    # print(gene)
    treat_coef=df_gene_total.loc[gene,"treat"]
    untreat_coef=df_gene_total.loc[gene,"untreat"]
    sns.kdeplot(treat_coef, color="blue", label="Better responders")
    sns.kdeplot(untreat_coef, color="red", label="Worse responders")
    plt.title(gene)
    plt.annotate(f"stats: {df_gene_total.loc[gene,'stats']}", xy=(0.4, 0.75), xycoords='axes fraction')
    plt.annotate(f"p-value: {df_gene_total.loc[gene,'pvalue']}", xy=(0.4, 0.8), xycoords='axes fraction')
    plt.legend()
    # plt.show()
    gene_save=gene.replace(" ","_")
    plt.savefig(f"{save_path_figures}{gene_save}_oracle_filtered_responders.pdf")
    plt.close()

# %%
df_nfatc2_treat=links_treat.filtered_links['CAF'].query("source=='NFATC2'")
df_nfatc2_treat.to_csv(f"{save_path_metadata}/nfatc2_treat_graph.csv",sep=",")

# %%
df_nfatc2_untreat=links_untreat.filtered_links['CAF'].query("source=='NFATC2'")
df_nfatc2_untreat.to_csv(f"{save_path_metadata}/nfatc2_untreat_graph.csv",sep=",")

# %%
df_nfatc2_extracted=dict()

# %%
nfact2_intersection=list(set(df_nfatc2_treat.target.unique()).intersection(set(df_nfatc2_untreat.target.unique())))
df_nfatc2_extracted['intersection']=nfact2_intersection

# %%
nfact2_intersection=list(set(df_nfatc2_treat.target.unique())-(set(df_nfatc2_untreat.target.unique())))
df_nfatc2_extracted['treat-untreat']=nfact2_intersection

# %%
nfact2_intersection=list(set(df_nfatc2_untreat.target.unique())-(set(df_nfatc2_treat.target.unique())))
df_nfatc2_extracted['untreat-treat']=nfact2_intersection

# %%
df_nfatc2_extracted

# %%
save_path_metadata

# %%
df = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in df_nfatc2_extracted.items()]))
df.to_csv(f"{save_path_metadata}/nfatc2_extracted.csv",sep=",",index=False)

#%%
df_nfatc2_treat=pd.read_csv(f"{save_path_metadata}/nfatc2_treat_graph.csv",sep=",",index_col=0)
df_nfatc2_untreat=pd.read_csv(f"{save_path_metadata}/nfatc2_untreat_graph.csv",sep=",",index_col=0)
# %%
dims=df_nfatc2_treat.shape

fig=plt.figure(figsize=(25,25))
plt.title("NFATC2-CAF-Treated",fontsize=24)
plt.annotate(f"number of edges: {dims[0]}", xy=(0.01, 0.95), xycoords='axes fraction',fontsize=16)
co.network_analysis.draw_network(df_nfatc2_treat,return_graph=False)
fig.savefig(f"{save_path_figures}/NFATC2-CAF-Network-Treated.pdf",bbox_inches='tight',dpi=300)
# plt.show()

# %%
dims=df_nfatc2_untreat.shape

fig=plt.figure(figsize=(25,25))
plt.title("NFATC2-CAF-Untreated",fontsize=24)
plt.annotate(f"number of edges: {dims[0]}", xy=(0.01, 0.95), xycoords='axes fraction',fontsize=30)
co.network_analysis.draw_network(df_nfatc2_untreat,return_graph=False)
fig.savefig(f"{save_path_figures}/NFATC2-CAF-Network-Untreated.pdf",bbox_inches='tight',dpi=300)
plt.show()

# %%
filtered_gene_list_clean=set(filtered_gene_list).intersection(set(oracle_treat.adata.var_names.to_list()))
filtered_gene_list_clean=list(sorted(filtered_gene_list_clean))

# %%
sc.pl.heatmap(oracle_treat.adata, var_names=filtered_gene_list_clean, groupby='level_2_annotation',cmap='inferno', use_raw=False,vmin=-2, vmax=2,swap_axes=True,dendrogram=True)

# %%
filtered_gene_list_clean=set(filtered_gene_list).intersection(set(oracle_untreat.adata.var_names.to_list()))
filtered_gene_list_clean=list(sorted(filtered_gene_list_clean))

# %%
sc.pl.heatmap(oracle_untreat.adata, var_names=filtered_gene_list_clean, groupby='level_2_annotation',cmap='inferno', use_raw=False,vmin=-2, vmax=2,swap_axes=True,dendrogram=True)

# %% [markdown]
# ## PLOT SHUFFLING RESULTS

# %% [markdown]
# ### TREATED
# %%
df_shuffle=pd.read_csv("/storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/reference_n_results/shuffled_scores/epithelial-CAF_treated_merged_score_shuffled_all.csv",sep=",",index_col=0)
df_shuffle.index.name="gene"
df_shuffle.head()

# %%
for cluster in df_shuffle['cluster'].unique():
    sns.kdeplot(df_shuffle.query("gene=='JUN' & cluster==@cluster")['eigenvector_centrality'], label=cluster)
plt.legend(loc='upper left',bbox_to_anchor=(1,1))
plt.title("JUN - Treated\nEigenvector centrality for 10000 shuffled networks")
plt.show()

for cluster in df_shuffle['cluster'].unique():
    sns.kdeplot(df_shuffle.query("cluster==@cluster")['eigenvector_centrality'], label=cluster)
plt.legend(loc='upper left',bbox_to_anchor=(1,1))
plt.title("Treated\nEigenvector centrality for 10000 shuffled networks")
plt.show()

# %%
df_stats=pd.DataFrame(columns=["cluster","stats","pvalue"])

for cluster in df_shuffle['cluster'].unique():
    dd=scipy.stats.ttest_1samp(links_treat.merged_score.query("cluster==@cluster")['eigenvector_centrality'],df_shuffle.query("cluster==@cluster")['eigenvector_centrality'].mean())
    df_stats=df_stats.append({"cluster":cluster,"stats":dd[0],"pvalue":dd[1]},ignore_index=True)
    
    fig=plt.figure(figsize=(6,5))
    sns.kdeplot(links_treat.merged_score.query("cluster==@cluster")['eigenvector_centrality'], label=f"Real")
    sns.kdeplot(df_shuffle.query("cluster==@cluster")['eigenvector_centrality'], label=f"10k_shuffled")
    plt.annotate(f"stats: {dd[0]}", xy=(0.4, 0.75), xycoords='axes fraction')
    plt.annotate(f"p-value: {dd[1]}", xy=(0.4, 0.8), xycoords='axes fraction')
    plt.legend(loc='upper right')
    sns.despine(fig=fig)
    plt.title(f"Treated - {cluster}\nEigenvector centrality for 10000 shuffled networks vs real")
    plt.show()

# df_stats
plt.close('all')

# %%
df_stats=pd.DataFrame(columns=["cluster","stats","pvalue"])

for cluster in df_shuffle['cluster'].unique():
    dd=scipy.stats.ks_2samp(links_treat.merged_score.query("cluster==@cluster")['eigenvector_centrality'],df_shuffle.query("cluster==@cluster")['eigenvector_centrality'])
    df_stats=df_stats.append({"cluster":cluster,"stats":dd[0],"pvalue":dd[1]},ignore_index=True)
    
    fig=plt.figure(figsize=(6,5))
    aa=sns.kdeplot(links_treat.merged_score.query("cluster==@cluster")['eigenvector_centrality'], label=f"Real")
    bb=sns.kdeplot(df_shuffle.query("cluster==@cluster")['eigenvector_centrality'], label=f"10k_shuffled")
    plt.annotate(f"stats: {dd[0]}", xy=(0.4, 0.75), xycoords='axes fraction')
    plt.annotate(f"p-value: {dd[1]}", xy=(0.4, 0.8), xycoords='axes fraction')
    plt.legend(loc='upper right')
    sns.despine(fig=fig)
    plt.grid(False)
    plt.title(f"Treated - {cluster}\nEigenvector centrality for 10000 shuffled networks vs real")
    plt.savefig(f"/storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/manuscript_figures/regev_treated_{cluster}_shuffled_ks.pdf",bbox_inches='tight',dpi=300)
    plt.show()

# df_stats
plt.close('all')

# %% [markdown]
# ### UNTREATED

# %%
df_shuffle=pd.read_csv("/storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/regev/reference_n_results/shuffled_scores/epithelial-CAF_untreated_merged_score_shuffled_all.csv",sep=",",index_col=0)
df_shuffle.index.name="gene"
df_shuffle.head()

# %%
for cluster in df_shuffle['cluster'].unique():
    sns.kdeplot(df_shuffle.query("gene=='JUN' & cluster==@cluster")['eigenvector_centrality'], label=cluster)
plt.legend(loc='upper left',bbox_to_anchor=(1,1))
plt.title("JUN - Untreated\nEigenvector centrality for 10000 shuffled networks")
plt.show()

for cluster in df_shuffle['cluster'].unique():
    sns.kdeplot(df_shuffle.query("cluster==@cluster")['eigenvector_centrality'], label=cluster)
plt.legend(loc='upper left',bbox_to_anchor=(1,1))
plt.title("Untreated\nEigenvector centrality for 10000 shuffled networks")
plt.show()

# %%
df_stats=pd.DataFrame(columns=["cluster","stats","pvalue"])

for cluster in df_shuffle['cluster'].unique():
    dd=scipy.stats.ttest_1samp(links_untreat.merged_score.query("cluster==@cluster")['eigenvector_centrality'],df_shuffle.query("cluster==@cluster")['eigenvector_centrality'].mean())
    df_stats=df_stats.append({"cluster":cluster,"stats":dd[0],"pvalue":dd[1]},ignore_index=True)
    
    fig=plt.figure(figsize=(6,5))
    sns.kdeplot(links_untreat.merged_score.query("cluster==@cluster")['eigenvector_centrality'], label=f"Real")
    sns.kdeplot(df_shuffle.query("cluster==@cluster")['eigenvector_centrality'], label=f"10k_shuffled")
    plt.annotate(f"stats: {dd[0]}", xy=(0.4, 0.75), xycoords='axes fraction')
    plt.annotate(f"p-value: {dd[1]}", xy=(0.4, 0.8), xycoords='axes fraction')
    plt.legend(loc='upper right')
    plt.title(f"Untreated - {cluster}\nEigenvector centrality for 10000 shuffled networks vs real")
    plt.show()

# df_stats
plt.close('all')

# %%
df_stats=pd.DataFrame(columns=["cluster","stats","pvalue"])

for cluster in df_shuffle['cluster'].unique():
    dd=scipy.stats.ks_2samp(links_untreat.merged_score.query("cluster==@cluster")['eigenvector_centrality'],df_shuffle.query("cluster==@cluster")['eigenvector_centrality'])
    df_stats=df_stats.append({"cluster":cluster,"stats":dd[0],"pvalue":dd[1]},ignore_index=True)
    
    fig=plt.figure(figsize=(6,5))
    aa=sns.kdeplot(links_untreat.merged_score.query("cluster==@cluster")['eigenvector_centrality'], label=f"Real")
    bb=sns.kdeplot(df_shuffle.query("cluster==@cluster")['eigenvector_centrality'], label=f"10k_shuffled")
    plt.annotate(f"stats: {dd[0]}", xy=(0.4, 0.75), xycoords='axes fraction')
    plt.annotate(f"p-value: {dd[1]}", xy=(0.4, 0.8), xycoords='axes fraction')
    plt.legend(loc='upper right')
    sns.despine(fig=fig)
    plt.grid(False)
    plt.title(f"Untreated - {cluster}\nEigenvector centrality for 10000 shuffled networks vs real")
    plt.savefig(f"/storage/sahuarea/samuele_cancellieri/PDAC_complete_analysis/data/manuscript_figures/regev_untreated_{cluster}_shuffled_ks.pdf",bbox_inches='tight',dpi=300)
    plt.show()

df_stats
plt.close('all')

# %% [markdown]
# ## PERTURBATION

# %% [markdown]
# ### PERTURBATION MERGED DATASET (TREATvUNTREAT)

# %%
oracle,links

# %%
links.filter_links(p=0.001,weight="coef_abs")
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True,GRN_unit='cluster')

# %%
filtered_gene_list_dict=dict()
for gene in filtered_gene_list:
    if gene in oracle_treat.adata.var_names:
        filtered_gene_list_dict[gene]=float(0)
    else:
        print(f"{gene} not in TFdict")

# %%
sorted(oracle.cluster_specific_TFdict['Poor response']['HMGA2'])

# %%
oracle.cluster_specific_TFdict['Poor response']

# %%
oracle.cluster_specific_TFdict['Poor response'].items()

# %%
oracle.cluster_specific_TFdict.keys()

# %%
# Explode the dictionary keys for each value and create a dataframe
exploded_data = []

for key_og in oracle.cluster_specific_TFdict.keys():
    for key, values in oracle.cluster_specific_TFdict[key_og].items():
        for value in values:
            exploded_data.append([key, value, key_og])

df_exploded = pd.DataFrame(exploded_data, columns=['source', 'target', 'cluster'])
df_exploded.to_csv(f"{save_path_metadata}/regev_malignant_response_cluster_specific_TFdict.csv",sep=",",index=False)

# %%
# Enter perturbation conditions to simulate signal propagation after the perturbation.
perturb_dict=dict()
perturb_genes=['ZBTB16', 'BCOR']
for gene in perturb_genes:
    perturb_dict[gene]=float(0)
oracle.simulate_shift(perturb_condition=perturb_dict, n_propagation=3, GRN_unit='cluster')

# %%
sc.pl.umap(oracle.adata, color=['level_2_annotation','response',*perturb_genes],cmap='inferno',vcenter=0,show=True, use_raw=False)
sc.pl.umap(oracle.adata, color=['level_2_annotation','response',*perturb_genes], layer='simulation_input',cmap='inferno',vcenter=0,show=True, use_raw=False)
sc.pl.umap(oracle.adata, color=['level_2_annotation','response',*perturb_genes], layer='simulated_count',cmap='inferno',vcenter=0,show=True, use_raw=False)

# %%
oracle.adata

# %%
df=sc.get.obs_df(oracle.adata,keys=['response',*perturb_genes],layer='simulation_input')
fig,ax=plt.subplots(1,1,figsize=(5,5))
sns.boxplot(data=df.melt(id_vars=['response'],value_vars=perturb_genes),x='response',y='value',hue='variable',ax=ax,hue_order=sorted(perturb_genes))
plt.xticks(rotation=90)

# %%
df=sc.get.obs_df(oracle.adata,keys=['response',*perturb_genes],layer='simulated_count')
fig,ax=plt.subplots(1,1,figsize=(5,5))
sns.boxplot(data=df.melt(id_vars=['response'],value_vars=perturb_genes),x='response',y='value',hue='variable',ax=ax,hue_order=sorted(perturb_genes))
plt.xticks(rotation=90)

# %%
oracle.adata.obs['response'].value_counts()

# %%
sc.tl.rank_genes_groups(oracle.adata, groupby='response', method='wilcoxon',use_raw=False,layer='simulated_count',groups=['Moderate response'],reference='Poor response')

# %%
sc.pl.rank_genes_groups_matrixplot(oracle.adata,dendrogram=False,standard_scale='var',n_genes=10)

# %%
df_simulated_count=sc.get.rank_genes_groups_df(oracle.adata,group='Moderate response',pval_cutoff=0.05)
df_simulated_count

# %%
sc.tl.rank_genes_groups(oracle.adata, groupby='response', method='wilcoxon',use_raw=False,layer='simulation_input',groups=['Moderate response'],reference='Poor response')

# %%
sc.pl.rank_genes_groups_matrixplot(oracle.adata,dendrogram=False,standard_scale='var',n_genes=10)

# %%
df_simulation_input=sc.get.rank_genes_groups_df(oracle.adata,group='Moderate response',pval_cutoff=0.05)
df_simulation_input

# %%
# Merge the two dataframes on the 'names' column
df_merged = pd.merge(df_simulated_count[['names', 'logfoldchanges']], df_simulation_input[['names', 'logfoldchanges']], on='names', suffixes=('_simulated', '_input'))

# Calculate the absolute difference in log fold changes
df_merged['logfoldchange_diff'] = abs(df_merged['logfoldchanges_simulated']) - abs(df_merged['logfoldchanges_input'])

# Sort the dataframe by the absolute difference in descending order
df_merged_sorted = df_merged.sort_values(by='logfoldchange_diff', ascending=False)

# Report the top differences
df_merged_sorted.head(50)

df_merged_sorted.to_csv(f"{save_path_metadata}/regev_malignant_response_simulated_vs_input.csv",sep=",",index=False)

# %%
save_path_metadata

# %% [markdown]
# ### PERTURBATION RESPONSE (BETTERvWORSE)

# %%
links.filter_links(p=0.001,weight="coef_abs")
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True, GRN_unit='cluster')

# %%
filtered_gene_list_dict=dict()
for gene in filtered_gene_list:
    if gene in oracle.adata.var_names:
        filtered_gene_list_dict[gene]=float(0)
    else:
        print(f"{gene} not in TFdict")

# %%
# Enter perturbation conditions to simulate signal propagation after the perturbation.
perturb_dict=dict()
perturb_genes=["BCOR","NFATC1","SPI1","ZBTB16"]
for gene in perturb_genes:
    perturb_dict[gene]=float(0)
oracle.simulate_shift(perturb_condition=perturb_dict, n_propagation=3, GRN_unit='cluster')

# %%
oracle.adata.layers

# %%
sc.get.obs_df(oracle.adata,keys=['BCOR'],layer="imputed_count").max()

# %%
sc.get.obs_df(oracle.adata,keys=['BCOR'],layer="simulation_input").max()

# %%
sc.get.obs_df(oracle.adata,keys=['BCOR'],layer="simulated_count").max()

# %%
filename

# %%
print("-".join(perturb_genes))
perturbed_genes="-".join(perturb_genes)

# %%
oracle.to_hdf5(f"{save_path_networks}/regev_epithelial(malignant)_perturbed_{perturbed_genes}.celloracle.oracle")

# %%
oracle.adata

# %% [markdown]
# #### HOTSPOT ANALYSIS

# %%
## hotspot already imported at top

# %%
score_dict=dict()
adata_cc=oracle.adata.copy()
sc.pp.filter_cells(adata_cc, min_counts=1)
sc.pp.filter_genes(adata_cc, min_cells=1)
mod_score=uti.myHotspot.hotspot_calculation(adata=adata_cc, counts='processed', score_dict=score_dict, save_paths=save_folders,filename='caf-epithelial(malignant-non_malignant)', perturbed_genes=perturbed_genes)

# %%
with open(f'regev_caf-epithelial(malignant-non_malignant)_hotspot_score_dict.pkl', 'wb') as fp:
    pickle.dump(score_dict, fp)
    print(f'dictionary saved successfully to {fp}')

# %%
mod_score.to_csv(f"{save_path_metadata}/regev_caf-epithelial(malignant-non_malignant)_mod_score.csv",sep=",")

# %% [markdown]
# ### PERTURBATION TREATED

# %%
links_treat.filter_links(p=0.001,weight="coef_abs")
oracle_treat.get_cluster_specific_TFdict_from_Links(links_object=links_treat)
oracle_treat.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

# %%
filtered_gene_list_dict=dict()
for gene in filtered_gene_list:
    if gene in oracle_treat.adata.var_names:
        filtered_gene_list_dict[gene]=float(0)
    else:
        print(f"{gene} not in TFdict")

# %%
# Enter perturbation conditions to simulate signal propagation after the perturbation.
perturb_dict=dict()
perturb_genes=["MECOM"]
for gene in perturb_genes:
    perturb_dict[gene]=float(0)
oracle_treat.simulate_shift(perturb_condition=perturb_dict, n_propagation=3)

# %%
oracle_treat.adata.layers

# %%
perturbed_genes="-".join(perturb_genes)
oracle_treat.to_hdf5(f"{save_path_networks}/regev_{filename_treat}_perturbated_{perturbed_genes}.celloracle.oracle")

# %%
def hotspot_calculation(adata, score_dict, counts, filename, layer=None,model='danb',latent_obsm='X_pca',umi_counts='total_counts'):
    # perturbed_genes="-".join(perturb_genes)
    adata.X=csc_matrix(adata.X)
    adata_hs = hs.Hotspot(adata, layer_key=layer, model=model, latent_obsm_key=latent_obsm, umi_counts_obs_key=umi_counts)
    adata_hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
    hs_results = adata_hs.compute_autocorrelations(jobs=48)
    hs_genes = hs_results.loc[hs_results.FDR < 0.05].index # Select genes
    local_correlations = adata_hs.compute_local_correlations(hs_genes, jobs=48) # jobs for parallelization
    modules = adata_hs.create_modules(min_gene_threshold=30, core_only=True, fdr_threshold=0.05)
    adata_hs.plot_local_correlations(yticklabels=False)
    plt.title(f"{counts}",size=10)
    plt.savefig(f"{save_path_figures}/{filename}_{perturbed_genes}_{counts}_local_correlations.pdf")
    plt.close()
    module_scores = adata_hs.calculate_module_scores()
    for col in module_scores.columns.tolist():
        adata.obs[f'mod_{counts}_{col}']=module_scores[col]
        sc.pl.umap(adata, color=['level_2_annotation',f'mod_{counts}_{col}'], ncols=2, wspace= 0.5, size=100, save=f"_{filename}_{perturbed_genes}_{counts}_mod_{col}.pdf", show=False)
    score_dict[f'{counts}']=adata_hs.results.join(adata_hs.modules)
    return module_scores

# %%
score_dict=dict()
# perturbed_genes="-".join(perturb_genes)
for counts in ["counts","normalized_count","imputed_count","simulated_count"]:
    print(counts.upper())
    # substitute the counts with the counts of interest (normalized, imputed, simulated) and run PCA and UMAP after basic filtering
    adata=oracle.adata.copy()
    adata.X=adata.layers[counts]
    sc.pp.filter_cells(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.pca(adata,svd_solver="arpack")
    sc.pl.pca(adata,color="level_2_annotation",save=f"_{filename}_{perturbed_genes}_{counts}.pdf",show=False,title=f"{filename}_{perturbed_genes}_{counts}")
    sc.tl.umap(adata)
    sc.pl.umap(adata,color='level_2_annotation',save=f"_{filename}_{perturbed_genes}_{counts}.pdf",show=False,title=f"{filename}_{perturbed_genes}_{counts}")
    mod_scores=hotspot_calculation(adata, score_dict, counts, filename=filename, model='danb',latent_obsm='X_pca',umi_counts='total_counts')
    mod_scores.to_csv(f"{save_path_metadata}/{filename}_{perturbed_genes}_{counts}_hotspot_scores.csv",sep=",")

with open(f'{filename}_{perturbed_genes}_hotspot_score_dict.pkl', 'wb') as fp:
    pickle.dump(score_dict, fp)
    print(f'dictionary saved successfully to {fp}')

# %% [markdown]
# ### PERTURBATION UNTREATED

# %%
links_untreat.filter_links(p=0.001,weight="coef_abs")
oracle_untreat.get_cluster_specific_TFdict_from_Links(links_object=links_untreat)
oracle_untreat.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

# %%
filtered_gene_list_dict=dict()
for gene in filtered_gene_list:
    if gene in oracle_untreat.adata.var_names:
        filtered_gene_list_dict[gene]=float(0)
    else:
        print(f"{gene} not in TFdict")

# %%
# Enter perturbation conditions to simulate signal propagation after the perturbation.
perturb_dict=dict()
perturb_genes=["MECOM"]
for gene in perturb_genes:
    perturb_dict[gene]=float(0)
oracle_untreat.simulate_shift(perturb_condition=perturb_dict, n_propagation=3)

# %%
oracle_untreat.adata.layers

# %%
perturbed_genes="-".join(perturb_genes)
oracle_untreat.to_hdf5(f"{save_path_networks}/regev_{filename_untreat}_perturbated_{perturbed_genes}.celloracle.oracle")

# %% [markdown]
# ### CELL STATE TRANSITION

# %%
print(f"running with {'-'.join(perturb_genes)}")

# Get transition probability
oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True, 
                                sampled_fraction=1)

# Calculate embedding 
oracle.calculate_embedding_shift(sigma_corr=0.05)

# %%
oracle_treat.adata

# %%
oracle_treat.adata.layers['delta_X']

# %%
fig, ax = plt.subplots(1, 1, figsize=(5, 5))

scale = 19
# Show quiver plot
oracle.plot_quiver(scale=scale,ax=ax)
# fig.savefig(f"{save_path_figures}/quiver_{filename_treat}_{'-'.join(perturb_genes)}.pdf",dpi=300, bbox_inches='tight')
# fig.show()
# %%
fig, ax = plt.subplots(1, 3,  figsize=[15, 6])

scale = 19
# Show quiver plot
oracle.plot_quiver(scale=scale, ax=ax[0])
ax[0].set_title(f"{'-'.join(perturb_genes)} simulation vector")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_quiver_random(scale=scale, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")

sc.pl.draw_graph(oracle.adata, color=[oracle.cluster_column_name], layer="counts", use_raw=False, cmap="inferno",ax=ax[2])

fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_quiver.pdf",dpi=300, bbox_inches='tight')
fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_quiver.svg",dpi=300, bbox_inches='tight')
# plt.subplots_adjust(wspace=)
plt.show()

# %%
n_grid = 40
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)

# %%
oracle.suggest_mass_thresholds(n_suggestion=20)

# %%
min_mass = 0.001
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)

# %%
perturb_genes

# %%
gois="-".join(perturb_genes)

fig, ax = plt.subplots(1, 2,  figsize=[13, 6])

scale_simulation = 0.02
# Show quiver plot
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0],show_background=True,s=5)
ax[0].set_title(f"Simulated cell identity shift vector: {gois} KO")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1],show_background=True,s=5)
ax[1].set_title(f"Randomized simulation vector")

fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_simulation_flow.pdf",dpi=300, bbox_inches='tight')
fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_simulation_flow.svg",dpi=300, bbox_inches='tight')
plt.show()

# %%
# Plot vector field with cell cluster
fig, axs = plt.subplots(1,2,figsize=[12, 8])

oracle.plot_cluster_whole(ax=axs[0], s=2)
oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=axs[0], show_background=False)
sc.pl.draw_graph(oracle.adata, color=[oracle.cluster_column_name],ax=axs[1],show=False)
# oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax, show_background=False)
fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_cluster_with_grid.pdf",dpi=300, bbox_inches='tight')
fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_cluster_with_grid.svg",dpi=300, bbox_inches='tight')

# %%
# Plot vector field with cell cluster
fig, ax = plt.subplots(figsize=[8, 8])

oracle.plot_cluster_whole(ax=ax, s=2)
# oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax, show_background=False)
fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_cluster_random_with_grid.pdf",dpi=300, bbox_inches='tight')
fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_cluster_random_with_grid.svg",dpi=300, bbox_inches='tight')

# %%
fig=sc.pl.umap(oracle.adata, color=[oracle.cluster_column_name],return_fig=True)
fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_umap.pdf",dpi=300, bbox_inches='tight')
fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_umap.svg",dpi=300, bbox_inches='tight')

# %%
fig=sc.pl.draw_graph(oracle.adata, color=[oracle.cluster_column_name],return_fig=True)
fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_draw-graph.pdf",dpi=300, bbox_inches='tight')
fig.savefig(f"{save_path_figures}/regev_{filename}_perturbed_{'-'.join(perturb_genes)}_draw-graph.svg",dpi=300, bbox_inches='tight')

# %% [markdown]
# ## PLOTS

# %% [markdown]
# ### RESETTING

# %%
uti.myPlotting.reset_plotting_settings()

# %%
plt.close('all')

# %%
fig=sc.pl.umap(oracle_treat.adata,color=['new_celltypes','NFATC2','MECOM'], ncols=1,return_fig=True)
fig.suptitle(f"{filename_treat}",fontsize=20)
plt.show()

fig=sc.pl.umap(oracle_untreat.adata,color=['new_celltypes','NFATC2','MECOM'], ncols=1,return_fig=True)
fig.suptitle(f"{filename_untreat}",fontsize=20)

# %% [markdown]
# ### EXTRACTION

# %%
save_path_matrices

# %%
oracle.adata.obs.level_2_annotation.value_counts()

# %%
links.merged_score

# %%
for cluster in links.cluster:
    dd=links.merged_score.query("cluster==@cluster")
    dd.sort_values('betweenness_centrality',ascending=False,inplace=True)
    print(dd)
    # dd.head(50).to_csv(f"{save_path_metadata}/regev_{filename}_{cluster.lower()}_top50_genes.csv",sep=",",index=True,index_label="gene")

# %%
dd=links.merged_score
dd.query("cluster=='Treated'").sort_values("betweenness_centrality",ascending=False).head(50).to_csv(f"{save_path_metadata}/regev_{filename.lower()}_treated_top50_genes.csv",sep=",",index=True,index_label="gene")
dd.query("cluster=='Untreated'").sort_values("betweenness_centrality",ascending=False).head(50).to_csv(f"{save_path_metadata}/regev_{filename.lower()}_untreated_top50_genes.csv",sep=",",index=True,index_label="gene")
# dd.to_csv(f"{save_path_metadata}/regev_caf-epithelial(malignant-non_malignant)_top100.csv",sep=",",index=True,index_label="gene")

# %%
with open(f"{save_path_metadata}/regev_caf-epithelial(malignant-non_malignant)_top100.txt","w") as ss:
    for elem in dd.index.unique().tolist():
        ss.write(f"{elem}\n")

# %% [markdown]
# ### TREATMENT AND RESPONSE

# %%
sc.pl.highest_expr_genes(oracle_treat.adata, n_top=20)
sc.pl.highest_expr_genes(oracle_untreat.adata, n_top=20)

# %%
sc.pl.highest_expr_genes(oracle.adata, n_top=20)

# %%
sc.pl.umap(oracle_treat.adata, color=['level_2_annotation'], ncols=1)
sc.pl.umap(oracle_untreat.adata, color=['level_2_annotation'], ncols=1)

# %% [markdown]
# ### SINGLE CELLTYPE TOP GENES

# %%
def plot_scores_as_rank(links, cluster,values, n_gene=50):
    """
    Pick up top n-th genes wich high-network scores and make plots.

    Args:
        links (Links object): See network_analisis.Links class for detail.
        cluster (str): Cluster nome to analyze.
        n_gene (int): Number of genes to plot. Default is 50.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """
    # values = ['degree_centrality_all',
    #               'degree_centrality_in', 'degree_centrality_out',
    #               'betweenness_centrality',  'eigenvector_centrality']
    for value in values:

        res = links.merged_score[links.merged_score.cluster == cluster]
        res = res[value].sort_values(ascending=False)
        res = res[:n_gene]

        # fig = plt.figure()

        # ppp=plt.scatter(res.values, range(len(res)))
        ppp=sns.scatterplot(x=res.values, y=range(len(res)))
        plt.yticks(range(len(res)), res.index.values)#, rotation=90)
        plt.xlabel(value)
        # plt.xlim((0,1))
        # plt.xticks(np.arange(0, 1.1, 0.5))
        plt.title(f"Top {n_gene} genes in {cluster}")
        plt.gca().invert_yaxis()
        # plt.subplots_adjust(left=0.5, right=0.85)

        # if not save is None:
        #     os.makedirs(save, exist_ok=True)
        #     path = os.path.join(save, f"ranked_values_in_{links.name}_{value}_{links.threshold_number}_in_{cluster}.{settings['save_figure_as']}")
        #     fig.savefig(path, transparent=True)
        # plt.show()
        # return fig
        return ppp

# %%
sns.reset_defaults()
plt.rcdefaults()

#%%
for cluster in links_treat.cluster:
    pplt=plot_scores_as_rank(links_treat,cluster=cluster,values=['eigenvector_centrality'], n_gene=20)
    # plt.savefig(f"{save_path_figures}/regev_{filename}_{cluster.lower()}_eigenvector_centrality.pdf",dpi=300, bbox_inches='tight')
    plt.show()
# %%
# import plotnine as pn
# Prepare data for plotnine
df_plot = links_treat.merged_score.copy()
df_plot.query("degree_out > 0", inplace=True)
df_plot = df_plot.sort_values(['cluster', 'eigenvector_centrality'], ascending=[True, False])
df_plot['rank'] = df_plot.groupby('cluster')['eigenvector_centrality'].rank(ascending=False, method='first')
top_n = 20
df_plot = df_plot[df_plot['rank'] <= top_n]
# df_plot.index.name = 'gene'
df_plot["gene"] = df_plot.index
df_plot.reset_index(drop=True, inplace=True)
df_plot = df_plot[df_plot['cluster'].isin(['ADM', 'CAF', 'Malignant', 'Ductal'])]

for cluster in df_plot['cluster'].unique():
    fig,ax=plt.subplots(figsize=(6, 6))
    
    df_plot_cluster = df_plot[df_plot['cluster'] == cluster]
    df_plot_cluster = df_plot_cluster.sort_values('eigenvector_centrality', ascending=False)
    sns.scatterplot(data=df_plot_cluster, x='eigenvector_centrality', y='gene',edgecolor="black", color="#0096FF", s=50, ax=ax)
    plt.title(f"Top {top_n} Eigenvector Centrality Treated\n{cluster}")
    plt.xlim(min(df_plot_cluster['eigenvector_centrality']) - 0.05, max(df_plot_cluster['eigenvector_centrality']) + 0.05)
    sns.despine()
    ax.grid(False)
    plt.tight_layout()
    fig.savefig(f"{save_path_figures}{filename_treat}_{cluster}_top{top_n}_eigenvector_centrality.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{save_path_figures}{filename_treat}_{cluster}_top{top_n}_eigenvector_centrality.svg", dpi=300, bbox_inches='tight')
    plt.show()
# %%
# Prepare data for plotnine
df_plot = links_untreat.merged_score.copy()
df_plot.query("degree_out > 0", inplace=True)
df_plot = df_plot.sort_values(['cluster', 'eigenvector_centrality'], ascending=[True, False])
df_plot['rank'] = df_plot.groupby('cluster')['eigenvector_centrality'].rank(ascending=False, method='first')
top_n = 20
df_plot = df_plot[df_plot['rank'] <= top_n]
# df_plot.index.name = 'gene'
df_plot["gene"] = df_plot.index
df_plot.reset_index(drop=True, inplace=True)
df_plot = df_plot[df_plot['cluster'].isin(['ADM', 'CAF', 'Malignant', 'Ductal'])]

for cluster in df_plot['cluster'].unique():
    fig,ax=plt.subplots(figsize=(6, 6))
    
    df_plot_cluster = df_plot[df_plot['cluster'] == cluster]
    df_plot_cluster = df_plot_cluster.sort_values('eigenvector_centrality', ascending=False)
    sns.scatterplot(data=df_plot_cluster, x='eigenvector_centrality', y='gene', edgecolor="black", color="#0096FF", s=50, ax=ax)
    plt.title(f"Top {top_n} Eigenvector Centrality Untreated\n{cluster}")
    plt.xlim(min(df_plot_cluster['eigenvector_centrality']) - 0.05, max(df_plot_cluster['eigenvector_centrality']) + 0.05)
    sns.despine()
    ax.grid(False)
    plt.tight_layout()
    fig.savefig(f"{save_path_figures}{filename_untreat}_{cluster}_top{top_n}_eigenvector_centrality.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{save_path_figures}{filename_untreat}_{cluster}_top{top_n}_eigenvector_centrality.svg", dpi=300, bbox_inches='tight')
    plt.show()
# %%
# links=links_treat

nrows=3
ncols=math.ceil(len(links.cluster)/nrows)
fig=plt.figure(figsize=(nrows*4, ncols*5))

measure="eigenvector_centrality"
# measure="betweenness_centrality"

for count,cluster in enumerate(links.cluster):
    fig.add_subplot(nrows, ncols, count+1)
    pplt=plot_scores_as_rank(links,cluster=cluster,values=[measure], n_gene=20)
    pplt.set_xlabel("")

fig.suptitle(f"{filename.upper()}\n{measure.capitalize()}", fontsize=12)
fig.subplots_adjust(top=0.96, hspace=0.25, wspace=0.5)
plt.tight_layout()
plt.show()
fig.savefig(f"{save_path_figures}{filename}_celltypes-{measure}.pdf", dpi=300)

# %% [markdown]
# ### CELLTYPES COMPARISON

# %%
def plot_score_comparison_2D(links, value, cluster1, cluster2, percentile=99, dot_color='black', annot_shifts=None, save=None, fillna_with_zero=True, plt_show=True):
    """
    Make a scatter plot that shows the relationship of a specific network score in two groups.

    Args:
        links (Links object): See network_analisis.Links class for detail.
        value (srt): The network score to be shown.
        cluster1 (str): Cluster nome to analyze. Network scores in the cluste1 are shown as x-axis.
        cluster2 (str): Cluster nome to analyze. Network scores in the cluste2 are shown as y-axis.
        percentile (float): Genes with a network score above the percentile will be shown with annotation. Default is 99.
        annot_shifts ((float, float)): Shift x and y cordinate for annotations.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """
    res = links.merged_score[links.merged_score.cluster.isin([cluster1, cluster2])][[value, "cluster"]]
    res = res.reset_index(drop=False)
    piv = pd.pivot_table(res, values=value, columns="cluster", index="index")
    if fillna_with_zero:
        piv = piv.fillna(0)
    else:
        piv = piv.fillna(piv.mean(axis=0))
    piv["sum"]=piv[cluster1]+piv[cluster2]
    
    c1_name=cluster1.replace(" ","_")
    c2_name=cluster2.replace(" ","_")
    # piv.sort_values(["sum"],ascending=False).to_csv(f"{save_path_metadata}{filename}_{value}_{c1_name}_vs_{c2_name}.csv")

    goi1 = piv[piv[cluster1] > np.percentile(piv[cluster1].values, percentile)].index
    goi2 = piv[piv[cluster2] > np.percentile(piv[cluster2].values, percentile)].index

    gois = np.union1d(goi1, goi2)
    # plt.close('all')
    # sns.reset_defaults()

    x, y = piv[cluster1], piv[cluster2]
    ppp=sns.scatterplot(x=x, y=y, markers="o",s=20)
    # plt.title(f"{value}")

    if annot_shifts is None:
        x_shift, y_shift = (x.max() - x.min())*0.03, (y.max() - y.min())*0.03
    else:
        x_shift, y_shift = annot_shifts
    texts=list()
    for goi in gois:
        x, y = piv.loc[goi, cluster1], piv.loc[goi, cluster2]
        texts.append(plt.text(x, y, goi, size=10))
        # _plot_goi(x, y, goi, {}, scatter=False, x_shift=x_shift, y_shift=y_shift)
    # adjust_text.adjust_text(texts, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
    
    if plt_show:
        plt.show()
        return gois,ppp,texts
    else:
        return gois,ppp,texts

# %%
# all_genes=list()
# fig, axs = plt.subplots(1,2, figsize=(20, 7))
# rep_count=0
# filename=filename_treat.split("_")[0]+"_"+"responders"


measure="eigenvector_centrality"
for count,sett in enumerate(list(itertools.combinations(links_untreat.cluster, 2))):
    cols=len(sett)
    cols=1
    plt.close('all')
    sns.reset_defaults()
    fig=plt.figure(figsize=(7,6))
    cluster1=str(sett[0])
    cluster2=str(sett[1])
    # plot eigenvector centrality
    # fig.add_subplot(1,cols,1)
    genes_t,plot_t,texts_t=plot_score_comparison_2D(links_untreat, measure, cluster1=cluster1,cluster2=cluster2,dot_color='blue', percentile=99, annot_shifts=None, save=None, fillna_with_zero=True, plt_show=False)
    # genes_u,plot_u,texts_u=plot_score_comparison_2D(links_untreat, measure,cluster1=cluster1,cluster2=cluster2,dot_color='red', percentile=99, annot_shifts=None, save=None, fillna_with_zero=True, plt_show=False)
    # plt.legend(['Treated','Untreated'],loc='upper left', bbox_to_anchor=(1, 1))
    # plt.legend()
    # adjust_text.adjust_text(texts=texts_t+texts_u, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
    adjust_text.adjust_text(texts=texts_t, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
    # plot betweenness centrality
    # fig.add_subplot(1,cols,2)
    # genes_t,plot_t,texts_t=plot_score_comparison_2D(links_treat, "betweenness_centrality", cluster1=cluster1,cluster2=cluster2,dot_color='blue', percentile=99, annot_shifts=None, save=None, fillna_with_zero=True, plt_show=False)
    # genes_u,plot_u,texts_u=plot_score_comparison_2D(links_untreat, "betweenness_centrality", cluster1=cluster1,cluster2=cluster2, dot_color='red',percentile=99, annot_shifts=None, save=None, fillna_with_zero=True, plt_show=False)
    # adjust_text.adjust_text(texts=texts_t+texts_u, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
    # #fix cluster names to savefig
    # cluster1=cluster1.replace(" ","_").strip()
    # cluster2=cluster2.replace(" ","_").strip()
    sns.despine()
    plt.grid(False)
    plt.xticks([0,0.5,1])
    plt.yticks([0,0.5,1])
    plt.title("Eigenvector Centrality - Untreated")
    # plt.show()
    plt.savefig(f"{save_path_figures}/{filename}_untreated_{measure}_{cluster1}_vs_{cluster2}.pdf", dpi=300,bbox_inches='tight')
    # fig.savefig(f"{save_path_figures}/{filename}_{measure}_{cluster1}_vs_{cluster2}.svg", dpi=300,bbox_inches='tight')
    
    plt.show()

# %%
sc.pl.umap(oracle.adata, color=['level_2_annotation','SLC26A3','response'], ncols=1)

# %%
cluster_cell="Malignant"
centrality="eigenvector_centrality"
# centrality="betweenness_centrality"
score_caf_treat=links_treat.merged_score.query("cluster==@cluster_cell").sort_values(by=centrality,ascending=False)
score_caf_untreat=links_untreat.merged_score.query("cluster==@cluster_cell").sort_values(by=centrality,ascending=False)

merged_score_caf=pd.merge(score_caf_treat,score_caf_untreat,how='inner',left_index=True,right_index=True,suffixes=('_treat','_untreat'))
merged_score_caf = merged_score_caf.loc[merged_score_caf.index.intersection(tf_list)]
merged_score_caf

# %%
# sns.reset_defaults()
sns.set_theme(style="white",font_scale=1.1,palette="deep")
plt.figure(figsize=(6,6))
sns.scatterplot(markers='o', edgecolor="black", color="#0096FF",s=50,x=merged_score_caf[f"{centrality}_treat"],y=merged_score_caf[f"{centrality}_untreat"])
# plt.figure.fra
# Calculate the 99th percentile for eigenvector centrality in treated and untreated conditions
percentile_treat = np.percentile(merged_score_caf[f"{centrality}_treat"], 90)
percentile_untreat = np.percentile(merged_score_caf[f"{centrality}_untreat"], 90)

# Identify genes above the 99th percentile in either treated or untreated conditions
top_genes_treat = merged_score_caf[merged_score_caf[f"{centrality}_treat"] > percentile_treat].index
top_genes_untreat = merged_score_caf[merged_score_caf[f"{centrality}_untreat"] > percentile_untreat].index
top_genes = np.union1d(top_genes_treat, top_genes_untreat)

# Annotate the scatter plot with gene names
texts_ll=list()
for gene in top_genes:
    x = merged_score_caf.loc[gene, f"{centrality}_treat"]
    y = merged_score_caf.loc[gene, f"{centrality}_untreat"]
    texts_ll.append(plt.text(x, y, gene, fontsize=10))

adjust_text.adjust_text(texts=texts_ll, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))

# plt.xlabel('Treated')
# plt.xlabel(oracle_treat.adata.obs["merged_response"].unique()[0])
# plt.ylabel(oracle_untreat.adata.obs["merged_response"].unique()[0])
plt.xlabel('Treated')
plt.ylabel('Untreated')
plt.title(f'{centrality} - ({cluster_cell})')
plt.xlim(right=max(merged_score_caf[f"{centrality}_treat"])*1.1)
plt.ylim(top=max(merged_score_caf[f"{centrality}_untreat"])*1.1)
plt.xticks([0,0.5,1])
plt.yticks([0,0.5,1])
# Add diagonal line
plt.plot([0, max(merged_score_caf[f"{centrality}_treat"].max()*1.1, merged_score_caf[f"{centrality}_untreat"].max()*1.1)], 
         [0, max(merged_score_caf[f"{centrality}_treat"].max()*1.1, merged_score_caf[f"{centrality}_untreat"].max()*1.1)], 
         'r--', alpha=0.5, linewidth=1)
sns.despine()
plt.grid(False)
plt.savefig(f"{save_path_figures}/{filename}_{centrality}_{cluster_cell}.pdf", dpi=300,bbox_inches='tight')
plt.savefig(f"{save_path_figures}/{filename}_{centrality}_{cluster_cell}.svg", dpi=300,bbox_inches='tight')
# plt.savefig(f"{save_path_figures}/{filename}_eigenvector_centrality_{cluster_cell}.svg", dpi=300,bbox_inches='tight')
plt.show()
# %% [markdown]
# ## Save PDFs results
# %%
print(f"{filename}_")

# %%
# subprocess.run(f"pdfunite {save_path_figures}{filename}*betweenness_centrality*.pdf {filename}_all_betweenness_centrality_out.pdf", shell=True)
# subprocess.run(f"pdfunite {save_path_figures}{filename}*_umap_*.pdf {filename}_all_umap_out.pdf", shell=True)

# %%
df_good_malignant = links_treat.merged_score.query("cluster == 'Malignant'").sort_values("eigenvector_centrality", ascending=False)
df_bad_malignant  = links_untreat.merged_score.query("cluster == 'Malignant'").sort_values("eigenvector_centrality", ascending=False)

# ---- Extract TFs from good (MM/MR) and bad (U/PR) malignant sets ----
# Full ordered TF lists (all TFs appearing in malignant cluster for each condition)
good_malignant_tfs = [g for g in df_good_malignant.index if g in tf_list]
bad_malignant_tfs  = [g for g in df_bad_malignant.index  if g in tf_list]

# Optionally pick top N TFs by eigenvector centrality (set to None to disable)
top_n_tf = 50  # change as needed
if top_n_tf is not None:
    top_good_malignant_tfs = good_malignant_tfs[:top_n_tf]
    top_bad_malignant_tfs  = bad_malignant_tfs[:top_n_tf]
else:
    top_good_malignant_tfs = good_malignant_tfs
    top_bad_malignant_tfs  = bad_malignant_tfs

# Sets for comparison
good_tf_set = set(good_malignant_tfs)
bad_tf_set  = set(bad_malignant_tfs)

intersection_tfs = good_tf_set & bad_tf_set
good_unique_tfs  = good_tf_set - bad_tf_set
bad_unique_tfs   = bad_tf_set - good_tf_set

print("=== MALIGNANT TF COMPARISON (GOOD vs BAD) ===")
print(f"Total TFs (good): {len(good_tf_set)}")
print(f"Total TFs (bad): {len(bad_tf_set)}")
print(f"Intersection TFs: {len(intersection_tfs)}")
print(f"Good-unique TFs: {len(good_unique_tfs)}")
print(f"Bad-unique TFs: {len(bad_unique_tfs)}")

print("\nIntersection TFs (sorted):")
print(sorted(intersection_tfs))
print("\nGood-unique TFs (sorted):")
print(sorted(good_unique_tfs))
print("\nBad-unique TFs (sorted):")
print(sorted(bad_unique_tfs))

# Save results
tf_comp_dict = {
    "good_all": sorted(good_tf_set),
    "bad_all": sorted(bad_tf_set),
    "intersection": sorted(intersection_tfs),
    "good_unique": sorted(good_unique_tfs),
    "bad_unique": sorted(bad_unique_tfs),
    f"good_top{top_n_tf}": top_good_malignant_tfs,
    f"bad_top{top_n_tf}": top_bad_malignant_tfs,
}

df_tf_comp = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in tf_comp_dict.items()]))
df_tf_comp.to_csv(f"{save_path_metadata}/malignant_good_vs_bad_TF_lists.csv", index=False)

# (Optional) write simple text files
with open(f"{save_path_metadata}/malignant_good_TFs.txt", "w") as f:
    for g in good_malignant_tfs: f.write(f"{g}\n")
with open(f"{save_path_metadata}/malignant_bad_TFs.txt", "w") as f:
    for g in bad_malignant_tfs: f.write(f"{g}\n")
with open(f"{save_path_metadata}/malignant_intersection_TFs.txt", "w") as f:
    for g in sorted(intersection_tfs): f.write(f"{g}\n")
with open(f"{save_path_metadata}/malignant_good_unique_TFs.txt", "w") as f:
    for g in sorted(good_unique_tfs): f.write(f"{g}\n")
with open(f"{save_path_metadata}/malignant_bad_unique_TFs.txt", "w") as f:
    for g in sorted(bad_unique_tfs): f.write(f"{g}\n")
    # === Functional annotation & clustering of TFs (automatic online retrieval + text clustering) ===
    # This block:
    #  1. Collects malignant TF lists (good/bad/intersection) already computed above
    #  2. Queries MyGene.info for functional summaries and GO BP/MF terms
    #  3. Builds a TF -> text corpus (summary + GO terms)
    #  4. Vectorizes with TF-IDF, selects optimal k via silhouette, clusters (KMeans)
    #  5. Derives human-readable cluster labels from top TF-IDF terms
    #  6. Saves results and keyword summaries
# %%
def get_gene_annotations(gene_list, batch_size=100):
    """Query MyGene.info for functional annotations"""
    if not _MG_AVAILABLE:
        print("MyGene not available, using empty annotations")
        return {gene: {"summary": "", "go_bp": [], "go_mf": []} for gene in gene_list}
    
    mg = MyGeneInfo()
    annotations = {}
    
    # Process genes in batches to avoid API limits
    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i+batch_size]
        try:
            results = mg.querymany(batch, scopes='symbol', 
                                 fields='summary,go.BP.term,go.MF.term', 
                                 species='human')
            
            for result in results:
                if 'query' in result:
                    gene = result['query']
                    annotations[gene] = {
                        "summary": result.get('summary', ''),
                        "go_bp": [],
                        "go_mf": []
                    }
                    
                    # Extract GO terms
                    if 'go' in result:
                        if 'BP' in result['go']:
                            bp_terms = result['go']['BP']
                            if isinstance(bp_terms, list):
                                annotations[gene]["go_bp"] = [term.get('term', '') for term in bp_terms if 'term' in term]
                            elif isinstance(bp_terms, dict) and 'term' in bp_terms:
                                annotations[gene]["go_bp"] = [bp_terms['term']]
                        
                        if 'MF' in result['go']:
                            mf_terms = result['go']['MF']
                            if isinstance(mf_terms, list):
                                annotations[gene]["go_mf"] = [term.get('term', '') for term in mf_terms if 'term' in term]
                            elif isinstance(mf_terms, dict) and 'term' in mf_terms:
                                annotations[gene]["go_mf"] = [mf_terms['term']]
        except Exception as e:
            print(f"Error querying batch {i//batch_size + 1}: {e}")
            # Fill with empty annotations for failed batch
            for gene in batch:
                if gene not in annotations:
                    annotations[gene] = {"summary": "", "go_bp": [], "go_mf": []}
    
    return annotations
# %%
# === Functional annotation & clustering of TFs (automatic online retrieval + text clustering) ===
# This block:
#  1. Collects malignant TF lists (good/bad/intersection) already computed above
#  2. Queries MyGene.info for functional summaries and GO BP/MF terms
#  3. Builds a TF -> text corpus (summary + GO terms)
#  4. Vectorizes with TF-IDF, selects optimal k via silhouette, clusters (KMeans)
#  5. Derives human-readable cluster labels from top TF-IDF terms
#  6. Saves results and keyword summaries

print("=== FUNCTIONAL ANNOTATION & CLUSTERING OF TFs ===\n")

def get_gene_annotations(gene_list, batch_size=100):
    """Query MyGene.info for functional annotations using MSigDB Hallmark pathways"""
    if not _MG_AVAILABLE:
        print("MyGene not available, using empty annotations")
        return {gene: {"summary": "", "hallmark_pathways": []} for gene in gene_list}
    
    mg = MyGeneInfo()
    annotations = {}
    
    # Process genes in batches to avoid API limits
    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i+batch_size]
        try:
            results = mg.querymany(batch, scopes='symbol', 
                                 fields='summary,pathway.hallmark', 
                                 species='human')
            
            for result in results:
                if 'query' in result:
                    gene = result['query']
                    annotations[gene] = {
                        "summary": result.get('summary', ''),
                        "hallmark_pathways": []
                    }
                    
                    # Extract Hallmark pathway terms
                    if 'pathway' in result and 'hallmark' in result['pathway']:
                        hallmark_terms = result['pathway']['hallmark']
                        if isinstance(hallmark_terms, list):
                            annotations[gene]["hallmark_pathways"] = [term.get('name', '') for term in hallmark_terms if 'name' in term]
                        elif isinstance(hallmark_terms, dict) and 'name' in hallmark_terms:
                            annotations[gene]["hallmark_pathways"] = [hallmark_terms['name']]
        except Exception as e:
            print(f"Error querying batch {i//batch_size + 1}: {e}")
            # Fill with empty annotations for failed batch
            for gene in batch:
                if gene not in annotations:
                    annotations[gene] = {"summary": "", "hallmark_pathways": []}
    
    return annotations

def build_text_corpus(annotations):
    """Build text corpus for each gene from annotations"""
    corpus = {}
    for gene, annot in annotations.items():
        text_parts = []
        
        # Add summary
        if annot["summary"]:
            text_parts.append(annot["summary"])
        
        # Add Hallmark pathway terms
        if annot["hallmark_pathways"]:
            text_parts.extend(annot["hallmark_pathways"])
        
        # Join all text
        corpus[gene] = " ".join(text_parts)
    
    return corpus

def find_optimal_clusters(vectors, gene_names, max_k=10):
    """Find optimal number of clusters using silhouette analysis"""
    if len(gene_names) < 3:
        return 1, None, []
    
    silhouette_scores = []
    max_k = min(max_k, len(gene_names) - 1)
    
    for k in range(2, max_k + 1):
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(vectors)
        silhouette_avg = silhouette_score(vectors, cluster_labels)
        silhouette_scores.append((k, silhouette_avg))
    
    if not silhouette_scores:
        return 1, None, []
    
    # Find k with highest silhouette score
    optimal_k = max(silhouette_scores, key=lambda x: x[1])[0]
    
    # Fit final model
    kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
    final_labels = kmeans.fit_predict(vectors)
    
    return optimal_k, kmeans, final_labels

def get_cluster_keywords(vectorizer, kmeans, cluster_labels, gene_names, top_n=5):
    """Extract top keywords for each cluster"""
    cluster_keywords = {}
    feature_names = vectorizer.get_feature_names_out()
    
    for cluster_id in range(kmeans.n_clusters):
        # Get center of this cluster
        center = kmeans.cluster_centers_[cluster_id]
        
        # Get top features (keywords) for this cluster
        top_indices = center.argsort()[-top_n:][::-1]
        top_keywords = [feature_names[i] for i in top_indices]
        
        # Get genes in this cluster
        cluster_genes = [gene_names[i] for i in range(len(gene_names)) if cluster_labels[i] == cluster_id]
        
        cluster_keywords[cluster_id] = {
            'keywords': top_keywords,
            'genes': cluster_genes,
            'label': ' & '.join(top_keywords[:3])  # Human-readable label
        }
    
    return cluster_keywords

# Collect all unique TFs from good/bad/intersection lists
all_malignant_tfs = list(set(good_malignant_tfs + bad_malignant_tfs))
print(f"Total unique malignant TFs to annotate: {len(all_malignant_tfs)}")

# Get annotations for all TFs
print("Querying MyGene.info for MSigDB Hallmark pathway annotations...")
tf_annotations = get_gene_annotations(all_malignant_tfs)

# Build text corpus
print("Building text corpus from annotations...")
tf_corpus = build_text_corpus(tf_annotations)

# Filter out TFs with empty corpus
valid_tfs = [tf for tf, text in tf_corpus.items() if text.strip()]
valid_corpus = [tf_corpus[tf] for tf in valid_tfs]

print(f"TFs with valid annotations: {len(valid_tfs)}")

if len(valid_tfs) >= 3:
    # Vectorize text using TF-IDF
    print("Vectorizing text with TF-IDF...")
    vectorizer = TfidfVectorizer(max_features=500, stop_words='english', 
                               min_df=2, max_df=0.8, ngram_range=(1,2))
    tfidf_matrix = vectorizer.fit_transform(valid_corpus)
    
    # Find optimal clusters
    print("Finding optimal number of clusters...")
    optimal_k, kmeans_model, cluster_assignments = find_optimal_clusters(
        tfidf_matrix.toarray(), valid_tfs, max_k=min(10, len(valid_tfs)//3))
    
    print(f"Optimal number of clusters: {optimal_k}")
    
    if optimal_k > 1 and kmeans_model is not None:
        # Get cluster keywords and assignments
        cluster_info = get_cluster_keywords(vectorizer, kmeans_model, cluster_assignments, valid_tfs)
        
        # Create comprehensive results DataFrame
        tf_cluster_results = []
        for tf in all_malignant_tfs:
            if tf in valid_tfs:
                cluster_id = cluster_assignments[valid_tfs.index(tf)]
                cluster_label = cluster_info[cluster_id]['label']
                keywords = ', '.join(cluster_info[cluster_id]['keywords'])
            else:
                cluster_id = -1  # No annotation
                cluster_label = "No annotation"
                keywords = ""
            
            # Determine TF category
            tf_category = "Intersection"
            if tf in good_unique_tfs:
                tf_category = "Good unique"
            elif tf in bad_unique_tfs:
                tf_category = "Bad unique"
            
            tf_cluster_results.append({
                'TF': tf,
                'Category': tf_category,
                'Cluster_ID': cluster_id,
                'Cluster_Label': cluster_label,
                'Top_Keywords': keywords,
                'Summary': tf_annotations.get(tf, {}).get('summary', ''),
                'Hallmark_Pathways': '; '.join(tf_annotations.get(tf, {}).get('hallmark_pathways', []))
            })
        
        df_tf_clusters = pd.DataFrame(tf_cluster_results)
        
        # Save results
        df_tf_clusters.to_csv(f"{save_path_metadata}/malignant_TF_hallmark_functional_clusters.csv", index=False)
        
        # Print cluster summary
        print("\n=== HALLMARK PATHWAY FUNCTIONAL CLUSTER SUMMARY ===")
        for cluster_id, info in cluster_info.items():
            print(f"\nCluster {cluster_id}: {info['label']}")
            print(f"Keywords: {', '.join(info['keywords'])}")
            print(f"Genes ({len(info['genes'])}): {', '.join(info['genes'])}")
        
        # Create visualization
        fig, axes = plt.subplots(2, 2, figsize=(20, 16))
        
        # Plot 1: Cluster composition by TF category
        cluster_category_counts = df_tf_clusters.groupby(['Cluster_Label', 'Category']).size().unstack(fill_value=0)
        cluster_category_counts.plot(kind='bar', stacked=True, ax=axes[0,0], 
                                   color=['#9D4EDD', '#0096FF', '#fb3310'])
        axes[0,0].set_title('TF Category Distribution by Hallmark Functional Cluster')
        axes[0,0].set_ylabel('Number of TFs')
        axes[0,0].legend(title='TF Category')
        axes[0,0].tick_params(axis='x', rotation=45)
        
        # Plot 2: Cluster sizes
        cluster_sizes = df_tf_clusters['Cluster_Label'].value_counts()
    axes[0,1].bar(range(len(cluster_sizes)), cluster_sizes.values, 
             color=plt.cm.inferno(np.linspace(0, 1, len(cluster_sizes))))
        axes[0,1].set_xticks(range(len(cluster_sizes)))
        axes[0,1].set_xticklabels(cluster_sizes.index, rotation=45, ha='right')
        axes[0,1].set_title('Hallmark Functional Cluster Sizes')
        axes[0,1].set_ylabel('Number of TFs')
        
        # Plot 3: Category proportions
        category_counts = df_tf_clusters['Category'].value_counts()
        axes[1,0].pie(category_counts.values, labels=category_counts.index, autopct='%1.1f%%',
                     colors=['#9D4EDD', '#0096FF', '#fb3310'])
        axes[1,0].set_title('Overall TF Category Distribution')
        
        # Plot 4: Top keywords word cloud style visualization
        all_keywords = []
        for info in cluster_info.values():
            all_keywords.extend(info['keywords'])
        keyword_counts = Counter(all_keywords)
        
        if keyword_counts:
            top_keywords = dict(keyword_counts.most_common(15))
            axes[1,1].barh(range(len(top_keywords)), list(top_keywords.values()))
            axes[1,1].set_yticks(range(len(top_keywords)))
            axes[1,1].set_yticklabels(list(top_keywords.keys()))
            axes[1,1].set_title('Most Frequent Hallmark Keywords Across Clusters')
            axes[1,1].set_xlabel('Frequency')
        
        plt.tight_layout()
        plt.savefig(f"{save_path_figures}/malignant_TF_hallmark_functional_clustering_analysis.pdf", 
                   bbox_inches='tight', dpi=300)
        plt.savefig(f"{save_path_figures}/malignant_TF_hallmark_functional_clustering_analysis.svg", 
                   bbox_inches='tight', dpi=300)
        plt.show()
        
        # Create detailed cluster heatmap showing TF categories within each functional cluster
        cluster_tf_matrix = df_tf_clusters.pivot_table(
            index='TF', columns='Cluster_Label', values='Cluster_ID', 
            fill_value=-1, aggfunc='first')
        cluster_tf_matrix = cluster_tf_matrix.replace(-1, np.nan)
        
        # Add TF category as row colors
        tf_category_colors = {'Good unique': '#0096FF', 'Bad unique': '#fb3310', 'Intersection': '#9D4EDD'}
        row_colors = [tf_category_colors.get(df_tf_clusters[df_tf_clusters['TF']==tf]['Category'].iloc[0], 'gray') 
                     for tf in cluster_tf_matrix.index]
        
        # Create clustermap
    g = sns.clustermap(cluster_tf_matrix.notna().astype(int), 
              row_colors=row_colors, cmap='inferno',
              figsize=(15, 20), row_cluster=True, col_cluster=False,
              linewidths=0.5, cbar_kws={'label': 'Cluster Membership'})
        
        # Add legend for row colors
        handles = [Patch(facecolor=color, label=category) 
                  for category, color in tf_category_colors.items()]
        g.ax_heatmap.legend(handles=handles, title='TF Category', 
                          bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.savefig(f"{save_path_figures}/malignant_TF_hallmark_functional_cluster_membership_heatmap.pdf", 
                   bbox_inches='tight', dpi=300)
        plt.show()
        
        print(f"\n=== RESULTS SAVED ===")
        print(f"Detailed annotations: {save_path_metadata}/malignant_TF_hallmark_functional_clusters.csv")
        print(f"Visualization: {save_path_figures}/malignant_TF_hallmark_functional_clustering_analysis.pdf")
        
    else:
        print("Unable to perform meaningful clustering (insufficient data or clusters)")
else:
    print("Insufficient TFs with valid annotations for clustering")
# %%
print("=== FUNCTIONAL ANNOTATION & CLUSTERING OF TFs ===\n")

def get_gene_annotations(gene_list, batch_size=100):
    """Query MyGene.info for functional annotations"""
    if not _MG_AVAILABLE:
        print("MyGene not available, using empty annotations")
        return {gene: {"summary": "", "go_bp": [], "go_mf": []} for gene in gene_list}
    
    mg = MyGeneInfo()
    annotations = {}
    
    # Process genes in batches to avoid API limits
    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i+batch_size]
        try:
            results = mg.querymany(batch, scopes='symbol', 
                                 fields='summary,go.BP.term,go.MF.term', 
                                 species='human')
            
            for result in results:
                if 'query' in result:
                    gene = result['query']
                    annotations[gene] = {
                        "summary": result.get('summary', ''),
                        "go_bp": [],
                        "go_mf": []
                    }
                    
                    # Extract GO terms
                    if 'go' in result:
                        if 'BP' in result['go']:
                            bp_terms = result['go']['BP']
                            if isinstance(bp_terms, list):
                                annotations[gene]["go_bp"] = [term.get('term', '') for term in bp_terms if 'term' in term]
                            elif isinstance(bp_terms, dict) and 'term' in bp_terms:
                                annotations[gene]["go_bp"] = [bp_terms['term']]
                        
                        if 'MF' in result['go']:
                            mf_terms = result['go']['MF']
                            if isinstance(mf_terms, list):
                                annotations[gene]["go_mf"] = [term.get('term', '') for term in mf_terms if 'term' in term]
                            elif isinstance(mf_terms, dict) and 'term' in mf_terms:
                                annotations[gene]["go_mf"] = [mf_terms['term']]
        except Exception as e:
            print(f"Error querying batch {i//batch_size + 1}: {e}")
            # Fill with empty annotations for failed batch
            for gene in batch:
                if gene not in annotations:
                    annotations[gene] = {"summary": "", "go_bp": [], "go_mf": []}
    
    return annotations

def build_text_corpus(annotations):
    """Build text corpus for each gene from annotations"""
    corpus = {}
    for gene, annot in annotations.items():
        text_parts = []
        
        # Add summary
        if annot["summary"]:
            text_parts.append(annot["summary"])
        
        # Add GO terms
        if annot["go_bp"]:
            text_parts.extend(annot["go_bp"])
        if annot["go_mf"]:
            text_parts.extend(annot["go_mf"])
        
        # Join all text
        corpus[gene] = " ".join(text_parts)
    
    return corpus

def find_optimal_clusters(vectors, gene_names, max_k=10):
    """Find optimal number of clusters using silhouette analysis"""
    if len(gene_names) < 3:
        return 1, None, []
    
    silhouette_scores = []
    max_k = min(max_k, len(gene_names) - 1)
    
    for k in range(2, max_k + 1):
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(vectors)
        silhouette_avg = silhouette_score(vectors, cluster_labels)
        silhouette_scores.append((k, silhouette_avg))
    
    if not silhouette_scores:
        return 1, None, []
    
    # Find k with highest silhouette score
    optimal_k = max(silhouette_scores, key=lambda x: x[1])[0]
    
    # Fit final model
    kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
    final_labels = kmeans.fit_predict(vectors)
    
    return optimal_k, kmeans, final_labels

def get_cluster_keywords(vectorizer, kmeans, cluster_labels, gene_names, top_n=5):
    """Extract top keywords for each cluster"""
    cluster_keywords = {}
    feature_names = vectorizer.get_feature_names_out()
    
    for cluster_id in range(kmeans.n_clusters):
        # Get center of this cluster
        center = kmeans.cluster_centers_[cluster_id]
        
        # Get top features (keywords) for this cluster
        top_indices = center.argsort()[-top_n:][::-1]
        top_keywords = [feature_names[i] for i in top_indices]
        
        # Get genes in this cluster
        cluster_genes = [gene_names[i] for i in range(len(gene_names)) if cluster_labels[i] == cluster_id]
        
        cluster_keywords[cluster_id] = {
            'keywords': top_keywords,
            'genes': cluster_genes,
            'label': ' & '.join(top_keywords[:3])  # Human-readable label
        }
    
    return cluster_keywords

# Collect all unique TFs from good/bad/intersection lists
all_malignant_tfs = list(set(good_malignant_tfs + bad_malignant_tfs))
print(f"Total unique malignant TFs to annotate: {len(all_malignant_tfs)}")

# Get annotations for all TFs
print("Querying MyGene.info for functional annotations...")
tf_annotations = get_gene_annotations(all_malignant_tfs)

# Build text corpus
print("Building text corpus from annotations...")
tf_corpus = build_text_corpus(tf_annotations)

# Filter out TFs with empty corpus
valid_tfs = [tf for tf, text in tf_corpus.items() if text.strip()]
valid_corpus = [tf_corpus[tf] for tf in valid_tfs]

print(f"TFs with valid annotations: {len(valid_tfs)}")

if len(valid_tfs) >= 3:
    # Vectorize text using TF-IDF
    print("Vectorizing text with TF-IDF...")
    vectorizer = TfidfVectorizer(max_features=500, stop_words='english', 
                               min_df=2, max_df=0.8, ngram_range=(1,2))
    tfidf_matrix = vectorizer.fit_transform(valid_corpus)
    
    # Find optimal clusters
    print("Finding optimal number of clusters...")
    optimal_k, kmeans_model, cluster_assignments = find_optimal_clusters(
        tfidf_matrix.toarray(), valid_tfs, max_k=min(10, len(valid_tfs)//3))
    
    print(f"Optimal number of clusters: {optimal_k}")
    
    if optimal_k > 1 and kmeans_model is not None:
        # Get cluster keywords and assignments
        cluster_info = get_cluster_keywords(vectorizer, kmeans_model, cluster_assignments, valid_tfs)
        
        # Create comprehensive results DataFrame
        tf_cluster_results = []
        for tf in all_malignant_tfs:
            if tf in valid_tfs:
                cluster_id = cluster_assignments[valid_tfs.index(tf)]
                cluster_label = cluster_info[cluster_id]['label']
                keywords = ', '.join(cluster_info[cluster_id]['keywords'])
            else:
                cluster_id = -1  # No annotation
                cluster_label = "No annotation"
                keywords = ""
            
            # Determine TF category
            tf_category = "Intersection"
            if tf in good_unique_tfs:
                tf_category = "Good unique"
            elif tf in bad_unique_tfs:
                tf_category = "Bad unique"
            
            tf_cluster_results.append({
                'TF': tf,
                'Category': tf_category,
                'Cluster_ID': cluster_id,
                'Cluster_Label': cluster_label,
                'Top_Keywords': keywords,
                'Summary': tf_annotations.get(tf, {}).get('summary', ''),
                'GO_BP_Terms': '; '.join(tf_annotations.get(tf, {}).get('go_bp', [])),
                'GO_MF_Terms': '; '.join(tf_annotations.get(tf, {}).get('go_mf', []))
            })
        
        df_tf_clusters = pd.DataFrame(tf_cluster_results)
        
        # Save results
        df_tf_clusters.to_csv(f"{save_path_metadata}/malignant_TF_functional_clusters.csv", index=False)
        
        # Print cluster summary
        print("\n=== FUNCTIONAL CLUSTER SUMMARY ===")
        for cluster_id, info in cluster_info.items():
            print(f"\nCluster {cluster_id}: {info['label']}")
            print(f"Keywords: {', '.join(info['keywords'])}")
            print(f"Genes ({len(info['genes'])}): {', '.join(info['genes'])}")
        
        # Create visualization
        fig, axes = plt.subplots(2, 2, figsize=(20, 16))
        
        # Plot 1: Cluster composition by TF category
        cluster_category_counts = df_tf_clusters.groupby(['Cluster_Label', 'Category']).size().unstack(fill_value=0)
        cluster_category_counts.plot(kind='bar', stacked=True, ax=axes[0,0], 
                                   color=['#9D4EDD', '#0096FF', '#fb3310'])
        axes[0,0].set_title('TF Category Distribution by Functional Cluster')
        axes[0,0].set_ylabel('Number of TFs')
        axes[0,0].legend(title='TF Category')
        axes[0,0].tick_params(axis='x', rotation=45)
        
        # Plot 2: Cluster sizes
        cluster_sizes = df_tf_clusters['Cluster_Label'].value_counts()
        axes[0,1].bar(range(len(cluster_sizes)), cluster_sizes.values, 
                     color=plt.cm.tab10(np.linspace(0, 1, len(cluster_sizes))))
        axes[0,1].set_xticks(range(len(cluster_sizes)))
        axes[0,1].set_xticklabels(cluster_sizes.index, rotation=45, ha='right')
        axes[0,1].set_title('Functional Cluster Sizes')
        axes[0,1].set_ylabel('Number of TFs')
        
        # Plot 3: Category proportions
        category_counts = df_tf_clusters['Category'].value_counts()
        axes[1,0].pie(category_counts.values, labels=category_counts.index, autopct='%1.1f%%',
                     colors=['#9D4EDD', '#0096FF', '#fb3310'])
        axes[1,0].set_title('Overall TF Category Distribution')
        
        # Plot 4: Top keywords word cloud style visualization
        all_keywords = []
        for info in cluster_info.values():
            all_keywords.extend(info['keywords'])
        keyword_counts = Counter(all_keywords)
        
        if keyword_counts:
            top_keywords = dict(keyword_counts.most_common(15))
            axes[1,1].barh(range(len(top_keywords)), list(top_keywords.values()))
            axes[1,1].set_yticks(range(len(top_keywords)))
            axes[1,1].set_yticklabels(list(top_keywords.keys()))
            axes[1,1].set_title('Most Frequent Keywords Across Clusters')
            axes[1,1].set_xlabel('Frequency')
        
        plt.tight_layout()
        plt.savefig(f"{save_path_figures}/malignant_TF_functional_clustering_analysis.pdf", 
                   bbox_inches='tight', dpi=300)
        plt.savefig(f"{save_path_figures}/malignant_TF_functional_clustering_analysis.svg", 
                   bbox_inches='tight', dpi=300)
        plt.show()
        
        # Create detailed cluster heatmap showing TF categories within each functional cluster
        cluster_tf_matrix = df_tf_clusters.pivot_table(
            index='TF', columns='Cluster_Label', values='Cluster_ID', 
            fill_value=-1, aggfunc='first')
        cluster_tf_matrix = cluster_tf_matrix.replace(-1, np.nan)
        
        # Add TF category as row colors
        tf_category_colors = {'Good unique': '#0096FF', 'Bad unique': '#fb3310', 'Intersection': '#9D4EDD'}
        row_colors = [tf_category_colors.get(df_tf_clusters[df_tf_clusters['TF']==tf]['Category'].iloc[0], 'gray') 
                     for tf in cluster_tf_matrix.index]
        
        # Create clustermap
        g = sns.clustermap(cluster_tf_matrix.notna().astype(int), 
                          row_colors=row_colors, cmap='Blues',
                          figsize=(15, 20), row_cluster=True, col_cluster=False,
                          linewidths=0.5, cbar_kws={'label': 'Cluster Membership'})
        
        # Add legend for row colors
        handles = [Patch(facecolor=color, label=category) 
                  for category, color in tf_category_colors.items()]
        g.ax_heatmap.legend(handles=handles, title='TF Category', 
                          bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.savefig(f"{save_path_figures}/malignant_TF_functional_cluster_membership_heatmap.pdf", 
                   bbox_inches='tight', dpi=300)
        plt.show()
        
        print(f"\n=== RESULTS SAVED ===")
        print(f"Detailed annotations: {save_path_metadata}/malignant_TF_functional_clusters.csv")
        print(f"Visualization: {save_path_figures}/malignant_TF_functional_clustering_analysis.pdf")
        
    else:
        print("Unable to perform meaningful clustering (insufficient data or clusters)")
else:
    print("Insufficient TFs with valid annotations for clustering")
# %%
# (Continue with downstream analysis if needed)

# Find intersection (genes high in both)
intersection_genes = (good_tf_set & bad_tf_set).intersection(tf_list)

# Find genes unique to good response
good_unique_genes = (good_tf_set - bad_tf_set).intersection(tf_list)

# Find genes unique to bad response
bad_unique_genes = (bad_tf_set - good_tf_set).intersection(tf_list)

# Print results
# print(f"Top {top_n_genes} genes analysis:")
print(f"Good response genes: {len(top_good_genes)}")
print(f"Bad response genes: {len(top_bad_genes)}")
print(f"Intersection (common): {len(intersection_genes)}")
print(f"Good response unique: {len(good_unique_genes)}")
print(f"Bad response unique: {len(bad_unique_genes)}")

print("\n=== INTERSECTION GENES (Common to both) ===")
print(sorted(list(intersection_genes)))

print("\n=== GOOD RESPONSE UNIQUE GENES ===")
print(sorted(list(good_unique_genes)))

print("\n=== BAD RESPONSE UNIQUE GENES ===")
print(sorted(list(bad_unique_genes)))

# Save results to files
results_dict = {
    'intersection': sorted(list(intersection_genes)),
    'good_unique': sorted(list(good_unique_genes)), 
    'bad_unique': sorted(list(bad_unique_genes))
}

# Save as CSV
df_results = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in results_dict.items()]))
df_results.to_csv(f"{save_path_metadata}/malignant_genes_comparison_good_vs_bad.csv", index=False)

# Create visualization
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot 1: Venn diagram style bar plot
categories = ['Good Unique', 'Intersection', 'Bad Unique']
counts = [len(good_unique_genes), len(intersection_genes), len(bad_unique_genes)]
colors = ['#0096FF', '#9D4EDD', '#fb3310']

axes[0].bar(categories, counts, color=colors)
axes[0].set_title('TF Set Overlap Analysis')
axes[0].set_ylabel('Number of Genes')
for i, count in enumerate(counts):
    axes[0].text(i, count + 1, str(count), ha='center', va='bottom', fontweight='bold')

# Plot 2: Eigenvector centrality comparison for intersection genes
if intersection_genes:
    good_centrality = df_good_malignant.loc[list(intersection_genes), 'eigenvector_centrality']
    bad_centrality = df_bad_malignant.loc[list(intersection_genes), 'eigenvector_centrality']
    
    # Add gene name annotations for intersection genes
    if len(intersection_genes) > 0:
        texts = []
        for gene in intersection_genes:
            x = good_centrality.loc[gene]
            y = bad_centrality.loc[gene]
            texts.append(axes[1].text(x, y, gene, fontsize=8))
        
        # Adjust text positions to avoid overlap
        adjust_text.adjust_text(texts, arrowprops=dict(arrowstyle="->", color='r', lw=0.5), ax=axes[1])
    axes[1].scatter(good_centrality, bad_centrality, alpha=0.7, s=50, color='purple')
    axes[1].plot([0, 1], [0, 1], 'r--', alpha=0.5)
    axes[1].set_xlabel('Good Response Centrality')
    axes[1].set_ylabel('Bad Response Centrality')
    axes[1].set_title('Centrality Comparison\n(Intersection Genes)')

# Plot 3: Top 10 genes from each category
top_display = 10
# Filter to get only TFs from unique genes
df_good_malignant_tfs = df_good_malignant.loc[good_unique_genes]
df_bad_malignant_tfs = df_bad_malignant.loc[bad_unique_genes]

all_unique_tfs = df_good_malignant_tfs.nlargest(top_display,columns=["eigenvector_centrality"]).index.tolist() + df_bad_malignant_tfs.nlargest(top_display,columns=["eigenvector_centrality"]).index.tolist()
gene_categories = ['Good'] * df_good_malignant_tfs.head(top_display).shape[0] + ['Bad'] * df_bad_malignant_tfs.head(top_display).shape[0]

if all_unique_tfs:
    centralities = []
    for gene in all_unique_tfs:
        if gene in df_good_malignant.index:
            centralities.append(df_good_malignant.loc[gene, 'eigenvector_centrality'])
        else:
            centralities.append(df_bad_malignant.loc[gene, 'eigenvector_centrality'])
    
    y_pos = range(len(all_unique_tfs))
    bar_colors = ['#0096FF' if cat == 'Good' else '#fb3310' for cat in gene_categories]
    
    axes[2].barh(y_pos, centralities, color=bar_colors)
    axes[2].set_yticks(y_pos)
    axes[2].set_yticklabels(all_unique_tfs, fontsize=8)
    axes[2].set_xlabel('Eigenvector Centrality')
    axes[2].set_title(f'Top {top_display} Unique TFs')
else:
    all_unique_genes = list(good_unique_genes)[:top_display] + list(bad_unique_genes)[:top_display]
    gene_categories = ['Good'] * min(top_display, len(good_unique_genes)) + ['Bad'] * min(top_display, len(bad_unique_genes))

plt.tight_layout()
plt.savefig(f"{save_path_figures}/malignant_genes_comparison_analysis.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{save_path_figures}/malignant_genes_comparison_analysis.svg", dpi=300, bbox_inches='tight')
plt.show()

# Filter for TFs only
intersection_tfs = intersection_genes & set(tf_list)
good_unique_tfs = good_unique_genes & set(tf_list)
bad_unique_tfs = bad_unique_genes & set(tf_list)

print(f"\n=== TRANSCRIPTION FACTORS ONLY ===")
print(f"Intersection TFs: {len(intersection_tfs)}")
print(f"Good unique TFs: {len(good_unique_tfs)}")  
print(f"Bad unique TFs: {len(bad_unique_tfs)}")

print("\nIntersection TFs:", sorted(list(intersection_tfs)))
print("Good unique TFs:", sorted(list(good_unique_tfs)))
print("Bad unique TFs:", sorted(list(bad_unique_tfs)))
# %%
gene_source="HMGA2"
df_target_genes=links_treat.filtered_links["Malignant"].query("source == @gene_source")
good_links=df_target_genes["target"].unique().tolist()
good_graph=co.network_analysis.draw_network(df_target_genes,return_graph=True)
good_graph.size()
# %%
df_target_genes=links_untreat.filtered_links["Malignant"].query("source == @gene_source")
bad_links=df_target_genes["target"].unique().tolist()
bad_graph=co.network_analysis.draw_network(df_target_genes,return_graph=True)
bad_graph.size()
# %%
good_graph.remove_nodes_from(bad_links)
good_graph.size()
# %%
bad_graph.remove_nodes_from(good_links)
bad_graph.size()
# %%
# Draw the graph
fig,ax=plt.subplots(figsize=(6, 6))
seed=1234
pos = nx.spring_layout(good_graph, k=0.5, iterations=5, seed=seed)
nx.draw(good_graph, pos, with_labels=True, node_size=100, node_color="skyblue", font_size=10, font_weight="bold", edge_color="gray")
plt.title("Gene Interaction Network with Pathways in Cancer (Good Response)")
# %%
# Draw the graph
fig,ax=plt.subplots(figsize=(6, 6))
seed=1234
pos = nx.spring_layout(bad_graph, k=0.5, iterations=5, seed=seed)
nx.draw(bad_graph, pos, with_labels=True, node_size=400, node_color="skyblue", font_size=10, font_weight="bold", edge_color="gray")
plt.title("Gene Interaction Network with Pathways in Cancer (Bad Response)")
# %%
sc.pl.violin(adata_malignant_original,keys=[gene_source],groupby="response",jitter=0.1,rotation=90)
# %%
sc.pl.dotplot(adata_malignant_original,var_names=[gene_source],groupby="response",use_raw=False)
# %%
# Compare links between good and bad response groups
print("=== COMPARING GOOD vs BAD RESPONSE LINKS ===\n")

# Get malignant cluster links for both groups
cluster_name = "Malignant"
good_links = links_treat.filtered_links[cluster_name]
bad_links = links_untreat.filtered_links[cluster_name]

# Create sets of source-target pairs for comparison
good_edges = set(zip(good_links['source'], good_links['target']))
bad_edges = set(zip(bad_links['source'], bad_links['target']))

# Find common edges (intersection)
common_edges = good_edges & bad_edges

# Find unique edges in good response
good_unique_edges = good_edges - bad_edges

# Find unique edges in bad response
bad_unique_edges = bad_edges - good_edges

print(f"Cluster: {cluster_name}")
print(f"Good response total edges: {len(good_edges)}")
print(f"Bad response total edges: {len(bad_edges)}")
print(f"Common edges: {len(common_edges)}")
print(f"Good unique edges: {len(good_unique_edges)}")
print(f"Bad unique edges: {len(bad_unique_edges)}")

# Analyze source (TF) level differences
good_sources = set(good_links['source'])
bad_sources = set(bad_links['source'])

common_sources = good_sources & bad_sources
good_unique_sources = good_sources - bad_sources
bad_unique_sources = bad_sources - good_sources

print(f"\n=== SOURCE (TF) LEVEL COMPARISON ===")
print(f"Total TFs in good response: {len(good_sources)}")
print(f"Total TFs in bad response: {len(bad_sources)}")
print(f"Common TFs: {len(common_sources)}")
print(f"Good unique TFs: {len(good_unique_sources)}")
print(f"Bad unique TFs: {len(bad_unique_sources)}")

# Filter for TFs only
good_unique_tfs = good_unique_sources & set(tf_list)
bad_unique_tfs = bad_unique_sources & set(tf_list)
common_tfs = common_sources & set(tf_list)

print(f"\nFiltered for Transcription Factors only:")
print(f"Common TFs: {len(common_tfs)}")
print(f"Good unique TFs: {len(good_unique_tfs)}")
print(f"Bad unique TFs: {len(bad_unique_tfs)}")

print(f"\nCommon TFs: {sorted(list(common_tfs))}")
print(f"Good unique TFs: {sorted(list(good_unique_tfs))}")
print(f"Bad unique TFs: {sorted(list(bad_unique_tfs))}")

# Save edge comparison results
edge_comparison_results = {
    'common_edges': [f"{source}->{target}" for source, target in common_edges],
    'good_unique_edges': [f"{source}->{target}" for source, target in good_unique_edges],
    'bad_unique_edges': [f"{source}->{target}" for source, target in bad_unique_edges]
}

# Save TF comparison results
tf_comparison_results = {
    'common_tfs': sorted(list(common_tfs)),
    'good_unique_tfs': sorted(list(good_unique_tfs)),
    'bad_unique_tfs': sorted(list(bad_unique_tfs))
}

# Convert to DataFrame and save
df_edge_results = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in edge_comparison_results.items()]))
df_edge_results.to_csv(f"{save_path_metadata}/{cluster_name}_edge_comparison_good_vs_bad.csv", index=False)

df_tf_results = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in tf_comparison_results.items()]))
df_tf_results.to_csv(f"{save_path_metadata}/{cluster_name}_tf_comparison_good_vs_bad.csv", index=False)

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Edge comparison
edge_categories = ['Common', 'Good Unique', 'Bad Unique']
edge_counts = [len(common_edges), len(good_unique_edges), len(bad_unique_edges)]
edge_colors = ['#9D4EDD', '#0096FF', '#fb3310']

axes[0, 0].bar(edge_categories, edge_counts, color=edge_colors)
axes[0, 0].set_title('Edge Comparison (Source-Target Pairs)')
axes[0, 0].set_ylabel('Number of Edges')
for i, count in enumerate(edge_counts):
    axes[0, 0].text(i, count + max(edge_counts)*0.01, str(count), ha='center', va='bottom', fontweight='bold')

# Plot 2: TF comparison
tf_categories = ['Common TFs', 'Good Unique TFs', 'Bad Unique TFs']
tf_counts = [len(common_tfs), len(good_unique_tfs), len(bad_unique_tfs)]

axes[0, 1].bar(tf_categories, tf_counts, color=edge_colors)
axes[0, 1].set_title('Transcription Factor Comparison')
axes[0, 1].set_ylabel('Number of TFs')
for i, count in enumerate(tf_counts):
    axes[0, 1].text(i, count + max(tf_counts)*0.01, str(count), ha='center', va='bottom', fontweight='bold')

# Plot 3: Target degree distribution for common TFs
if common_tfs:
    common_tf_target_counts = []
    for tf in common_tfs:
        good_targets = len(good_links[good_links['source'] == tf])
        bad_targets = len(bad_links[bad_links['source'] == tf])
        common_tf_target_counts.extend([good_targets, bad_targets])
    
    axes[1, 0].hist([good_targets for tf in common_tfs for good_targets in [len(good_links[good_links['source'] == tf])]], 
                    bins=20, alpha=0.7, color='#0096FF', label='Good Response')
    axes[1, 0].hist([bad_targets for tf in common_tfs for bad_targets in [len(bad_links[bad_links['source'] == tf])]], 
                    bins=20, alpha=0.7, color='#fb3310', label='Bad Response')
    axes[1, 0].set_title('Target Count Distribution for Common TFs')
    axes[1, 0].set_xlabel('Number of Targets per TF')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].legend()

# Plot 4: Top unique TFs by degree
top_n = min(10, max(len(good_unique_tfs), len(bad_unique_tfs)))
if good_unique_tfs or bad_unique_tfs:
    # Get target counts for unique TFs
    good_unique_degrees = [(tf, len(good_links[good_links['source'] == tf])) for tf in good_unique_tfs]
    bad_unique_degrees = [(tf, len(bad_links[bad_links['source'] == tf])) for tf in bad_unique_tfs]
    
    # Sort by degree and take top N
    good_unique_degrees = sorted(good_unique_degrees, key=lambda x: x[1], reverse=True)[:top_n]
    bad_unique_degrees = sorted(bad_unique_degrees, key=lambda x: x[1], reverse=True)[:top_n]
    
    # Combine for plotting
    all_unique_tfs = [item[0] for item in good_unique_degrees] + [item[0] for item in bad_unique_degrees]
    all_degrees = [item[1] for item in good_unique_degrees] + [item[1] for item in bad_unique_degrees]
    all_colors = ['#0096FF'] * len(good_unique_degrees) + ['#fb3310'] * len(bad_unique_degrees)
    
    y_pos = range(len(all_unique_tfs))
    axes[1, 1].barh(y_pos, all_degrees, color=all_colors)
    axes[1, 1].set_yticks(y_pos)
    axes[1, 1].set_yticklabels(all_unique_tfs, fontsize=8)
    axes[1, 1].set_xlabel('Number of Target Genes')
    axes[1, 1].set_title(f'Top {top_n} Unique TFs by Degree')
    
    # Add legend
    legend_elements = [Patch(facecolor='#0096FF', label='Good Response'),
                      Patch(facecolor='#fb3310', label='Bad Response')]
    axes[1, 1].legend(handles=legend_elements)

plt.tight_layout()
plt.savefig(f"{save_path_figures}/{cluster_name}_links_comparison_analysis.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{save_path_figures}/{cluster_name}_links_comparison_analysis.svg", dpi=300, bbox_inches='tight')
plt.show()
# %%
## gseapy already imported at top

tf_list_from_malignant_response = {"Untreated": [], "Poor response": [],"Minimal response": [], "Moderate response": []}
for response in links_malignant.merged_score["cluster"].unique():
    mg_score_tmp = links_malignant.merged_score[links_malignant.merged_score["cluster"] == response]
    tf_list_from_malignant_response[response] = mg_score_tmp.index.unique().tolist()
    tf_list_from_malignant_response[response] = [tf for tf in tf_list_from_malignant_response[response] if tf in tf_list]

gene_sets = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "GO_Biological_Process_2025",
    "GO_Molecular_Function_2023",
    "GO_Molecular_Function_2025",
    "MSigDB_Hallmark_2020",
    "Reactome_2022",
    "Reactome_Pathways_2024"
]

# Define colors for each gene set
gene_set_colors = {
    "KEGG_2021_Human": "#1f77b4",
    "GO_Biological_Process_2023": "#ff7f0e",
    "GO_Biological_Process_2025": "#2ca02c", 
    "GO_Molecular_Function_2023": "#d62728",
    "GO_Molecular_Function_2025": "#9467bd",
    "MSigDB_Hallmark_2020": "#8c564b",
    "Reactome_2022": "#e377c2",
    "Reactome_Pathways_2024": "#7f7f7f"
}

res_whole_tf=pd.DataFrame()
for response in tf_list_from_malignant_response:
    res=gp.enrichr(gene_list=tf_list_from_malignant_response[response],
                   gene_sets=gene_sets,
                   cutoff=0.05)
    res.results.to_csv(f"{save_path_metadata}/enrichment_malignant_{response}_wholeTF.csv",sep=",",index=False)
    res_filtered = res.results[res.results["Adjusted P-value"] < 0.05]
    res_filtered["response"] = response
    res_whole_tf = pd.concat([res_whole_tf,res_filtered], ignore_index=True)

    barplot(
    res_filtered,
    column="Adjusted P-value",
    group="Gene_set",
    size=10,
    top_term=10,
    title=f"Enrichment whole TF Analysis",
    figsize=(10, 20),
    color=gene_set_colors,
    ofname=f"{save_path_figures}/enrichment_malignant_{response}_wholeTF.pdf"
    )
# %%
# Create TF classification based on enriched terms (not response)
tf_term_classification = {}
term_tf_dict = {}
gene_source="wholeTF"

# Process both good and poor response results if they exist
all_results = []
for response in res_whole_tf["response"].unique():
    res_data = res_whole_tf.query("response == @response")
    res_data = res_data[res_data["Adjusted P-value"] < 0.05]
    res_data=res_data.query("Gene_set == 'MSigDB_Hallmark_2020'")
    all_results.append((response, res_data))

# Extract TFs and their associated terms
for response_type, res_data in all_results:
    res_filtered = res_data[res_data["Adjusted P-value"] < 0.05]
    
    for _, row in res_filtered.iterrows():
        term = row['Term']
        genes = row['Genes'].split(';')
        
        # Filter for TFs only
        tf_genes = [gene for gene in genes if gene in tf_list]
        
        # Store TF classifications
        for tf in tf_genes:
            if tf not in tf_term_classification:
                tf_term_classification[tf] = set()
            tf_term_classification[tf].add(term)
        
        # Store terms and their TFs
        if term not in term_tf_dict:
            term_tf_dict[term] = set()
        term_tf_dict[term].update(tf_genes)

# Convert sets to lists for easier handling
for tf in tf_term_classification:
    tf_term_classification[tf] = list(tf_term_classification[tf])

for term in term_tf_dict:
    term_tf_dict[term] = list(term_tf_dict[term])

# Create a DataFrame for TF-Term associations
tf_term_data = []
for tf, terms in tf_term_classification.items():
    for term in terms:
        tf_term_data.append({
            'TF': tf,
            'Term': term,
            'Term_Category': term.split('_')[0] if '_' in term else term.split(' ')[0]
        })

df_tf_terms = pd.DataFrame(tf_term_data)

# Save the classification
df_tf_terms.to_csv(f"{save_path_metadata}/TF_term_classification_{gene_source}.csv", index=False)

# Print summary
print(f"=== TF CLASSIFICATION BY ENRICHED TERMS ===")
print(f"Total TFs classified: {len(tf_term_classification)}")
print(f"Total unique terms: {len(term_tf_dict)}")

# Show top terms by number of TFs
term_counts = df_tf_terms['Term'].value_counts()
print(f"\nTop 10 terms by TF count:")
print(term_counts.head(10))

# Create visualization of TF-term associations
# plt.figure(figsize=(15, 10))
# Create a binary matrix for heatmap
terms_to_plot = term_counts.head(15).index  # Top 15 terms
tfs_to_plot = df_tf_terms[df_tf_terms['Term'].isin(terms_to_plot)]['TF'].unique()

binary_matrix = pd.DataFrame(0, index=tfs_to_plot, columns=terms_to_plot)
for _, row in df_tf_terms.iterrows():
    if row['Term'] in terms_to_plot and row['TF'] in tfs_to_plot:
        binary_matrix.loc[row['TF'], row['Term']] = 1

sns.clustermap(binary_matrix, cmap='inferno', figsize=(15, 15), 
               row_cluster=True, col_cluster=False, 
               linewidths=0.5, square=False,cbar_pos=None)
plt.tight_layout()
plt.savefig(f"{save_path_figures}/TF_term_association_heatmap_{gene_source}.pdf", bbox_inches='tight', dpi=300)
plt.show()
# %%
gene_source = "HMGA2"  # Example gene source, can be changed as needed
# %%
# Extract links from gene source in both conditions
gene_source_links_good = links_good.filtered_links["Malignant"].query("source == @gene_source")
gene_source_links_poor = links_poor.filtered_links["Malignant"].query("source == @gene_source")

print(f"=== LINKS FOR {gene_source} ===")
print(f"Good response links: {len(gene_source_links_good)}")
print(f"Poor response links: {len(gene_source_links_poor)}")

# Get functional annotations for the target genes
if gene_source_links_good.shape[0] > 0 or gene_source_links_poor.shape[0] > 0:
    # Collect all unique target genes
    all_targets = set()
    if gene_source_links_good.shape[0] > 0:
        all_targets.update(gene_source_links_good['target'].unique())
    if gene_source_links_poor.shape[0] > 0:
        all_targets.update(gene_source_links_poor['target'].unique())
    
    all_targets = list(all_targets)
    print(f"Getting functional annotations for {len(all_targets)} target genes...")
    
    # Get gene annotations using the function defined earlier
    target_annotations = get_gene_annotations(all_targets)
    
    # Add functional annotations to the gene_source_links dataframes
    if gene_source_links_good.shape[0] > 0:
        gene_source_links_good['target_summary'] = gene_source_links_good['target'].map(
            lambda x: target_annotations.get(x, {}).get('summary', '')
        )
        gene_source_links_good['target_go_bp'] = gene_source_links_good['target'].map(
            lambda x: '; '.join(target_annotations.get(x, {}).get('go_bp', []))
        )
        gene_source_links_good['target_go_mf'] = gene_source_links_good['target'].map(
            lambda x: '; '.join(target_annotations.get(x, {}).get('go_mf', []))
        )
    
    if gene_source_links_poor.shape[0] > 0:
        gene_source_links_poor['target_summary'] = gene_source_links_poor['target'].map(
            lambda x: target_annotations.get(x, {}).get('summary', '')
        )
        gene_source_links_poor['target_go_bp'] = gene_source_links_poor['target'].map(
            lambda x: '; '.join(target_annotations.get(x, {}).get('go_bp', []))
        )
        gene_source_links_poor['target_go_mf'] = gene_source_links_poor['target'].map(
            lambda x: '; '.join(target_annotations.get(x, {}).get('go_mf', []))
        )

# Sort the dataframes by target for better organization
gene_source_links_good = gene_source_links_good.sort_values('target')
gene_source_links_poor = gene_source_links_poor.sort_values('target')
# Save the links
gene_source_links_good.to_csv(f"{save_path_metadata}/{gene_source}_links_good_response.csv", sep=",", index=False)
gene_source_links_poor.to_csv(f"{save_path_metadata}/{gene_source}_links_poor_response.csv", sep=",", index=False)
# Save Excel files alongside CSV files

# Save gene_source_links with annotations as Excel
with pd.ExcelWriter(f"{save_path_metadata}/{gene_source}_links_comparison.xlsx") as writer:
    gene_source_links_good.to_excel(writer, sheet_name='Good_Response_Links', index=False)
    gene_source_links_poor.to_excel(writer, sheet_name='Poor_Response_Links', index=False)
    
    # Create summary sheet
    summary_data = {
        'Metric': ['Good Response Targets', 'Poor Response Targets', 'Intersection Targets', 
                  'Good Unique Targets', 'Poor Unique Targets'],
        'Count': [len(good_targets), len(poor_targets), len(intersection_targets), 
                 len(good_unique_targets), len(poor_unique_targets)]
    }
    pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)
    
    # Save target lists
    pd.DataFrame({'Intersection_Targets': sorted(list(intersection_targets))}).to_excel(
        writer, sheet_name='Intersection_Targets', index=False)
    pd.DataFrame({'Good_Unique_Targets': sorted(list(good_unique_targets))}).to_excel(
        writer, sheet_name='Good_Unique_Targets', index=False)
    pd.DataFrame({'Poor_Unique_Targets': sorted(list(poor_unique_targets))}).to_excel(
        writer, sheet_name='Poor_Unique_Targets', index=False)

# Find target overlap and differences
good_targets = set(gene_source_links_good['target'].tolist())
poor_targets = set(gene_source_links_poor['target'].tolist())

intersection_targets = good_targets & poor_targets
good_unique_targets = good_targets - poor_targets
poor_unique_targets = poor_targets - good_targets

print(f"\nTarget comparison for {gene_source}:")
print(f"Good response targets: {len(good_targets)}")
print(f"Poor response targets: {len(poor_targets)}")
print(f"Intersection targets: {len(intersection_targets)}")
print(f"Good unique targets: {len(good_unique_targets)}")
print(f"Poor unique targets: {len(poor_unique_targets)}")

# Create comparison visualization
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Target overlap comparison
categories = ['Good Unique', 'Intersection', 'Poor Unique']
counts = [len(good_unique_targets), len(intersection_targets), len(poor_unique_targets)]
colors = ['#0096FF', '#9D4EDD', '#fb3310']

axes[0].bar(categories, counts, color=colors)
axes[0].set_title(f'{gene_source} Target Overlap')
axes[0].set_ylabel('Number of Targets')
for i, count in enumerate(counts):
    axes[0].text(i, count + 1, str(count), ha='center', va='bottom', fontweight='bold')

# Plot 2: Coefficient comparison for intersection targets
if intersection_targets:
    good_coefs = gene_source_links_good.set_index('target')['coef_abs']
    poor_coefs = gene_source_links_poor.set_index('target')['coef_abs']
    
    intersection_good_coefs = [good_coefs.get(target, 0) for target in intersection_targets]
    intersection_poor_coefs = [poor_coefs.get(target, 0) for target in intersection_targets]
    
    axes[1].scatter(intersection_good_coefs, intersection_poor_coefs, alpha=0.7, s=50)
    axes[1].plot([0, max(max(intersection_good_coefs), max(intersection_poor_coefs))], 
                [0, max(max(intersection_good_coefs), max(intersection_poor_coefs))], 'r--', alpha=0.5)
    axes[1].set_xlabel('Good Response Coefficient')
    axes[1].set_ylabel('Poor Response Coefficient')
    axes[1].set_title(f'{gene_source} Coefficient Comparison\n(Intersection Targets)')
    
    # Add text annotations for top targets
    top_targets = sorted(intersection_targets, 
                        key=lambda x: good_coefs.get(x, 0) + poor_coefs.get(x, 0), 
                        reverse=True)[:10]
    texts = []
    for target in top_targets:
        x = good_coefs.get(target, 0)
        y = poor_coefs.get(target, 0)
        texts.append(axes[1].text(x, y, target, fontsize=10))
    
    adjust_text.adjust_text(texts, arrowprops=dict(arrowstyle="->", color='r', lw=0.5), ax=axes[1])

sns.despine(ax=axes[0])
sns.despine(ax=axes[1])
plt.tight_layout()
plt.savefig(f"{save_path_figures}/{gene_source}_links_comparison.pdf", bbox_inches='tight', dpi=300)
plt.show()
# %%
gene_sets = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "MSigDB_Hallmark_2020",
    # "MSigDB_Oncogenic_Signatures",
    # "SynGO_2024",
    # "Tabula_Sapiens",
    # "CellMarker_2024",
    # "Rummagene_transcription_factors",
]
res_good_response=None
res_poor_response=None

plt.close("all")

try:
    # Run GSEA analysis
    res=gp.enrichr(gene_list=links_good.filtered_links["Malignant"].query("source == @gene_source")["target"].tolist(),
                gene_sets=gene_sets,
                # outdir=f"{save_path_metadata}/gsea_{filename}_treat",
                cutoff=0.05)
    res.results.to_csv(f"{save_path_metadata}/enrichment_{gene_source}_good_response.csv",sep=",",index=False)

    barplot(
        res.results,
        column="Adjusted P-value",
        group="Gene_set",
        size=10,
        top_term=10,
        title=f"{gene_source} - Enrichment Good Response Analysis",
        figsize=(10, 20),
        color={
            "KEGG_2021_Human": "#1f77b4",
            "GO_Biological_Process_2023": "#ff7f0e", 
            "GO_Molecular_Function_2023": "#2ca02c",
            "MSigDB_Hallmark_2020": "#d62728",
        },
        ofname=f"{save_path_figures}/enrichment_{gene_source}_good_response_barplot.pdf"
    )

    res_filtered = res.results[res.results["Adjusted P-value"] < 0.05]
    res_good_response = res_filtered.copy()
    res_filtered=res_filtered.query("Gene_set=='MSigDB_Hallmark_2020'")
    

    dotplot(res_filtered,
            ofname=f"{save_path_figures}/enrichment_{gene_source}_good_response_dotplot.pdf"
    )
except:
    print(f"Error in GSEA analysis for {gene_source} in treated links. Skipping...")
    pass

try:
    # Run GSEA analysis
    res=gp.enrichr(gene_list=links_poor.filtered_links["Malignant"].query("source == @gene_source")["target"].tolist(),
            gene_sets=gene_sets,
            # outdir=f"{save_path_metadata}/gsea_{filename}_treat",
            cutoff=0.05)
    res.results.to_csv(f"{save_path_metadata}/enrichment_{gene_source}_poor_response.csv",sep=",",index=False)

    barplot(
        res.results,
        column="Adjusted P-value",
        group="Gene_set",
        size=10,
        top_term=10,
        title=f"{gene_source} - Enrichment Poor Response Analysis",
        figsize=(10, 20),
        color={
            "KEGG_2021_Human": "#1f77b4",
            "GO_Biological_Process_2023": "#ff7f0e", 
            "GO_Molecular_Function_2023": "#2ca02c",
            "MSigDB_Hallmark_2020": "#d62728",
        },
        ofname=f"{save_path_figures}/enrichment_{gene_source}_poor_response_barplot.pdf"
    )

    res_filtered = res.results[res.results["Adjusted P-value"] < 0.05]
    res_poor_response = res_filtered.copy()
    res_filtered=res_filtered.query("Gene_set=='MSigDB_Hallmark_2020'")

    dotplot(res_filtered,
            ofname=f"{save_path_figures}/enrichment_{gene_source}_poor_response_dotplot.pdf"
    )
except:
    print(f"Error in GSEA analysis for {gene_source} in untreated links. Skipping...")
    pass
# %%
# Create network visualization with clustered enrichment annotations
plt.close('all')

# Extract target genes and enrichment terms from both good and bad response
good_targets = links_good.filtered_links["Malignant"].query("source == @gene_source")["target"].tolist()
bad_targets = links_poor.filtered_links["Malignant"].query("source == @gene_source")["target"].tolist()

# Get enrichment results for both conditions
good_enrichment_terms = {}
bad_enrichment_terms = {}

if res_good_response is not None:
    for _, row in res_good_response[res_good_response["Gene_set"] == "MSigDB_Hallmark_2020"].iterrows():
        term = row['Term'].replace("HALLMARK_", "").replace("_", " ")
        genes = row['Genes'].split(';')
        for gene in genes:
            if gene not in good_enrichment_terms:
                good_enrichment_terms[gene] = []
            good_enrichment_terms[gene].append(term)

if res_poor_response is not None:
    for _, row in res_poor_response[res_poor_response["Gene_set"] == "MSigDB_Hallmark_2020"].iterrows():
        term = row['Term'].replace("HALLMARK_", "").replace("_", " ")
        genes = row['Genes'].split(';')
        for gene in genes:
            if gene not in bad_enrichment_terms:
                bad_enrichment_terms[gene] = []
            bad_enrichment_terms[gene].append(term)

# Create network graphs
good_df = links_good.filtered_links["Malignant"].query("source == @gene_source")
bad_df = links_poor.filtered_links["Malignant"].query("source == @gene_source")

# Build NetworkX graphs
good_graph = nx.from_pandas_edgelist(good_df, source='source', target='target', 
                                   edge_attr=['coef_mean'], create_using=nx.DiGraph())
bad_graph = nx.from_pandas_edgelist(bad_df, source='source', target='target', 
                                  edge_attr=['coef_mean'], create_using=nx.DiGraph())

# Create unified annotation dictionary
all_enrichment_terms = {}
all_enrichment_terms.update(good_enrichment_terms)
for gene, terms in bad_enrichment_terms.items():
    if gene not in all_enrichment_terms:
        all_enrichment_terms[gene] = terms
    else:
        all_enrichment_terms[gene].extend(terms)

# Remove duplicates from combined annotations
for gene in all_enrichment_terms:
    all_enrichment_terms[gene] = list(set(all_enrichment_terms[gene]))

# Create pathway groups based on most frequent terms
pathway_groups = {}
for gene, terms in all_enrichment_terms.items():
    if terms:
        # Use the first (most significant) term as the primary pathway
        primary_pathway = terms[0]
        if primary_pathway not in pathway_groups:
            pathway_groups[primary_pathway] = []
        pathway_groups[primary_pathway].append(gene)

# Create color mapping for pathways
import matplotlib.cm as cm
pathway_colors = {}
color_palette = cm.get_cmap('tab20')
for i, pathway in enumerate(pathway_groups.keys()):
    pathway_colors[pathway] = color_palette(i / len(pathway_groups))

# Function to create clustered layout based on pathway annotations
def create_pathway_layout(graph, pathway_groups, gene_source):
    pos = {}
    
    # Position the source gene at the center
    pos[gene_source] = (0, 0)
    
    # Calculate positions for each pathway group
    num_groups = len(pathway_groups)
    angle_step = 2 * math.pi / num_groups if num_groups > 0 else 0
    
    for i, (pathway, genes) in enumerate(pathway_groups.items()):
        # Calculate group center position
        angle = i * angle_step
        group_center_x = 3 * math.cos(angle)
        group_center_y = 3 * math.sin(angle)
        
        # Position genes within the group using a circular layout
        genes_in_graph = [g for g in genes if g in graph.nodes()]
        if len(genes_in_graph) > 0:
            sub_angle_step = 2 * math.pi / len(genes_in_graph) if len(genes_in_graph) > 1 else 0
            radius = 0.8  # Radius for genes within each group
            
            for j, gene in enumerate(genes_in_graph):
                sub_angle = j * sub_angle_step
                pos[gene] = (
                    group_center_x + radius * math.cos(sub_angle),
                    group_center_y + radius * math.sin(sub_angle)
                )
    
    # Position remaining genes (without pathway annotation) in outer ring
    remaining_genes = [g for g in graph.nodes() if g not in pos]
    if remaining_genes:
        outer_radius = 5
        outer_angle_step = 2 * math.pi / len(remaining_genes)
        for i, gene in enumerate(remaining_genes):
            angle = i * outer_angle_step
            pos[gene] = (
                outer_radius * math.cos(angle),
                outer_radius * math.sin(angle)
            )
    
    return pos

# Create visualization for good response network
# fig, axes = plt.subplots(1, 2, figsize=(24, 12))
fig, ax = plt.subplots(1, 1, figsize=(12, 12))

# Good response network
pos_good = create_pathway_layout(good_graph, pathway_groups, gene_source)

# Create node colors based on pathway annotation
good_node_colors = []
good_node_labels = {}
for node in good_graph.nodes():
    if node == gene_source:
        good_node_colors.append('red')
        good_node_labels[node] = node
    else:
        # Find pathway for this node
        node_pathway = None
        for pathway, genes in pathway_groups.items():
            if node in genes:
                node_pathway = pathway
                break
        
        if node_pathway:
            good_node_colors.append(pathway_colors[node_pathway])
            good_node_labels[node] = node
        else:
            good_node_colors.append('lightgray')
            good_node_labels[node] = node

# Draw good response network
nx.draw_networkx_nodes(good_graph, pos_good, node_color=good_node_colors, 
                      node_size=300, alpha=0.8, ax=ax)
# Draw edges with color based on coef_mean: red if negative, green otherwise
edge_colors = []
for u, v, data in good_graph.edges(data=True):
    coef = data.get('coef_mean', 0)
    edge_colors.append('red' if coef < 0 else 'green')

nx.draw_networkx_edges(
    good_graph, pos_good, 
    edge_color=edge_colors, alpha=0.5, 
    arrows=True, arrowsize=10, ax=ax
)
nx.draw_networkx_labels(good_graph, pos_good, good_node_labels, font_size=8, 
                       font_weight='bold', ax=ax)

ax.set_title(f"{gene_source} Network - Good Response\n(Genes clustered by Hallmark pathways)", 
                  fontsize=14, fontweight='bold')
ax.axis('off')

# Create legend for pathways
legend_elements = []
legend_elements.append(Patch(facecolor='red', label=f'{gene_source} (Source TF)'))
for pathway, color in pathway_colors.items():
    if len(pathway) > 30:  # Truncate long pathway names
        pathway = pathway[:30] + "..."
    legend_elements.append(Patch(facecolor=color, label=pathway))
legend_elements.append(Patch(facecolor='lightgray', label='No annotation'))

plt.figlegend(handles=legend_elements, loc='center',
              bbox_to_anchor=(0.5, 0.01), 
              ncol=3, fontsize=10)

fig.savefig(f"{save_path_figures}/{gene_source}_clustered_networks_by_pathways_good.pdf", 
            bbox_inches='tight', dpi=300)
fig.savefig(f"{save_path_figures}/{gene_source}_clustered_networks_by_pathways_good.svg", 
            bbox_inches='tight', dpi=300)
fig.savefig(f"{save_path_figures}/{gene_source}_clustered_networks_by_pathways_good.tiff", 
            bbox_inches='tight', dpi=300)
plt.show()
# %% Bad Response Network
fig, ax = plt.subplots(1, 1, figsize=(12, 12))
# Bad response network
pos_bad = create_pathway_layout(bad_graph, pathway_groups, gene_source)

# Create node colors for bad response
bad_node_colors = []
bad_node_labels = {}
for node in bad_graph.nodes():
    if node == gene_source:
        bad_node_colors.append('red')
        bad_node_labels[node] = node
    else:
        # Find pathway for this node
        node_pathway = None
        for pathway, genes in pathway_groups.items():
            if node in genes:
                node_pathway = pathway
                break
        
        if node_pathway:
            bad_node_colors.append(pathway_colors[node_pathway])
            bad_node_labels[node] = node
        else:
            bad_node_colors.append('lightgray')
            bad_node_labels[node] = node

# Draw bad response network
nx.draw_networkx_nodes(bad_graph, pos_bad, node_color=bad_node_colors, 
                      node_size=300, alpha=0.8, ax=ax)
# Draw edges with color based on coef_mean: red if negative, green otherwise
edge_colors = []
for u, v, data in bad_graph.edges(data=True):
    coef = data.get('coef_mean', 0)
    edge_colors.append('red' if coef < 0 else 'green')

nx.draw_networkx_edges(
    bad_graph, pos_bad, 
    edge_color=edge_colors, alpha=0.5, 
    arrows=True, arrowsize=10, ax=ax
)
nx.draw_networkx_labels(bad_graph, pos_bad, bad_node_labels, font_size=8, 
                       font_weight='bold', ax=ax)

ax.set_title(f"{gene_source} Network - Poor Response\n(Genes clustered by Hallmark pathways)", 
              fontsize=14, fontweight='bold')
ax.axis('off')

# Create legend for pathways
legend_elements = []
legend_elements.append(Patch(facecolor='red', label=f'{gene_source} (Source TF)'))
for pathway, color in pathway_colors.items():
    if len(pathway) > 30:  # Truncate long pathway names
        pathway = pathway[:30] + "..."
    legend_elements.append(Patch(facecolor=color, label=pathway))
legend_elements.append(Patch(facecolor='lightgray', label='No annotation'))

plt.figlegend(handles=legend_elements, loc='center', 
              bbox_to_anchor=(0.5, 0.01), 
              ncol=3, fontsize=10)

fig.savefig(f"{save_path_figures}/{gene_source}_clustered_networks_by_pathways_bad.pdf", 
            bbox_inches='tight', dpi=300)
fig.savefig(f"{save_path_figures}/{gene_source}_clustered_networks_by_pathways_bad.svg", 
            bbox_inches='tight', dpi=300)
fig.savefig(f"{save_path_figures}/{gene_source}_clustered_networks_by_pathways_bad.tiff", 
            bbox_inches='tight', dpi=300)
plt.show()

# plt.tight_layout()
# plt.savefig(f"{save_path_figures}/{gene_source}_clustered_networks_by_pathways.pdf", 
#             bbox_inches='tight', dpi=300)
# plt.savefig(f"{save_path_figures}/{gene_source}_clustered_networks_by_pathways.svg", 
#             bbox_inches='tight', dpi=300)
# plt.show()

# Create summary statistics
print(f"=== NETWORK CLUSTERING SUMMARY FOR {gene_source} ===")
print(f"Total pathways identified: {len(pathway_groups)}")
print(f"Good response network: {good_graph.number_of_nodes()} nodes, {good_graph.number_of_edges()} edges")
print(f"Bad response network: {bad_graph.number_of_nodes()} nodes, {bad_graph.number_of_edges()} edges")

print(f"\n=== PATHWAY GROUPS ===")
for pathway, genes in pathway_groups.items():
    print(f"{pathway}: {len(genes)} genes")
    print(f"  Genes: {', '.join(genes[:10])}{'...' if len(genes) > 10 else ''}")

# Save pathway grouping results
pathway_df = []
for pathway, genes in pathway_groups.items():
    for gene in genes:
        good_response_target = gene in good_targets
        bad_response_target = gene in bad_targets
        pathway_df.append({
            'Gene': gene,
            'Pathway': pathway,
            'Good_Response_Target': good_response_target,
            'Bad_Response_Target': bad_response_target,
            'Both_Responses': good_response_target and bad_response_target
        })

df_pathway_results = pd.DataFrame(pathway_df)
df_pathway_results.to_csv(f"{save_path_metadata}/{gene_source}_pathway_clustering_results.csv", index=False)

print(f"\nResults saved to: {save_path_metadata}/{gene_source}_pathway_clustering_results.csv")
# %%
store_scores_good=list()
if res_good_response is None:
    print("Skipping score calculation and plotting due to missing enrichment results.")
else:
    rr_tmp=res_good_response.query("Gene_set == 'MSigDB_Hallmark_2020'")
    for term in rr_tmp['Term'].unique():
        rr_tmp_term = rr_tmp.query("Term == @term")
        genes = rr_tmp_term['Genes'].iloc[0].split(';')
        term = term.replace(" ", "_")
        store_scores_good.append(term)
        sc.tl.score_genes(adata_filtered, gene_list=genes, random_state=42, score_name=f"{gene_source}_{term}_good_response", use_raw=False)
# %%
order=['Untreated', 'Poor response', 'Minimal response', 'Moderate response']
# order_pairs = list(itertools.combinations(order, 2))
order_pairs=[('Untreated', 'Poor response'),
            ('Untreated', 'Minimal response'),
            ('Untreated', 'Moderate response'),
            ('Poor response', 'Moderate response'),
            ('Minimal response', 'Moderate response')]
# statistical test for annotation
ks2=StatTest(scipy.stats.ks_2samp, test_long_name="KS2sample",test_short_name="KS2")

# %%
# Calculate number of subplots needed
n_scores = len(store_scores_good)
ncols = 4  # Number of columns
nrows = math.ceil(n_scores / ncols)
palette_response = {
    "Untreated": "#55A868",
    "Poor response": "#8172B3",
    "Minimal response": "#4C72B0",
    "Moderate response": "#E32116",
}

# Create subplots
fig, axes = plt.subplots(nrows, ncols, figsize=(12, 5*nrows))

# Flatten axes array if more than one row
if nrows > 1:
    axes = axes.flatten()
elif nrows == 1:
    axes = [axes] if ncols == 1 else axes

# Plot each score in its own subplot
for idx, term in enumerate(sorted(store_scores_good)):
    ax = axes[idx] if n_scores > 1 else axes
    sns.boxplot(data=adata_filtered.obs, x='response', y=f"{gene_source}_{term}_good_response",
                order=order,
                ax=ax, palette=palette_response, showfliers=True, fill=False)
    annotator = Annotator(ax, order_pairs, data=adata_filtered.obs, x='response', y=f"{gene_source}_{term}_good_response", order=order)
    annotator.configure(test=ks2, text_format='star', loc='inside')
    annotator.apply_and_annotate()
    ax.grid(False)
    ax.hlines(0, *ax.get_xlim(), color='gray', linestyle='--')
    sns.despine(ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_ylabel("")
    ax.set_title(f"{term.replace('_', ' ')}")

# Hide empty subplots if any
if n_scores < len(axes):
    for idx in range(n_scores, len(axes)):
        axes[idx].set_visible(False)

fig.suptitle(f"{gene_source} - Enrichment Good Response Scores - Hallmark2020", fontsize=16, y=1.02)
fig.tight_layout()
fig.savefig(f"{save_path_figures}/enrichment_{gene_source}_good_response_scores.pdf", 
            bbox_inches='tight', dpi=300)
fig.savefig(f"{save_path_figures}/enrichment_{gene_source}_good_response_scores.svg", 
            bbox_inches='tight', dpi=300)
fig.show()
# %%
store_scores_poor=list()
if res_poor_response is None:
    print("Skipping score calculation and plotting due to missing enrichment results.")
else:
    rr_tmp=res_poor_response.query("Gene_set == 'MSigDB_Hallmark_2020'")
    for term in rr_tmp['Term'].unique():
        genes = rr_tmp[rr_tmp['Term'] == term]['Genes'].iloc[0].split(';')
        term = term.replace(" ", "_")
        store_scores_poor.append(term)
        sc.tl.score_genes(adata_filtered, gene_list=genes, random_state=42, score_name=f"{gene_source}_{term}_poor_response", use_raw=False)
# %%
# Calculate number of subplots needed
n_scores = len(store_scores_poor)
ncols = 4  # Number of columns
nrows = math.ceil(n_scores / ncols)

# Create subplots
fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 5*nrows))

# Flatten axes array if more than one row
if nrows > 1:
    axes = axes.flatten()
elif nrows == 1:
    axes = [axes] if ncols == 1 else axes

# Plot each score in its own subplot
for idx, term in enumerate(sorted(store_scores_poor)):
    ax = axes[idx] if n_scores > 1 else axes
    sns.boxplot(data=adata_filtered.obs, x='response', y=f"{gene_source}_{term}_poor_response",
                order=order,
                ax=ax, palette=palette_response, showfliers=True, fill=False)
    annotator = Annotator(ax, order_pairs, data=adata_filtered.obs, x='response', y=f"{gene_source}_{term}_poor_response", order=order)
    annotator.configure(test=ks2, text_format='star', loc='inside')
    annotator.apply_and_annotate()
    ax.grid(False)
    ax.hlines(0, *ax.get_xlim(), color='gray', linestyle='--')
    sns.despine(ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_ylabel("")
    ax.set_title(f"{term.replace('_', ' ')}")

# Hide empty subplots if any
if n_scores < len(axes):
    for idx in range(n_scores, len(axes)):
        axes[idx].set_visible(False)

fig.suptitle(f"{gene_source} - Enrichment Poor Response Scores - Hallmark2020", fontsize=16, y=1.02)
fig.tight_layout()
fig.savefig(f"{save_path_figures}/enrichment_{gene_source}_poor_response_scores.pdf", 
            bbox_inches='tight', dpi=300)
fig.savefig(f"{save_path_figures}/enrichment_{gene_source}_poor_response_scores.svg", 
            bbox_inches='tight', dpi=300)
fig.show()
# %%
# Combine good and bad response enrichment scores for paired comparison
combined_scores = sorted(list(set(store_scores_good) & set(store_scores_poor)))  # Get common enrichment terms

if combined_scores:
    # Calculate number of subplots needed
    n_scores = len(combined_scores)
    ncols = 4  # Number of columns
    nrows = math.ceil(n_scores / ncols)
    
    # Create subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 5*nrows))
    
    # Flatten axes array if more than one row
    if nrows > 1:
        axes = axes.flatten()
    elif nrows == 1:
        axes = [axes] if ncols == 1 else axes
    
    # Plot each paired score comparison
    for idx, term in enumerate(combined_scores):
        ax = axes[idx] if n_scores > 1 else axes
        
        # Prepare data for plotting - combine both good and bad response scores
        plot_data = []
        for response in order:
            # Get good response scores
            if f"{gene_source}_{term}_good_response" in adata_filtered.obs.columns:
                good_scores = adata_filtered.obs[adata_filtered.obs['response'] == response][f"{gene_source}_{term}_good_response"]
                for score in good_scores:
                    plot_data.append({
                        'response': response,
                        'score': score,
                        'source_type': 'Good Response Targets'
                    })
            
            # Get bad response scores  
            if f"{gene_source}_{term}_poor_response" in adata_filtered.obs.columns:
                bad_scores = adata_filtered.obs[adata_filtered.obs['response'] == response][f"{gene_source}_{term}_poor_response"]
                for score in bad_scores:
                    plot_data.append({
                        'response': response,
                        'score': score,
                        'source_type': 'Bad Response Targets'
                    })
        
        # Convert to DataFrame
        df_plot = pd.DataFrame(plot_data)

        # Alternative: Create pairs within each response group (Good vs Bad targets)
        response_pairs_within_group = []
        for response in order:
            response_pairs_within_group.append(((response, 'Good Response Targets'), 
                                               (response, 'Bad Response Targets')))

        # Create paired boxplot
        sns.boxplot(data=df_plot, x='response', y='score', hue='source_type',
                    order=order, ax=ax, 
                    palette={'Good Response Targets': '#0096FF', 'Bad Response Targets': '#fb3310'},
                    showfliers=True, fill=False,legend=False)
        annotator = Annotator(ax, response_pairs_within_group, data=df_plot,
                              x="response",y="score",hue="source_type", order=order)
        annotator.configure(test=ks2, text_format='star', loc='inside')
        annotator.apply_and_annotate()
        
        ax.grid(False)
        ax.hlines(0, *ax.get_xlim(), color='gray', linestyle='--')
        sns.despine(ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_ylabel("Gene Score")
        ax.set_xlabel("")
        ax.set_title(f"{term.replace('_', ' ')}")
    
    # Hide empty subplots if any
    if n_scores < len(axes):
        for idx in range(n_scores, len(axes)):
            axes[idx].set_visible(False)
    
    fig.suptitle(f"{gene_source} - Paired Enrichment Scores Comparison\nGood vs Bad Response Targets", 
                 fontsize=16, y=1.02)
    fig.tight_layout()
    fig.savefig(f"{save_path_figures}/enrichment_{gene_source}_paired_comparison_scores.pdf", 
                bbox_inches='tight', dpi=300)
    fig.savefig(f"{save_path_figures}/enrichment_{gene_source}_paired_comparison_scores.svg", 
                bbox_inches='tight', dpi=300)
    fig.show()
else:
    print("No common enrichment terms found between good and bad response groups.")
    
# %%
# Get unique scores from good and bad that are not in common
unique_good_scores = sorted(list(set(store_scores_good) - set(store_scores_poor)))
unique_bad_scores = sorted(list(set(store_scores_poor) - set(store_scores_good)))

print(f"Unique scores in good response: {len(unique_good_scores)}")
print(f"Unique scores in bad response: {len(unique_bad_scores)}")

# Plot unique good response scores
if unique_good_scores:
    # Calculate number of subplots needed for unique good scores
    n_scores = len(unique_good_scores)
    ncols = 4  # Number of columns
    nrows = math.ceil(n_scores / ncols)
    
    # Create subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 6*nrows))
    
    # Flatten axes array if more than one row
    if nrows > 1:
        axes = axes.flatten()
    elif nrows == 1:
        axes = [axes] if ncols == 1 else axes
    
    # Plot each score in its own subplot
    for idx, term in enumerate(unique_good_scores):
        ax = axes[idx] if n_scores > 1 else axes
        sns.boxplot(data=adata_filtered.obs, x='response', y=f"{gene_source}_{term}_good_response",
                    order=order,
                    ax=ax, palette=palette_response, showfliers=True, fill=False)
        annotator = Annotator(ax, order_pairs, data=adata_filtered.obs, x='response', 
                            y=f"{gene_source}_{term}_good_response", order=order)
        annotator.configure(test=ks2, text_format='star', loc='inside')
        annotator.apply_and_annotate()
        ax.grid(False)
        ax.hlines(0, *ax.get_xlim(), color='gray', linestyle='--')
        sns.despine(ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_ylabel("")
        ax.set_title(f"{term.replace('_', ' ')}")
    
    # Hide empty subplots if any
    if n_scores < len(axes):
        for idx in range(n_scores, len(axes)):
            axes[idx].set_visible(False)
    
    plt.suptitle(f"{gene_source} - Unique Good Response Scores - Hallmark2020", fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(f"{save_path_figures}/enrichment_{gene_source}_unique_good_response_scores.pdf", 
                bbox_inches='tight', dpi=300)
    plt.show()

# Plot unique bad response scores
if unique_bad_scores:
    # Calculate number of subplots needed for unique bad scores
    n_scores = len(unique_bad_scores)
    ncols = 4 if n_scores > 4 else n_scores  # Number of columns
    nrows = math.ceil(n_scores / ncols)
    
    # Create subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 6*nrows))
    
    # Flatten axes array if more than one row
    if nrows > 1:
        axes = axes.flatten()
    elif nrows == 1:
        axes = [axes] if ncols == 1 else axes
    
    # Plot each score in its own subplot
    for idx, term in enumerate(unique_bad_scores):
        ax = axes[idx] if n_scores > 1 else axes
        sns.boxplot(data=adata_filtered.obs, x='response', y=f"{gene_source}_{term}_poor_response",
                    order=order,
                    ax=ax, palette=palette_response, showfliers=True, fill=False)
        annotator = Annotator(ax, order_pairs, data=adata_filtered.obs, x='response', 
                            y=f"{gene_source}_{term}_poor_response", order=order)
        annotator.configure(test=ks2, text_format='star', loc='inside')
        annotator.apply_and_annotate()
        ax.grid(False)
        ax.hlines(0, *ax.get_xlim(), color='gray', linestyle='--')
        sns.despine(ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_ylabel("")
        ax.set_title(f"{term.replace('_', ' ')}")
    
    # Hide empty subplots if any
    if n_scores < len(axes):
        for idx in range(n_scores, len(axes)):
            axes[idx].set_visible(False)
    
    plt.suptitle(f"{gene_source} - Unique Bad Response Scores - Hallmark2020", fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(f"{save_path_figures}/enrichment_{gene_source}_unique_bad_response_scores.pdf", 
                bbox_inches='tight', dpi=300)
    plt.show()

# Print summary of what we found
print(f"\n=== UNIQUE SCORE SUMMARY ===")
print(f"Unique to good response: {unique_good_scores}")
print(f"Unique to bad response: {unique_bad_scores}")
# %%
# Get target genes for gene_source in good and bad response
good_targets = set(links_treat.filtered_links["Malignant"].query("source == @gene_source")["target"].tolist())
bad_targets = set(links_untreat.filtered_links["Malignant"].query("source == @gene_source")["target"].tolist())

# Find intersection and unique targets
intersection_targets = good_targets & bad_targets
good_unique_targets = good_targets - bad_targets
bad_unique_targets = bad_targets - good_targets

print(f"=== TARGET COMPARISON FOR {gene_source} ===")
print(f"Good response targets: {len(good_targets)}")
print(f"Bad response targets: {len(bad_targets)}")
print(f"Intersection targets: {len(intersection_targets)}")
print(f"Good unique targets: {len(good_unique_targets)}")
print(f"Bad unique targets: {len(bad_unique_targets)}")

print(f"\nIntersection targets: {sorted(list(intersection_targets))}")
print(f"Good unique targets: {sorted(list(good_unique_targets))}")
print(f"Bad unique targets: {sorted(list(bad_unique_targets))}")

# Save results
target_comparison_results = {
    'intersection_targets': sorted(list(intersection_targets)),
    'good_unique_targets': sorted(list(good_unique_targets)),
    'bad_unique_targets': sorted(list(bad_unique_targets))
}

df_target_results = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in target_comparison_results.items()]))
df_target_results.to_csv(f"{save_path_metadata}/{gene_source}_target_comparison_good_vs_bad.csv", index=False)

# Create visualization
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot 1: Bar plot of target counts
categories = ['Good Unique', 'Intersection', 'Bad Unique']
counts = [len(good_unique_targets), len(intersection_targets), len(bad_unique_targets)]
colors = ['#0096FF', '#9D4EDD', '#fb3310']

axes[0].bar(categories, counts, color=colors)
axes[0].set_title(f'{gene_source} Target Gene Overlap')
axes[0].set_ylabel('Number of Target Genes')
for i, count in enumerate(counts):
    axes[0].text(i, count + 1, str(count), ha='center', va='bottom', fontweight='bold')

# Plot 2: Venn diagram style representation
venn = venn2([good_targets, bad_targets], ('Good Response', 'Bad Response'), ax=axes[1])
axes[1].set_title(f'{gene_source} Target Overlap (Venn)')

# Plot 3: Top targets by coefficient for each category
if intersection_targets or good_unique_targets or bad_unique_targets:
    # Get coefficient information for visualization
    good_coef = links_treat.filtered_links["Malignant"].query("source == @gene_source").set_index('target')['coef_abs']
    bad_coef = links_untreat.filtered_links["Malignant"].query("source == @gene_source").set_index('target')['coef_abs']
    
    # Prepare data for plotting top targets
    plot_targets = []
    plot_coefs = []
    plot_categories = []
    
    # Add top intersection targets
    if intersection_targets:
        intersection_coefs = [(target, (good_coef.get(target, 0) + bad_coef.get(target, 0))/2) 
                             for target in intersection_targets]
        intersection_coefs.sort(key=lambda x: x[1], reverse=True)
        for target, coef in intersection_coefs[:5]:  # Top 5
            plot_targets.append(target)
            plot_coefs.append(coef)
            plot_categories.append('Intersection')
    
    # Add top good unique targets
    if good_unique_targets:
        good_unique_coefs = [(target, good_coef.get(target, 0)) for target in good_unique_targets]
        good_unique_coefs.sort(key=lambda x: x[1], reverse=True)
        for target, coef in good_unique_coefs[:5]:  # Top 5
            plot_targets.append(target)
            plot_coefs.append(coef)
            plot_categories.append('Good Unique')
    
    # Add top bad unique targets
    if bad_unique_targets:
        bad_unique_coefs = [(target, bad_coef.get(target, 0)) for target in bad_unique_targets]
        bad_unique_coefs.sort(key=lambda x: x[1], reverse=True)
        for target, coef in bad_unique_coefs[:5]:  # Top 5
            plot_targets.append(target)
            plot_coefs.append(coef)
            plot_categories.append('Bad Unique')
    
    # Create horizontal bar plot
    y_pos = range(len(plot_targets))
    bar_colors = ['#9D4EDD' if cat == 'Intersection' else '#0096FF' if cat == 'Good Unique' else '#fb3310' 
                 for cat in plot_categories]
    
    axes[2].barh(y_pos, plot_coefs, color=bar_colors, alpha=0.8)
    axes[2].set_yticks(y_pos)
    axes[2].set_yticklabels(plot_targets, fontsize=8)
    axes[2].set_xlabel('Coefficient (abs)')
    axes[2].set_title(f'Top {gene_source} Targets by Category')
    
    # Add legend
    legend_elements = [Patch(facecolor='#9D4EDD', label='Intersection'),
                      Patch(facecolor='#0096FF', label='Good Unique'),
                      Patch(facecolor='#fb3310', label='Bad Unique')]
    axes[2].legend(handles=legend_elements, loc='lower right')

plt.tight_layout()
plt.savefig(f"{save_path_figures}/{gene_source}_target_comparison_analysis.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{save_path_figures}/{gene_source}_target_comparison_analysis.svg", dpi=300, bbox_inches='tight')
plt.show()
# %%
df_exp=sc.get.obs_df(adata_filtered, keys=["response", "PPARG", "HMGA2"],use_raw=False)
df_exp
# %%
df_exp_melt = df_exp.melt(id_vars=["response"], var_name="TF", value_name="expression")
df_exp_melt
# %%
fig,ax=plt.subplots(figsize=(10,6))
sns.violinplot(data=df_exp_melt, x="TF", y="expression",palette=palette_response, hue="response",hue_order=order,ax=ax,fill=False)
sns.move_legend(ax, "upper right")
sns.despine(ax=ax)
ax.grid(False)
ax.set_title("PPARG and HMGA2 Expression by Response Type", fontsize=16)

fig.savefig(f"{save_path_figures}/PPARG-HMGA2_TF_expression_by_response.pdf", dpi=300, bbox_inches='tight')
fig.tight_layout()
fig.show()
# %%
# Extract top 1000 TFs for each cluster from treated and untreated merged_score
top_n = 50

# Get unique clusters that exist in both datasets
common_clusters = set(links_treat.merged_score['cluster'].unique()) & set(links_untreat.merged_score['cluster'].unique())
# Remove specific clusters from common_clusters
clusters_to_remove = {'CAF', 'myCAF', 'Ductal (atypical)',"Acinar"}
common_clusters = common_clusters - clusters_to_remove
# Initialize list to store data for heatmap
heatmap_data_list = []

for cluster in common_clusters:
    # Get top TFs for treated
    df_treat_cluster = links_treat.merged_score[
        (links_treat.merged_score['cluster'] == cluster) & 
        (links_treat.merged_score.index.isin(tf_list))
    ].query('eigenvector_centrality>0.2')
    
    # Get top TFs for untreated
    df_untreat_cluster = links_untreat.merged_score[
        (links_untreat.merged_score['cluster'] == cluster) & 
        (links_untreat.merged_score.index.isin(tf_list))
    ].query('eigenvector_centrality>0.2')

    # Get union of top TFs from both conditions
    all_top_tfs = set(df_treat_cluster.index) | set(df_untreat_cluster.index)
    
    for tf in all_top_tfs:
        # Get eigenvector centrality values (0 if TF not in top list)
        treat_centrality = df_treat_cluster.loc[tf, 'eigenvector_centrality'] if tf in df_treat_cluster.index else 0
        untreat_centrality = df_untreat_cluster.loc[tf, 'eigenvector_centrality'] if tf in df_untreat_cluster.index else 0
        
        # Calculate difference (treat + (-1 * untreat))
        centrality_difference = treat_centrality + (-1 * untreat_centrality)
        
        heatmap_data_list.append({
            'TF': tf,
            'Cluster': cluster,
            'Treated_Centrality': treat_centrality,
            'Untreated_Centrality': untreat_centrality,
            'Centrality_Difference': centrality_difference
        })

# Create DataFrame and pivot for heatmap
df_heatmap = pd.DataFrame(heatmap_data_list)  # Filter out rows with zero difference
df_heatmap.to_csv(f"{save_path_metadata}/regev_{filename}_top{top_n}_TFs_centrality_difference.csv", index=False)
heatmap_matrix = df_heatmap.pivot(index='TF', columns='Cluster', values='Centrality_Difference')
heatmap_matrix_larg=heatmap_matrix.nlargest(top_n, 'Malignant')  # Get top 20 TFs based on centrality difference
heatmap_matrix_small=heatmap_matrix.nsmallest(top_n, 'Malignant')  # Get top 20 TFs based on centrality difference

# Combine both DataFrames
heatmap_matrix = pd.concat([heatmap_matrix_larg, heatmap_matrix_small])
heatmap_matrix = heatmap_matrix[~heatmap_matrix.index.duplicated(keep='first')]
# Fill NaN values with 0
heatmap_matrix = heatmap_matrix.fillna(0)

# Create clustermap
# fig,ax = plt.subplots(figsize=(15, 20))
plt.close('all')
sns.set(font_scale=1.2)
sns.clustermap(heatmap_matrix, square=True, cmap='RdBu_r', annot=False,
               row_cluster=True, col_cluster=True, figsize=(20, 20),
               linewidths=1.5,metric='cosine',vmin=-1, vmax=1)
# plt.suptitle("Top TFs Eigenvector Centrality Difference\nGood Response - Bad Response", fontsize=20)
plt.tight_layout()
plt.savefig(f"{save_path_figures}/top{top_n}_TFs_centrality_difference_heatmap.pdf", bbox_inches='tight', dpi=300)
plt.savefig(f"{save_path_figures}/top{top_n}_TFs_centrality_difference_heatmap.svg", bbox_inches='tight', dpi=300)
plt.show()

# %%
# Extract top 1000 TFs for each cluster from treated and untreated merged_score
top_n = 30

# Get unique clusters that exist in both datasets
common_clusters = set(links_treat.merged_score['cluster'].unique()) & set(links_untreat.merged_score['cluster'].unique())
# Remove specific clusters from common_clusters
clusters_to_remove = {'CAF', 'myCAF', 'Ductal (atypical)',"Acinar"}
common_clusters = common_clusters - clusters_to_remove
# Initialize list to store data for heatmap
heatmap_data_list = []

for cluster in common_clusters:
    # Get top TFs for treated
    df_treat_cluster = links_treat.merged_score[
        (links_treat.merged_score['cluster'] == cluster) & 
        (links_treat.merged_score.index.isin(tf_list))
    ].query('eigenvector_centrality>0.2')
    
    # Get top TFs for untreated
    df_untreat_cluster = links_untreat.merged_score[
        (links_untreat.merged_score['cluster'] == cluster) & 
        (links_untreat.merged_score.index.isin(tf_list))
    ].query('eigenvector_centrality>0.2')

    # Get intersection of top TFs from both conditions
    all_top_tfs = set(df_treat_cluster.index).intersection(set(df_untreat_cluster.index))

    for tf in all_top_tfs:
        # Get eigenvector centrality values (0 if TF not in top list)
        treat_centrality = df_treat_cluster.loc[tf, 'eigenvector_centrality'] if tf in df_treat_cluster.index else 0
        untreat_centrality = df_untreat_cluster.loc[tf, 'eigenvector_centrality'] if tf in df_untreat_cluster.index else 0
        
        # Calculate difference (treat + (-1 * untreat))
        centrality_difference = treat_centrality + (-1 * untreat_centrality)
        
        heatmap_data_list.append({
            'TF': tf,
            'Cluster': cluster,
            'Treated_Centrality': treat_centrality,
            'Untreated_Centrality': untreat_centrality,
            'Centrality_Difference': centrality_difference
        })

# Create DataFrame and pivot for heatmap
df_heatmap = pd.DataFrame(heatmap_data_list)  # Filter out rows with zero difference
heatmap_matrix = df_heatmap.pivot(index='TF', columns='Cluster', values='Centrality_Difference')
# Fill NaN values with 0
heatmap_matrix = heatmap_matrix.fillna(0)

# Create clustermap
# fig,ax = plt.subplots(figsize=(15, 20))
plt.close('all')
sns.set(font_scale=2)
sns.clustermap(heatmap_matrix, square=True, cmap='RdBu_r', annot=False,
               row_cluster=True, col_cluster=True, figsize=(20, 20),
               linewidths=2,metric='cosine',vmin=-1,vmax=1)
# plt.suptitle("Top TFs Eigenvector Centrality Difference\nGood Response - Bad Response", fontsize=20)
plt.tight_layout()
plt.savefig(f"{save_path_figures}/top{top_n}_TFs_centrality_difference_intersection_heatmap.pdf", bbox_inches='tight', dpi=300)
plt.savefig(f"{save_path_figures}/top{top_n}_TFs_centrality_difference_intersection_heatmap.svg", bbox_inches='tight', dpi=300)
plt.show()

# %%
filename="regev_malignant_response"
oracle_malignant = co.load_hdf5(f"{save_path_networks}/{filename}.celloracle.oracle")
# oracle_malignant
links_malignant = co.load_hdf5(f"{save_path_networks}/{filename}.celloracle.links")
# %%
adata_malignant = oracle_malignant.adata.copy()
adata_malignant.obs['response'] = pd.Categorical(adata_malignant.obs['response'], categories=["Untreated", "Poor response",  "Minimal response", "Moderate response"],ordered=True)
# %%
# top_n_per_cluster = 100

# Filter for TFs only and get top TFs for each cluster
cluster_tf_dict = {}
for cluster in links_malignant.merged_score['cluster'].unique():
    cluster_tfs = links_malignant.merged_score[
        (links_malignant.merged_score['cluster'] == cluster) & 
        (links_malignant.merged_score.index.isin(tf_list))
    ]
    cluster_tf_dict[cluster] = set(cluster_tfs.index)

# Find common TFs across all clusters
all_clusters = list(cluster_tf_dict.keys())
common_tfs_all = set.intersection(*cluster_tf_dict.values()) if cluster_tf_dict else set()

print(f"=== COMMON TFs ACROSS ALL CLUSTERS ===")
print(f"Number of common TFs: {len(common_tfs_all)}")
print(f"Common TFs: {sorted(list(common_tfs_all))}")

# Group clusters by response type
moderate_minimal_clusters = [cluster for cluster in all_clusters if 'Moderate response' in cluster or 'Minimal response' in cluster]
poor_untreated_clusters = [cluster for cluster in all_clusters if 'Poor response' in cluster or 'Untreated' in cluster]

print(f"\nModerate/Minimal clusters: {moderate_minimal_clusters}")
print(f"Poor/Untreated clusters: {poor_untreated_clusters}")

# Find TFs specific to moderate+minimal response clusters
if moderate_minimal_clusters:
    moderate_minimal_tfs = set.intersection(*[cluster_tf_dict[cluster] for cluster in moderate_minimal_clusters])
else:
    moderate_minimal_tfs = set()

# Find TFs specific to poor+untreated response clusters  
if poor_untreated_clusters:
    poor_untreated_tfs = set.intersection(*[cluster_tf_dict[cluster] for cluster in poor_untreated_clusters])
else:
    poor_untreated_tfs = set()

# Find TFs that appear only in moderate+minimal clusters (not in poor+untreated)
moderate_minimal_only_tfs = moderate_minimal_tfs - poor_untreated_tfs

# Find TFs that appear only in poor+untreated clusters (not in moderate+minimal)
poor_untreated_only_tfs = poor_untreated_tfs - moderate_minimal_tfs

print(f"\n=== TFs SPECIFIC TO RESPONSE GROUPS ===")
print(f"TFs common in Moderate/Minimal clusters: {len(moderate_minimal_tfs)}")
print(f"Moderate/Minimal TFs: {sorted(list(moderate_minimal_tfs))}")

print(f"\nTFs common in Poor/Untreated clusters: {len(poor_untreated_tfs)}")
print(f"Poor/Untreated TFs: {sorted(list(poor_untreated_tfs))}")

print(f"\nTFs ONLY in Moderate/Minimal (not in Poor/Untreated): {len(moderate_minimal_only_tfs)}")
print(f"Moderate/Minimal ONLY TFs: {sorted(list(moderate_minimal_only_tfs))}")

print(f"\nTFs ONLY in Poor/Untreated (not in Moderate/Minimal): {len(poor_untreated_only_tfs)}")
print(f"Poor/Untreated ONLY TFs: {sorted(list(poor_untreated_only_tfs))}")

# Save results
results_dict = {
    'common_all_clusters': sorted(list(common_tfs_all)),
    'moderate_minimal_common': sorted(list(moderate_minimal_tfs)),
    'poor_untreated_common': sorted(list(poor_untreated_tfs)),
    'moderate_minimal_only': sorted(list(moderate_minimal_only_tfs)),
    'poor_untreated_only': sorted(list(poor_untreated_only_tfs))
}

df_results = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in results_dict.items()]))
df_results.to_csv(f"{save_path_metadata}/malignant_response_TF_analysis.csv", index=False)

# Create visualization for top 10 TFs by eigenvector centrality in each response group
plt.close('all')
fig, axes = plt.subplots(2, 2, figsize=(20, 15))

# Get top 10 TFs for each response group by eigenvector centrality
response_groups = ['Moderate response', 'Minimal response', 'Poor response', 'Untreated']
colors = ['#E32116', '#4C72B0', '#8172B3', '#55A868']

for idx, (response_group, color) in enumerate(zip(response_groups, colors)):
    row = idx // 2
    col = idx % 2
    ax = axes[row, col]
    
    # Get TFs for this response group
    response_tfs = links_malignant.merged_score[
        (links_malignant.merged_score['cluster'] == response_group) & 
        (links_malignant.merged_score.index.isin(tf_list))
    ].sort_values('eigenvector_centrality', ascending=False)
    
    # Get top 10 TFs
    top_10_tfs = response_tfs.head(10)
    
    # Create horizontal bar plot
    y_positions = range(len(top_10_tfs))
    ax.barh(y_positions, top_10_tfs['eigenvector_centrality'], 
            color=color, alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # Customize the plot
    ax.set_yticks(y_positions)
    ax.set_yticklabels(top_10_tfs.index, fontsize=10)
    ax.set_xlabel('Eigenvector Centrality', fontsize=12)
    ax.set_title(f'Top 10 TFs - {response_group}', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    ax.set_xlim(0, 1)

    sns.despine(ax=ax)
    ax.grid(False)

    # Add value labels on bars
    for i, (tf, centrality) in enumerate(zip(top_10_tfs.index, top_10_tfs['eigenvector_centrality'])):
        ax.text(centrality + 0.01, i, f'{centrality:.3f}', 
                va='center', ha='left', fontsize=9)

plt.tight_layout()
plt.savefig(f"{save_path_figures}/top10_TFs_eigenvector_centrality_by_response.pdf", 
            bbox_inches='tight', dpi=300)
plt.savefig(f"{save_path_figures}/top10_TFs_eigenvector_centrality_by_response.svg", 
            bbox_inches='tight', dpi=300)
plt.show()

# Create a combined comparison plot
fig, ax = plt.subplots(figsize=(15, 12))

# Collect all top 10 TFs and their data
all_top_tfs_data = []
for response_group, color in zip(response_groups, colors):
    response_tfs = links_malignant.merged_score[
        (links_malignant.merged_score['cluster'] == response_group) & 
        (links_malignant.merged_score.index.isin(tf_list))
    ].sort_values('eigenvector_centrality', ascending=False)
    
    top_10_tfs = response_tfs.head(10)
    for tf, centrality in zip(top_10_tfs.index, top_10_tfs['eigenvector_centrality']):
        all_top_tfs_data.append({
            'TF': tf,
            'Response': response_group,
            'Eigenvector_Centrality': centrality,
            'Color': color
        })

# Create DataFrame and sort by centrality
df_all_top = pd.DataFrame(all_top_tfs_data)
df_all_top = df_all_top.sort_values(by=["TF",'Eigenvector_Centrality'], ascending=[False, True])

# Create the plot
y_positions = range(len(df_all_top))
bars = ax.barh(y_positions, df_all_top['Eigenvector_Centrality'], 
               color=df_all_top['Color'], alpha=0.8, edgecolor='black', linewidth=0.5)

# Count occurrences of each TF and create horizontal lines to separate groups
tf_counts = df_all_top['TF'].value_counts()
current_position = 0
previous_tf = None

for i, tf in enumerate(df_all_top['TF']):
    # Add horizontal line when TF changes
    if previous_tf is not None and tf != previous_tf:
        ax.axhline(y=i-0.5, color='black', linewidth=0.8, alpha=0.6)
    
    # Add text annotation with count for the first occurrence of each TF
    if previous_tf != tf:
        count = tf_counts[tf]
        ax.text(max(df_all_top['Eigenvector_Centrality']) * 1.02, i, 
                f'n={count}', va='center', ha='left', fontsize=8, 
                bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.7))
    
    previous_tf = tf

# Customize the plot
ax.set_yticks(y_positions)
# Create labels that show TF name only once per group
tf_labels = []
tf_groups = {}
current_pos = 0

# First pass: identify positions for each TF
for i, tf in enumerate(df_all_top['TF']):
    if tf not in tf_groups:
        tf_groups[tf] = []
    tf_groups[tf].append(i)

# Second pass: create labels with TF name at center position
for i, tf in enumerate(df_all_top['TF']):
    positions = tf_groups[tf]
    if len(positions) % 2 == 1:  # odd number - use middle position
        center_pos = positions[len(positions) // 2]
    else:  # even number - use center between two middle positions
        center_pos = (positions[len(positions) // 2 - 1] + positions[len(positions) // 2]) / 2
    
    # Only add label if current position is closest to center
    if i == min(positions, key=lambda x: abs(x - center_pos)):
        tf_labels.append(tf)
    else:
        tf_labels.append('')

ax.set_yticklabels(tf_labels, fontsize=10)
# ax.set_yticklabels([f"{tf} ({resp})" for tf, resp in 
#                     zip(df_all_top['TF'], df_all_top['Response'])], fontsize=9)
ax.set_xlabel('Eigenvector Centrality', fontsize=12)
ax.set_title('Top 10 TFs by Eigenvector Centrality - All Response Groups', 
             fontsize=14, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

# Add legend
import matplotlib.patches as mpatches

legend_elements = [mpatches.Patch(facecolor=color, label=response) 
                   for response, color in zip(response_groups, colors)]
ax.legend(handles=legend_elements, loc='lower right')

plt.grid(False)
sns.despine()
plt.tight_layout()
plt.savefig(f"{save_path_figures}/all_top10_TFs_eigenvector_centrality_combined.pdf", 
            bbox_inches='tight', dpi=300)
plt.savefig(f"{save_path_figures}/all_top10_TFs_eigenvector_centrality_combined.svg", 
            bbox_inches='tight', dpi=300)
plt.show()
# %%
# Extract expression data for the unique TFs from df_all_top
unique_tfs = df_all_top["TF"].unique().tolist()

sc.pl.matrixplot(adata_malignant_original, var_names=sorted(unique_tfs),standard_scale='var',
                 groupby='response', use_raw=False, cmap='RdBu_r', dendrogram=False, 
                 swap_axes=True, save="all_top10_TFs_expression_matrix.pdf")
# %%
# Perform enrichment analysis for different TF groups
print("=== PERFORMING ENRICHMENT ANALYSIS FOR TF GROUPS ===\n")

gene_sets = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "MSigDB_Hallmark_2020",
    "Reactome_2022"
]

# Define colors for each gene set
gene_set_colors = {
    "KEGG_2021_Human": "#1f77b4",
    "GO_Biological_Process_2023": "#ff7f0e",
    "GO_Molecular_Function_2023": "#2ca02c",
    "MSigDB_Hallmark_2020": "#d62728",
    "Reactome_2022": "#8c564b"
}

# Dictionary to store all enrichment results
enrichment_results = {}

# 1. Enrichment for common TFs across all clusters
if common_tfs_all:
    print(f"1. Enriching common TFs across all clusters ({len(common_tfs_all)} TFs)...")
    try:
        res_common = gp.enrichr(
            gene_list=sorted(list(common_tfs_all)),
            gene_sets=gene_sets,
            cutoff=0.05
        )
        
        res_common.results.to_csv(f"{save_path_metadata}/enrichment_common_tfs_all_clusters.csv", sep=",", index=False)
        enrichment_results['common_all'] = res_common.results
        
        # Create barplot for common TFs
        res_filtered = res_common.results[res_common.results["Adjusted P-value"] < 0.05]
        if not res_filtered.empty:
            barplot(
                res_filtered,
                column="Adjusted P-value",
                group="Gene_set",
                size=10,
                top_term=15,
                title="Enrichment Analysis - Common TFs Across All Clusters",
                figsize=(12, 15),
                color=gene_set_colors,
                ofname=f"{save_path_figures}/enrichment_common_tfs_all_clusters.pdf"
            )
        print(f"   Found {len(res_filtered)} significant terms")
    except Exception as e:
        print(f"   Error in common TFs enrichment: {e}")

# 2. Enrichment for moderate/minimal response TFs
if moderate_minimal_tfs:
    print(f"\n2. Enriching moderate/minimal response TFs ({len(moderate_minimal_tfs)} TFs)...")
    try:
        res_moderate_minimal = gp.enrichr(
            gene_list=sorted(list(moderate_minimal_tfs)),
            gene_sets=gene_sets,
            cutoff=0.05
        )
        
        res_moderate_minimal.results.to_csv(f"{save_path_metadata}/enrichment_moderate_minimal_tfs.csv", sep=",", index=False)
        enrichment_results['moderate_minimal'] = res_moderate_minimal.results
        
        # Create barplot for moderate/minimal TFs
        res_filtered = res_moderate_minimal.results[res_moderate_minimal.results["Adjusted P-value"] < 0.05]
        if not res_filtered.empty:
            barplot(
                res_filtered,
                column="Adjusted P-value",
                group="Gene_set",
                size=10,
                top_term=15,
                title="Enrichment Analysis - Moderate/Minimal Response TFs",
                figsize=(12, 15),
                color=gene_set_colors,
                ofname=f"{save_path_figures}/enrichment_moderate_minimal_tfs.pdf"
            )
        print(f"   Found {len(res_filtered)} significant terms")
    except Exception as e:
        print(f"   Error in moderate/minimal TFs enrichment: {e}")

# 3. Enrichment for poor/untreated response TFs  
if poor_untreated_tfs:
    print(f"\n3. Enriching poor/untreated response TFs ({len(poor_untreated_tfs)} TFs)...")
    try:
        res_poor_untreated = gp.enrichr(
            gene_list=sorted(list(poor_untreated_tfs)),
            gene_sets=gene_sets,
            cutoff=0.05
        )
        
        res_poor_untreated.results.to_csv(f"{save_path_metadata}/enrichment_poor_untreated_tfs.csv", sep=",", index=False)
        enrichment_results['poor_untreated'] = res_poor_untreated.results
        
        # Create barplot for poor/untreated TFs
        res_filtered = res_poor_untreated.results[res_poor_untreated.results["Adjusted P-value"] < 0.05]
        if not res_filtered.empty:
            barplot(
                res_filtered,
                column="Adjusted P-value",
                group="Gene_set",
                size=10,
                top_term=15,
                title="Enrichment Analysis - Poor/Untreated Response TFs",
                figsize=(12, 15),
                color=gene_set_colors,
                ofname=f"{save_path_figures}/enrichment_poor_untreated_tfs.pdf"
            )
        print(f"   Found {len(res_filtered)} significant terms")
    except Exception as e:
        print(f"   Error in poor/untreated TFs enrichment: {e}")

# 4. Enrichment for moderate/minimal ONLY TFs
if moderate_minimal_only_tfs:
    print(f"\n4. Enriching moderate/minimal ONLY TFs ({len(moderate_minimal_only_tfs)} TFs)...")
    try:
        res_moderate_minimal_only = gp.enrichr(
            gene_list=sorted(list(moderate_minimal_only_tfs)),
            gene_sets=gene_sets,
            cutoff=0.05
        )
        
        res_moderate_minimal_only.results.to_csv(f"{save_path_metadata}/enrichment_moderate_minimal_only_tfs.csv", sep=",", index=False)
        enrichment_results['moderate_minimal_only'] = res_moderate_minimal_only.results
        
        # Create barplot for moderate/minimal ONLY TFs
        res_filtered = res_moderate_minimal_only.results[res_moderate_minimal_only.results["Adjusted P-value"] < 0.05]
        if not res_filtered.empty:
            barplot(
                res_filtered,
                column="Adjusted P-value",
                group="Gene_set",
                size=10,
                top_term=15,
                title="Enrichment Analysis - Moderate/Minimal ONLY TFs",
                figsize=(12, 15),
                color=gene_set_colors,
                ofname=f"{save_path_figures}/enrichment_moderate_minimal_only_tfs.pdf"
            )
        print(f"   Found {len(res_filtered)} significant terms")
    except Exception as e:
        print(f"   Error in moderate/minimal ONLY TFs enrichment: {e}")

# 5. Enrichment for poor/untreated ONLY TFs
if poor_untreated_only_tfs:
    print(f"\n5. Enriching poor/untreated ONLY TFs ({len(poor_untreated_only_tfs)} TFs)...")
    try:
        res_poor_untreated_only = gp.enrichr(
            gene_list=sorted(list(poor_untreated_only_tfs)),
            gene_sets=gene_sets,
            cutoff=0.05
        )
        
        res_poor_untreated_only.results.to_csv(f"{save_path_metadata}/enrichment_poor_untreated_only_tfs.csv", sep=",", index=False)
        enrichment_results['poor_untreated_only'] = res_poor_untreated_only.results
        
        # Create barplot for poor/untreated ONLY TFs
        res_filtered = res_poor_untreated_only.results[res_poor_untreated_only.results["Adjusted P-value"] < 0.05]
        if not res_filtered.empty:
            barplot(
                res_filtered,
                column="Adjusted P-value",
                group="Gene_set",
                size=10,
                top_term=15,
                title="Enrichment Analysis - Poor/Untreated ONLY TFs",
                figsize=(12, 15),
                color=gene_set_colors,
                ofname=f"{save_path_figures}/enrichment_poor_untreated_only_tfs.pdf"
            )
        print(f"   Found {len(res_filtered)} significant terms")
    except Exception as e:
        print(f"   Error in poor/untreated ONLY TFs enrichment: {e}")

# 6. Enrichment for good vs bad paired comparison (from earlier analysis)
if 'good_unique_tfs' in locals() and 'bad_unique_tfs' in locals():
    print(f"\n6. Enriching good response unique TFs ({len(good_unique_tfs)} TFs)...")
    try:
        res_good_unique = gp.enrichr(
            gene_list=sorted(list(good_unique_tfs)),
            gene_sets=gene_sets,
            cutoff=0.05
        )
        
        res_good_unique.results.to_csv(f"{save_path_metadata}/enrichment_good_unique_tfs.csv", sep=",", index=False)
        enrichment_results['good_unique'] = res_good_unique.results
        
        # Create barplot for good unique TFs
        res_filtered = res_good_unique.results[res_good_unique.results["Adjusted P-value"] < 0.05]
        if not res_filtered.empty:
            barplot(
                res_filtered,
                column="Adjusted P-value",
                group="Gene_set",
                size=10,
                top_term=15,
                title="Enrichment Analysis - Good Response Unique TFs",
                figsize=(12, 15),
                color=gene_set_colors,
                ofname=f"{save_path_figures}/enrichment_good_unique_tfs.pdf"
            )
        print(f"   Found {len(res_filtered)} significant terms")
    except Exception as e:
        print(f"   Error in good unique TFs enrichment: {e}")

    print(f"\n7. Enriching bad response unique TFs ({len(bad_unique_tfs)} TFs)...")
    try:
        res_bad_unique = gp.enrichr(
            gene_list=sorted(list(bad_unique_tfs)),
            gene_sets=gene_sets,
            cutoff=0.05
        )
        
        res_bad_unique.results.to_csv(f"{save_path_metadata}/enrichment_bad_unique_tfs.csv", sep=",", index=False)
        enrichment_results['bad_unique'] = res_bad_unique.results
        
        # Create barplot for bad unique TFs
        res_filtered = res_bad_unique.results[res_bad_unique.results["Adjusted P-value"] < 0.05]
        if not res_filtered.empty:
            barplot(
                res_filtered,
                column="Adjusted P-value",
                group="Gene_set",
                size=10,
                top_term=15,
                title="Enrichment Analysis - Bad Response Unique TFs",
                figsize=(12, 15),
                color=gene_set_colors,
                ofname=f"{save_path_figures}/enrichment_bad_unique_tfs.pdf"
            )
        print(f"   Found {len(res_filtered)} significant terms")
    except Exception as e:
        print(f"   Error in bad unique TFs enrichment: {e}")

    # Enrichment for intersection TFs (good and bad common)
    if 'intersection_tfs' in locals() and intersection_tfs:
        print(f"\n8. Enriching intersection TFs ({len(intersection_tfs)} TFs)...")
        try:
            res_intersection = gp.enrichr(
                gene_list=sorted(list(intersection_tfs)),
                gene_sets=gene_sets,
                cutoff=0.05
            )
            
            res_intersection.results.to_csv(f"{save_path_metadata}/enrichment_intersection_tfs.csv", sep=",", index=False)
            enrichment_results['intersection'] = res_intersection.results
            
            # Create barplot for intersection TFs
            res_filtered = res_intersection.results[res_intersection.results["Adjusted P-value"] < 0.05]
            if not res_filtered.empty:
                barplot(
                    res_filtered,
                    column="Adjusted P-value",
                    group="Gene_set",
                    size=10,
                    top_term=15,
                    title="Enrichment Analysis - Intersection TFs (Good & Bad Common)",
                    figsize=(12, 15),
                    color=gene_set_colors,
                    ofname=f"{save_path_figures}/enrichment_intersection_tfs.pdf"
                )
            print(f"   Found {len(res_filtered)} significant terms")
        except Exception as e:
            print(f"   Error in intersection TFs enrichment: {e}")

# Create a comprehensive comparison plot of all enrichment results
print("\n=== CREATING COMPREHENSIVE ENRICHMENT COMPARISON ===")

# Combine all results with group labels
all_enrichment_data = []
group_labels = {
    'common_all': 'Common All Clusters',
    'moderate_minimal': 'Moderate/Minimal',
    'poor_untreated': 'Poor/Untreated',
    'moderate_minimal_only': 'Moderate/Minimal Only',
    'poor_untreated_only': 'Poor/Untreated Only',
    'good_unique': 'Good Unique',
    'bad_unique': 'Bad Unique',
    'intersection': 'Good/Bad Intersection'
}

for group_key, group_label in group_labels.items():
    if group_key in enrichment_results:
        df_group = enrichment_results[group_key].copy()
        df_group = df_group[df_group["Adjusted P-value"] < 0.05]
        df_group = df_group[df_group["Gene_set"] == "MSigDB_Hallmark_2020"]  # Focus on Hallmark pathways
        df_group['TF_Group'] = group_label
        all_enrichment_data.append(df_group)

if all_enrichment_data:
    # Combine all enrichment data
    df_all_enrichment = pd.concat(all_enrichment_data, ignore_index=True)
    
    # Save combined results
    df_all_enrichment.to_csv(f"{save_path_metadata}/all_tf_groups_enrichment_comparison.csv", sep=",", index=False)
    
    # Create a comparison dotplot
    if not df_all_enrichment.empty:
        # Get top 20 most significant terms across all groups
        top_terms = df_all_enrichment.nsmallest(20, 'Adjusted P-value')['Term'].unique()
        df_plot = df_all_enrichment[df_all_enrichment['Term'].isin(top_terms)]
        
        # Create pivot table for heatmap
        heatmap_data = df_plot.pivot_table(
            index='Term', 
            columns='TF_Group', 
            values='Adjusted P-value',
            fill_value=1.0
        )
        
        # Transform p-values for better visualization (-log10)
        heatmap_data_log = -np.log10(heatmap_data)
        
        # Create clustermap
        plt.figure(figsize=(15, 12))
        sns.clustermap(
            heatmap_data_log,
            cmap='Reds',
            annot=False,
            fmt='.2f',
            figsize=(20, 15),
            row_cluster=True,
            col_cluster=True,
            linewidths=0.5,
            cbar_kws={'label': '-log10(Adjusted P-value)'}
        )
        plt.tight_layout()
        plt.savefig(f"{save_path_figures}/tf_groups_enrichment_comparison_heatmap.pdf", 
                    bbox_inches='tight', dpi=300)
        plt.show()
        
        print(f"Created comprehensive enrichment comparison with {len(df_all_enrichment)} significant terms")
    else:
        print("No significant enrichment terms found for comparison")
else:
    print("No enrichment results available for comparison")

print("\n=== ENRICHMENT ANALYSIS COMPLETED ===")
print("Results saved in:")
for group_key, group_label in group_labels.items():
    if group_key in enrichment_results:
        print(f"  - {group_label}: {save_path_metadata}/enrichment_{group_key}_tfs.csv")
# %%
# Extract targets for each TF group and perform enrichment analysis
print("=== EXTRACTING TARGETS FOR TF GROUPS AND PERFORMING ENRICHMENT ===\n")

# Define the TF groups we want to analyze
tf_groups_dict = {}
if 'common_tfs_all' in locals() and common_tfs_all:
    tf_groups_dict['common_all'] = common_tfs_all
if 'moderate_minimal_tfs' in locals() and moderate_minimal_tfs:
    tf_groups_dict['moderate_minimal'] = moderate_minimal_tfs
if 'poor_untreated_tfs' in locals() and poor_untreated_tfs:
    tf_groups_dict['poor_untreated'] = poor_untreated_tfs
if 'good_unique_tfs' in locals() and good_unique_tfs:
    tf_groups_dict['good_unique'] = good_unique_tfs
if 'bad_unique_tfs' in locals() and bad_unique_tfs:
    tf_groups_dict['bad_unique'] = bad_unique_tfs
if 'intersection_tfs' in locals() and intersection_tfs:
    tf_groups_dict['intersection'] = intersection_tfs

# Gene sets for enrichment analysis
target_gene_sets = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "MSigDB_Hallmark_2020",
    "Reactome_2022"
]

# Define colors for each gene set
target_gene_set_colors = {
    "KEGG_2021_Human": "#1f77b4",
    "GO_Biological_Process_2023": "#ff7f0e",
    "GO_Molecular_Function_2023": "#2ca02c",
    "MSigDB_Hallmark_2020": "#d62728",
    "Reactome_2022": "#8c564b"
}

# Dictionary to store all target enrichment results
target_enrichment_results = {}

# Function to extract targets for a list of TFs from both response groups
def extract_targets_for_tfs(tf_list, links_good, links_bad, cluster_name="Malignant"):
    all_targets = set()
    
    for tf in tf_list:
        # Extract targets from good response
        if cluster_name in links_good.filtered_links:
            good_targets = links_good.filtered_links[cluster_name].query("source == @tf")['target'].tolist()
            all_targets.update(good_targets)
        
        # Extract targets from bad response
        if cluster_name in links_bad.filtered_links:
            bad_targets = links_bad.filtered_links[cluster_name].query("source == @tf")['target'].tolist()
            all_targets.update(bad_targets)
    
    return list(all_targets)

# Process each TF group
for group_name, tf_set in tf_groups_dict.items():
    print(f"\nProcessing TF group: {group_name} ({len(tf_set)} TFs)")
    
    # Extract targets for all TFs in this group
    targets = extract_targets_for_tfs(list(tf_set), links_treat, links_untreat)
    
    print(f"  Found {len(targets)} unique targets")
    
    # Save target list
    targets_df = pd.DataFrame({'Target_Genes': sorted(targets)})
    targets_df.to_csv(f"{save_path_metadata}/targets_{group_name}_tfs.csv", index=False)
    
    # Perform enrichment analysis on targets
    if len(targets) >= 3:  # Minimum genes for enrichment
        try:
            res_targets = gp.enrichr(
                gene_list=sorted(targets),
                gene_sets=target_gene_sets,
                cutoff=0.05
            )
            
            # Save full results
            res_targets.results.to_csv(f"{save_path_metadata}/enrichment_targets_{group_name}_tfs.csv", sep=",", index=False)
            target_enrichment_results[group_name] = res_targets.results
            
            # Filter significant results
            res_filtered = res_targets.results[res_targets.results["Adjusted P-value"] < 0.05]
            
            if not res_filtered.empty:
                print(f"    Found {len(res_filtered)} significant enrichment terms")
                
                # Create barplot for targets
                barplot(
                    res_filtered,
                    column="Adjusted P-value",
                    group="Gene_set",
                    size=10,
                    top_term=15,
                    title=f"Target Enrichment Analysis - {group_name.replace('_', ' ').title()} TFs",
                    figsize=(12, 15),
                    color=target_gene_set_colors,
                    ofname=f"{save_path_figures}/enrichment_targets_{group_name}_tfs.pdf"
                )
                
                # Create dotplot for Hallmark pathways
                res_hallmark = res_filtered[res_filtered["Gene_set"] == "MSigDB_Hallmark_2020"]
                if not res_hallmark.empty:
                    dotplot(res_hallmark,
                            title=f"Hallmark Pathways - {group_name.replace('_', ' ').title()} TF Targets",
                            ofname=f"{save_path_figures}/enrichment_targets_hallmark_{group_name}_tfs.pdf")
            else:
                print(f"    No significant enrichment terms found")
                
        except Exception as e:
            print(f"    Error in target enrichment for {group_name}: {e}")
    else:
        print(f"    Too few targets ({len(targets)}) for enrichment analysis")

# Create comparison analysis of target enrichments
print("\n=== CREATING TARGET ENRICHMENT COMPARISON ===")

# Combine all target enrichment results for comparison
all_target_enrichment_data = []
target_group_labels = {
    'common_all': 'Common All Targets',
    'moderate_minimal': 'Moderate/Minimal Targets', 
    'poor_untreated': 'Poor/Untreated Targets',
    'good_unique': 'Good Unique Targets',
    'bad_unique': 'Bad Unique Targets',
    'intersection': 'Intersection Targets'
}

for group_key, group_label in target_group_labels.items():
    if group_key in target_enrichment_results:
        df_group = target_enrichment_results[group_key].copy()
        df_group = df_group[df_group["Adjusted P-value"] < 0.05]
        df_group = df_group[df_group["Gene_set"] == "MSigDB_Hallmark_2020"]  # Focus on Hallmark
        df_group['Target_Group'] = group_label
        all_target_enrichment_data.append(df_group)

if all_target_enrichment_data:
    # Combine all target enrichment data
    df_all_target_enrichment = pd.concat(all_target_enrichment_data, ignore_index=True)
    
    # Save combined results
    df_all_target_enrichment.to_csv(f"{save_path_metadata}/all_tf_target_groups_enrichment_comparison.csv", sep=",", index=False)
    
    if not df_all_target_enrichment.empty:
        # Get top terms for comparison
        top_terms = df_all_target_enrichment.nsmallest(30, 'Adjusted P-value')['Term'].unique()
        df_plot = df_all_target_enrichment[df_all_target_enrichment['Term'].isin(top_terms)]
        
        # Create pivot table for heatmap
        heatmap_data = df_plot.pivot_table(
            index='Term', 
            columns='Target_Group', 
            values='Adjusted P-value',
            fill_value=1.0
        )
        
        # Transform p-values for better visualization (-log10)
        heatmap_data_log = -np.log10(heatmap_data)
        
        # Create clustermap
        plt.figure(figsize=(15, 12))
        sns.clustermap(
            heatmap_data_log,
            cmap='Blues',
            annot=False,
            fmt='.2f',
            figsize=(18, 12),
            row_cluster=True,
            col_cluster=True,
            linewidths=0.5,
            cbar_kws={'label': '-log10(Adjusted P-value)'}
        )
        plt.tight_layout()
        plt.savefig(f"{save_path_figures}/tf_target_groups_enrichment_comparison_heatmap.pdf", 
                    bbox_inches='tight', dpi=300)
        plt.show()
        
        print(f"Created target enrichment comparison with {len(df_all_target_enrichment)} significant terms")
        
        # Create a summary table of top pathways per group
        summary_data = []
        for group in heatmap_data.columns:
            group_data = heatmap_data[group].dropna().nsmallest(5)  # Top 5 most significant
            for pathway, pvalue in group_data.items():
                summary_data.append({
                    'Target_Group': group,
                    'Pathway': pathway.replace('HALLMARK_', '').replace('_', ' '),
                    'Adj_P_value': pvalue,
                    'Neg_Log_P': -np.log10(pvalue)
                })
        
        df_summary = pd.DataFrame(summary_data)
        df_summary.to_csv(f"{save_path_metadata}/tf_target_groups_top_pathways_summary.csv", index=False)
        
        # Create a summary barplot
        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        axes = axes.flatten()
        
        for idx, group in enumerate(heatmap_data.columns[:6]):  # Plot up to 6 groups
            if idx < len(axes):
                ax = axes[idx]
                group_data = heatmap_data[group].dropna().nsmallest(10)
                
                if not group_data.empty:
                    y_pos = range(len(group_data))
                    pathway_labels = [p.replace('HALLMARK_', '').replace('_', ' ')[:30] + ('...' if len(p) > 30 else '') 
                                    for p in group_data.index]
                    
                    bars = ax.barh(y_pos, -np.log10(group_data.values), color='steelblue', alpha=0.7)
                    ax.set_yticks(y_pos)
                    ax.set_yticklabels(pathway_labels, fontsize=8)
                    ax.set_xlabel('-log10(Adj P-value)', fontsize=10)
                    ax.set_title(f'{group}', fontsize=12, fontweight='bold')
                    ax.grid(axis='x', alpha=0.3)
                    
                    # Add value labels
                    for i, (pathway, pvalue) in enumerate(group_data.items()):
                        ax.text(-np.log10(pvalue) + 0.1, i, f'{pvalue:.2e}', 
                               va='center', ha='left', fontsize=7)
        
        # Hide empty subplots
        for idx in range(len(heatmap_data.columns), len(axes)):
            axes[idx].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(f"{save_path_figures}/tf_target_groups_top_pathways_summary.pdf", 
                    bbox_inches='tight', dpi=300)
        plt.show()
        
    else:
        print("No significant target enrichment terms found for comparison")
else:
    print("No target enrichment results available for comparison")

print(f"\n=== TARGET ENRICHMENT ANALYSIS COMPLETED ===")
print("Target lists and enrichment results saved in:")
for group_name in tf_groups_dict.keys():
    print(f"  - Targets: {save_path_metadata}/targets_{group_name}_tfs.csv")
    if group_name in target_enrichment_results:
        print(f"  - Enrichment: {save_path_metadata}/enrichment_targets_{group_name}_tfs.csv")
# %%
# Extract targets from links_malignant for the TF groups and perform enrichment analysis
print("=== EXTRACTING TARGETS FROM MALIGNANT LINKS FOR TF GROUPS ===\n")

# Define the clusters we're interested in based on links_malignant
malignant_clusters = links_malignant.merged_score['cluster'].unique()
print(f"Available malignant clusters: {malignant_clusters}")

# Group clusters by response type for malignant analysis
poor_untreated_clusters_malignant = [c for c in malignant_clusters if c in ['Poor response', 'Untreated']]
minimal_moderate_clusters_malignant = [c for c in malignant_clusters if c in ['Minimal response', 'Moderate response']]

print(f"Poor/Untreated clusters: {poor_untreated_clusters_malignant}")
print(f"Minimal/Moderate clusters: {minimal_moderate_clusters_malignant}")

# Function to extract targets for TFs from malignant links
def extract_malignant_targets_for_tfs(tf_list, links_malignant, clusters_list):
    """Extract targets for TFs from specified clusters in malignant links"""
    all_targets = set()
    
    for cluster in clusters_list:
        if cluster in links_malignant.filtered_links:
            cluster_links = links_malignant.filtered_links[cluster]
            for tf in tf_list:
                tf_targets = cluster_links.query("source == @tf")['target'].tolist()
                all_targets.update(tf_targets)
    
    return list(all_targets)

# Dictionary to store malignant target results
malignant_target_results = {}

# 1. Extract targets for TFs from poor/untreated response clusters
if poor_untreated_tfs and poor_untreated_clusters_malignant:
    print(f"\n1. Extracting targets for poor/untreated TFs from malignant clusters...")
    poor_untreated_targets_malignant = extract_malignant_targets_for_tfs(
        list(poor_untreated_tfs), links_malignant, poor_untreated_clusters_malignant
    )
    malignant_target_results['poor_untreated'] = poor_untreated_targets_malignant
    print(f"   Found {len(poor_untreated_targets_malignant)} targets")

# 2. Extract targets for TFs from minimal/moderate response clusters
if minimal_moderate_clusters_malignant and 'moderate_minimal_tfs' in locals():
    print(f"\n2. Extracting targets for minimal/moderate TFs from malignant clusters...")
    minimal_moderate_targets_malignant = extract_malignant_targets_for_tfs(
        list(moderate_minimal_tfs), links_malignant, minimal_moderate_clusters_malignant
    )
    malignant_target_results['minimal_moderate'] = minimal_moderate_targets_malignant
    print(f"   Found {len(minimal_moderate_targets_malignant)} targets")

# 3. Extract targets for common TFs across all malignant clusters
if common_tfs_all:
    print(f"\n3. Extracting targets for common TFs from all malignant clusters...")
    common_targets_malignant = extract_malignant_targets_for_tfs(
        list(common_tfs_all), links_malignant, list(malignant_clusters)
    )
    malignant_target_results['common_all'] = common_targets_malignant
    print(f"   Found {len(common_targets_malignant)} targets")

# 4. Extract targets for good/bad unique TFs (from earlier binary analysis)
if 'good_unique_tfs' in locals() and good_unique_tfs:
    print(f"\n4. Extracting targets for good response unique TFs...")
    good_unique_targets_malignant = extract_malignant_targets_for_tfs(
        list(good_unique_tfs), links_malignant, minimal_moderate_clusters_malignant
    )
    malignant_target_results['good_unique'] = good_unique_targets_malignant
    print(f"   Found {len(good_unique_targets_malignant)} targets")

if 'bad_unique_tfs' in locals() and bad_unique_tfs:
    print(f"\n5. Extracting targets for bad response unique TFs...")
    bad_unique_targets_malignant = extract_malignant_targets_for_tfs(
        list(bad_unique_tfs), links_malignant, poor_untreated_clusters_malignant
    )
    malignant_target_results['bad_unique'] = bad_unique_targets_malignant
    print(f"   Found {len(bad_unique_targets_malignant)} targets")

# Compare targets between poor/untreated and minimal/moderate groups
if 'poor_untreated' in malignant_target_results and 'minimal_moderate' in malignant_target_results:
    poor_untreated_targets_set = set(malignant_target_results['poor_untreated'])
    minimal_moderate_targets_set = set(malignant_target_results['minimal_moderate'])
    
    # Find common and unique targets
    common_targets_malignant_groups = poor_untreated_targets_set & minimal_moderate_targets_set
    poor_untreated_unique_targets = poor_untreated_targets_set - minimal_moderate_targets_set
    minimal_moderate_unique_targets = minimal_moderate_targets_set - poor_untreated_targets_set
    
    print(f"\n=== TARGET COMPARISON BETWEEN RESPONSE GROUPS ===")
    print(f"Poor/Untreated targets: {len(poor_untreated_targets_set)}")
    print(f"Minimal/Moderate targets: {len(minimal_moderate_targets_set)}")
    print(f"Common targets: {len(common_targets_malignant_groups)}")
    print(f"Poor/Untreated unique targets: {len(poor_untreated_unique_targets)}")
    print(f"Minimal/Moderate unique targets: {len(minimal_moderate_unique_targets)}")
    
    # Add these to malignant_target_results for enrichment
    malignant_target_results['common_response_groups'] = list(common_targets_malignant_groups)
    malignant_target_results['poor_untreated_unique'] = list(poor_untreated_unique_targets)
    malignant_target_results['minimal_moderate_unique'] = list(minimal_moderate_unique_targets)

# Save target lists
for group_name, targets in malignant_target_results.items():
    targets_df = pd.DataFrame({'Target_Genes': sorted(targets)})
    targets_df.to_csv(f"{save_path_metadata}/malignant_targets_{group_name}.csv", index=False)
    print(f"Saved {len(targets)} targets for {group_name}")

# Gene sets for enrichment analysis
malignant_gene_sets = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023", 
    "MSigDB_Hallmark_2020",
    "Reactome_2022"
]

malignant_gene_set_colors = {
    "KEGG_2021_Human": "#1f77b4",
    "GO_Biological_Process_2023": "#ff7f0e",
    "GO_Molecular_Function_2023": "#2ca02c",
    "MSigDB_Hallmark_2020": "#d62728",
    "Reactome_2022": "#8c564b"
}

# Perform enrichment analysis on malignant targets
malignant_enrichment_results = {}

print(f"\n=== PERFORMING ENRICHMENT ANALYSIS ON MALIGNANT TARGETS ===")

for group_name, targets in malignant_target_results.items():
    if len(targets) >= 5:  # Minimum genes for meaningful enrichment
        print(f"\nEnriching targets for {group_name} ({len(targets)} targets)...")
        
        try:
            res_malignant = gp.enrichr(
                gene_list=sorted(targets),
                gene_sets=malignant_gene_sets,
                cutoff=0.05
            )
            
            # Save results
            res_malignant.results.to_csv(
                f"{save_path_metadata}/malignant_enrichment_targets_{group_name}.csv", 
                sep=",", index=False
            )
            malignant_enrichment_results[group_name] = res_malignant.results
            
            # Filter significant results
            res_filtered = res_malignant.results[res_malignant.results["Adjusted P-value"] < 0.05]
            
            if not res_filtered.empty:
                print(f"  Found {len(res_filtered)} significant terms")
                
                # Create barplot
                barplot(
                    res_filtered,
                    column="Adjusted P-value",
                    group="Gene_set",
                    size=10,
                    top_term=15,
                    title=f"Malignant Target Enrichment - {group_name.replace('_', ' ').title()}",
                    figsize=(12, 15),
                    color=malignant_gene_set_colors,
                    ofname=f"{save_path_figures}/malignant_enrichment_targets_{group_name}.pdf"
                )
                
                # Create dotplot for Hallmark pathways
                res_hallmark = res_filtered[res_filtered["Gene_set"] == "MSigDB_Hallmark_2020"]
                if not res_hallmark.empty:
                    dotplot(res_hallmark,
                           title=f"Hallmark - {group_name.replace('_', ' ').title()} Malignant Targets",
                           ofname=f"{save_path_figures}/malignant_hallmark_targets_{group_name}.pdf")
            else:
                print(f"  No significant enrichment found")
                
        except Exception as e:
            print(f"  Error in enrichment for {group_name}: {e}")
    else:
        print(f"\nSkipping {group_name}: too few targets ({len(targets)})")

# Create comparative analysis of malignant target enrichments
if malignant_enrichment_results:
    print(f"\n=== CREATING MALIGNANT TARGET ENRICHMENT COMPARISON ===")
    
    # Combine results focusing on Hallmark pathways
    all_malignant_enrichment = []
    
    for group_name, results_df in malignant_enrichment_results.items():
        df_group = results_df.copy()
        df_group = df_group[df_group["Adjusted P-value"] < 0.05]
        df_group = df_group[df_group["Gene_set"] == "MSigDB_Hallmark_2020"]
        df_group['Target_Group'] = group_name.replace('_', ' ').title()
        all_malignant_enrichment.append(df_group)
    
    if all_malignant_enrichment:
        df_combined_malignant = pd.concat(all_malignant_enrichment, ignore_index=True)
        df_combined_malignant.to_csv(
            f"{save_path_metadata}/malignant_all_target_enrichment_comparison.csv", 
            sep=",", index=False
        )
        
        # Create comparison heatmap
        if not df_combined_malignant.empty:
            # Get top pathways
            top_pathways = df_combined_malignant.nsmallest(25, 'Adjusted P-value')['Term'].unique()
            df_heatmap = df_combined_malignant[df_combined_malignant['Term'].isin(top_pathways)]
            
            # Create pivot table
            heatmap_matrix = df_heatmap.pivot_table(
                index='Term',
                columns='Target_Group', 
                values='Adjusted P-value',
                fill_value=1.0
            )
            
            # Transform to -log10 for better visualization
            heatmap_matrix_log = -np.log10(heatmap_matrix)
            
            # Create clustermap
            plt.figure(figsize=(12, 10))
            sns.clustermap(
                heatmap_matrix_log,
                cmap='Reds',
                annot=False,
                figsize=(15, 12),
                row_cluster=True,
                col_cluster=True,
                linewidths=0.5,
                cbar_kws={'label': '-log10(Adjusted P-value)'}
            )
            plt.tight_layout()
            plt.savefig(f"{save_path_figures}/malignant_target_enrichment_comparison_heatmap.pdf",
                       bbox_inches='tight', dpi=300)
            plt.show()
            
            print(f"Created comparison heatmap with {len(df_combined_malignant)} significant terms")
            
            # Create individual comparison plots for key contrasts
            key_comparisons = [
                ('Poor Untreated', 'Minimal Moderate'),
                ('Poor Untreated Unique', 'Minimal Moderate Unique'),
                ('Common All', 'Common Response Groups')
            ]
            
            for comp1, comp2 in key_comparisons:
                if comp1 in heatmap_matrix.columns and comp2 in heatmap_matrix.columns:
                    fig, ax = plt.subplots(figsize=(6, 6))
                    
                    # Get data for both groups
                    comp1_data = heatmap_matrix[comp1].dropna()
                    comp2_data = heatmap_matrix[comp2].dropna()
                    
                    # Find common pathways
                    common_pathways = set(comp1_data.index) & set(comp2_data.index)
                    
                    if common_pathways:
                        common_pathways = list(common_pathways)
                        x_vals = [-np.log10(comp1_data.loc[common_pathways])]
                        y_vals = [-np.log10(comp2_data.loc[common_pathways])]
                        
                        ax.scatter(x_vals, y_vals, alpha=0.7, s=50)
                        
                        # Add pathway labels for most significant points
                        for pathway in common_pathways[:10]:  # Top 10 most significant
                            x = -np.log10(comp1_data.loc[pathway])
                            y = -np.log10(comp2_data.loc[pathway])
                            texts = []
                            texts.append(ax.text(x, y, pathway.replace('HALLMARK_', '').replace('_', ' ')[:20], fontsize=8, alpha=0.8))
                            
                        # Use adjusttext to fix overlapping labels
                        adjust_text.adjust_text(texts=texts, ax=ax, arrowprops=dict(arrowstyle='->', color='r', lw=0.5))
                        # ax.plot([0, max(x_vals)], [0, max(x_vals)], 'r--', alpha=0.5)
                        ax.set_xlabel(f'-log10(Adj P-value) {comp1}')
                        ax.set_ylabel(f'-log10(Adj P-value) {comp2}')
                        ax.set_title(f'Pathway Enrichment Comparison\n{comp1} vs {comp2}')
                        ax.grid(False)
                        sns.despine(ax=ax)
                        
                        plt.tight_layout()
                        plt.savefig(f"{save_path_figures}/malignant_pathway_comparison_{comp1.replace(' ', '_').lower()}_vs_{comp2.replace(' ', '_').lower()}.pdf",
                                   bbox_inches='tight', dpi=300)
                        plt.show()

print(f"\n=== MALIGNANT TARGET ENRICHMENT ANALYSIS COMPLETED ===")
print("Results saved:")
for group_name in malignant_target_results.keys():
    print(f"  - Targets: {save_path_metadata}/malignant_targets_{group_name}.csv")
    if group_name in malignant_enrichment_results:
        print(f"  - Enrichment: {save_path_metadata}/malignant_enrichment_targets_{group_name}.csv")
# %%
genes=heatmap_matrix.index.unique().tolist()
# Extract genes and metadata from adata_epithelial
gene_expression_data = sc.get.obs_df(adata_epithelial, keys=genes, use_raw=False)
metadata = adata_epithelial.obs[['merged_response', 'Level 2 Annotation']].copy()

# Combine expression data with metadata
combined_data = pd.concat([metadata, gene_expression_data], axis=1)

# Group by response and cell type, then calculate mean expression
grouped_data = combined_data.groupby(['merged_response', 'Level 2 Annotation'])[genes].mean()

# Reset index to make response and cell type regular columns
grouped_data = grouped_data.reset_index()

# Create a unique identifier for each group
grouped_data['group_id'] = grouped_data['merged_response'].astype(str) + '_' + grouped_data['Level 2 Annotation'].astype(str)

# Set group_id as index and transpose for clustermap
heatmap_data = grouped_data.set_index('group_id')[genes].T

# Standardize the expression values (z-score normalization)
heatmap_data_standardized = heatmap_data.apply(zscore, axis=1)

# Create clustermap
plt.figure(figsize=(8, 8))
sns.set(font_scale=1.0)
clustermap = sns.clustermap(heatmap_data_standardized, 
                          cmap='RdBu_r', 
                          center=0,
                          figsize=(12, 15),
                          row_cluster=True, 
                          col_cluster=True,
                          linewidths=0.5,
                          cbar_kws={'label': 'Standardized Expression (Z-score)'})

plt.suptitle('Gene Expression Clustermap by Response and Cell Type', y=0.98, fontsize=14)
plt.tight_layout()
plt.savefig(f"{save_path_figures}/gene_expression_clustermap_response_celltype.pdf", bbox_inches='tight', dpi=300)
plt.show()
# %%
sc.tl.rank_genes_groups(adata_malignant, 'merged_response', method='wilcoxon', use_raw=False)
# %%
df_rank_malignant=sc.get.rank_genes_groups_df(adata_malignant, group=None,pval_cutoff=0.05)
# %%
df_rank_malignant_tf_filtered = df_rank_malignant[df_rank_malignant['names'].isin(tf_list)]
df_rank_malignant_tf_filtered
# %%
# Get top 10 TFs for each response group and create heatmap
top_n_tfs = 10

# Get top TFs for each group based on absolute log fold change
df_mm_mr = df_rank_malignant_tf_filtered[df_rank_malignant_tf_filtered['group'] == 'MM/MR'].nlargest(top_n_tfs, 'logfoldchanges')
df_u_pr = df_rank_malignant_tf_filtered[df_rank_malignant_tf_filtered['group'] == 'U/PR'].nlargest(top_n_tfs, 'logfoldchanges')

# Combine the top TFs from both groups
top_tfs_combined = pd.concat([df_mm_mr, df_u_pr])

# Pivot the data to create a matrix for heatmap
heatmap_data = top_tfs_combined.pivot(index='names', columns='group', values='logfoldchanges')

# Fill NaN values with 0
heatmap_data = heatmap_data.fillna(0)

# Create heatmap
plt.figure(figsize=(6, 8))
sns.set(font_scale=1.2)
sns.clustermap(heatmap_data, annot=True, cmap='RdBu_r', center=0, 
            fmt='.2f', linewidths=0.5)
# plt.title(f'Top {top_n_tfs} TFs Log Fold Changes by Response Group')
plt.xlabel('Response Group')
plt.ylabel('Transcription Factors')
plt.tight_layout()
# plt.savefig(f"{save_path_figures}/top{top_n_tfs}_TFs_logfoldchange_heatmap.pdf", bbox_inches='tight', dpi=300)
plt.show()

# %%
sc.pl.heatmap(adata_malignant, var_names=top_tfs_combined['names'].unique(),swap_axes=True
              ,use_raw=False,groupby='merged_response', cmap='RdBu_r', standard_scale='var',
              dendrogram=False, figsize=(7, 9))
# %%
