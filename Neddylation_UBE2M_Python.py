#!/usr/bin/env python
# coding: utf-8

# # 安装环境相关的包

# In[1]:


#pip install jupyter_contrib_nbextensions && jupyter contrib nbextensions install
#安装变量显示器Variable Inspector


# In[2]:


#pip install scanpy -i https://pypi.tuna.tsinghua.edu.cn/simple


# In[ ]:


#%pip uninstall matplotlib -y
#%pip uninstall matplotlib-inline -y
#%pip install matplotlib
#%pip install matplotlib-inline


# In[1]:


pip uninstall infercnvpy -y


# In[3]:


pip install infercnvpy==0.4.3 #0.4.3能成功


# In[4]:


#pip install ipywidgets


# In[1]:


#pip install infercnvpy -i https://pypi.tuna.tsinghua.edu.cn/simple some-package


# In[6]:


#$pip install git+https://github.com/icbi-lab/infercnvpy.git@main
#没成功


# In[5]:


#pip install ipykernel


# In[1]:


#pip install diopy


# In[21]:


pip install numpy==1.23.0


# In[5]:


pip install pandas==1.5.3


# In[3]:


pip install  anndata==0.10.5


# In[6]:


pip install pyarrow


# In[58]:


pip install openpyxl


# In[1]:


get_ipython().system('pip list')


# # 导入环境的包 ##

# In[1]:


import os 
print(os.path.abspath('.'))
#看文件保存位置


# In[2]:


import scanpy as sc
#导入scanpy包


# In[56]:


import pandas as pd
#导入pandas 数据框 包


# In[5]:


import loompy


# In[2]:


import diopy
#这是R和python数据转换的包


# In[3]:


import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import warnings

warnings.simplefilter("ignore")

sc.settings.set_figure_params(figsize=(5, 5))


# In[31]:


sc.settings.set_figure_params(figsize=(12, 12))


# In[4]:


import infercnvpy as cnv


# # 一些常用代码

# In[272]:


##保存obs/metadata数据
import pandas as pd
adata.obs.to_excel("Epithelial_T_nedd.xlsx",engine='openpyxl')


# In[276]:


#读取metadata数据
metadata_allcell = pd.read_csv('metadata_allcell.csv',index_col=0)
metadata_allcell


# In[69]:


##文件保存路径/位置
import os 
print(os.path.abspath('.'))


# In[ ]:


##更改文件保存路径/位置
import os
os.chdir('/home/wling_32/jupyter_home/NGmetadata/cNMF')


# In[203]:


##保存文件
sc.write('Epithelial_T_afterumap.h5ad', adata)


# In[222]:


##读取文件
adata=sc.read('/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T_afterumap.h5ad')
#allcell=sc.read('/home/wling_32/NGmetadata/Rawdata/allcell.h5ad')


# In[269]:


##删除某一列
#adata.obs = adata.obs.drop(columns=["Sox_aucell"]) 
#adata.obs = adata.obs.drop(columns=['cNMF_Cluster1',"cNMF_Cluster2","cNMF_Cluster3","cNMF_Cluster4","cNMF_Cluster5","cNMF_Cluster6","cNMF_Cluster7","cNMF_Cluster"]) #删除anndata的特定列


# In[23]:


#提取anndata的umap数据
UMAP1=adata.obsm["X_umap"][:,0]
UMAP2=adata.obsm["X_umap"][:,1]
X_umap=pd.DataFrame({'UMAP1': UMAP1, 'UMAP2': UMAP2})
X_umap
X_umap.to_csv('Epithelial_T_umap.csv', index=False)


# In[5]:


#adata  = object
adata  = diopy.input.read_h5(file = '/home/wling_32/jupyter_home/NGmetadata/sce_Epithelial.h5')
#这是另外一种方法dior 转换R和Python文件


# In[4]:


#print(sc.__file__)


# In[ ]:


#data=sc.read('/home/wling_32/NGmetadata/Rawdata/allcell.h5ad')
#导入allcell的文件，已经经过seuratdisk转为h5ad格式方便python读取


# In[223]:


print(allcell)


# In[ ]:


import matplotlib.pyplot as plt
#sc.pp.filter_cells(data, min_genes=200)
#sc.pp.filter_genes(data, min_cells=3)
sc.pl.highest_expr_genes(data, n_top=20)
#plt.savefig("1.pdf")


# In[ ]:


# 运行代码部分


# In[ ]:


sc.pl.umap(data, color=['UBE2M'])


# In[ ]:


sc.pl.umap(allcell, color='cell.type', legend_loc='on data', title='', frameon=False, save='.pdf')


# In[ ]:


category_counts = data.obs['cell.type'].value_counts()
print(category_counts)


# In[ ]:


data_epi = data[data.obs["cell.type"] == "Epithelial"]


# In[ ]:


print(data_epi)


# In[ ]:


sc.write('Epithelial.h5ad', data_epi)


# In[3]:


import scanpy as sc
Epithelial=sc.read('/home/wling_32/jupyter_home/NGmetadata/Epithelial.h5ad')
print(Epithelial)


# In[4]:


Epithelial.obs['sample.origin'].value_counts()


# In[1]:


#import os
#os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"


# In[11]:


import sys
print(sys.version)


# # InferCNV代码运行

# 
# 

# In[ ]:


adata  = diopy.input.read_h5(file = '/home/wling_32/jupyter_home/NGmetadata/sce_Epithelial.h5')


# In[6]:


adata.obs['sample.origin'].value_counts()


# In[7]:


Epithelial = adata


# In[8]:


print(Epithelial)


# In[9]:


Epithelial.obs['Sample_type'] = Epithelial.obs['sample.origin'].apply(
    lambda x: 'Normal cell' if x in ['Normal'] else 'Tumor cell'
)


# In[12]:


print(Epithelial)
Epithelial.obs['Sample_type'].value_counts()


# In[13]:


Epithelial.obs['Sample_type'].unique()


# In[16]:


gene_annotations = pd.read_csv('geneannotation.csv')
gene_annotations.head()


# In[17]:


# 将染色体位置信息添加到 adata.var
Epithelial.var['chromosome'] = ''  # 初始化染色体列
Epithelial.var['start'] = 0   # 初始化起始位置列
Epithelial.var['end'] = 0     # 初始化结束位置列


# In[18]:


# 匹配基因ID并添加染色体位置信息
for gene_id in Epithelial.var_names:
    if gene_id in gene_annotations['gene_symbol'].values:
        # 找到匹配的基因并获取其染色体位置信息
        chromosome, start, end = gene_annotations.loc[gene_annotations['gene_symbol'] == gene_id, ['chromosome', 'start', 'end']].values[0]
        Epithelial.var.loc[gene_id, 'chromosome'] = chromosome
        Epithelial.var.loc[gene_id, 'start'] = start
        Epithelial.var.loc[gene_id, 'end'] = end


# In[22]:


#去掉NA   这部没做
#Epithelial = Epithelial[: , Epithelial.var.chromosome.notna()]
#Epithelial.var


# In[23]:


#run infercnvpy
cnv.tl.infercnv(
    Epithelial,
    reference_key="Sample_type",
    reference_cat=[
        "Normal cell",
    ],
    window_size=250,
)


# In[65]:


cnv.pl.chromosome_heatmap(Epithelial, groupby="Sample_type",figsize=(16, 6),save="TvsN_cnv_heatmap.pdf")


# In[34]:


#Clustering by CNV profiles and identifying tumor cells
cnv.tl.pca(Epithelial)
cnv.pp.neighbors(Epithelial)
cnv.tl.leiden(Epithelial)

cnv.pl.chromosome_heatmap(Epithelial, groupby="cnv_leiden", dendrogram=True,save="cnvcluster_heatmap.pdf")


# In[66]:


cnv.pl.chromosome_heatmap(Epithelial, groupby="cnv_leiden", dendrogram=True,figsize=(16, 30),save="cnvcluster_heatmap.pdf")


# In[35]:


cnv.tl.umap(Epithelial)
cnv.tl.cnv_score(Epithelial)


# In[74]:


#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
#ax4.axis("off")
cnv.pl.umap(
    Epithelial,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    legend_fontsize= 'x-large',
    #ax=ax1,  #指定在哪一行,
    edges_width=1,
    show=True,
    save="cnvcluster_umap.pdf",
)


# In[39]:


cnv.pl.umap(Epithelial, 
            color="cnv_score", 
            #ax=ax2, 
            show=True,
            save="cnvscore_umap.pdf")


# In[40]:


cnv.pl.umap(Epithelial, 
            color="Sample_type", 
            #ax=ax3,
            save="sampletype_umap.pdf")  #前面这三幅图主要是让我们通过肉眼观察这个infercnvcluster从而判断哪些属于恶性细胞


# In[49]:


Epithelial


# In[44]:


Epithelial.obs['cnv_leiden'].value_counts()


# In[60]:


#这里先保存Epithelial文件
sc.write('Epithelial.h5ad', Epithelial)


# In[59]:


#导出anndata数据
#diopy.output.write_h5(Epithelial, file = 'sce_Epithelial_2.h5')##报错 Error: setting an array element with a sequence.
#sc.write('Epithelial.h5ad', Epithelial)  在R里无法导入 算了
#干脆导入metadata也就是obs吧
Epithelial.obs.to_excel("epithelial_obs.xlsx",engine='openpyxl')


# In[76]:


#方法1  通过前面这三幅图，让我们肉眼观察这个infercnvcluster，从而判断哪些属于恶性细胞   --------用了 选0 1 25作为tumor中的normal
Epithelial.obs["cnv_status"] = "Tumor"
Epithelial.obs.loc[
    Epithelial.obs["cnv_leiden"].isin(["0", "1", "25"]) & Epithelial.obs["Sample_type"].isin(["Tumor cell"]), "cnv_status"
] = "Normal"
Epithelial.obs.loc[
    Epithelial.obs["Sample_type"].isin(["Normal cell"]), "cnv_status"
] = "Normal"
Epithelial.obs["cnv_status"].value_counts()


# In[ ]:


#方法2  找normal细胞的cnvscore中，没用
Epithelial.obs["cnv_status"] = "normal"
Epithelial.obs.loc[adata.obs["cnv_score"]>0.008, "cnv_status"] = (
    "tumor"
)


# In[77]:


cnv.pl.chromosome_heatmap(Epithelial, groupby="cnv_status",figsize=(16, 6),save="cnv_status_heatmap.pdf")


# In[78]:


sc.write('Epithelial.h5ad', Epithelial)


# # cNMF

# ## 提取cnv_status为Tumor的anndata

# In[22]:


import os 
print(os.path.abspath('.'))
#看文件保存位置


# In[21]:


import os
os.chdir('/home/wling_32/jupyter_home/NGmetadata/cNMF')
##更改文件保存位置


# In[88]:


Epithelial.obs["cnv_status"].value_counts()


# In[89]:


Epithelial_T = Epithelial[Epithelial.obs["cnv_status"] == "Tumor"]


# In[90]:


print(Epithelial_T)


# In[91]:


sc.write('Epithelial_T.h5ad', Epithelial_T)


# ## 安装以及加载包

# In[2]:


pip install scanpy


# In[1]:


pip install cnmf


# In[ ]:


pip install importlib-metadata


# In[148]:


pip install openpyxl #---如果pandas导不excel文件出来，要装openpyxl


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')

import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib.pyplot as plt
from IPython.display import Image
import scanpy as sc
from cnmf import cNMF
np.random.seed(14)


# ## 初步运行cNMF

# In[56]:


adata=sc.read('/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T.h5ad')


# In[57]:


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, min_counts=200)  
sc.pp.filter_genes(adata, min_cells=3)


# In[150]:


adata


# In[58]:


## plot log10 # counts per cell
plt.hist(np.log10(adata.obs['n_counts']), bins=100)
_ = plt.xlabel('log10 Counts Per cell')
_ = plt.ylabel('# Cells')


# In[153]:


numiter=20 # Number of NMF replicates. Set this to a larger value ~200 for real data. We set this to a relatively low value here for illustration at a faster speed
numhvgenes=2000 ## Number of over-dispersed genes to use for running the actual factorizations

## Results will be saved to [output_directory]/[run_name] which in this example is example_PBMC/cNMF/pbmc_cNMF
output_directory = '/home/wling_32/jupyter_home/NGmetadata/cNMF'
if not os.path.exists(output_directory):
    os.mkdir(output_directory)
run_name = 'Epithelial_T_cNMF'

## Specify the Ks to use as a space separated list in this case "5 6 7 8 9 10"
K = ' '.join([str(i) for i in range(5,11)])

## To speed this up, you can run it for only K=7-8 with the option below
#K = ' '.join([str(i) for i in range(7,9)])

seed = 14 ## Specify a seed pseudorandom number generation for reproducibility

## Path to the filtered counts dataset we output previously
countfn = '/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T.h5ad'


# In[154]:


## Initialize the cnmf object that will be used to run analyses
cnmf_obj = cNMF(output_dir=output_directory, name=run_name)


# In[13]:


## Prepare the data, I.e. subset to 2000 high-variance genes, and variance normalize
cnmf_obj.prepare(counts_fn=countfn, components=np.arange(5,11), n_iter=20, seed=14, num_highvar_genes=2000)


# In[14]:


## Specify that the jobs are being distributed over a single worker (total_workers=1) and then launch that worker
cnmf_obj.factorize(worker_i=0, total_workers=1)


# In[20]:


## Using GNU parallel
## This took 4 minutes in our testing
numworkers = 4
factorize_cmd = 'nohup parallel python ../cnmf.py factorize --output-dir /home/wling_32/jupyter_home/NGmetadata/cNMF/ --name Epithelial_T_cNMF --worker-index {} ::: 0 1 2 3'
print('Factorize command to simultaneously run factorization over %d cores using GNU parallel:\n%s' % (numworkers, factorize_cmd))
#!{factorize_cmd}


# In[16]:


cnmf_obj.combine()


# In[155]:


cnmf_obj.k_selection_plot(close_fig=False)  #7是最稳定的


# In[21]:


print('This saves the corresponding figure to the following file: %s' % cnmf_obj.paths['k_selection_plot'])


# In[22]:


kselect_plot_cmd = 'cnmf k_selection_plot --output-dir /home/wling_32/jupyter_home/NGmetadata/cNMF/ --name Epithelial_T_cNMF'
print('K selection plot command: %s' % kselect_plot_cmd)
#!{kselect_plot_cmd}


# ## 筛选k值与密度阈值

# In[157]:


#下一步计算给定 K 选择的共识解决方案。首先在没有任何异常值过滤的情况下运行它以查看它的样子。 
#将密度阈值设置为 >= 2.00（两个单位向量之间的最大可能距离）可确保不会过滤任何内容
selected_K = 7
density_threshold = 2.00


# In[24]:


cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)


# In[158]:


density_threshold = 0.15


# In[26]:


cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)


# In[31]:


get_ipython().system(' ls /home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T_cNMF')

#我们最感兴趣的是密度阈值为 0.15 的用法和gene_spectra_score文件：

#.Epithelial_T_cNMF.gene_spectra_score.k_7.dt_0_15.txt
#.Epithelial_T_cNMF.usages.k_7.dt_0_15.consensus.txt


# ## 加载 AnnData 对象，TPT 对其进行归一化，均值和方差归一化每个基因，运行 PCA 并运行 UMAP

# In[3]:


adata=sc.read('/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T.h5ad')


# In[4]:


adata


# In[5]:


## Obtain high variance genes that were used for cNMF as these were saved to a text file
hvgs = open('/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T_cNMF/Epithelial_T_cNMF.overdispersed_genes.txt').read().split('\n')


# In[6]:


sc.pp.normalize_per_cell(adata, counts_per_cell_after=10**4) ## TPT normalization


# In[7]:


## Set log-normalized data to the raw attribute of the AnnData object to make it easy to plot expression levels of individual genes.
## This does not log normalize the actual AnnData data matrix
adata.raw = sc.pp.log1p(adata.copy(), copy=True)


# In[8]:


## Subset out only the high-variance genes

adata = adata[:,hvgs]


# In[9]:


## Mean and variance normalize the genes

sc.pp.scale(adata)


# In[10]:


## Run PCA

sc.pp.pca(adata)


# In[11]:


## Make a scree plot to determine number of PCs to use for UMAP

sc.pl.pca_variance_ratio(adata, log=True)


# In[9]:


#pip install importlib_metadata==5.0.0--转成python3.10后就不用装了


# In[12]:


## Construct the nearest neighbor graph for UMAP

sc.pp.neighbors(adata, n_neighbors=50, n_pcs=15)


# In[13]:


## Run UMAP   这一步一直有问题-No module named 'importlib.metadata'  我转成python3.10内核
sc.tl.umap(adata)


# In[389]:


#sc.tl.leiden(adata) 这一步相当于聚类，咱们也不用聚类 都是cNMF聚类的


# In[35]:


sc.write('Epithelial_T_afterumap.h5ad', adata)


# ## 生成规范化的文件（主要是7种聚类中各自的topgene）

# In[159]:


usage_norm, gep_scores, gep_tpm, topgenes = cnmf_obj.load_results(K=selected_K, density_threshold=density_threshold)
usage_norm.columns = ['cNMF_Cluster%d' % i for i in usage_norm.columns]


# In[44]:


usage_file = cnmf_obj.paths['consensus_usages__txt'] % (selected_K, '0_15')
print(usage_file)


# In[45]:


gene_scores_file = cnmf_obj.paths['gene_spectra_score__txt'] % (selected_K, '0_15')
print(gene_scores_file)


# In[46]:


gene_tpm_file = cnmf_obj.paths['gene_spectra_tpm__txt'] % (selected_K, '0_15')
print(gene_tpm_file)


# In[160]:


top_genes = []
ngenes = 100
for gep in gep_scores.columns:
    top_genes.append(list(gep_scores.sort_values(by=gep, ascending=False).index[:ngenes]))
    
top_genes = pd.DataFrame(top_genes, index=gep_scores.columns).T
top_genes


# In[48]:


topgenes.head(20)


# In[51]:


pip install openpyxl


# In[52]:


#保存topgenes
import pandas as pd
topgenes.to_excel("topgenes.xlsx",engine='openpyxl')


# In[196]:


##删除某一列
#adata.obs = adata.obs.drop(columns=["index\t1\t2\t3\t4\t5\t6\t7"]) 
#adata.obs = adata.obs.drop(columns=['cNMF_Cluster1',"cNMF_Cluster2","cNMF_Cluster3","cNMF_Cluster4","cNMF_Cluster5","cNMF_Cluster6","cNMF_Cluster7","cNMF_Cluster"]) #删除anndata的特定列


# In[199]:


adata


# ## 绘制cNMF Umap图以及定义cNMF亚群

# In[198]:


#最普通版--  添加Epithelial_T_cNMF.usages.k_7.dt_0_15.consensus信息到adata中
#adata=sc.read('/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T_afterumap.h5ad')
new_columns = ['cNMF_Cluster1',"cNMF_Cluster2","cNMF_Cluster3","cNMF_Cluster4","cNMF_Cluster5","cNMF_Cluster6","cNMF_Cluster7"]
usage_norm = pd.read_csv('/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T_cNMF/Epithelial_T_cNMF.usages.k_7.dt_0_15.consensus.txt',sep='\t',index_col=0,)
usage_norm.columns = new_columns
usage_norm

adata.obs = pd.merge(left=adata.obs, right=usage_norm, how='left', left_index=True, right_index=True)


# In[200]:


columns_list = ['cNMF_Cluster1', 'cNMF_Cluster2', 'cNMF_Cluster3', 'cNMF_Cluster4', 'cNMF_Cluster5', 'cNMF_Cluster6', 'cNMF_Cluster7']
adata.obs['cNMF_Cluster'] = adata.obs[columns_list].idxmax(axis=1)


# In[ ]:


adata
sc.write('Epithelial_T_afterumap.h5ad', adata)


# In[204]:


sc.pl.umap(adata, color=usage_norm.columns,
           use_raw=True, ncols=3, vmin=0, vmax=1,save="cNMF_Cluster_umap.pdf")


# In[205]:


sc.pl.umap(adata, color=["cNMF_Cluster"],
           use_raw=True, 
           palette='Set2',frameon=False,
           save="cNMF_Cluster_umap2.pdf")
#sc.pl.umap(adata, color='cNMF_Cluster', add_outline=True, #legend_loc='on data',
#               legend_fontsize=12, legend_fontoutline=2,frameon=False,
#               title='cNMF_Cluster', palette='Set1')


# In[40]:


#sc.pl.embedding_density(adata,color=usage_norm.columns,color_map='gnuplot2',add_outline=True)


# In[40]:


#pip install omicverse 


# In[2]:


#pip install torch==1.13.1+cu117 torchvision==0.14.1+cu117 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu117


# In[4]:


import torch
print("torch版本：",torch.__version__)


# In[5]:


#pip install torch-geometric


# In[4]:


#高级版 omicverse
import scanpy as sc
import omicverse as ov
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import patheffects

from matplotlib import gridspec
import matplotlib.pyplot as plt


# In[2]:


adata=sc.read('/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T_afterumap.h5ad')


# In[143]:


sc.write('Epithelial_T_afterumap.h5ad', adata)


# In[206]:


adata


# In[10]:


#从这一步开始到下面3行其实都不用
#output_directory = '/home/wling_32/jupyter_home/NGmetadata/cNMF'
#if not os.path.exists(output_directory):
#    os.mkdir(output_directory)
#run_name = 'Epithelial_T_cNMF'


# In[11]:


#cnmf_obj = cNMF(output_dir=output_directory, name=run_name)


# In[12]:


#selected_K = 7
#density_threshold = 0.15


# In[13]:


#usage_norm, gep_scores, gep_tpm, topgenes = cnmf_obj.load_results(K=selected_K, density_threshold=density_threshold)
#usage_norm.columns = ['cNMF_Cluster%d' % i for i in usage_norm.columns]


# In[168]:


result_dict = cnmf_obj.load_results(K=selected_K, density_threshold=density_threshold)


# In[207]:


result_dict


# In[ ]:


#adata=sc.read('/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T_afterumap.h5ad')


# In[39]:


result_dict[3].columns  ##这里只能用数字表示哪个数组


# In[47]:


#cnmf_obj.get_results(adata,result_dict)


# In[56]:


#ov.pl.embedding(adata, basis='X_umap',color=result_dict[0].columns,
#           use_raw=False, ncols=3, vmin=0, vmax=1,frameon='small')


# In[214]:


fig, ax = plt.subplots(figsize=(4,4))
ov.pl.embedding(
    adata,
    basis="X_umap",
    color=['cNMF_Cluster7'],
    frameon='small',
    title="cNMF_Cluster7",
    #legend_loc='on data',
    legend_fontsize=14,
    legend_fontoutline=2,
    #size=10,
    ax=ax,
    #legend_loc=True, 
    add_outline=False, 
    #add_outline=True,
    outline_color='black',
    outline_width=1,
    show=True,palette='Set1',
    save="cNMF_Cluster_umapcluster7.pdf"
)


# In[33]:


get_ipython().run_line_magic('matplotlib', 'inline')

import matplotlib.pyplot as plt

# 绘制热图
plt.figure()
# 这里添加你的热图绘制代码
cnmf_obj.consensus(k=selected_K, density_threshold=density_threshold, show_clustering=True, close_clustergram_fig=False)
# 保存为PDF格式
plt.savefig('heatmap.pdf', format='pdf')  ###——————这种绘制的不太行 无法在AI中编辑


# In[215]:


fig, ax = plt.subplots(figsize=(4,4))
ov.pl.embedding(
    adata,
    basis="X_umap",
    color=['cNMF_Cluster'],
    frameon='small',
    title="cNMF_Cluster",
    #legend_loc='on data',
    legend_fontsize=14,
    legend_fontoutline=2,
    #size=10,
    ax=ax,
    #legend_loc=True, 
    add_outline=False, 
    #add_outline=True,
    outline_color='black',
    outline_width=1,
    show=True,palette=["#F3AE63","#D55640","#DBBD99","#65BA8E","#E2BECB","#73558B","#7BBDDB"],
    save="cNMF_Cluster_umap3.pdf"
)


# In[76]:


##绘制cNMF_cluster的markergenge--气泡图
plot_genes=[]
for i in result_dict[3].columns:
    plot_genes+=result_dict[3][i][:3].values.reshape(-1).tolist()


# In[119]:


from matplotlib.colors import LinearSegmentedColormap


# In[217]:


colors = [(0, '#BDB6AE'),  # 灰色
          (0.5, '#FFFFFF'),  # 白色
          (1, '#A10923')]  # 红色  #B33A52 #AD1F37 #A10923 #DB092C #D62E0D


# In[218]:


cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)


# In[123]:


new_colors=sns.color_palette("PuOr_r", 10)[1:8]
sns.palplot(sns.color_palette("PuOr_r", 10)[1:9])


# In[220]:


sc.pl.dotplot(adata,plot_genes,
              "cNMF_Cluster", 
              dendrogram=False,
              standard_scale='var',
              cmap=cmap,
              save="cNMF_cluster的markergenge3_1.pdf")


# In[51]:


ax = sc.pl.tracksplot(adata, plot_genes, groupby='cNMF_Cluster', dendrogram=False)--不用


# In[67]:


ax = sc.pl.heatmap(adata, plot_genes, groupby='cNMF_Cluster', 
                   #layer='scaled',
                   vmin=-2, vmax=2, cmap='RdBu_r', dendrogram=False, swap_axes=True, figsize=(11,4))--不用


# ## 绘制cNMF FR图

# In[23]:


import os
os.chdir('/home/wling_32/jupyter_home/NGmetadata/Neddylation_AUcell')


# In[ ]:


sc.tl.draw_graph(adata) #很长时间 已保存于拟时序分析PAGA过程中


# In[24]:


ov.pl.embedding(adata,
                basis='X_draw_graph_fr',
               color='Neddylation_aucell',
               frameon='small',
               show=True,
               save="T_Neddylation_AUCell_FR.pdf",)


# In[26]:


import os
os.chdir('/home/wling_32/jupyter_home/NGmetadata/cNMF')


# In[28]:


import matplotlib.pyplot as plt


# In[25]:


ov.pl.embedding(adata,
                basis='X_draw_graph_fr',
               color='cNMF_Cluster',
               frameon='small',
               show=True,
               save="T_cNMF_Cluster_FR.pdf",)


# In[35]:


fig, ax = plt.subplots(figsize=(4,4))
ov.pl.embedding(
    adata,
    basis="X_draw_graph_fr",
    color=['cNMF_Cluster7'],
    frameon='small',
    title="cNMF_Cluster7",
    #legend_loc='on data',
    legend_fontsize=14,
    legend_fontoutline=2,
    #size=10,
    ax=ax,
    #legend_loc=True, 
    add_outline=False, 
    #add_outline=True,
    outline_color='black',
    outline_width=1,
    show=True,palette='Set1',
    save="cNMF_Cluster_FRcluster7.pdf"
)


# # 添加Neddylation评分

# In[ ]:


##我使用R的allcell数据，构建了AUCell评分，然后添加到python中的allcell的obs中使用，
#因为R的allcell数据的基因是初始的3w多，python里的只有高辨基因2000，不够AUCell识别neddylation通路基因


# In[6]:


##文件保存路径/位置
import os 
print(os.path.abspath('.'))


# In[230]:


##更改文件保存路径/位置
import os
os.chdir('/home/wling_32/jupyter_home/NGmetadata/Neddylation_AUcell')


# ## 绘制allcell的Umap图

# In[16]:


import omicverse as ov  ##用的omicverse，使用omicverse手册
import scanpy as sc
import scvelo as scv
import pandas as pd


# In[ ]:


allcell=sc.read('/home/wling_32/NGmetadata/Rawdata/allcell.h5ad')


# In[313]:


allcell


# In[315]:


from matplotlib import patheffects
fig, ax = plt.subplots(figsize=(4,4))
ov.pl.embedding(allcell,
    basis="X_umap",
    color=['cell.type'],title='',
                   show=False, legend_loc=None, add_outline=False, 
                   frameon='small',legend_fontoutline=2,ax=ax
                 )

ov.pl.embedding_adjust(
    allcell,
    basis="X_umap",
    groupby='cell.type',
    ax=ax,
    adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
    text_kwargs=dict(fontsize=12 ,weight='bold',
                     path_effects=[patheffects.withStroke(linewidth=2, foreground='w')] ),
)
plt.savefig("allcell_umap.pdf")


# ## allcell打分

# In[ ]:


#读取metadata数据
metadata_allcell = pd.read_csv('metadata_allcell.csv',index_col=0)
metadata_allcell


# In[279]:


print(allcell.obs.index)  #看行名是否相同


# In[281]:


allcell.obs["Neddylation_aucell"]=metadata_allcell["AUC"]


# In[283]:


ov.pl.embedding(allcell,
                basis='X_umap',
               color='Neddylation_aucell',
               frameon='small',
               show=True,
               save="Neddylation_AUCell.pdf",)


# In[312]:


sc.pl.embedding(allcell,
                basis='umap',
          color=["Neddylation_aucell"])


# In[310]:


fig, ax = plt.subplots(figsize=(8,3))
ov.pl.bardotplot(allcell,groupby='cell.type',color='Neddylation_aucell',figsize=(8.3),
           ax=ax,
          ylabel='Expression',
           bar_kwargs={'alpha':0.5,'linewidth':2,'width':0.6,'capsize':4},
           scatter_kwargs={'alpha':0.01,'s':1,'marker':'o'})

#ov.pl.add_palue(ax,line_x1=3,line_x2=4,line_y=0.1,
#          text_y=0.02,
#          text='$p={}$'.format(round(0.001,3)),
#          fontsize=11,fontcolor='#000000',
#             horizontalalignment='center',)
plt.savefig('Neddylation_aucell_barplot.png')


# In[303]:


import pandas as pd
import seaborn as sns
#sns.set_style('white')

ov.pl.single_group_boxplot(allcell,groupby='cell.type',
             color='Neddylation_aucell',
             type_color_dict=dict(zip(pd.Categorical(allcell.obs['cell.type']).categories, ['#1f77b4',
  '#ff7f0e',
  '#279e68',
  '#d62728',
  '#aa40fc',
  '#8c564b',
  '#e377c2',
  '#b5bd61',
  '#17becf',
  '#aec7e8',
  '#ffbb78'])),
             x_ticks_plot=True,
             figsize=(5,2),
             kruskal_test=True,
             ylabel='Neddylation_aucell',
             legend_plot=False,
             bbox_to_anchor=(1,1),
             title='Expression',
             scatter_kwargs={'alpha':0.8,'s':5,'marker':'o'},
             point_number=15,
             sort=False,
             )
plt.grid(False)
plt.xticks(rotation=90,fontsize=12)
plt.savefig("Allcell_Neddylation_boxplot.pdf")


# In[291]:


#allcell.uns


# ## 恶性上皮打分

# In[5]:


adata=sc.read('/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T_afterumap.h5ad')


# In[18]:


UMAP1=adata.obsm["X_umap"][:,0]
UMAP2=adata.obsm["X_umap"][:,1]
X_umap=pd.DataFrame({'UMAP1': UMAP1, 'UMAP2': UMAP2})
X_umap


# In[ ]:





# In[317]:


metadata_allcell = pd.read_csv('metadata_allcell.csv',index_col=0)
metadata_allcell


# In[321]:


filtered_metadata_allcell = metadata_allcell[metadata_allcell.index.isin(adata.obs.index)]
filtered_metadata_allcell


# In[322]:


print(adata.obs.index)  


# In[323]:


adata.obs["Neddylation_aucell"]=filtered_metadata_allcell["AUC"]


# In[324]:


ov.pl.embedding(adata,
                basis='X_umap',
               color='Neddylation_aucell',
               frameon='small',
               show=True,
               save="T_Neddylation_AUCell.pdf",)


# In[325]:


import pandas as pd
import seaborn as sns
#sns.set_style('white')

ov.pl.single_group_boxplot(adata,groupby='cNMF_Cluster',
             color='Neddylation_aucell',
             type_color_dict=dict(zip(pd.Categorical(adata.obs['cNMF_Cluster']).categories, ["#F3AE63","#D55640","#DBBD99","#65BA8E","#E2BECB","#73558B","#7BBDDB"])),
             x_ticks_plot=True,
             figsize=(5,2),
             kruskal_test=True,
             ylabel='Neddylation_aucell',
             legend_plot=False,
             bbox_to_anchor=(1,1),
             title='Expression',
             scatter_kwargs={'alpha':0.8,'s':5,'marker':'o'},
             point_number=15,
             sort=False,
             )
plt.grid(False)
plt.xticks(rotation=90,fontsize=12)
plt.savefig("Tumor_Neddylation_boxplot.pdf")


# In[328]:


sc.write('/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T_afterumap.h5ad', adata)


# ## 正常和恶性上皮的打分

# In[331]:


adata.obs["Sample_type"].value_counts()


# In[334]:


Epithelial=sc.read('/home/wling_32/jupyter_home/NGmetadata/Epithelial.h5ad')


# In[336]:


Epithelial.obs["Sample_type"].value_counts()


# In[337]:


filtered_metadata_allcell = metadata_allcell[metadata_allcell.index.isin(Epithelial.obs.index)]
filtered_metadata_allcell


# In[338]:


print(Epithelial.obs.index)  


# In[341]:


print(Epithelial.obs.index == filtered_metadata_allcell.index)


# In[342]:


# 比较两个DataFrame中的列
result = Epithelial.obs.index == filtered_metadata_allcell.index

# 检查是否有任何False值
if result.any() == False:
    print("存在不相等的元素")
else:
    print("所有元素都相等")


# In[343]:


Epithelial.obs["Neddylation_aucell"]=filtered_metadata_allcell["AUC"]


# In[344]:


Epithelial


# In[372]:


import pandas as pd
Epithelial.obs.to_excel("/home/wling_32/jupyter_home/NGmetadata/epithelial_obs.xlsx",engine='openpyxl')


# In[ ]:


##这部分去R绘制箱线图或者云雨图了


# # 拟时序分析（PAGA、Cytotrace2）

# In[81]:


##文件保存路径/位置
import os 
print(os.path.abspath('.'))


# In[45]:


##更改文件保存路径/位置
import os
os.chdir('/home/wling_32/jupyter_home/NGmetadata/Pseudotime')


# ## Cytotrace2

# ### R版本（跑完导入到python）

# In[59]:


adata=sc.read('/home/wling_32/jupyter_home/NGmetadata/Pseudotime/Epithelial_T_afterumap.h5ad')


# In[47]:


cytotrace_metadata = pd.read_csv('cytotrace_metadata.csv',index_col=0)
cytotrace_metadata


# In[49]:


print(adata.obs.index==cytotrace_metadata.index)


# In[52]:


new_columns_df = cytotrace_metadata[["CytoTRACE2_Score","CytoTRACE2_Potency","CytoTRACE2_Relative","preKNN_CytoTRACE2_Score","preKNN_CytoTRACE2_Potency"]]
new_columns_df


# In[60]:


adata.obs = adata.obs.join(new_columns_df)


# In[61]:


adata.obs["CytoTRACE2_Potency"].value_counts()


# In[62]:


adata


# In[70]:


sc.write('Epithelial_T_afterumap.h5ad', adata)#此部分加了X_draw_graph_fr 和 CytoTRACE2数据


# In[72]:


fig, ax = plt.subplots(figsize=(4,4))
ov.pl.embedding(
    adata,
    basis="X_draw_graph_fr",
    color=['CytoTRACE2_Relative'],
    frameon='small',
    title="CytoTRACE2_Relative",
    #legend_loc='on data',
    legend_fontsize=14,
    legend_fontoutline=2,
    #size=10,
    ax=ax,
    #legend_loc=True, 
    add_outline=False, 
    #add_outline=True,
    outline_color='black',
    outline_width=1,
    show=True,palette='Set1',
    save="CytoTRACE2_RelativeScore_FR.pdf"
)


# In[68]:


fig, ax = plt.subplots(figsize=(4,4))
ov.pl.embedding(
    adata,
    basis="X_draw_graph_fr",
    color=['CytoTRACE2_Potency'],
    frameon='small',
    title="CytoTRACE2_Potency",
    #legend_loc='on data',
    legend_fontsize=14,
    legend_fontoutline=2,
    #size=10,
    ax=ax,
    #legend_loc=True, 
    add_outline=False, 
    #add_outline=True,
    outline_color='black',
    outline_width=1,
    show=True,palette=["#CFC7C7","#CC1D1D","#CC7878"],
    save="CytoTRACE2_Potency_FR.pdf"
)


# ## PAGA

# ### 运行PAGA及初步出图

# In[4]:


adata=sc.read('/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T_afterumap.h5ad')


# In[22]:


adata


# In[42]:


import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns 


# In[6]:


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80)


# In[7]:


results_file = '/home/wling_32/jupyter_home/NGmetadata/cNMF/Pseudotime' 


# In[8]:


#sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='cNMF_Cluster', legend_loc='on data')#时间长


# In[9]:


sc.write('Epithelial_T_afterumap.h5ad', adata) ##这里保存的版本多了X_draw_graph的内容


# In[74]:


sc.tl.paga(adata, groups='cNMF_Cluster')


# In[76]:


adata


# In[75]:


sc.pl.paga(adata, color=['cNMF_Cluster', 'Neddylation_aucell'])


# In[41]:


#sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['cNMF_Cluster','UBE2M'], 
                 #legend_loc='on data'
                )


# In[77]:


#设置根节点 然后运行PAGA的拟时序
adata.uns['iroot'] = np.flatnonzero(adata.obs['cNMF_Cluster']  == 'cNMF_Cluster2')[0]
sc.tl.dpt(adata)


# In[79]:


sc.pl.draw_graph(adata, color=['cNMF_Cluster', 'CytoTRACE2_Score',"Neddylation_aucell"], 
                 #legend_loc='on data'
                ) #这图不行

