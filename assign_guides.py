import crispat
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import os


platforms = ['BD']
#samples = ['vivo2','vivo3','vivo4','vivo5']
samples = ['vivo1']
methods = ['TE']

for pf in platforms:
  for samp in samples:
    for method in methods:
      guide_nam = f"{pf}_{samp}"
      print(guide_nam)
      input_file = f"../Seurat/vivo1/gRNA_counts_{samp}.csv"
      outdir = guide_nam
      # create h5ad
      crispat.create_anndata_from_csv(input_file, "./")
      # rename file
      ad_file = f"gRNA_counts_{guide_nam}.h5ad"
      os.rename('gRNA_counts.h5ad', ad_file)
      # check Anndata
      adata = ad.read_h5ad(ad_file)
      # adata.obs['batch'] = adata.obs.index.to_series().str.split('A').str[0].astype(int)
      # adata.write_h5ad(ad_file)
      # Assign gRNAs
      crispat.ga_ratio(ad_file, [0.3], outdir + "/ratios/", UMI_threshold=3)  # ratio thresholds
      crispat.ga_gauss(ad_file, outdir + "/gauss/", inference='vi', UMI_threshold=3)  # Gaussian mixture model (Cell Ranger approach)
      #crispat.ga_poisson_gauss(ad_file, outdir + "/poisson_gauss/", UMI_threshold=3)  # Poisson-Gaussian mm (Replogle et al. Cell 2022)
      crispat.ga_2beta(ad_file, outdir + "/2-BetaMM/", UMI_threshold=3)  # 2-Beta mixture model
      # compare assignments
      #ga_dict={'BetaMM':['2-BetaMM',None],'Gauss': ['gauss', None], 'Poisson-Gauss': ['poisson_gauss', None],'Ratio_30': ['ratios', 0.3]}
      ga_dict={'BetaMM':['2-BetaMM',None],'Gauss': ['gauss', None], 'Ratio_30': ['ratios', 0.3]}
      perts = crispat.load_assignments(ga_dict, guide_nam + '/')
      crispat.plot_n_assigned_cells(perts)
      plt.savefig(f'plot_n_assigned_cells_{guide_nam}.png', dpi=300, bbox_inches='tight')
      single_perts = perts.groupby(['method', 'cell']).filter(lambda x: len(x) == 1)
      single_perts = single_perts.astype({'cell': 'str'})
      crispat.plot_intersection_heatmap(single_perts)
      plt.savefig(f'heatmap_intersection_{guide_nam}.png', dpi=300, bbox_inches='tight')


nam="vivo1Only"
for method in methods:
    print(nam)
    ga_dict={'BetaMM':['2-BetaMM',None],'Gauss': ['gauss', None],'Ratio_30': ['ratios', 0.3]}
    perts1 = crispat.load_assignments(ga_dict, 'BD_vivo1/')
    #perts2 = crispat.load_assignments(ga_dict, 'BD_vivo2/')
    #perts3 = crispat.load_assignments(ga_dict, 'BD_vivo3/')
    #perts4 = crispat.load_assignments(ga_dict, 'BD_vivo4/')
    #perts5 = crispat.load_assignments(ga_dict, 'BD_vivo5/')
    perts1.method = "vivo1_"+perts1.method
    #perts2.method = "vivo2_"+perts2.method
    #perts3.method = "vivo3_"+perts3.method
    #perts4.method = "vivo4_"+perts4.method
    #perts5.method = "vivo5_"+perts5.method
    perts = perts1
    #perts = pd.concat([perts2, perts3, perts4, perts5], ignore_index=True)
    perts = perts[perts['method'].str.endswith("Gauss")]
    perts.to_csv(f"combined_assignments_{nam}_Gauss_only.tsv", sep='\t', index=False)
    crispat.plot_n_assigned_cells(perts)
    plt.savefig(f'in_vivo_n_assigned_cells_{nam}_Gauss_only.png', dpi=300, bbox_inches='tight')


