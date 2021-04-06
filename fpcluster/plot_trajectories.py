import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import numpy as np

def plot_trajectories(ds, figsize=(16,6), add_colorbar=False, ax=None, add_labels=False):
    if ax==None:
        ax = plt.gca(projection=ccrs.PlateCarree())
    
    # fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection':ccrs.PlateCarree()})
    # ax.coastlines()
    for cluster in ds.nclust:
        xclust = ds.xclust.sel(nclust=cluster.values)
        yclust = ds.yclust.sel(nclust=cluster.values)
        zclust = ds.zclust.sel(nclust=cluster.values)
        p =np.array([ds.xclust.sel(nclust=cluster.values) ,ds.yclust.sel(nclust=cluster.values)]).T.reshape(-1,1,2)


        segments = np.concatenate([p[:-1], p[1:]], axis=1)
        lc = LineCollection(segments)
        lc.set_array(zclust.values)
        lc.set_linewidth(2)

        line = ax.add_collection(lc)
        if add_labels:
            ax.text(xclust[14].values ,yclust[14].values, '{:.1f} %'.format(ds.fclust.sel(nclust=cluster.values).values),
                verticalalignment='center', horizontalalignment='center')
    if add_colorbar:
        fig.colorbar(line,pad=0.02, label='Meters above sea level [m]')
    return line
