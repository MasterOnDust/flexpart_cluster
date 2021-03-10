import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import numpy as np

def plot_trajectories(ds, figsize=(16,6)):
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection':ccrs.PlateCarree()})
    ax.coastlines()
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
        ax.text(xclust[round(len(ds.xclust)/2)].values ,yclust[round(len(xclust)/2)].values, '{:.1f} %'.format(ds.fclust.sel(nclust=cluster.values).values),
                verticalalignment='center', horizontalalignment='center')
    # ax.legend()
    fig.colorbar(line,pad=0.02, label='Meters above sea level [m]')