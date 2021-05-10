import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import cartopy
import numpy as np

def plot_cluster_centriod(ds, figsize=(16,6), add_colorbar=False, ax=None, add_labels=False):
    """
    DESCRIPTION:
    ===========
        Takes in a dataset containing the cluster and trajectories and plots them. 

    """


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

def plot_cluster(lons,lats,height,labels, 
                    ax=None, alpha=.4,
                    add_colorbar=False, ncols=3, 
                    extent=[60,120,20,60],figsize=(15,15), vmin=0,vmax=7000):
    n_trajec = len(labels)
    clusts = np.unique(labels)
    n_clust = len(clusts)
    if isinstance(ax, np.ndarray):
        ax = ax.ravel()
    else:
        nrows = int(np.ceil(n_clust/ncols)) 
        # nrows += n_clust % ncols
        fig,ax = plt.subplots(nrows=nrows,ncols=ncols,
                    subplot_kw={'projection':ccrs.PlateCarree()}, figsize=figsize)
        ax = ax.ravel()
        for ax_i in ax: 
            ax_i.set_extent(extent)
            ax_i.add_feature(cartopy.feature.BORDERS, linestyle=':')
            ax_i.scatter(lons[0,0],  lats[0,0], marker='*', s=50, zorder=1200, color='black')
    for x,y,z,n in zip(lons.T, lats.T, height.T,labels):
        p = np.array([x,y]).T.reshape(-1,1,2)
    
        segments = np.concatenate([p[:-1], p[1:]], axis=1)
        lc = LineCollection(segments,alpha=.4)
        lc.set_array(z)
        lc.set_clim(vmin,vmax)
        line = ax[n].add_collection(lc)

    for i in range(n_clust):
        ax[i].text(0.9,0.1, np.count_nonzero(labels==i),transform=ax[i].transAxes,  horizontalalignment='center')
    
    if add_colorbar:
        fig = plt.gcf()
        fig.colorbar(line, ax=ax.tolist())

