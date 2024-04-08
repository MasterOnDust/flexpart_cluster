import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import cartopy
import numpy as np
from matplotlib import collections
try:
    from sklearn.metrics import DistanceMetric
except ImportError:
    print('Adaptive_Kmeans does not work without sklearn')
import shapely

def plot_center_trajectory(lons, lats, height, 
                            add_colorbar=False, ax=None,p0=None,
                             add_labels=False, vmin=100, vmax=7000, 
                             cbar_label='Meters above sea level [m]', extent=[70,30,60,120],
                          receptor_marker_color='black', join_trajec_to_receptor=False,
                          cmap=None, norm=None):
    """
    DESCRIPITON:
    ===========
        Plot the cluster centriod trajectories:

    ARGS:
    =====
        lons : longtidues of center trajectoru
        lats : latitudes of center trajectory
        height: height of center trajectory

    KWARGS:
    =======
        ax : cartopy geoAxes
        add_colorbar : boolean, whether to include colorbar or not
        p0: starting loctation of trajectory [lon0,lat0]
        vmin,vmax : determine limits of colorbar
        cbar_label : what label the colorbar should have 

    """


    if ax==None:
        ax = plt.gca(projection=ccrs.PlateCarree())
        ax.set_extent(extent)
        ax.add_feature(cartopy.feature.BORDERS, linestyle=':')


    
    if isinstance(p0,list) or isinstance(p0, np.ndarray):

        ax.scatter(p0[0],  p0[1], marker='*', s=75, zorder=1200, facecolor=receptor_marker_color, linewidth=.8, edgecolor='black')
    if join_trajec_to_receptor:
        lons_array = np.zeros(len(lons)+1)
        lats_array = np.zeros(len(lats)+1)
        height_array = np.zeros(len(height)+1)
        lons_array[0] = p0[0]
        lats_array[0] = p0[1]
        lons_array[1:] = lons[:]
        lats_array[1:] = lats[:]
        height_array[0] = height[0]
        height_array[1:] = height[:]
        
    else:
        lons_array = lons
        lats_array = lats
        height_array = height
        
    # print(lons_array[0],lons[0])
    p = np.array((lons_array,lats_array)).T.reshape(-1,1,2)
    segments = np.concatenate([p[:-1], p[1:]], axis=1)
    lc = LineCollection(segments)
    lc.set_array(height_array)
    lc.set_linewidth(2)
    if norm:
        lc.set_norm(norm)
    else:
        lc.set_clim(vmin=vmin, vmax=vmax)
    if cmap:
        lc.set_cmap(cmap)
    line = ax.add_collection(lc)

    if add_colorbar:
        fig=plt.gcf()
        fig.colorbar(line,pad=0.02, label=cbar_label)
    return line

def plot_cluster_centriods(lons, lats, height, 
                            fclust=None,add_colorbar=False, ax=None,p0=None,
                             add_labels=False, vmin=100, vmax=7000, 
                             cbar_label='Meters above sea level [m]', extent=[70,30,60,120],
                           receptor_marker_color='black'):
    """
    DESCRIPITON:
    ===========
        Plot the cluster centriod trajectories:

    ARGS:
    =====
        lons : longtidues of centriod trajectories
        lats : latitudes of centriod trajectories
        height: height of centriod trajectories,


    KWARGS:
    =======
        flcust : percentage of each trajectory belonging to each cluster
        ax : cartopy geoAxes
        add_colorbar : boolean, whether to include colorbar or not
        p0: starting loctation of trajectory [lon0,lat0]
        vmin,vmax : determine limits of colorbar
        cbar_label : what label the colorbar should have 


    """
    
    if ax==None:
        ax = plt.gca(projection=ccrs.PlateCarree())
        ax.set_extent(extent)
        ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
    
    if isinstance(p0,list) or isinstance(p0, np.ndarray):

        ax.scatter(p0[0],  p0[1], marker='*', s=50, zorder=1200, color=receptor_marker_color)

    for cluster in range(lons.shape[0]):
        p = np.array((lons[cluster,:],lats[cluster,:])).T.reshape(-1,1,2)
        segments = np.concatenate([p[:-1], p[1:]], axis=1)
        lc = LineCollection(segments)
        lc.set_array(height[cluster,:])
        lc.set_linewidth(2)
        lc.set_clim(vmin=vmin, vmax=vmax)
        line = ax.add_collection(lc)
        if isinstance(fclust, np.ndarray):
            ax.text(lons[cluster,14] ,lats[cluster,14], '{:.1f} %'.format(fclust[cluster]),
                verticalalignment='center', horizontalalignment='center')
    if add_colorbar:
        fig=plt.gcf()
        fig.colorbar(line,pad=0.02, label=cbar_label)
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





def get_distance_circles(lons,lats, xcenter,ycenter, method='std', weight=None):
    """
    DESCRIPTION
    ===========
        Calcluate the distance circle, for each point along the centriod trajectory. 
        The radius of the distance circle could either be the standard deviation 
        of the distance from the center trajectory, or mean distance from center 
        trajectory. If the weight argument is provided weighted mean and standard 
        deviation is used.
        
    ARGUMENTS
    =========
        lons : numpy.array (n_trajectory_points, n_trajectories) longitudes of all trajectories
        lats : numpy.array (n_trajectory_points, n_trajectories) latitudes of all trajectories
        xcenter : numpy.array (n_trajectory_points) longitudes of center trajectory
        ycenter : numpy.array (n_trajectory_points) lattitudes of center trajectory
    
    KWARGS
    ======
        method : std/mean, determine what kind of distance circle is calculated.
        weight : calculated weighted standard deviation / mean
        
    OUTPUT
    ======
    
        list of shapely circles 
    
    """
    
    ct_dist = dist_from_center(lons,lats,xcenter,ycenter)*1000

    m_dist = np.average(ct_dist,weights=weight, axis=1)

    if method == 'std':
        variance = np.zeros((lons.shape[0]))
        for i in range(lons.shape[0]):
            variance[i] = np.average((ct_dist[i,:]-m_dist[i])**2, weights=weight)
        ct_dist = np.sqrt(variance)
    elif method == 'mean':
        ct_dist = m_dist
    else:
        raise(ValueError('{} is not a valid method, use either mean or std'.format(method)))

    geoms = []
    for i in range(lons.shape[0]):

        cp_0 = cartopy.geodesic.Geodesic().circle(xcenter[i],ycenter[i], ct_dist[i], 200)
        geoms.append(shapely.geometry.Polygon(cp_0))
    
    return geoms

def dist_from_center(lons,lats,xcenter,ycenter):
    """
    DESCRIPTION
    ===========
        Calculates the distance from the center
    
    
    
    """
    dist_metric = DistanceMetric.get_metric('haversine')
    dists = np.zeros_like(lons)
    center = np.stack((ycenter,xcenter),axis=-1)
    center = np.deg2rad(center)
    for i in range(lons.shape[1]):
        temp_t = np.stack((lats[:,i],lons[:,i]),axis=-1)
#         dists[:,i] = paired_distances(center,temp_t)
        temp_t = np.deg2rad(temp_t)
        dists[:,i] = np.diag(dist_metric.pairwise(center, temp_t))*6731
    
    return dists
