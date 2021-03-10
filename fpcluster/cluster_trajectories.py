# Author Ove Haugvaldstadstad

from .trajectory_clustering import clustering
import xarray as xr

def cluster_trajectories(lons, lats, heights, nclust=5, max_iters=100):
    """
    DESCRIPTION:
    ============
        Cluster trajectories accoring to Dorling et al. (1992) and is based on the 
        fortran clustering routine in FLEXPART (Stohl et al. 2006), but modified to 
        cluster trajectories rather than particle postions. The criterion used to 
        assign each trajectory to a cluster is to minimize the distrance to the centroid 
        cluster.  
    
    USAGE:
    ======
        ds = cluster_trajectories(lons,lats,heights, nclust=5,max_iters=100)
        
        input:
            lons : array shape(btime, n) btime, number of time step, n number
                   of trajectories.
            lats : array shape(btime, n) btime, number of time step, n number
                   of trajectories.
            heights : array shape(btime, n) btime, number of time step, n number
                   of trajectories.
            nclust : number of clusters to cluster the trajectories.
            max_iters : maxium number interations  
        
        return:
            ds : xarray Dataset. 
    
    """
    
    x_clust, y_clust, z_clust, fclust,rms ,rmsclust, zrms, niters = clustering(lons,
                                                                  lats,heights,ncluster=nclust, niters=max_iters)
    ds =xr.Dataset()
    ds = ds.assign_coords(nclust=[i for i in range(1,nclust+1)])
    ds['fclust'] = xr.DataArray(fclust, dims='nclust', attrs={
        'long_name':'percentage of each trajectories belonging to a cluster'})
    ds['xclust'] = xr.DataArray(x_clust, dims=['btime','nclust'], attrs={
        'long_name' : 'x coordinate of centroid cluster trajectory',
        'units' : 'longitude'
    })
    ds['yclust'] = xr.DataArray(y_clust, dims=['btime','nclust'], attrs={
        'long_name' : 'y coordinate of centroid cluster trajectory',
        'units' : 'latitude'
    })
    ds['zclust'] = xr.DataArray(z_clust, dims=['btime','nclust'], attrs={
        'long_name' : 'cluster height',
        'units' : 'm'
    })
    ds['rms'] = rms
    ds['rms'].attrs['long_name'] = 'total horizontal rms distance after clustering '
    ds['rmsclust']=xr.DataArray(rmsclust, dims='nclust', attrs={
        'long_name' : 'horizontal rms distance for each individual cluster',
        'units' : 'm'
    })
    ds['n_iterations'] = niters
    ds['n_iterations'].attrs['long_name'] = 'number interations done'
    ds['zrms'] = zrms
    ds['zrms'].attrs['long_name'] = 'total vertical rms distance after clustering' 
    ds['zrms'].attrs['units'] = 'm'
    ds.attrs['title'] = 'Centriod cluster trajectories'
    ds.attrs['source'] =  ("""Dorling, S.R., Davies, T.D. and Pierce, C.E. (1992): 
    Cluster analysis: a technique for estimating the synoptic meteorological 
    controls on air and precipitation chemistry - method and applications.   
    Atmospheric Environment 26A, 2575-2581.""")
    ds['n_iterations'] = niters
    

    return ds