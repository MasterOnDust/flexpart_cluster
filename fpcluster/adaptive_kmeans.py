
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter,AutoMinorLocator
try:
    from sklearn.metrics import DistanceMetric
    from sklearn.cluster import KMeans
except ImportError:
    print('Adaptive_Kmeans does not work without sklearn')
class Adaptive_KMeans():

    def __init__(self, lons,lats,height, distance_metric='euclidean'):
        """
        DESCRIPTION:
        ===========
            Adaptive K-means clustering algorithm which starts with kmax clusters and 
            then merges the closest cluster until there are only two clusters left, 
            after each merger the score is calculated and there is a large percentage
            change in score when moving from n to n-1 clusters the n cluster is retrained 



        """

        self.lons = lons
        self.lats = lats
        self.height = height

        # if distance_metric == 'euclidean':
        #     dist = euclidean_distances(lons.T,lats.T)
        #     # distances are mirrored, only keep upper triangle values
        #     dist = dist[np.triu_indices(dist.shape[0])]
        #     # sklearn require data to be of shape 
        #     self.dist = dist.reshape((-1,1))

        trajecs = np.stack([lons.T,lats.T], axis=-1)
        self.traj = trajecs.reshape((trajecs.shape[0],trajecs.shape[1]*2))

        # self.traj = traj

    


    def cluster(self,n_init=20,kmax=30, weight=None):
        """
        DESCRIPTION:
        ============
            Adaptive k-means clustering that start with kmax clusters. And for 
            each iteration the two closest cluster centriods are merged and new
            cluster centriod for the merged cluster is calculated. The new centriods
            are then feed into the KMeans algorithm again. 
        
        attributes: scores, relative change in score between current and previous cluster.

        """
        self.labels = np.zeros((kmax-1,self.traj.shape[0])) 
        self.scores = []
        self.KMeans = {}
        self.kmax = kmax
        self.distortions = []
        dist_metric = DistanceMetric.get_metric('haversine')
        score_change = 0
        nums = np.arange(2,kmax+1,1)[::-1]
        for num in nums:
            if num==kmax:
                km = KMeans(n_clusters=num, n_init=n_init).fit(self.traj, sample_weight=weight)
                score_change = km.inertia_
                self.distortions.append(km.inertia_)
                km_prev = km
            else:
                km = KMeans(n_clusters=num,init=cluster_centriod,n_init=1).fit(self.traj, sample_weight=weight)
                score_change = ((km.inertia_ -km_prev.inertia_)/km_prev.inertia_)*100
                km_prev= km
                self.distortions.append(km.inertia_)
                self.scores.append(score_change)
            cluster_centriod = km.cluster_centers_
            s = cluster_centriod.shape
            cluster_centriod = cluster_centriod.reshape((s[0],int(s[1]/2),2))[:,:,[1,0]]
            nclust = np.arange(0,cluster_centriod.shape[0])
            perm = combinations(nclust,2)
            dist_min=10e10
            # Find the two clusters that are closest to each other, This also need update
            for pair in perm:
                # dist = haversine_distances(np.radians(cluster_centriod[pair[0]]),
                # np.radians(cluster_centriod[pair[1]])).sum()* 6371000/1000
                dist = np.diag(dist_metric.pairwise(np.radians(cluster_centriod[pair[0]]),
                np.radians(cluster_centriod[pair[1]]))*6731).sum()
                if dist < dist_min:
                    closest_centriods = pair
                    dist_min = dist
            clust = self.traj[np.where((km.labels_==pair[0]) | (km.labels_==pair[1])),:].squeeze()
            clust = clust.reshape((clust.shape[0],int(clust.shape[1]/2),2))
            # calculate center of mass of new merged trajectory.
            merged_clust = center_of_mass_trajectory(clust[:,:,0].T,clust[:,:,1].T)
            
            merged_clust = np.array(merged_clust).reshape((cluster_centriod.shape[1],2))
            nclust = np.delete(nclust,[pair])
            new_cluster_centriod = np.zeros((num-1,int(s[1]/2),2))
            new_cluster_centriod[0:-1,:,:] =cluster_centriod[nclust,:,:]
            new_cluster_centriod[-1,:,:]=merged_clust
            cluster_centriod = new_cluster_centriod.reshape((num-1,s[1]))

            # Store data
            self.KMeans[num] = km
    def score_plot(self):
        ax = plt.gca()
        x = np.arange(1,self.kmax+1)
        ax.scatter(x[1:-1],self.scores[::-1], color='black')
        ax.set_xlabel('Number of clusters')
        ax.set_ylabel('Relative change in TSS (%)')
        ax.set_xticks(x[::3])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_locator(  AutoMinorLocator(3))

    def get_cluster_trajecs(self,k):
        """
        DESCRIPTION:
        ============
            Return the trajectories belowing to the KMeans cluster analysis 
            of k clusters. 

        """
        trajecs_cluster = {}
        c = self.KMeans[k]
        labels = c.labels_
        num_cluster = np.unique(labels)
        n_trajecs = len(labels)
        for l in num_cluster:
            clabels = np.where(labels==l)
            lons = self.lons[:,clabels].squeeze()
            height = self.height[:,clabels].squeeze()
            lats = self.lats[:,clabels].squeeze()
            trajecs_cluster[f'lats{l}'] = lats
            trajecs_cluster[f'lons{l}'] = lons
            trajecs_cluster[f'height{l}'] = height
            trajecs_cluster[f'fclust{l}'] = (np.count_nonzero(labels==l)/n_trajecs)*100
        return trajecs_cluster


    def cluster_centriods(self, k):
        """
        DESCRIPTION:
        ============
            Returns the x and y center position for each cluster aswell as 
            the average height. 

            return: xclust, yclust, zclust, flcust

        """

        c = self.KMeans[k]
        ccenters = c.cluster_centers_
        cshape = ccenters.shape
        ccenters = ccenters.reshape((cshape[0],int(cshape[1]/2),2))
        cshape = ccenters.shape
        cheight = np.zeros((cshape[0],cshape[1]))
        labels = c.labels_
        fclust = np.zeros(cshape[0])
        for i in range(k):
            n_cluster_trajec = np.count_nonzero(labels==i)
            heights=self.height[:,np.where(labels==i)].squeeze()
            cheight[i,:] = heights.mean(axis=1)
            fclust[i] = (n_cluster_trajec/len(labels))*100 

        return ccenters[:,:,0],ccenters[:,:,1] ,cheight, fclust
        
                 



def center_of_mass_trajectory(xl,yl, weights=None):
    """Calculates the center of mass trajectory for all trajectories"""

    xll = np.deg2rad(xl)
    yll = np.deg2rad(yl)
    x = np.cos(yll)*np.sin(xll)
    y = -1*np.cos(yll)*np.cos(xll)
    z = np.sin(yll)
    xav = np.average(x,axis=1, weights=weights)
    yav = np.average(y,axis=1, weights=weights)
    zav = np.average(z,axis=1, weights=weights)
    
    xcenter = np.arctan2(xav,-1.*yav)
    ycenter = np.arctan2(zav,np.sqrt(xav*xav+yav*yav))

    xcenter = np.rad2deg(xcenter)
    ycenter = np.rad2deg(ycenter)

    return xcenter,ycenter
