# Authour Ove Haugvaldstad

import numpy as np
import pandas as pd


def read_trajectories(paths, kind='drydep'):
    """
    DESCRIPTION
    ===========
        Takes paths to flexpart trajectories.txt output files and 
        read them in. Wether the trajectories are from a wetdep or drydep/conentration 
        flexpart simulation has to be specified since, the first time step from the wet-
        deposition might not be useable since many trajectories get terminated as a 
        consequence of the drydepostion scheme. 
        
        There are few assumptions on the trajectories.txt files:
         (1) this function has been made for reading in backwards trajectories 
         (2) trajectories.txt only contain trajectories from a single location.
             since the trajectories. 
         (3) Only return the centriod postions of the partices, at each timestep for each
             release.
             
    USAGE:
    ======
        lons, lats, heights  = read_trajectories(paths, kind='drydep')
        
        input:
        
            paths, list of paths to trajectories.txt output files
            kind, (optional) valid arguments: drydep, conc, wetdep, if wetdep the first time 
                step is skipped
        
        output:
            
            lons : centriod x postion for the trajectories array shape (btime+1, n) btime is timesteps
                   and n is the numnber of trajectories. The first value is set to the release
                   starting point .
            lats : centriod y postion for the trajectories array shape (btime+1, n) btime is timesteps
                   and n is the numnber of trajectories. The first value is set to the release 
                   starting point.
            height : mean height of the centriod trajectory array shape (btime,n).
    """
    
    trajectories = []
    for path in paths:
        with open(path, 'r') as f:
            f.readline()
            f.readline()
            numpoint=int(f.readline().strip())
            df = load_trajectories(path, skip_rows=3*numpoint+3)
            locations = df[1]
            df=df[0]
        trajectories.append(df)
    lon0, lat0 =  float(locations.iloc[0,1]), float(locations.iloc[0,2])
    df = pd.concat(trajectories)
    if kind == 'drydep' or kind == 'conc':
        btime_steps= df.time.loc[df.time < 0].unique()
    else: 
        btime_steps= df.time.loc[df.time < 0].unique()[1:]

    print()
    numpoints = len(df.loc[df.time == btime_steps[5]])
    n_btime_steps = len(btime_steps)
    lons = np.zeros((n_btime_steps+1,numpoints))
    lats = np.zeros((n_btime_steps+1,numpoints))
    heights = np.zeros((n_btime_steps,numpoints))
    lons[0,:] = lon0
    lats[0,:] = lat0

    for i,t in enumerate(btime_steps):
        lons[i+1,:] = df.loc[df.time==t][['lon']].T
        lats[i+1,:] = df.loc[df.time==t][['lat']].T
        heights[i,:] = df.loc[df.time==t][['height']].T
 
    return lons, lats, heights


def load_trajectories(path, nclusters=5, skip_rows=24):
    cluster_list = []
    cluster_names = ['xcluster', 'ycluster', 'zcluster', 'fcluster',
    'rmscluster']
    for i in range(nclusters):
        for cn in cluster_names:
            cluster_list.append(cn + '(' +str(i)+ ')')

    with open(path,'r') as trajecFile:
        header = trajecFile.readline().split(' ')
        trajecFile.readline()
        nlocs = int(trajecFile.readline().strip())
        lines = [next(trajecFile) for i in range(nlocs*3)]
        locs = [line.strip().split(' ')[0] for i,line in enumerate(lines[2::3])]
        locs_dict = {i+1:line.strip() for i,line in enumerate(lines[2::3])}
            



    intial_data = [line1.strip() +  line2.strip() for line1, line2 in zip(lines[::3],lines[1::3])]

    inital_data = [line.split() for line in intial_data]

    df_head = pd.DataFrame(np.array(inital_data)[:,1:4])
    df_head = df_head.rename(columns={1:'time', 2: 'lon', 3:'lat'})

    df_head['locations'] = locs
        
    s_time = pd.to_datetime(header[0] + header[1])

    cols = ['location', 'time', 'lon', 'lat',
         'height', 'mean topography',
         'mean mixing height', 'mean tropopause height', 'mean PV index',
         'rms distance', 'rms', 'zrms distance', 'zrms',
         'fraction mixing layer', 'fraction PV<2pvu',
         'fraction in troposphere'] + cluster_list
    
    
    df = pd.read_csv(path, sep='\s+',
                    skiprows=lambda x: x <skip_rows, names=cols)
    
    for key, location in locs_dict.items():
        df.loc[df.loc[:,'location']==key, 'location'] = location
    
    df.loc[:,['location']] = df.loc[:,['location']].astype('category')
#     time_p_rel = s_time + pd.to_timedelta(sec_p_rel, unit = 's')
    df['start time'] = s_time
    return df, df_head
