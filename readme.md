# FLEXPART clustering

Python wrapper of the fortan trajectory clustering routine in FLEXPART (Stohl et al. 2005).

## Usage:

```
from fpcluster.clustering import clustering

print(clustering.__doc__)

>> Wrapper for ``clustering``.

Parameters
----------
xl : input rank-1 array('f') with bounds (n)
yl : input rank-1 array('f') with bounds (n)
zl : input rank-1 array('f') with bounds (n)

Other Parameters
----------------
n : input int, optional
    Default: len(xl)
ncluster : input int, optional
    Default: 5

Returns
-------
xclust : rank-1 array('f') with bounds (ncluster)
yclust : rank-1 array('f') with bounds (ncluster)
zclust : rank-1 array('f') with bounds (ncluster)
fclust : rank-1 array('f') with bounds (ncluster)
rms : float
rmsclust : rank-1 array('f') with bounds (ncluster)
zrms : float

```

## Installation 

` git clone https://github.com/Ovewh/flexpart_cluster.git`

Then install the package using pip:

`pip install flexpart_cluster `