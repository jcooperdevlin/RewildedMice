import numpy as np
import fitsne
import sys

###ftsne

c = np.loadtxt(sys.argv[1], usecols=range(1,17), dtype=np.float)
Y = fitsne.FItSNE(c)
np.savetxt(sys.argv[2], Y, delimiter='\t')

#
### UMAP
import umap

reducer = umap.UMAP()

#umap combo
embedding = reducer.fit_transform(c)
np.savetxt(sys.argv[3], embedding, delimiter='\t')

