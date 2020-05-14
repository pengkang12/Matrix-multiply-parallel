import h5py
import numpy as np

for i in range(500, 1000, 1000):
    img = np.random.random((i, i))

    f1 = h5py.File(str(i)+'_nocomp.h5', 'w')
    f1.create_dataset('img', data=img)
    f1.close()
   
    f2 = h5py.File(str(i)+'_complevel_9.h5', 'w')
    f2.create_dataset('img', data=img, compression='gzip', compression_opts=9)
    f2.close()

    f3 = h5py.File(str(i)+'_complevel_0.h5', 'w')
    f3.create_dataset('img', data=img, compression='gzip', compression_opts=0)
    f3.close()

