# Analysis Code

## Installation
To install the analysis modules in your anaconda environment, first navigate to this folder and make sure you're in the desired conda environment. My environment is 'analysis', so I would type the following command in the terminal.

```console
[user@klone-login01 analysis]$ conda activate analysis
```

Then install the packages by typing:

Mac:
```console
[user@klone-login01 analysis]$ pip install -e .
[user@klone-login01 analysis]$ pip install .
```

Hyak:
```console
[user@klone-login01 analysis]$ pip install --user -e .
[user@klone-login01 analysis]$ pip install --user .
```

## Computing correlations

```python
import correlations
import numpy as np
sort_path = '/run/user/1000/gvfs/smb-share:server=128.95.10.183,share=data/data/sorted/20230228C'
sorter_name='kilosort2'
file_names = ['data001','data010']
target_ids = [1,2,3]
spike_times, cell_ids = correlations.get_spike_times(sort_path,sorter_name,file_names)
ccg, lags = correlations.compute_crosscorrelations(spike_times, cell_ids, target_ids)

import hdf5storage
mdic = {'ccg':ccg, 'lags':lags, 'cell_ids':cell_ids, 'target_ids':target_ids}
filepath = '/data/data/sorted/20230228C_chunk2_xcorr.mat'
hdf5storage.savemat( filepath, mdic, format=7.3, matlab_compatible=True, compress=False )

```
