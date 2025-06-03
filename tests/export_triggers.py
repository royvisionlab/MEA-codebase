import numpy as np
from typing import Tuple, Dict, Union, List, Set, Any
import pickle

def save_trigger_pickle_file(pickle_file_path: str,
                            array_id: int,
                            n_samples: int,
                            trigger_times: np.ndarray,
                            neuron_spike_time_offset: int = 0) -> None:
    trigger_data_dict = {
        'trigger_times': trigger_times,
        'array_id': array_id,
        'n_samples': n_samples,
        'neuron_spike_time_offset': neuron_spike_time_offset
    }

    with open(pickle_file_path, 'wb') as pfile:
        pickle.dump(trigger_data_dict, pfile)