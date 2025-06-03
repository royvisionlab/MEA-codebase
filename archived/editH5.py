
import h5py
import numpy as np


# filepath = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220909C_test.h5'
filepath = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220909C.h5'


blocks = [
    '/experiment-a3a2b3b4-d67e-4e76-bb82-cd19944840a1/epochGroups/epochGroup-1f4ee741-eb0e-48f6-a843-aceae4fe2887/epochBlocks/manookinlab.protocols.FlashedSpatialNoise-beee9c4a-18a0-4579-8791-80638867b463',
    '/experiment-a3a2b3b4-d67e-4e76-bb82-cd19944840a1/epochGroups/epochGroup-e91913cf-a96f-4a08-b426-348a623476a2/epochBlocks/manookinlab.protocols.NaturalImageFlashPlusNoise-d73eaee7-ec11-43c4-abd7-cb6b09d042f7',
    '/experiment-a3a2b3b4-d67e-4e76-bb82-cd19944840a1/epochGroups/epochGroup-e91913cf-a96f-4a08-b426-348a623476a2/epochBlocks/manookinlab.protocols.NaturalImageFlashPlusNoise-6e4b1d95-fad6-4c65-8943-c012212fa27a',
    '/experiment-a3a2b3b4-d67e-4e76-bb82-cd19944840a1/epochGroups/epochGroup-95bfc121-3e50-4c8f-8b21-eb92fb7ab84a/epochBlocks/manookinlab.protocols.FlashedSpatialNoise-1483c85c-ef4c-4503-92ae-7f4832b05d8b'
]
new_names = [
    b'20220909C\\data026.bin',
    b'20220909C\\data027.bin',
    b'20220909C\\data028.bin',
    b'20220909C\\data029.bin'
]
# new_names = [
#     "20220909C\\data026.bin",
#     "20220909C\\data027.bin",
#     "20220909C\\data028.bin",
#     "20220909C\\data029.bin"
# ]

# file = h5py.File (filepath)

with h5py.File(filepath, 'r+') as file:
    for count, block in enumerate(blocks):
        new_name = new_names[count]
        file[block+'/properties/'].attrs['dataFileName'] = np.string_(new_name)
        # file[block+'/properties/'].attrs.modify('dataFileName', np.char.encode(new_name, "utf-8"))
        for k in file[block+'/epochs/'].keys():
            epoch_path = block + '/epochs/' + k + '/protocolParameters'
            file[epoch_path].attrs['dataFileName'] = np.string_(new_name)

with h5py.File(filepath, 'r+') as file:
    for count, block in enumerate(blocks):
        new_name = new_names[count]
        file[block+'/properties/'].attrs['dataFileName'] = new_name
        for k in file[block+'/epochs/'].keys():
            epoch_path = block + '/epochs/' + k + '/protocolParameters'
            file[epoch_path].attrs['dataFileName'] = new_name

with h5py.File(filepath) as file:
    for count, block in enumerate(blocks):
        new_name = new_names[count]
        print( file[block+'/properties/'].attrs['dataFileName'] )
        print(new_name)
        for k in file[block+'/epochs/'].keys():
            epoch_path = block + '/epochs/' + k + '/protocolParameters'
            # file[epoch_path].attrs['dataFileName'] = new_name
            # print(file[epoch_path].attrs['dataFileName'])

# data000 - > data026
new_name = b'20220909C\\data026.bin'
block = '/experiment-a3a2b3b4-d67e-4e76-bb82-cd19944840a1/epochGroups/epochGroup-1f4ee741-eb0e-48f6-a843-aceae4fe2887/epochBlocks/manookinlab.protocols.FlashedSpatialNoise-beee9c4a-18a0-4579-8791-80638867b463'

file[block+'/properties/'].attrs['dataFileName'] = new_name
# file[block+'/properties/'].attrs.modify('dataFileName',b'20220909C\\data026.bin')

for k in file[block+'/epochs/'].keys():
    epoch_path = block + '/epochs/' + k + '/protocolParameters'
    # file[epoch_path].attrs['dataFileName'] = new_name
    print(file[epoch_path].attrs['dataFileName'])




# for k in file[block+'/properties/'].attrs.keys():
#     print(k)

# file.attrs.modify('dataFileName',new_name)

# file.close()
