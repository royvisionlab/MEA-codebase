import hdf5storage
import numpy as np
import struct
import os

param_keys =  ['ID', 'classID', 'x0', 'y0', 'SigmaX', 'SigmaY', 'Theta', 'gAmp', 'contourX', 'contourY', 'contourArea', 'simpleContourX', 'simpleContourY', 'simpleContourArea', 'RedTimeCourse', 'GreenTimeCourse', 'BlueTimeCourse', 't1', 't2', 't3', 'a1', 'a2', 'a3', 'n1', 'n2', 'n3', 'tOffset', 'dot', 'dot2', 'srm', 'rl', 'amp1', 'amp2', 'amp3', 'blueness', 'RedVTimeCourse', 'GreenVTimeCourse', 'BlueVTimeCourse', 't1V', 't2V', 't3V', 'a1V', 'a2V', 'a3V', 'n1V', 'n2V', 'n3V', 'tVOffset', 'dotV', 'dot2V', 'srmV', 'rlV', 'amp1V', 'amp2V', 'amp3V', 'bluenessV', 'Auto', 'acfBinning', 'nSpikes', 'acfMean', 'acfRMS']
param_values = ['Double', 'String', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'DoubleArray', 'DoubleArray', 'Double', 'DoubleArray', 'DoubleArray', 'Double', 'DoubleArray', 'DoubleArray', 'DoubleArray', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'DoubleArray', 'DoubleArray', 'DoubleArray', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'DoubleArray', 'Double', 'Double', 'Double', 'Double']
# param_keys = ['ID', 'classID', 'x0', 'y0', 'SigmaX', 'SigmaY', 'Theta', 'gAmp', 'contourX', 'contourY', 'contourArea', 'simpleContourX', 'simpleContourY', 'simpleContourArea', 'RedTimeCourse', 'GreenTimeCourse', 'BlueTimeCourse', 't1', 't2', 't3', 'a1', 'a2', 'a3', 'n1', 'n2', 'n3', 'tOffset', 'dot', 'dot2', 'srm', 'rl', 'amp1', 'amp2', 'amp3', 'blueness', 'RedVTimeCourse', 'GreenVTimeCourse', 'BlueVTimeCourse', 't1V', 't2V', 't3V', 'a1V', 'a2V', 'a3V', 'n1V', 'n2V', 'n3V', 'tVOffset', 'dotV', 'dot2V', 'srmV', 'rlV', 'amp1V', 'amp2V', 'amp3V', 'bluenessV', 'Auto', 'acfBinning', 'nSpikes', 'acfMean', 'acfRMS', 'contamination']
# param_values = ['Double', 'String', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'DoubleArray', 'DoubleArray', 'Double', 'DoubleArray', 'DoubleArray', 'Double', 'DoubleArray', 'DoubleArray', 'DoubleArray', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'DoubleArray', 'DoubleArray', 'DoubleArray', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'Double', 'DoubleArray', 'Double', 'Double', 'Double', 'Double', 'Double']

# data_sizes = [10, 9, 10, 10, 10, 10, 10, 10, 6, 6, 10, 6, 6, 10, 496, 496, 496, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 6, 6, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 1608, 10, 10, 10, 10, 10]
data_sizes = [10, 9, 10, 10, 10, 10, 10, 10, 6, 6, 10, 6, 6, 10, 408, 408, 408, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 6, 6, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 1608, 10, 10, 10, 10]

paramNames = dict()
for i, key in enumerate(param_keys):
    paramNames[key] = param_values[i]

# paramNames = {'ID': 'Double', 'classID': 'String', 'x0': 'Double', 'y0': 'Double', 'SigmaX': 'Double', 'SigmaY': 'Double',
#     'Theta': 'Double', 'gAmp': 'Double', 'contourX': 'DoubleArray', 'contourY': 'DoubleArray', 'contourArea': 'Double',
#     'simpleContourX': 'DoubleArray', 'simpleContourY':  'DoubleArray', 'simpleContourArea': 'Double', 'EIx0': 'Double',
#     'EIy0': 'Double', 'EISigmaX': 'Double', 'EISigmaY': 'Double', 'EITheta': 'Double', 'EIgAmp': 'Double', 'RedTimeCourse': 'DoubleArray', 
#     'GreenTimeCourse': 'DoubleArray', 'BlueTimeCourse': 'DoubleArray', 't1': 'Double', 't2': 'Double', 't3': 'Double', 
#     'a1': 'Double', 'a2': 'Double', 'a3': 'Double', 'n1': 'Double', 'n2': 'Double', 'n3': 'Double', 'tOffset': 'Double', 
#     'dot': 'Double', 'dot2': 'Double', 'srm': 'Double', 'rl': 'Double', 'amp1': 'Double', 'amp2': 'Double', 'amp3': 'Double',
#     'blueness': 'Double', 'RedVTimeCourse': 'DoubleArray', 'GreenVTimeCourse': 'DoubleArray', 'BlueVTimeCourse': 'DoubleArray',
#     't1V': 'Double', 't2V': 'Double', 't3V': 'Double', 'a1V': 'Double', 'a2V': 'Double', 'a3V': 'Double', 'n1V': 'Double',
#     'n2V': 'Double', 'n3V': 'Double', 'tVOffset': 'Double', 'dotV': 'Double', 'dot2V': 'Double', 'srmV': 'Double', 'rlV': 'Double',
#     'amp1V': 'Double', 'amp2V': 'Double', 'amp3V': 'Double', 'bluenessV': 'Double', 'Auto': 'DoubleArray', 'acfBinning': 'Double',
#     'nSpikes': 'Double', 'acfMean': 'Double', 'acfRMS': 'Double', 'contamination': 'Double'}

# data_sizes = [10,9,10,10,10,10,10,10,6,6,10,6,6,10,496,496,496,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,6,6,6,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,1608,10,10,10,10]

filepath = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220823C/noise/noise.params'
mpath = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220823C/20220823C_yass_sta.mat'
mdic = hdf5storage.loadmat(mpath)
timecourse_matrix = mdic['timecourse_matrix']
gauss_params = mdic['gauss_params']
acf = mdic['acf']
isi = mdic['isi']
sta = mdic['sta']
cluster_id = mdic['cluster_id']
spike_count = mdic['spike_count']

x0 = gauss_params[:, 2]
y0 = sta.shape[2] - gauss_params[:, 3]
sigma_x = gauss_params[:, 4]
sigma_y = gauss_params[:, 5]
theta = gauss_params[:, 1]
# spike_count = np.sum(isi, axis=1)

# timecourse_matrix = timecourse_matrix[:,:50,:]
# acf = acf[:,:200]

from vision_utils import ParamsWriter, STAWriter
wr = ParamsWriter(filepath=filepath, cluster_id=cluster_id)
wr.write(timecourse_matrix[:,::-1,:], acf, spike_count[:,0], x0, y0, sigma_x, sigma_y, theta, isi_binning=1.0)

# Test the STAWriter class.
sw = STAWriter(filepath='/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220823C/noise/noise.sta')
sw.write(sta, cluster_id, stixel_size=30.0, frame_refresh=1000/120)

# Old stuff

# fp = open(filepath, 'wb')

# nParams = len(param_keys)
# n_cells = timecourse_matrix.shape[0]
# maxNeurons = 10000
# # HEADER_SIZE = 1500
# cellNameString = [0, 0, 0, 3, 65, 108, 108]

# header = struct.pack('>III', nParams, n_cells, maxNeurons)
# fp.write(header)

# for param in paramNames:
#     for i in range(3):
#         fp.write(struct.pack('>b', 0))
#     l = len(param)
#     fp.write(struct.pack('>b', l))
#     for k in range(l):
#         fp.write(struct.pack('>s', bytes(param[k], 'utf-8')))
#     for i in range(3):
#         fp.write(struct.pack('>b', 0))
#     value = paramNames[param]
#     l = len(value)
#     fp.write(struct.pack('>b', l))
#     for k in range(l):
#         fp.write(struct.pack('>s', bytes(value[k],'utf-8')))



# # locs = [2481307]
# # for i, value in enumerate(data_sizes):
# #     locs.append(locs[i] + data_sizes[i])
# new_locs = []
# new_locs.append(2441282)
# count = 0
# for cell in range(n_cells):
#     for param in range(nParams):
#         new_locs.append(new_locs[count] + data_sizes[param])
#         count += 1
# new_locs = new_locs[:-1]

# # new_locs = np.zeros((nParams*n_cells+1,),dtype=int)
# # new_locs[0] = 2481307 #2441282
# # for i in range(nParams):
# #     new_locs[i+1] = new_locs[i] + data_sizes[i]

# # for cellID in range(n_cells-1):
# #     for i in range(nParams):
# #         new_locs[(cellID+1)*nParams + i + 1] = new_locs[(cellID+1)*nParams + i] + data_sizes[i]

# # new_locs = new_locs[:-1]

# # Write the locations.
# fp.seek(nParams*21-1, os.SEEK_SET)
# for loc in new_locs:
#     fp.write(struct.pack('>I', loc-2))

# end_of_fid = new_locs[-1]+8
# num = end_of_fid - fp.tell()
# for i in range(num):
#     fp.write(struct.pack('>b', 0))

# # Add the cells.
# param_count = 0
# for cellID in range(n_cells):
#     for i in range(nParams):
#         fp.seek(new_locs[param_count]-2, os.SEEK_SET)
#         # fp.seek(new_locs[i + cellID*nParams], os.SEEK_SET)
#         param_count += 1
#         param_key = list(paramNames.keys())[i]
#         param_value = list(paramNames.values())[i]
#         if param_value == 'String':
#             fp.write(struct.pack('>h', 327))
#             for v in cellNameString:
#                 fp.write(struct.pack('>b', v))
#         elif list(paramNames.values())[i] == 'DoubleArray':
#             is_array = False
#             if (param_key == 'RedTimeCourse' or 
#                 param_key == 'GreenTimeCourse' or
#                 param_key == 'BlueTimeCourse'):
#                 fp.write(struct.pack('>h', 319))
#                 fp.write(struct.pack('>h', 404))
#                 fp.write(struct.pack('>h', 0))
#                 fp.write(struct.pack('>h', 50))
#                 is_array = True
#             elif param_key == 'Auto':
#                 fp.write(struct.pack('>h', 319))
#                 fp.write(struct.pack('>h', 1604))
#                 fp.write(struct.pack('>h', 0))
#                 fp.write(struct.pack('>h', 200))
#                 is_array = True
#             else:
#                 fp.write(struct.pack('>h', 260))
#                 fp.write(struct.pack('>h', 0))
#                 fp.write(struct.pack('>h', 0))
            
#             if is_array:
#                 if param_key == 'RedTimeCourse':
#                     my_array = timecourse_matrix[cellID, :50, 0]
#                 elif param_key == 'GreenTimeCourse':
#                     my_array = timecourse_matrix[cellID, :50, 1]
#                 elif param_key == 'BlueTimeCourse':
#                     my_array = timecourse_matrix[cellID, :50, 2]
#                 else:
#                     my_array = acf[cellID, :]
#                 for val in my_array:
#                     fp.write(struct.pack('>d', val))
#         else:
#             if param_key == 'ID':
#                 k = cellID + 1
#             elif param_key == 'x0':
#                 k = np.max([1.0, gauss_params[cellID, 2]])
#             elif param_key == 'y0':
#                 k = sta.shape[2] - gauss_params[cellID, 3]
#                 k = np.max([1.0, k])
#             elif param_key == 'SigmaX':
#                 k = np.max([1.0, gauss_params[cellID, 4]])
#             elif param_key == 'SigmaY':
#                 k = np.max([1.0, gauss_params[cellID, 5]])
#             elif param_key == 'Theta':
#                 k = gauss_params[cellID, 1]
#             elif param_key == 'nSpikes':
#                 k = np.sum(isi[cellID, :])
#             elif param_key == 'tOffset':
#                 k = 0.0
#             elif param_key == 'acfBinning':
#                 k = 1.0
#             elif param_key == 'acfMean':
#                 # k = 50.872
#                 k = np.mean(acf[cellID, :])
#             elif param_key == 'acfRMS':
#                 k = 43.1999
#             elif (param_key == 'blueness' or param_key == 'amp1' or param_key == 'amp2' or param_key == 'amp3'):
#                 k = 0.0
#             else:
#                 k = 0.0
            
#             fp.write(struct.pack('>h', 200))
#             fp.write(struct.pack('>d', k))

# fp.close()


# import codecs

# dummy_params_file = '/Users/michaelmanookin/Documents/GitRepos/Manookin-Lab/data999.old_params'
# f = open(dummy_params_file, 'rb')

# f.seek(0, os.SEEK_SET)
# nParams = struct.unpack('>I', f.read(4))[0]
# ncells = struct.unpack('>I', f.read(4))[0]
# struct.unpack('>I', f.read(4))

# l_keys = []
# l_values = []
# i_keys = []
# i_values = []

# byte_count = 3
# for i in range(nParams):
#     for kk in range(3):
#         f.read(1)
#         byte_count += 1
#     # Key
#     l = struct.unpack('>b', f.read(1))[0]
#     byte_count += 1
#     l_keys.append(l)
#     string = ''
#     for kk in range(l):
#         b = struct.unpack('>s', f.read(1))[0]
#         byte_count += 1
#         string += codecs.decode(b)
#     i_keys.append(string)
#     for kk in range(3):
#         f.read(1)
#         byte_count += 1
#     # Value
#     l = struct.unpack('>b', f.read(1))[0]
#     byte_count += 1
#     l_values.append(l)
#     string = ''
#     for kk in range(l):
#         b = struct.unpack('>s', f.read(1))[0]
#         byte_count += 1
#         string += codecs.decode(b)
#     i_values.append(string)

# locs = []
# for i in range(nParams*n_cells):
#     locs.append(struct.unpack('>I', f.read(4))[0])

# # f.seek(1307, os.SEEK_SET)

# f.seek(locs[1]+2, os.SEEK_SET)
# foo = []
# for i in range(7):
#     foo.append(struct.unpack('>b', f.read(1))[0])

# # d_locs = []
# # for i in range(len(locs)-1):
# #     d_locs.append(locs[i+1] - locs[i])

# nums = []
# for i in range(nParams+1):
#     nums.append(struct.unpack('>I', f.read(4))[0])

# diffs = []
# for i in range(nParams):
#     diffs.append(nums[i+1] - nums[i])

# nums = []
# while f.tell() < 1500:
#     nums.append(struct.unpack('>I', f.read(4))[0])




# string = ''
# my_s = []
# for i in range(1500):
#     b = struct.unpack('>s', f.read(1))[0]
#     string += codecs.decode(b)
#     my_s.append(struct.unpack('>b', f.read(1))[0])

# string = ''.join(map(chr, my_s))


