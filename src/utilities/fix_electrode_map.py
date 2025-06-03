import visionloader as vl
import visionwriter as vw
import os
import argparse

# Define the mapping from the original electrode order to the new order.
idx_map = [3,   2,   1,   0,   7,   6,   5,   4,  11,  10,   9,   8,  15,  14,  13,  12,  19,  18,  17,  16,  23,  22,  21,  20,  27,  26,  25,  24,  31,  30,  29,  28,  35,  34,  33,  32,  39,  38,  37,  36,  43,  42,  41,  40,  47,  46,  45,  44,  51,  50,  49,  48,  55,  54,  53,  52,  59,  58,  57,  56,  63,  62,  61,  60,  67,  66,  65,  64,  71,  70,  69,  68,  75,  74,  73,  72,  79,  78,  77,  76,  83,  82,  81,  80,  87,  86,  85,  84,  91,  90,  89,  88,  95,  94,  93,  92,  99,  98,  97,  96, 103, 102, 101, 100, 107, 106, 105, 104, 111, 110, 109, 108, 115, 114, 113, 112, 119, 118, 117, 116, 123, 122, 121, 120, 127, 126, 125, 124, 135, 134, 133, 132, 131, 130, 129, 128, 143, 142, 141, 140, 139, 138, 137, 136, 151, 150, 149, 148, 147, 146, 145, 144, 159, 158, 157, 156, 155, 154, 153, 152, 167, 166, 165, 164, 163, 162, 161, 160, 175, 174, 173, 172, 171, 170, 169, 168, 183, 182, 181, 180, 179, 178, 177, 176, 191, 190, 189, 188, 187, 186, 185, 184, 199, 198, 197, 196, 195, 194, 193, 192, 207, 206, 205, 204, 203, 202, 201, 200, 215, 214, 213, 212, 211, 210, 209, 208, 223, 222, 221, 220, 219, 218, 217, 216, 231, 230, 229, 228, 227, 226, 225, 224, 239, 238, 237, 236, 235, 234, 233, 232, 247, 246, 245, 244, 243, 242, 241, 240, 255, 254, 253, 252, 251, 250, 249, 248, 259, 258, 257, 256, 263, 262, 261, 260, 267, 266, 265, 264, 271, 270, 269, 268, 275, 274, 273, 272, 279, 278, 277, 276, 283, 282, 281, 280, 287, 286, 285, 284, 291, 290, 289, 288, 295, 294, 293, 292, 299, 298, 297, 296, 303, 302, 301, 300, 307, 306, 305, 304, 311, 310, 309, 308, 315, 314, 313, 312, 319, 318, 317, 316, 323, 322, 321, 320, 327, 326, 325, 324, 331, 330, 329, 328, 335, 334, 333, 332, 339, 338, 337, 336, 343, 342, 341, 340, 347, 346, 345, 344, 351, 350, 349, 348, 355, 354, 353, 352, 359, 358, 357, 356, 363, 362, 361, 360, 367, 366, 365, 364, 371, 370, 369, 368, 375, 374, 373, 372, 379, 378, 377, 376, 383, 382, 381, 380, 391, 390, 389, 388, 387, 386, 385, 384, 399, 398, 397, 396, 395, 394, 393, 392, 407, 406, 405, 404, 403, 402, 401, 400, 415, 414, 413, 412, 411, 410, 409, 408, 423, 422, 421, 420, 419, 418, 417, 416, 431, 430, 429, 428, 427, 426, 425, 424, 439, 438, 437, 436, 435, 434, 433, 432, 447, 446, 445, 444, 443, 442, 441, 440, 455, 454, 453, 452, 451, 450, 449, 448, 463, 462, 461, 460, 459, 458, 457, 456, 471, 470, 469, 468, 467, 466, 465, 464, 479, 478, 477, 476, 475, 474, 473, 472, 487, 486, 485, 484, 483, 482, 481, 480, 495, 494, 493, 492, 491, 490, 489, 488, 503, 502, 501, 500, 499, 498, 497, 496, 511, 510, 509, 508, 507, 506, 505, 504]

if __name__ == '__main__':
    '''
    Average and save all of the EI files that were combined for YASS spike sorting.
    Example: 
        python fix_electrode_map.py /gscratch/retina/data/sorted/20220531C -a yass -f data005 data006
    '''

    parser = argparse.ArgumentParser(description='Average EI files..')
    parser.add_argument('rootdir', type=str, help='Path to sorted spikes and .ei files')
    parser.add_argument('-f','--file', nargs='+', default=None, type=str, help='name of data file(s) to analyze i.e. data003 (default: None)')
    parser.add_argument('-a','--algorithm', default='yass', type=str, help='sorting algorithm (default: yass)')

    args = parser.parse_args()

    if type(args.file) is str:
        file_names = [args.file]
    else:
        file_names = args.file

    for count, data_file in enumerate(file_names):
        ei_path = os.path.join(args.rootdir, data_file, args.algorithm)
        eir = vl.EIReader(ei_path, data_file)

        writeable_ei_by_cell_id = dict()

        electrode_map = eir.get_electrode_map()

        if count == 0:
            array_id = eir.array_id
        eis = eir.get_all_eis_by_cell_id()
        eir.close()

        for ei_cell_id, readable_ei in eis.items():
            ei_matrix = readable_ei.ei
            ei_error = readable_ei.ei_error
            n_spikes = readable_ei.n_spikes
            ei_matrix = ei_matrix[idx_map,:] # reorder the electrodes
            ei_error = ei_error[idx_map,:] # reorder the electrodes
            if count == 0:
                left_samples = readable_ei.nl_points
                right_samples = readable_ei.nr_points
                count += 1
            writeable_ei = vw.WriteableEIData(ei_matrix, ei_error, n_spikes)
            writeable_ei_by_cell_id[ei_cell_id] = writeable_ei

        # Write the updated EI file.
        ei_writer = vw.EIWriter(ei_path, data_file, left_samples, right_samples, array_id, True)
        ei_writer.write_eis_by_cell_id(writeable_ei_by_cell_id)
    
