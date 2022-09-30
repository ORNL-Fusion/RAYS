import numpy as np
import math

# Parse in_list into segments of length len(ref_list).  Subtract rev_list
# from segments and calculate 2-norm. Use numpy arrrays

# Convert list to numpy array.  Split list into sublists of length len_sub_list
def split_flat_list(in_list, len_sub_list):

# check that len(in_list) is a multiple of len_sub_list
    if(float(int(len(in_list)/len_sub_list)) != len(in_list)/len_sub_list):
        message = 'split_flat_list: len(in_list not divoisible by len_sub_list)'
        print(message)
        raise

    d1 = int(len(in_list)/len_sub_list)
    x = np.array(in_list)
    return x.reshape(d1, len_sub_list)

#_________________________________________________________________________________________________
# Takes a flat list of numbers consisting of a concatenation of vectors each of length 
# len_sub_list.  Converts in_list to numpy array and splits into individual vector sub arrays.
# Takes a slice of each sub array with given offset and subtracts off a ref_array of the
# same length as the slice.  Calculates the 2-norm of the difference array.  Returns a
# numpy array containing the difference 2-norm for each vector in in_list
# 
def deviation_norm_sub_subarray(in_list, len_sub_list, ref_list, offset, norm = 1.):
    p_array = split_flat_list(in_list, len_sub_list)
    ref_array = np.array(ref_list)
    sub_list = []
    for sub_array in p_array:
        sub_sub_array = sub_array[offset : offset + len(ref_array)]
        diff = sub_sub_array - ref_array
        x = math.sqrt(np.sum(np.square(diff)))/norm
        sub_list.append(x)
    return np.array(sub_list)   


#_________________________________________________________________________________________________

if __name__ == '__main__':

    in_list = [1.,2.,3.,4.,5.,6.,7.,8.]
    len_sub_list = 4
    print(split_flat_list(in_list, len_sub_list))
    
    print(' ')
    ref_list = [1.,1.]
    offset = 2
    dev = deviation_norm_sub_subarray(in_list, len_sub_list, ref_list, offset, 2.)
    print(dev)