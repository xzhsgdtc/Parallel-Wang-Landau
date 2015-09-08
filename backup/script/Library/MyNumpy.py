############################################
###
#
#   This file contais methods used to process
#   ndarray in numpy
#
#   @Author: Jerry Shi
#   @Date: May 7th, 2013
#
###
############################################


import numpy as np
import sys

#################
# @description:
#   This method try to extract certain columns of given ndarray
#   and return the "slimed" ndarray
# @input:
#   1. ndarray
#   2. column indexes
# @return:
#   1.slimed ndarray
#
def slim(content, index):
    size = len(index)
    if max(index) > content.shape[1]:
        print "Error, MyNumpy.slim(...), index (%d) out of range!" % max(index)
        sys.exit(1)
    length = content.shape[0]
    result = np.reshape(content[:,index[0]], (length, 1))
    for i in range(1,size):
        result = np.append(result, np.reshape(content[:,index[i]], (length, 1)), axis = 1)
    return result


