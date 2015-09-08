############################################
###
#
#   This file contains codes for derivatives
#
#   @Author: Jerry Shi
#   @Date: May 23rd, 2013
#
###
############################################




import numpy as np


###
# @Method: five-point stencil
# @Input: data as numpy ndarray
#   Columns:    [0]x    [1]y
# @Output: first derivative and second derivative
#   Columns:    [0]x    [1]y    [2]dy_dx    [4]d2y_dx2
def derivatives(data, h=1):
    num = data.shape[0]
    result = np.zeros((num-4*h,4))
    dx = data[h,0] - data[0,0]
    for i in range(2*h, num-2*h):
        result[i-2*h,:] = [  data[i,0]   ,
                    data[i,1]   ,
                    (-data[i+2*h,1] + 8*data[i+h,1] -8*data[i-h,1]+ data[i-2*h,1])/(12*dx)  ,
                    (-data[i+2*h,1] + 16*data[i+h,1]-30*data[i,1] +16*data[i-h,1]- data[i-2*h,1])/(12*np.square(dx))  ]
    return result

