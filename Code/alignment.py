from numpy import array, zeros, full, argmin, inf, ndim
from scipy.spatial.distance import cdist
from math import isinf
from dtw import *
import numpy as np
import local_dtw

class TrajectoryAligner(object):

    # def __init__(self):
    #     self.w = w
    #     self.warp = warp

    def align(self,matrix,keep_internals=False,open_end=False,open_begin=False,step_pattern="symmetric2"):
        print("THE CLASS IS CANCELED - MOVE TO CALLING THE STATIC METHOD")
        return dtw(np.ascontiguousarray(matrix),keep_internals=keep_internals,open_end=open_end,open_begin=open_begin,step_pattern=step_pattern)

def align(x,y,distance_matrix,how,threshold,step_pattern="symmetric2",input_type="distance",window_type="none",window_args={}):
    assert how in ["local", "global","open_begin_end","open_begin"]
    if input_type == "distance":
        if how=="local":
            ii,jj = local_dtw.local_align(distance_matrix,threshold)
            result = AlignmentResult("local",ii=ii,jj=jj)
            return result
        if how=="open_begin_end":
            # dtw_object = dtw(np.ascontiguousarray(distance_matrix), step_pattern=rabinerJuangStepPattern(4,"c",True), open_end=True, open_begin=True)
            dtw_object = dtw(np.ascontiguousarray(distance_matrix), step_pattern="asymmetric", open_end=True, open_begin=True)
            ii, jj = dtw_object.index1, dtw_object.index2
            # print(ii,jj)
            result = AlignmentResult("open_begin_end", ii=ii, jj=jj, DTW=dtw_object)
            return result

        if how=="open_begin":
            dtw_object = dtw(np.ascontiguousarray(distance_matrix), step_pattern="asymmetric", open_end=False, open_begin=True)
            ii, jj = dtw_object.index1, dtw_object.index2
            # print(ii,jj)
            result = AlignmentResult("open_begin", ii=ii, jj=jj, DTW=dtw_object)
            return result


        if how=="global":
            return dtw(np.ascontiguousarray(distance_matrix),step_pattern=step_pattern,window_type=window_type,window_args=window_args)
    elif input_type == "measurements":
        if how=="local":
            print("NOT YET IMPLEMENTED")
            # ii,jj = local_dtw.local_align(distance_matrix,threshold)
            # result = AlignmentResult("local",ii=ii,jj=jj)
            # return result
        if how=="open_begin_end":
            dtw_object = dtw(x,y, step_pattern="asymmetric", open_end=True, open_begin=True, keep_internals=True)
            ii, jj = dtw_object.index1, dtw_object.index2
            # print(ii,jj)
            result = AlignmentResult("open_begin_end", ii=ii, jj=jj, DTW=dtw_object)
            # dtw_object.plot(type="twoway", offset=-2)
            return result
        if how=="global":
            print("NOT YET IMPLEMENTED")





class AlignmentResult(object):
    def __init__(self,how,ii=None,jj=None,DTW=None):
        self.dtw_object=DTW
        self.how=how
        if how=="global":
            self.normalized_distance = DTW.normalizedDistance
        self.ii = ii
        self.jj = jj

    def distance_on_matrix(self,M):
        index1,index2 = self.dtw_object.index1,self.dtw_object.index2
        n,m = len(index1),len(index2)
        assert n==m
        distance=0
        for i in range(n):
            a=index1[i]
            b=index2[i]
            distance+=M.iloc[a,b]

        return distance/(n+m)


# def calculate_alignment_metrices():

def foo():
    print("boo")
