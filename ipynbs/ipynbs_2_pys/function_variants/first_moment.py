#!/usr/bin/env python
# coding: utf-8



def first_moment(cube,arr,oned=False):
    zeroth = np.nansum(cube,axis=0)
    
    
    if oned:
        num = np.nansum(cube * arr)
    else:
        num = np.nansum(cube * arr[:,np.newaxis,np.newaxis],axis=0)
        
    return num/zeroth, np.log10(zeroth*1.3*1.97e20)


