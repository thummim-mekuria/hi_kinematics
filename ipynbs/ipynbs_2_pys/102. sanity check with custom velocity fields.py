#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd

from astropy.modeling import models,fitting
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c


x,y=np.meshgrid(np.arange(0,10,0.1),np.arange(0,100,10))


x.shape


xx,yy = np.meshgrid(np.arange(100),np.arange(10))


init = models.Planar2D(1.,1.,0.)
fttr = fitting.LMLSQFitter()


xfit = fttr(init,xx,yy,x)


xfit


yfit = fttr(init,xx,yy,y)


yfit


deg10 = x*np.cos(10*u.deg) + y*np.sin(10*u.deg)


plt.imshow(deg10,origin='lower')
plt.plot(np.arange(100),np.arange(100)*np.cos(87*u.deg),color='red')


deg10_fit = fttr(init,xx,yy,deg10)


deg10_fit


np.arctan(deg10_fit.slope_y/deg10_fit.slope_x).to(u.deg)


initin = models.Linear1D(1.,0.)
ftrtr = fitting.LinearLSQFitter()


1.96961551*0.05


0.0984807755/300


ftrtr(initin,np.ravel(xx*300),np.ravel(deg10))


ftrtr(initin,np.ravel(yy),np.ravel(deg10))


plt.scatter(np.ravel(xx),np.ravel(deg10))


plt.scatter(np.ravel(yy),np.ravel(deg10))
plt.plot(np.ravel(yy),np.ravel(yy)*1.736,color='red')


#double check that the Omega is being recovered correctly for angled velocity field

#calculate the omega,theta for both the co and hi


#see if i can recover nias Omega and theta if i do the dv/dl & dv/db on pixels 
#and divide by d -> ig kind of the straightforward application of the formula from goodman


#also add in noise (np.random.random()) and see if the slopes 
#from fitting plane vs fitting x & y individually are still consistent


plt.imshow(y)
