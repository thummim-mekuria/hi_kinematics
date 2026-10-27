#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd


import matplotlib.pyplot as plt


np.where(names=='97.0-10.0-35.0')


lbv[1304700]


len(names)


1305800/2845800


names[1305800]


lbv.transpose()[2][855800:905800]


len(np.arange(180,-180,-1))*len(np.arange(-50,52,2))*len(np.arange(-100,101,1.3))


2772000/1e6


lbv = np.array(np.meshgrid(np.arange(180,-180,-1),np.arange(-50,52,2),np.arange(-100,101,1.3))).T.reshape(-1,3)


names = ['-'.join(list([str(round(lbv[j][i],2)) for i in range(3)])) for j in range(len(lbv))]


names[57000]


df = pd.DataFrame(lbv,columns=['lon','lat','vel'])


df['name'] = names


df['vunc'] = np.full(len(df),0.6)


df['pfar'] = np.full(len(df),0.5)


df = df[['name','lon','lat','vel','vunc','pfar']]


df.to_csv('sources_info.inp',index=False)


names[::100]
