#===================================================================================
import pandas as pd
#===================================================================================


#===================================================================================
sources_inp = pd.read_csv('sources_info.inp',delim_whitespace=True,skiprows=1)
summary_prt = pd.read_csv('summary.prt',delim_whitespace=True,skiprows=2)
#-----------------------------------------------------------------------------------
sm = summary_prt.loc[1:] #get rid of second !commented out header
sm = sm[['!','Long.','Lat.','+/-']] #only choose the columns we want
sm.columns = ['lon','lat','vel','dist'] #add column names
#-----------------------------------------------------------------------------------
#recast everything into floats from 'object'
sm['lon'] = sm['lon'].astype('float64')
sm['lat'] = sm['lat'].astype('float64')
sm['vel'] = sm['vel'].astype('float64')
sm['dist'] = sm['dist'].astype('float64')
#-----------------------------------------------------------------------------------
sources_inp.columns = ['_idx','lon','lat','vel','pfar','vunc'] #rename column names
src = sources_inp[['_idx','lon','lat','vel']] #then only choose the relevant ones
#-----------------------------------------------------------------------------------
src['_idx'] = src['_idx'].astype('int') #recast idx's as int
src['lon'] = round(src['lon'],2) #round so same sigfigs as sm
src['lat'] = round(src['lat'],2)
src['vel'] = round(src['vel'],1)
#-----------------------------------------------------------------------------------
mg = src.merge(sm,on=['lon','lat','vel'])
#===================================================================================


#===================================================================================
ppv_cat = pd.read_csv('lvs_cat.csv')
#-----------------------------------------------------------------------------------
big_leafy_cat = ppv_cat.merge(mg,on='_idx') 
#the biggg catalog -> has dend_idx, d, the stuff calculated in ppv_catalog
#-----------------------------------------------------------------------------------
rsnbl_leafy_cat = big_leafy_cat[['_idx','lon','lat','vel','dist','area_exact']]
rsnbl_leafy_cat.columns = \
['_idx','lon<deg>','lat<deg>','vel<km/s>','dist<kpc>','area_sky<deg2>']
#-----------------------------------------------------------------------------------
rsnbl_leafy_cat['area_sky<deg2>'] = round(rsnbl_leafy_cat['area_sky<deg2>'],2)
#===================================================================================


#===================================================================================
rsnbl_leafy_cat.to_csv('src_sm_ppv_merged.csv',index=False)
#===================================================================================








