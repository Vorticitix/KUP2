#%%
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import sys
from scipy.stats import ttest_ind
import xarray as xr
path = '/storage/climatestor/PleioCEP/doensen/data/'
path2 = '/storage/climatestor/PleioCEP/doensen/data/extracted/'
#%%

loc = 'eur'
file_cyc = 'fort_36_total_{}_full.txt'.format(loc)
file_T = 'T850/T850_0005_3601_anom_full_{}.nc'.format(loc)
file_PREC = 'PREC/PREC_anom_{}.nc'.format(loc)
file_NAOI = 'try_PSL/NAOI_full.nc'

T = xr.open_dataset(path2+file_T).squeeze()
prec = xr.open_dataset(path2+file_PREC).squeeze()
prec['PRECT'] = prec.PRECT*1e3*86400
NAOI = xr.open_dataset(path2+file_NAOI).squeeze()


df=pd.read_csv(path+file_cyc).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
#df_rcp=pd.read_csv(path+file_rcp).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
# med = pd.read_csv(path+file_med).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
# rcp = pd.read_csv(path+file_rcp).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
# rcp = rcp.where(rcp['year']>=3515).dropna()
# df = pd.concat([med,rcp])
df = df.where(df['year']>=5).dropna(how='all')
df = df.where((df['month']>=12) | (df['month']<=2)).dropna(how='all')


T_sel = T.sel(time=T.time[T.time.dt.year>5]).T
T_sel = T_sel.sel(time=T_sel.time[(T_sel.time.dt.month>=12)|(T_sel.time.dt.month<=2)])

_, index = np.unique(prec.time, return_index=True)
prec = prec.isel(time=index)
prec_sel = prec.sel(time=prec.time[((prec.time.dt.month>=12)|(prec.time.dt.month<=2))&(prec.time.dt.year>5)]).PRECT

NAOI_sel = NAOI.sel(time=NAOI.time[NAOI.time.dt.year>5])
NAOI_sel = NAOI_sel.sel(time=NAOI_sel.time[(NAOI_sel.time.dt.month>=12)|(NAOI_sel.time.dt.month<=2)])
NAOI_sel = NAOI_sel.reindex(time=T_sel.time.values).PSL
#df['yearseason']=(df['year'].astype(int)-1502).astype(str)+'_'+df['month'].map(seasons_dic)
df['yearmonth']=np.around(df['year'].values + df['month'].values/12 -1/12,3)
#df['cen'] = (df['year']-1502)//100





#df_groupby= df.groupby('yearseason')
df_gb_max = df[['zrad','zdep','precmean','preccmean','preclmean','wsmean','iic']].groupby('iic').max()
df_gb_min = df[['slp','iic','agetot']].groupby('iic').min()
df_gb_yas = df[['year','yearmonth','iic']].groupby('iic').first()
df_gb_iic = pd.concat([df_gb_max,df_gb_min,df_gb_yas],axis=1)
df_gb_iic = df_gb_iic.where((df_gb_iic['slp']!=0)
                            &(df_gb_iic['zrad']!=0)
                            &(df_gb_iic['zdep']!=0)).dropna(how='all')
agetot = df_gb_iic.groupby('yearmonth').mean()[['agetot','year']]
df_gb_iic = df_gb_iic.drop(columns='agetot')
#%%
#Load in list of strongest eruptions
erup = pd.read_csv(path+'cyclone_volcano/eruptions.txt')
#Sort eruption chronologically and test wether there is an eruption 10 years before an eruption and 5 years after an eruption
erup_sort = erup.sort_values('time')
bool_arr = np.full(erup.time.values.shape,True)
for i in np.arange(len(erup.time)):
    if i < len(erup_sort.time.values)-1:
        diff_lower = np.abs(erup_sort.time.values[i-1]-erup_sort.time.values[i])
        diff_higher = np.abs(erup_sort.time.values[i+1]-erup_sort.time.values[i])
        if (diff_lower<10)|(diff_higher<5):
            bool_arr[i]=False
    else:
        bool_arr[i]=False

#Only keep eruptions that do not overlap other eruptions
erup_sep_sort = erup_sort[bool_arr==True]
erup_sep = erup_sep_sort.sort_values('colmass',ascending=False)
#erup_sep = erup_sep[(erup_sep['lat']<20)&(erup_sep['lat']>-20)]
#Select 30 strongest eruptions
strong = erup_sep['time'][:30]
# %%
less = []; more1 =  []; more2 =  []; more3 =  []; more4 =  []; more5 =  []
less_xr = []; more1_xr =  []; more2_xr =  []; more3_xr =  []; more4_xr =  []; more5_xr =  []
less_at = []; more1_at =  []; more2_at =  []; more3_at =  []; more4_at =  []; more5_at =  []


vars = ['precmean','wsmean','zdep','agetot','T','PREC','NAOI']
for year in strong:
    dic = dict.fromkeys([1, 2, 3, 4, 5])
    intyear = int(year)
    yearsless = np.arange(intyear - 5,intyear)
    yearsmore = np.arange(intyear+1,intyear+6)
    assert(intyear+5==yearsmore[-1])
    assert(intyear-5==yearsless[0])
    dfless = df_gb_iic[np.isin(df_gb_iic['year'],yearsless)]
    T_less = T_sel[np.isin(T_sel.time.dt.year,yearsless)]
    prec_less = prec_sel[np.isin(prec_sel.time.dt.year,yearsless)]
    NAOI_less = NAOI_sel[np.isin(NAOI_sel.time.dt.year,yearsless)]
    agetot_less = agetot[np.isin(agetot['year'],yearsless)]
    mean = dfless.mean()
    mean_T = T_less.mean()
    mean_P = prec_less.mean()
    mean_NAOI =NAOI_less.mean()
    mean_agetot = agetot_less.mean()
    mean_agetot['year'] = 0
    dfless_norm = dfless - mean
    dfless_norm.loc[:,'year']=dfless['year']
    T_less_norm = T_less - float(mean_T)
    prec_less_norm = prec_less - float(mean_P)
    NAOI_less_norm = NAOI_less - float(mean_NAOI)
    agetot_less_norm = agetot_less - mean_agetot
    agetot_less_norm = agetot_less_norm['agetot']
    dfmore = df_gb_iic[np.isin(df_gb_iic['year'],yearsmore)]
    T_more = T_sel[np.isin(T_sel.time.dt.year,yearsmore)]
    prec_more = prec_sel[np.isin(prec_sel.time.dt.year,yearsmore)]
    NAOI_more = NAOI_sel[np.isin(NAOI_sel.time.dt.year,yearsmore)]
    agetot_more = agetot[np.isin(agetot['year'],yearsmore)]
    dfmore_norm = dfmore - mean
    dfmore_norm.loc[:,'year']=dfmore['year']
    T_more_norm = T_more - float(mean_T)
    prec_more_norm = prec_more - float(mean_P)
    NAOI_more_norm = NAOI_more - float(mean_NAOI)
    agetot_more_norm = agetot_more - mean_agetot
    for i,yearmore in enumerate(yearsmore):
        T_more_year = T_more_norm[T_more_norm.time.dt.year==yearmore]
        prec_more_year = prec_more_norm[prec_more_norm.time.dt.year==yearmore]
        NAOI_more_year = NAOI_more_norm[NAOI_more_norm.time.dt.year==yearmore]
        dfmore_norm_year = dfmore_norm[dfmore_norm['year']==yearmore]
        agetot_norm_year = agetot_more_norm[agetot_more_norm['year']==yearmore]['agetot']
        dic[i+1]=[dfmore_norm_year,agetot_norm_year,T_more_year,prec_more_year,NAOI_more_year]
        assert((T_more_year.time.dt.year==yearmore).all())
        assert((prec_more_year.time.dt.year==yearmore).all())
        assert((NAOI_more_year.time.dt.year==yearmore).all())
        assert((dfmore_norm_year['year']==yearmore).all())
        
    fig,axz= plt.subplots(2,4,figsize=(16,9))
    for i,var in enumerate(vars):        
        if i < 3:
            var_less = dfless_norm[var]
            var_more1 = dic[1][0][var]
            var_more2 = dic[2][0][var]
            var_more3 = dic[3][0][var]
            var_more4 = dic[4][0][var]
            var_more5 = dic[5][0][var]
        elif i >= 3:
            if i == 3:
                var_less = agetot_less_norm
            elif i == 4:
                var_less = T_less_norm
            elif i ==5:
                var_less = prec_less_norm
            elif i == 6:
                var_less = NAOI_less_norm
            var_more1 = dic[1][i-2]
            var_more2 = dic[2][i-2]
            var_more3 = dic[3][i-2]
            var_more4 = dic[4][i-2]
            var_more5 = dic[5][i-2]
        
        #p = ttest_ind(var_less,var_more,equal_var=False)[1]
        data = [var_less,var_more1,var_more2,var_more3,var_more4,var_more5]
        ax = axz.flat[i]
        ax.boxplot(data)
        ax.grid()
        ax.set_xticklabels(['Before Eruption','1','2','3','4','5'])
        ax.set_title('{}'.format(var))
    fig.suptitle('{}'.format(intyear))
    #fig.savefig('/storage/climatestor/PleioCEP/doensen/data/cyclone_volcano/figs/volc_all_{}.png'.format(intyear),dpi=300)
    plt.close(fig)
    less.append(dfless_norm)
    more1.append(dic[1][0])
    more2.append(dic[2][0])
    more3.append(dic[3][0])
    more4.append(dic[4][0])
    more5.append(dic[5][0])
    less_xr.append([T_less_norm,prec_less_norm,NAOI_less_norm])
    more1_xr.append(dic[1][2:5])
    more2_xr.append(dic[2][2:5])
    more3_xr.append(dic[3][2:5])
    more4_xr.append(dic[4][2:5])
    more5_xr.append(dic[5][2:5])
    less_at.append(agetot_less_norm)
    more1_at.append(dic[1][1])
    more2_at.append(dic[2][1])
    more3_at.append(dic[3][1])
    more4_at.append(dic[4][1])
    more5_at.append(dic[5][1])
    
df_less_total = pd.concat(less)
df_more_total1 = pd.concat(more1)
df_more_total2 = pd.concat(more2)
df_more_total3 = pd.concat(more3)
df_more_total4 = pd.concat(more4)
df_more_total5 = pd.concat(more5)
less_xr_total = [xr.concat([row[i] for row in less_xr],dim='time') for i in range(3)]
more_xr_total1 = [xr.concat([row[i] for row in more1_xr],dim='time') for i in range(3)]
more_xr_total2 = [xr.concat([row[i] for row in more2_xr],dim='time') for i in range(3)]
more_xr_total3 = [xr.concat([row[i] for row in more3_xr],dim='time') for i in range(3)]
more_xr_total4 = [xr.concat([row[i] for row in more4_xr],dim='time') for i in range(3)]
more_xr_total5 = [xr.concat([row[i] for row in more5_xr],dim='time') for i in range(3)]
at_less_total = pd.concat(less_at)
at_more_total1 = pd.concat(more1_at)
at_more_total2 = pd.concat(more2_at)
at_more_total3 = pd.concat(more3_at)
at_more_total4 = pd.concat(more4_at)
at_more_total5 = pd.concat(more5_at)
#%%
fig,axz= plt.subplots(2,4,figsize=(16,9))
for i,var in enumerate(vars):
    if i < 3:
        var_less = df_less_total[var]
        var_more1 = df_more_total1[var]
        var_more2 = df_more_total2[var]
        var_more3 = df_more_total3[var]
        var_more4 = df_more_total4[var]
        var_more5 = df_more_total5[var]
    elif i == 3:
        var_less = at_less_total
        var_more1 = at_more_total1 
        var_more2 = at_more_total2 
        var_more3 = at_more_total3 
        var_more4 = at_more_total4 
        var_more5 = at_more_total5 
    elif i >= 4:
        var_less = less_xr_total[i-4] 
        var_more1 = more_xr_total1[i-4]
        var_more2 = more_xr_total2[i-4]
        var_more3 = more_xr_total3[i-4]
        var_more4 = more_xr_total4[i-4]
        var_more5 = more_xr_total5[i-4]
    p1 = ttest_ind(var_less,var_more1,equal_var=False)[1]
    p2 = ttest_ind(var_less,var_more2,equal_var=False)[1]
    p3 = ttest_ind(var_less,var_more3,equal_var=False)[1]
    p4 = ttest_ind(var_less,var_more4,equal_var=False)[1]
    p5 = ttest_ind(var_less,var_more5,equal_var=False)[1]
    data = [var_less,var_more1,var_more2,var_more3,var_more4,var_more5]
    ax = axz.flat[i]
    ax.boxplot(data)
    ax.grid()
    ax.set_xticklabels(['Before \nEruption','1\np={:.2f}'.format(p1),
                        '2\np={:.2f}'.format(p2),
                        '3\np={:.2f}'.format(p3),
                        '4\np={:.2f}'.format(p4),
                        '5\np={:.2f}'.format(p5)])
    ax.set_title('{}'.format(var))
fig.tight_layout()
fig.savefig('/storage/climatestor/PleioCEP/doensen/data/cyclone_volcano/figs/volc_all_{}.png'.format(loc),dpi=300)
# %%
