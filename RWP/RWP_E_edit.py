"""
-----------------------------------------------------------------------------------------
--> Purpose:
	Compute the following Rossby Wave Packet amplitude and phase locally in space and time. 

--> Data requirements:
	Meridional wind speed (anomaly) as a function of time, latitude, and longitude at an isobaric level of the upper-troposphere, gridded at a latitude-longitude grid of at least 6-hourly and 2°x2° resolution.

--> Author: Georgios Fragkoulidis (gfragkou@uni-mainz.de)

--> Last update: 15 July 2022

--> References:
	1) Fragkoulidis, G., Wirth, V., Bossmann, P., and Fink, A. H., 2018: Linking Northern Hemisphere temperature extremes to Rossby wave packets, Quarterly Journal of the Royal Meteorological Society, doi:10.1002/qj.3228
	2) Fragkoulidis, G. and Wirth, V., 2020: Local Rossby wave packet amplitude, phase speed, and group velocity: Seasonal variability and their role in temperature extremes, Journal of Climate, doi:10.1175/JCLI-D-19-0377.1  (hereafter, FW20)  
	3) Fragkoulidis, G., 2022: Decadal variability in extratropical Rossby wave packet amplitude, phase, and phase speed, Weather Clim. Dynam. Discuss, doi:10.5194/wcd-2022-28   
-----------------------------------------------------------------------------------------
"""
#%%
###################################################
######## IMPORT MODULES ###########################
###################################################
import time
import datetime as dt
from netCDF4 import Dataset, date2num
from numpy import mean, sin, cos, arctan2, pi, zeros, zeros_like, max, min, nan, arange, imag, real, array, meshgrid, linspace
from multiprocess import Pool, cpu_count
from math import sqrt, floor, ceil
from scipy.fft import fft, ifft
from scipy import signal
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import xarray as xr
import sys


###################################################
########### FUNCTIONS #############################

#%%
def wledit(y, lat, lmin, lmax):
	"""
	############################################################################################################################
	- Wavelength restriction: Take a 1-D function and return it restricted to the (lmin,lmax) wavelength range. 
	############################################################################################################################
	- INPUT:
		* y: 1-D function (typically V)
		* lat: latitude 
		* lmin, lmax: wavelength range
	############################################################################################################################
	"""
	# Wavenumbers corresponding to the given wavelength limits (= latitude circle length / wavelength)
	wn_min = 360*111.319491667*cos(lat*pi/180)/lmax
	wn_max = 360*111.319491667*cos(lat*pi/180)/lmin     
	# Make an array of 0.1 resolution with limits that depend on restriction range and extend a bit further.       
	spectrum = [round(i,1) for i in arange(round(wn_min,1)-1.5,round(wn_max,1)+1.6,0.1).tolist()] 
	# Keep the integers in this array
	spectrum_integers = [i for i in arange(ceil(min(spectrum)),floor(max(spectrum))+1,1).tolist() if i > 0] 
	# Return a Tukey window, also known as a tapered cosine window. For alpha=0.3, 30% of the window will have a lower than 1 weight.
	window = signal.tukey(len(spectrum),alpha=0.3)   
	
	# Remove the mean before the FFT.        
	y = y - mean(y)
	# FFT  
	ffty = fft(y)
	# Go through the spectrum and put the corresponding Tukey window weights in a masking list
	mask = zeros(len(ffty))
	for k in range(0,len(spectrum_integers)):
		# The spectrum index of the integer wavenumber k:
		j = spectrum.index(spectrum_integers[k])
		# Mask contains the Tukey weights of the contributing wavenumbers. These weights depend on the range restriction we choose and the parameter alpha. The Tukey window element that corresponds to the integer wavenumber k is window[j]:
		mask[spectrum_integers[k]] = window[j] 
	    
	# Mask according to the previously defined weights
	mask_ffty = ffty*mask
	# IFFT
	y_edit = ifft(mask_ffty)
	# Since the ignored negative frequencies would contribute the same as the positive ones:
	y_edit = 2*y_edit.real 
	return y_edit
        
def wledit_lowpass(y, lat, lmin):
	"""
	Take a 1-D function and discard wavelengths below a certain threshold. 
	-y: V series to edit
	-lat: Latitude
	-lmin: minimum wavelength
	"""
	wn_max = (360*cos(lat*pi/180)*111.319491667)/lmin    
	spectrum = [round(i,1) for i in arange(0,round(wn_max,1)+1.6,0.1).tolist()] # Make an array of 0.1 resolution with limits that depend on our restriction range and extend a bit further.
	spectrum_integers = list(range(0,floor(max(spectrum))+1,1))  #[arange(0,floor(max(spectrum))+1,1).tolist()] # Take the sub-array that only contains integers
	#print(ceil(min(spectrum)),spectrum_integers)
	window = signal.tukey(len(spectrum),alpha=0.3)   # Return a Tukey window, also known as a tapered cosine window. For alpha=0.3, 30% of the window will have a lower than 1 weight.
	#print(window)
	# Change the first non-1 values to 1 so that we have a Tukey window that has a tapered edge only on the high wavenumbers.
	for i in range(len(window)):
		if window[i]<1:
			window[i]=1
		elif window[i]==1:
			break
	ffty = fft(y)
	mask = zeros(len(ffty))
	for k in range(0,len(spectrum_integers)):
		j = spectrum.index(spectrum_integers[k])
		#print(k,spectrum_integers[k],spectrum.index(spectrum_integers[k]),j,window[j])
		mask[spectrum_integers[k]] = window[j] # the mask is an array that contains the Tukey weights of the contributing wavenumbers. These weights depend on the range restriction we choose and the parameter alpha.
	mask_ffty = ffty*mask
	y_edit = ifft(mask_ffty)
	y_edit  = 2*y_edit.real # Since the ignored negative frequencies would contribute the same as the positive ones
	y_edit = y_edit - ffty.real[0]/len(ffty) # so that the zero frequency is not counted twice
	return y_edit

def hilbert(y,N,envelope=True):
	"""
	############################################################################################################################
	- Analytic signal and Envelope calculation using the Hilbert transform technique (Marple, 1999)
	############################################################################################################################
	- INPUT:
		* y: 1-D function for which we want the envelope
		* N: Size of function 
	- OUTPUT:
		* Either the analytic signal or the envelope
	############################################################################################################################
	"""
	# Check whether the signal has even or odd length
	if N%2==0:
		a = int(N/2)
	else:
		if N>1:
			a = int((N-1)/2)
		else:
			a = 0
	# FFT of y
	z = fft(y)
	# Zero-out the negative frequencies
	z[a+1:N] = 0
	# Double the positive frequencies except from the 0th and (N/2)th ones
	z = 2*z
	z[0] = z[0]/2
	if N%2==0: # for the even-length case, we also have the Nyquist frequency in the spectrum. This is shared between the positive and negative frequencies so we need to keep it once (see Marple 1999). For odd lengths, there is no Nyquist frequency in the spectrum.
		z[a] = z[a]/2
	# Inverse FFT to get the analytic signal
	an_sig = ifft(z)
	# Envelope
	env = abs(an_sig)
	if envelope:
		output = env
	else:
		output = an_sig

	return output
             
def circumference(lat):
		'''
		Circumference at a given latitude in an ellipsoid Earth = Circumference at Equator * cos(lat) = 360 * width of 1 degree at the Equator (111.319491667km) or 2*pi*Earth's radius at Equator (6378137m) * cos(lat).
		'''
		a = 6378137 # radius at the Equator
		b = 6356752.3142 # radius at the poles
		ellip = 1/sqrt(1-((a**2-b**2)/a**2)*(sin(lat*pi/180))**2)# factor that accounts for the ellipsoid shape of Earth (not perfect circle)
		lat_circum = 2 * pi * 6378137 * cos(lat*pi/180) * ellip
		return lat_circum 

def my_colormap(reversed=False):	
	C = array([
	[252,  254,  255],[250,  253,  255],[247,  252,  254],[244,  251,  254],[242,  250,  254],[239,  249,  254],[236,  248,  253],\
	[234,  247,  253],[231,  246,  253],[229,  245,  253],[226,  244,  253],[223,  243,  252],[221,  242,  252],[218,  241,  252],\
	[215,  240,  252],[213,  239,  252],[210, 238,  251],[207, 237,  251],[205, 236,  251],[202,  235,  251],[199,  234,  250],\
	[197,  233,  250],[194,  232,  250],[191,  231,  250],[189,  230,  250],[186,  229,  249],[181,  227,  249],[178,  226,  249],\
	[176,  225,  249],[170,  223,  248],[168,  222,  248],[165,  221,  248],[162,  220,  247],[157,  218,  247],[155,  216,  246],\
	[152,  214,  245],[150,  212,  243],[148,  210,  242],[143,  206,  240],[141,  204,  238],[139,  202,  237],[136,  200,  236],\
	[134,  197,  235],[132,  195,  234],[129,  193,  232],[127,  191,  231],[125,  189,  230],[123,  187,  229],[120,  185,  228],\
	[118,  183,  226],[116,  181,  225],[111,  177,  223],[109,  175,  221],[106,  173,  220],[104,  171,  219],[102,  169,  218],\
	[100,  167,  217],[97 , 165 , 215 ],[95 , 163 , 214 ],[93 , 160 , 213 ],[90 , 158 , 212 ],[88 , 156 , 211 ],[86 , 154 , 209 ],\
	[79 , 148 , 206 ],[77 , 146 , 204 ],[72 , 142 , 202 ],[72 , 143 , 198 ],[72 , 144 , 195 ],[72 , 145 , 191 ],[72 , 146 , 188 ],\
	[72 , 147 , 184 ],[72 , 148 , 181 ],[72 , 149 , 177 ],[72 , 150 , 173 ],[72 , 151 , 170 ],[72 , 153 , 166 ],[72 , 154 , 163 ],\
	[72 , 155 , 159 ],[72 , 156 , 156 ],[72 , 157 , 152 ],[72 , 158 , 148 ],[72 , 159 , 145 ],[72 , 160 , 141 ],[72 , 161 , 138 ],\
	[73 , 162 , 134 ],[73 , 163 , 131 ],[73 , 164 , 127 ],[73 , 165 , 124 ],[73 , 166 , 120 ],[73 , 167 , 116 ],[73 , 168 , 113 ],\
	[73 , 169 , 109 ],[73 , 170 , 106 ],[73 , 172 , 102 ],[73 , 173 ,  99 ],[73 , 174 ,  95 ],[73 , 175 ,  91 ],[73 , 176 ,  88 ],\
	[73 , 177 ,  84 ],[73 , 178 ,  81 ],[73 , 179 ,  77 ],[73 , 181 ,  70 ],[78 , 182 ,  71 ],[83 , 184 ,  71 ],[87 , 185 ,  72 ],\
	[92 , 187 ,  72 ],[97 , 188 ,  73 ],[106,  191,   74],[111,  192,   75],[116,  193,   75],[121,  195,   76],[126,  196,   77],\
	[130,  198,   77],[135,  199,   78],[140,  200,   78],[145,  202,   79],[150,  203,   80],[154,  204,   80],[159,  206,   81],\
	[164,  207,   81],[169,  209,   82],[173,  210,   82],[178,  211,   83],[183,  213,   84],[188,  214,   84],[193,  215,   85],\
	[197,  217,   85],[202,  218,   86],[207,  220,   87],[212,  221,   87],[217,  222,   88],[221,  224,   88],[226,  225,   89],\
	[231,  226,   90],[240,  229,   91],[245,  231,   91],[250,  229,   91],[250,  225,   89],[250,  222,   88],[249,  218,   86],\
	[249,  212,   84],[249,  208,   82],[249,  205,   81],[249,  201,   80],[249 , 198,   78],[249 , 195,   77],[248 , 191,   75],\
	[248 , 188,   74],[248 , 181,   71],[248 , 178,   70],[248 , 174,   69],[248 , 171,   67],[247 , 167,   66],[247 , 164,   64],\
	[247 , 160,   63],[247 , 157,   62],[247 , 154,   60],[247 , 150,   59],[247 , 147,   58],[246 , 143,   56],[246 , 140,   55],\
	[246 , 137,   53],[246 , 133,   52],[246 , 130,   51],[246 , 126,   49],[246 , 123,   48],[246 , 120,   47],[245 , 113,   44],\
	[245 , 106,   41],[244 , 104,   41],[243 , 102,   41],[242 , 100,   41],[241 ,  98,   41],[240 ,  96,   41],[239 ,  94,   41],\
	[239 ,  92,   41],[238 ,  90,   41],[237 ,  88,   41],[236 ,  86,   41],[235 ,  84,   41],[234 ,  82,   41],[233 ,  80,   41],\
	[232 ,  78,   41],[231 ,  76,   41],[230 ,  74,   41],[229 ,  72,   41],[228 ,  70,   41],[228 ,  67,   40],[227 ,  65,   40],\
	[226 ,  63,   40],[225 ,  61,   40],[224 ,  59,   40],[223 ,  57,   40],[222 ,  55,   40],[221 ,  53,   40],[220 ,  51,   40],\
	[219 ,  49,   40],[218 ,  47,   40],[217 ,  45,   40],[217 ,  43,   40],[216 ,  41,   40],[215 ,  39,   40],[214 ,  37,   40],\
	[213 ,  35,   40],[211 ,  31,   40],[209 ,  31,   40],[207 ,  30,   39],[206 ,  30,   39],[204 ,  30,   38],[202 ,  30,   38],\
	[200 ,  29,   38],[199 ,  29,   37],[197 ,  29,   37],[195 ,  29,   36],[193 ,  28,   36],[192 ,  28,   36],[190 ,  28,   35],\
	[188 ,  27,   35],[186 ,  27,   34],[185 ,  27,   34],[183 ,  27,   34],[181 ,  26,   33],[179 ,  26,   33],[178 ,  26,   32],\
	[176 ,  26,   32],[174 ,  25,   31],[172 ,  25,   31],[171 ,  25,   31],[169 ,  25,   30],[167 ,  24,   30],[165 ,  24,   29],\
	[164 ,  24,   29],[162 ,  23,   29],[160 ,  23,   28],[158 ,  23,   28],[157 ,  23,   27],[155 ,  22,   27],[153 ,  22,   27],\
	[151 ,  22,   26],[150 ,  22,   26],[146 ,  21,   25]])
	if reversed: C = C[::-1]
	colmap = mpl.colors.ListedColormap(C/255.0)
	return colmap

def loop1(k): 
	for i in range(noLats):
		v[k,i,:] = wledit(v[k,i,:],lats[i],lmin,lmax)
	return v[k,:,:]
	
def loop2(k):
	for j in range(noLons):
		for i in range(0,hwm): # latitudes close to the North Pole
			filt = sum(v[k,0:i+bm+1,j]*Hm[hwm-i:]) # the built in sum is faster than numpy's sum here
			v_latavg[k,i,j] = filt/sum(Hm[hwm-i:])			
		for i in range(hwm,noLats-hwm): 
			filt = sum(v[k,i-bm:i+bm+1,j]*Hm[:])   
			v_latavg[k,i,j] = filt/sum(Hm)   
		for i in range(noLats-hwm,noLats): # latitudes close to the South Pole
			sp = hwm-(noLats-i)+1   
			filt = sum(v[k,i-bm:,j]*Hm[:-sp])
			v_latavg[k,i,j] = filt/sum(Hm[:-sp])
	return v_latavg[k,:,:]
	
def loop3(k): 
	for i in range(noLats):
		env[k,i,:] = hilbert(v[k,i,:], noLons, envelope=True)
		env[k,i,:] = wledit_lowpass(env[k,i,:],lats[i],2*lmin)
	return env[k,:,:]
	
# def loop4(k): 
# 	for i in range(noLats):
# 		analsig = hilbert(v[k,i,:], noLons, envelope=False)
# 		hilb_phase[k,i,:] = arctan2(imag(analsig),real(analsig))
# 	return hilb_phase[k,:,:]

# def loop5(k):
# 	phase_label_nh_k = phase_label_nh[k,:,:]
# 	hilb_phase_nh_k = hilb_phase_nh[k,:,:]
# 	phase_label_nh_k[hilb_phase_nh_k > 0] = 1                               # NH ridges
# 	phase_label_nh_k[(hilb_phase_nh_k <= 0) & (hilb_phase_nh_k > -4)] = -1  # NH troughs
# 	# I now set grid points outside RWPs to "maskval"
# 	for i in range(0,eqlat):
# 		for j in range(-1,noLons-1):
# 			if not (env[k-1,i,j] >= env_thres and env[k,i,j] >= env_thres and env[k+1,i,j] >= env_thres and env[k,i,j-1] >= env_thres and env[k,i,j+1] >= env_thres): # if not all these are true, then phase index is maskval
# 				phase_label_nh_k[i,j] = maskval
# 	return phase_label_nh_k[:,:]

# def loop6(k):
# 	phase_label_sh_k = phase_label_sh[k,:,:]
# 	hilb_phase_sh_k = hilb_phase_sh[k,:,:]
# 	phase_label_sh_k[hilb_phase_sh_k >= 0] = -1                             # SH troughs
# 	phase_label_sh_k[(hilb_phase_sh_k < 0) & (hilb_phase_sh_k > -4)] = 1# SH ridges
# 	# I now set grid points outside RWPs to "maskval"
# 	for n,i in enumerate(range(eqlat,noLats)):
# 		for j in range(-1,noLons-1):
# 			if not (env[k-1,i,j] >= env_thres and env[k,i,j] >= env_thres and env[k+1,i,j] >= env_thres and env[k,i,j-1] >= env_thres and env[k,i,j+1] >= env_thres): # if not all these are true, then phase index is maskval
# 				phase_label_sh_k[n,j] = maskval
# 	return phase_label_sh_k[:,:]



###################################################
########### MAIN PROGRAM ##########################
###################################################

path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_era5_1981_2010/'
file = 'ERA5_V300_all.nc'
ds = xr.open_dataset(path+file).squeeze().drop_vars('Z')

for year in arange(1981,2011):
    try:
        ds_new = ds.copy(deep=True)
    except FileNotFoundError:
        continue
    print(year)
    v = ds.V.values
    start_time = time.time()
    cores      = cpu_count()-2  # How many CPU cores will be used for the parallel computations. cpu_count() gives the total number of (virtual) cores in the machine.
    
    ###################################################
    ######## SET UP ###################################
    ###################################################
    env_thres   = 15             # Compute RWP properties only where the RWP amplitude exceeds this threshold (in m/s)
    lmin        = 2000           # Minimum wavelength for zonal filtering (in km)
    lmax        = 10000          # Maximum wavelength for zonal filtering (in km)
    maskval     = nan            # Value given to grid points outside RWP objects (where Cp and Cg are not defined).
    save_files  = True
    
    # Select date (yy:year, mm:month, dd:day, tt:hour)
    # yy=year; mm=9; dd=22; tt=0
    
    # Load exemplary dataset: ERA5 Meridional wind at 300hPa for September 2018 (6-hourly, 2x2, anomaly from climatology) 
    # filein = Dataset('era51_mars_v300_6hourly_anom_from_smoothed04_clim_2018_09.nc', mode='r', format='NETCDF4')	
    # varname = 'v'; levname = 'level'; timname = 'time'
    
    # time_data  = filein.variables[timname][:] # Reading in the time values
    # time_units = filein.variables[timname].units    # Reading in the time units
    # cal_temps  = filein.variables[timname].calendar # Calendar in use (proleptic_gregorian)
    # presslev   = filein.variables[levname][:]
    lons  = ds['lon'].values
    lats  = ds['lat'].values
    noLats = len(lats) # number of latitudes
    noLons = len(lons) # number of longitudes
    res = abs(int(lons[1]-lons[0]))  # Horizontal resolution (in degrees) (assuming same resolution in lat and lon; here, 2 degrees)
    # timeres = abs(time_data[1]-time_data[0])  # Temporal resolution (in hours)
    timeres = 6
    eqlat = int((noLats-1)/2)  # Equator latitude index
    
    # mydate = dt.datetime(yy,mm,dd,tt)
    # time_value = date2num(mydate,units=time_units,calendar=cal_temps)
    # t_list = list(time_data)               
    # tim = t_list.index(time_value)
    # tt1 = tim-1 # Compute 3 subsequent timesteps (required for phase/group speed)
    # tt2 = tim+2
    
    # Read meridional wind data
    # v = filein.variables[varname][tt1:tt2,0,:,:]
    noTimes = v.shape[0] # number of time steps
    calc_times1 = range(noTimes) # how many snapshots do I compute
    calc_times2 = range(1,noTimes-1)
    
    hann_width = int(15/res)    # Hann window width (in grid points) for the meridional filtering. The higher this width the stronger the smoothing.
    hthres     = int(20/res)    # parameter L0 in FW20 (20 degrees longitude; 10 grid points in the case of this 2x2 dataset)
    hthres_y   = int(10/res)    # parameter N0 in FW20 (10 degrees latitude; 5 grid points in the case of this 2x2 dataset)
    zwin       = hthres         # window extension for finding the h-length sum(El) maximum
    zwin_y     = hthres_y       # window extension for finding the h-length sum(En) maximum
    
    cc2 = 1/(2*timeres*60*60)
    latid = zeros(noLats)
    deg_len = zeros(noLats)
    cc6 = zeros(noLats)
    for j in range(noLats):
    	latid[j] = 90-res*j
    	deg_len[j] = circumference(latid[j])/360 # the length of one degree longitude at a given latitude circle
    	cc6[j] = 1/(2*res*deg_len[j])
    deg_len_y = circumference(0)/360   # the length of one degree latitude
    cc6y = 1/(2*res*deg_len_y)
    
    ###################################################
    ######## COMPUTATIONS #############################
    ###################################################
    
    ###########################################################################################
    # Spatial (zonal and meridional) filtering. The zonal filtering sometimes introduces discontinuities in the meridional direction, so here we apply a weak smoothing to reduce such instances.
    print('Zonal filtering')
    # with Pool(cores) as p: v[0:noTimes,:,:] = p.map(loop1, calc_times1)
    with Pool(cores) as p: v[:,:,:] = p.map(loop1, calc_times1)
    
    print('Meridional filtering')
    xm = arange(0, hann_width, 1)
    bm = int((len(xm)-1)/2)
    Hm = 0.5-0.5*cos(2*pi*xm/(hann_width-1)) # Hann window function 
    hwm = int(hann_width/2)
    v_latavg = zeros_like(v)
    # with Pool(cores) as p: v[0:noTimes,:,:] = p.map(loop2, calc_times1)
    with Pool(cores) as p: v[:,:,:] = p.map(loop2, calc_times1)
    del v_latavg
    
    print('Envelope')
    env = zeros_like(v)
    # with Pool(cores) as p: env[0:noTimes,:,:] = p.map(loop3, calc_times1)
    with Pool(cores) as p: env[:,:,:] = p.map(loop3, calc_times1)
    
    ds_new.V[:,:,:]=env[:,:,:]
    ds_new.to_netcdf(path+'RWP_E_{:04d}.nc'.format(year))

# print('Phase')
# hilb_phase = zeros_like(v)
# with Pool(cores) as p: hilb_phase[0:noTimes,:,:] = p.map(loop4, calc_times1)

# # Mask possible unphysical values
# hilb_phase[(hilb_phase<-pi) | (hilb_phase>pi)] = maskval

# print('Phase labelling')
# hilb_phase_nh = hilb_phase[:,:eqlat,:]
# hilb_phase_sh = hilb_phase[:,eqlat:,:] # this one includes the equator as well. If a RWP passes over the equator (rare), then what is a trough in NH (southerlies downstream of northerlies) will be a ridge in the SH, so the label will abruptly change from -1 in NH to 1 in SH.
# phase_label_nh = zeros_like(hilb_phase_nh)
# phase_label_sh = zeros_like(hilb_phase_sh)
# with Pool(cores) as p: phase_label_nh[1:noTimes-1,:,:] = p.map(loop5, calc_times2)
# with Pool(cores) as p: phase_label_sh[1:noTimes-1,:,:] = p.map(loop6, calc_times2)

# ################
# print("----- Computations lasted %s minutes -----" % ((time.time() - start_time)/60))
# ################


# #### Plot
# levels_v = arange(-40,41,2.5); ticks_v = levels_v[::4]
# levels_env = arange(15,40.5,2.5); ticks_env = levels_env[::2]
# levels_phase = linspace(-pi,pi,81); ticks_phase = levels_phase[::10]; ticklabels_phase = ['-π','-3π/4','-π/2','-π/4','0','+π/4','+π/2','+3π/4','+π']
# levels_phase_label = [-1,0,1]; ticks_phase_label = [-0.5,0.5]; ticklabels_phase_label = ['Trough (-1)','Ridge (+1)']

# my_cmap0 = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue',colors =['darkblue','royalblue','white','indianred','darkred'],N=len(levels_v))
# my_cmap1 = mcolors.LinearSegmentedColormap.from_list(name='bwr',colors =['darkblue','royalblue','white','white','indianred','darkred'],N=80)
# my_cmap2 = plt.cm.get_cmap('ocean_r')
# my_cmap4 = mcolors.ListedColormap(['#5a86ad', 'goldenrod'])
# my_cmap5 = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue',colors =['darkgreen','lightsteelblue','#5a86ad','lightsteelblue','rebeccapurple','gold','darkgoldenrod','gold','darkgreen'],N=80)
# props1 = dict(boxstyle='round', facecolor='gold', alpha=1)

# # Fields to plot:
# v_raw_plot = filein.variables['v'][tim,0,:,:]
# v_filt_plot = v[1,:,:]
# env_plot = env[1,:,:]
# hilb_phase_plot = hilb_phase_nh[1,:,:]
# phase_label_plot = phase_label_nh[1,:,:]

# # Shift the data and longitudes according to the chosen projection. 
# eastern_bound = 90 # Eastern map longitude boundary
# v_raw_plot, lons1 = shiftgrid(eastern_bound, v_raw_plot, lons, start=False) 
# v_filt_plot, lons1 = shiftgrid(eastern_bound, v_filt_plot, lons, start=False) 
# env_plot, lons1 = shiftgrid(eastern_bound, env_plot, lons, start=False) 
# hilb_phase_plot, lons1 = shiftgrid(eastern_bound, hilb_phase_plot, lons, start=False)
# phase_label_plot, lons1 = shiftgrid(eastern_bound, phase_label_plot, lons, start=False) 

# # Concatenate first longitude to the end, such that no gap appears in the map
# v_raw_plot, lons2 = addcyclic(v_raw_plot, lons1)
# v_filt_plot, lons2 = addcyclic(v_filt_plot, lons1)
# env_plot, lons2 = addcyclic(env_plot, lons1)
# hilb_phase_plot, lons2 = addcyclic(hilb_phase_plot, lons1)
# phase_label_plot, lons2 = addcyclic(phase_label_plot, lons1)

# LON, LAT = meshgrid(lons2, lats)
# LONnh, LATnh = meshgrid(lons2, lats[:eqlat])

# def mybasemap(myax,var,var_levels,var_ticks,mycmap,cbar_label,cbar_ticklabels=None,phase=False):
# 	m1 = Basemap(projection='cyl',llcrnrlon=eastern_bound-360,llcrnrlat=15, urcrnrlon=eastern_bound, urcrnrlat=85, resolution='c')
# 	if phase: 
# 		cs = m1.pcolormesh(LONnh, LATnh, var, shading='nearest', latlon='true', cmap=mycmap, norm=mcolors.BoundaryNorm(var_levels, ncolors=mycmap.N, clip=True))
# 	else: 
# 		cs = m1.contourf(LON, LAT, var, var_levels, latlon='true', extend='both', cmap=mycmap)
# 	m1.drawparallels((0,20,40,60,80), labels=[1,0,0,0], linewidth=0.3, fontsize=7, textcolor='gray')
# 	m1.drawmeridians((0,60,120,180,240,300), labels=[0,0,1,0], linewidth=0.3, fontsize=7, textcolor='gray')
# 	m1.drawcoastlines(linewidth=0.7,color='gray')
# 	cbar = m1.colorbar(cs, location='bottom', size="12%", pad="20%", ticks=var_ticks)
# 	cbar.set_label(cbar_label, fontsize=14)
# 	if cbar_ticklabels==None: cbar_ticklabels=var_ticks
# 	cbar.ax.set_xticklabels(cbar_ticklabels)
# 	cbar.ax.tick_params(labelsize=12)  

# fig = plt.figure(1, figsize=(11,15))
# plt.suptitle('%s/%s/%s - %s0UTC' %(dd,mm,yy,tt),x=0.52,y=0.92, fontsize=15, weight='bold')
# ax1 = plt.subplot(5,1,1)
# mybasemap(ax1,v_raw_plot,levels_v,ticks_v,my_cmap1,r"Meridional wind anomaly (v') at 300 hPa [$ms^{-1}$]")
# ax1.text(0.01, 0.95, "raw", transform=ax1.transAxes, fontsize=15, verticalalignment='top', bbox=props1)
# ax2 = plt.subplot(5,1,2)
# mybasemap(ax2,v_filt_plot,levels_v,ticks_v,my_cmap1,r"Meridional wind anomaly (v') at 300 hPa [$ms^{-1}$]")
# ax2.text(0.01, 0.95, "filtered", transform=ax2.transAxes, fontsize=15, verticalalignment='top', bbox=props1)
# ax3 = plt.subplot(5,1,3)
# mybasemap(ax3,hilb_phase_plot,levels_phase,ticks_phase,my_cmap5,"Phase (Φ) at 300 hPa (rad)",cbar_ticklabels=ticklabels_phase,phase=True)
# ax4 = plt.subplot(5,1,4)
# mybasemap(ax4,phase_label_plot,levels_phase_label,ticks_phase_label,my_cmap4,"Phase index at 300 hPa",cbar_ticklabels=ticklabels_phase_label,phase=True)
# ax5 = plt.subplot(5,1,5)
# mybasemap(ax5,env_plot,levels_env,ticks_env,my_cmap2,r"Amplitude (E) at 300 hPa [$ms^{-1}$]")
# plt.subplots_adjust(hspace=0.5)
# plt.savefig('rwp_E-phase_%d-%d-%d.jpg' %(dd,mm,yy), bbox_inches='tight', dpi=400)
#filein.close()

################
print("----- Plotting finished -----")
################

# %%
