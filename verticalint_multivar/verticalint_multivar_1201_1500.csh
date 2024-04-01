#!/bin/bash 

#set yy = 3380
#set eyy = 3514
#set eyy = 3381



cd ../../data/extracted


#while ( $yy < $eyy )


#if ( $yy < 1000 ) then
#set yy = "0"$yy 
#endif
#echo $yy
for yy in {1201..1500};do
echo ${yy#0}
if [[ ${yy#0} -le 1246 ]]; then
	datpath="/storage/climatestor/PleioCEP/kim/cesm1.2_data/archive/TRANS.3501BP/atm/hist"
elif [[ ${yy#0} -gt 1246 ]] && [[ ${yy#0} -le 2699 ]]; then
	datpath="/storage/climatestor/PleioCEP/kim/cesm1.2_data/archive/TRANS.3501BP_1/atm/hist"
elif [[ ${yy#0} -gt 2699 ]] && [[ ${yy#0} -le 3351 ]]; then
	datpath="/storage/climatestor/PleioCEP/kim/cesm1.2_data/archive/TRANS.3501BP_2/atm/hist"
elif [[ ${yy#0} -gt 3352 ]]; then
	datpath="/storage/climatestor/PleioCEP/kim/cesm1.2_data/archive/TRANS.3501BP_3/atm/hist"
fi

/bin/rm Interpolation_Z_hsigma_ps_$yy.ncl


cat > Interpolation_Z_hsigma_ps_$yy.ncl <<EOF
;*********************************************
; vert_2_500.ncl / This script interpolates sigma to pressure
; coordinates in outputs the geopotential height (for all plevs chosen)
;*********************************************
load "/usr/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "/usr/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;*********************************************
begin
;*************************************************
; read in data
;*************************************************

;******only years
;******year has to be given with the function call: 'for a in { }; do ncl
;******year=... Interpolation_... ;done'

path1="$datpath/TRANS.3501BP.cam.h1.$yy-01-01-00000.nc"
;path1="/storage/climatestor/PleioCEP/doensen/data/extracted/TRANS.3501BP.cam.h1.????-01-01-00000_sel.nc

in1 = addfile(path1,"r")
intVar_geo = in1->Z3 ;*****intVar is the Variable which is going to be interpolated
psl = in1->PSL
precc = in1->PRECC
precl = in1->PRECL
intVar_u = in1->U
intVar_v = in1->V
;in2 = addfile("PS-00.nc","r")
ps = in1->PS

;in3 = addfile(path1,"r")
hbcofa = in1->hyam
hbcofb = in1->hybm

;***************************************************
; interpolate to pressure levels
;***************************************************

;----------------------------------------------------------------------
;   define other arguments required by vinth2p
;----------------------------------------------------------------------

;   type of interpolation: 1 = linear, 2 = log, 3 = loglog
     interp = 2

;   is extrapolation desired if data is outside the range of PS
;    extrap = False
     extrap = True
plevs =(/ 1000.0 /)

plevs!0     = "plevs"
plevs&plevs =  plevs
plevs@long_name = "Pressure"
plevs@unit = "hPa"

p0 = 1000 ;set reference pressure
ii = 1 ;unused constant but needs to be provided

intVar_PS = vinth2p(intVar_geo,hbcofa,hbcofb,plevs,ps,interp,p0,ii,extrap)
intVar_PS = intVar_PS(:,:,::-1,:)
plevs = plevs*100
intVar_PS&plevs = plevs
intVar_PS&plevs@unit = "Pa"

plevs =(/ 850.0 /)

plevs!0     = "plevs"
plevs&plevs =  plevs
plevs@long_name = "Pressure"
plevs@unit = "hPa"

p0 = 1000 ;set reference pressure
ii = 1 ;unused constant but needs to be provided

intVar_U = vinth2p(intVar_u,hbcofa,hbcofb,plevs,ps,interp,p0,ii,extrap)
intVar_V = vinth2p(intVar_v,hbcofa,hbcofb,plevs,ps,interp,p0,ii,extrap)
intVar_WS = wind_speed(intVar_U,intVar_V)
intVar_WS = intVar_WS(:,:,::-1,:)
plevs = plevs*100
intVar_WS&plevs = plevs
intVar_WS&plevs@unit = "Pa"
;system("echo saving")

;setfileoption("nc","Format","NetCDF4Classic")
;fileout=getenv("FZ")

fileout="TRANS.3501BP.cam.h1.$yy-01-01.sel.nc"
;fileout="test.nc"
fout=addfile(fileout,"c")

;fileout="test_int.nc"
;fout=addfile(fileout,"c")

precc = precc(:,::-1,:)
precl = precl(:,::-1,:)
psl = psl(:,::-1,:)
prect = precc + precl

fout->Z3 = intVar_PS      ; write into new file 
fout->PSL = psl
fout->PRECT = prect
fout->PRECC = precc
fout->PRECL = precl
fout->WS = intVar_WS
fout = rm_single_dims(fout)
system("echo new file for GPH")        ; print path and new file to screen as confirmation

end
EOF

ncl < Interpolation_Z_hsigma_ps_$yy.ncl

done
exit
