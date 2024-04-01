#!/bin/csh

set yy = 3380
#set eyy = 3514
set eyy = 3381



cd ../data/extracted


while ( $yy < $eyy )


if ( $yy < 1000 ) then
set yy = "0"$yy 
endif
echo $yy
/bin/rm Interpolation_Z_hsigma_ps.ncl


cat > Interpolation_Z_hsigma_ps.ncl <<EOF
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

path1="/storage/climatestor/PleioCEP/kim/cesm1.2_data/archive/TRANS.3501BP_3/atm/hist/TRANS.3501BP.cam.h1.$yy-01-01-00000.nc"
;path1="/storage/climatestor/PleioCEP/doensen/data/extracted/TRANS.3501BP.cam.h1.????-01-01-00000_sel.nc

in1 = addfile(path1,"r")
intVar = in1->Z3 ;*****intVar is the Variable which is going to be interpolated
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

intVar_PS = vinth2p(intVar,hbcofa,hbcofb,plevs,ps,interp,1000,1,extrap)
intVar_PS = intVar_PS(:,:,::-1,:)
plevs = plevs*100
intVar_PS&plevs = plevs
intVar_PS&plevs@unit = "Pa"

;system("echo saving")

;setfileoption("nc","Format","NetCDF4Classic")
;fileout=getenv("FZ")

fileout="TRANS.3501BP.cam.h1.$yy-01-01.z1000.nc"
fout=addfile(fileout,"c")

;fileout="test_int.nc"
;fout=addfile(fileout,"c")

fout->Z3 = intVar_PS      ; write into new file 
system("echo new file for GPH")        ; print path and new file to screen as confirmation

end
EOF

ncl < Interpolation_Z_hsigma_ps.ncl

set yy = `/usr/bin/expr $yy + 1`
end

exit
