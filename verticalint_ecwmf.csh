#!/bin/csh




set yy = 3481
#set eyy = 3514
set eyy = 3482



cd ../data/extracted


while ( $yy < $eyy )


if ( $yy < 1000 ) then
set yy = "0"$yy 
endif
echo $yy



/bin/rm Interpolation_Z_hsigma_ps.ncl


cat > Interpolation_Z_hsigma_ps.ncl <<EOF
;*********************************************
; vert_2_500.ncl / This script interpolates sigma to pressure coordinates in outputs the geopotential height (for all plevs chosen)
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
;******year has to be given with the function call: 'for a in {  }; do ncl
;******year=... Interpolation_... ;done'

path1="/storage/climatestor/PleioCEP/kim/cesm1.2_data/archive/TRANS.3501BP_3/atm/hist/TRANS.3501BP.cam.h1.$yy-01-01-00000.nc"

in = addfile(path1,"r")


;----------------------------------------------------------------------
 ;   read needed variables from file
 ;----------------------------------------------------------------------

 Z3   = in->Z3      ; select variable to be converted
 P0mb = 1000.
 hyam = in->hyam    ; get a coefficiants
 hybm = in->hybm    ; get b coefficiants
 PS   = in->PS      ; get pressure
 TBOT = in->TBOT    ; get temperature at lowest layer (closest to surface)
 dims = dimsizes(Z3)
 nlevs= dims(1)
 PHIS = Z3(:,nlevs-1,:,:)*9.81   ; get geopotential [m^2/s^2] at the bottom (lowest layer)

 ;----------------------------------------------------------------------
 ;   define other arguments required by vinth2p
 ;----------------------------------------------------------------------

 ;   type of interpolation: 1 = linear, 2 = log, 3 = loglog
 interp = 2

 ;   is extrapolation desired if data is outside the range of PS
 ;   extrap = False
 extrap = True

;  A scalar integer indicating which variable to interpolate: 1 = temperature, -1 = geopotential height, 0 = all others. 

varflg = -1

;   create an array of desired pressure levels:
plevs =(/ 1000.0 /)
			
plevs!0     = "plevs"
plevs&plevs =  plevs
plevs@long_name = "Pressure"
plevs@unit = "hPa"
  
;intVar_PS = vinth2p(intVar,hbcofa,hbcofb,plevs,ps,interp,1000,1,extrap)
;plevs = plevs*100
;intVar_PS&plevs = plevs
;intVar_PS&plevs@unit = "Pa"

intVar_PS = vinth2p_ecmwf(Z3,hyam,hybm,plevs,PS,interp,P0mb,1,extrap,varflg,TBOT,PHIS) 

intVar_PS!0 = "time"
intVar_PS!1 = "lev"
intVar_PS!2 = "lat"
intVar_PS!3 = "lon"
intVar_PS&time = in->time
intVar_PS&lev  = plevs
intVar_PS&lat  = in->lat
intVar_PS&lon  = in->lon
intVar_PS@units     = "m" 
intVar_PS@long_name = "Geopotential Height (above sea level)"

;system("echo saving")	
;setfileoption("nc","Format","NetCDF4Classic")
;fileout=getenv("FZ")

;fileout="TRANS.3501BP.cam.h1.$yy-01-01.z10000_ecmwf.nc"
fileout="test_false.nc"

system("rm " + fileout)        ; remove any pre-existing file


fout=addfile(fileout,"c")

fout->Z3 = intVar_PS(:,:,::-1,:) 	 ; write into new file
system("echo new file for GPH")  ; print path and new file to screen as confirmation

end
EOF



ncl < Interpolation_Z_hsigma_ps.ncl

#mv test2.nc BPRD_trans.cam2.h1.$yy"-"01-01.z1000.nc
#cdo -f srv copy  test2.nc BPRD_trans.cam2.h1.$yy"-"01-01.z1000.srv


set yy = `/usr/bin/expr $yy + 1`
#set zeiger = `/usr/bin/expr $zeiger + 100000`
#echo $yy $zeiger
end


exit
