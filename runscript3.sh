#!/bin/csh

mkdir /storage/climatestor/PleioCEP/doensen/data/cyclone_1205_1804

cd /storage/climatestor/PleioCEP/doensen/data/cyclone_1205_1804

#/bin/rm *
set y = 1575
#set eyy = 859
set eyy = 1584

set datpath = "/storage/climatestor/PleioCEP/doensen/data/extracted/"



cdo -f srv copy /storage/climatestor/PleioCEP/doensen/data/extracted/oro.nc oro.srv




while ( $y < $eyy )

set zeiger = `/usr/bin/expr $y - 5`
set zeiger = `/usr/bin/expr $zeiger / 10`
set zeiger = `/usr/bin/expr $zeiger \* 300000`
echo $zeiger

set arr = (0 1 2 3 4 5 6 7 8 9)
set i = 1
set ei = 11
while ($i < $ei)
set yy = `/usr/bin/expr $y + $i - 1`
set yy = `/usr/bin/printf "%04d" $yy`
set arr[$i]=$yy
echo ${arr[$i]}
set i = `/usr/bin/expr $i + 1`
end


#echo $i
#cdo ymonsub /alphadata04/blumer/no_Backup/Data/Postprocessed/BPRD_trans_1/Z1000/DJF/Z1000gp_DJF_${yy}_`expr $yy + 9`.srv /data0cd2/raible/sandro/mean/djfmean_1980_2009.srv tmp.srv
/bin/rm xatrs_CCSM_0.9x1.25.f
/bin/rm extract_slp_prec.f90
/bin/rm xatrs_CCSM_0.9x1.25.exe
/bin/rm extract_slp_prec.exe
#/bin/rm tmpslp.nc tmpprec.nc
#/bin/rm tmpslp2.nc tmppre2c.nc tmpz1000.srv



#if ($y == 1375) then
#cdo -O mergetime $datpath/TRANS.3501BP.cam.h1.${arr[1]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[2]}-01-01.sel_alt.nc  $datpath/TRANS.3501BP.cam.h1.${arr[3]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[4]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[5]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[7]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[8]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[9]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[10]}-01-01.sel_alt.nc TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"01-01.sel.nc
#
#elif ($y == 1575) then
#cdo -O mergetime $datpath/TRANS.3501BP.cam.h1.${arr[1]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[3]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[4]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[5]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[6]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[7]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[8]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[9]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[10]}-01-01.sel_alt.nc TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"01-01.sel.nc
#
#else
#cdo -O mergetime $datpath/TRANS.3501BP.cam.h1.${arr[1]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[2]}-01-01.sel_alt.nc  $datpath/TRANS.3501BP.cam.h1.${arr[3]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[4]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[5]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[6]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[7]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[8]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[9]}-01-01.sel_alt.nc $datpath/TRANS.3501BP.cam.h1.${arr[10]}-01-01.sel_alt.nc TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"01-01.sel.nc
#fi 
echo $y


#if ($y == 1375) then
#cdo -O mergetime $datpath/TRANS.3501BP.cam.h0.${arr[1]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[2]}-Qtot.nc  $datpath/TRANS.3501BP.cam.h0.${arr[3]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[4]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[5]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[7]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[8]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[9]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[10]}-Qtot.nc TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"Qtot.nc
#
#elif ($y == 1585) then
#cdo -O mergetime $datpath/TRANS.3501BP.cam.h0.${arr[1]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[3]}-Qtot.nc  $datpath/TRANS.3501BP.cam.h0.${arr[4]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[5]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[6]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[7]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[8]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[9]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[10]}-Qtot.nc TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"Qtot.nc
#
#
#else
#cdo -O mergetime $datpath/TRANS.3501BP.cam.h0.${arr[1]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[2]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[3]}-Qtot.nc  $datpath/TRANS.3501BP.cam.h0.${arr[4]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[5]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[6]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[7]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[8]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[9]}-Qtot.nc $datpath/TRANS.3501BP.cam.h0.${arr[10]}-Qtot.nc TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"Qtot.nc
#
echo "${arr[1]}"_"${arr[10]}"
cdo -f srv copy -selvar,Z3 TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"01-01.sel.nc  tmpz1000.srv
cdo -f srv copy -selvar,PRECT TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"01-01.sel.nc  PRECT.srv
cdo -f srv copy -selvar,PSL TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"01-01.sel.nc PSL.srv
cdo -f srv copy -selvar,WS TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"01-01.sel.nc WS.srv
cdo -f srv copy -selvar,PRECL TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"01-01.sel.nc PRECL.srv
cdo -f srv copy -selvar,PRECC TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"01-01.sel.nc PRECC.srv
cdo -f srv copy -invertlat -selvar,Q TRANS.3501BP.cam.h1.${arr[1]}"_"${arr[10]}"_"Qtot.nc Q.srv

cat > xatrs_CCSM_0.9x1.25.f<<EOF
!*********************************************************************
!          analysis of cyclones and anticyclones
!                  file xatrs.f
!*********************************************************************
!                Richard Blender
!     Meteorologisches Institut, Universitaet Hamburg
!         Bundesstr. 55, D-20146 Hamburg, FRG
!       phone 49.40.4123-5070, fax 49.40.4123-5066
!            email: blender at dkrz.d400.dedo
!                   May 1996
!           f90 version and radius Jan 2007
!----------------------------------------------------------------------
!
!                brief description
!
!     this program analyses a sequence of geopotential or vorticity
!     fields to detect low/high pressure systems or vorticity extrema
!     and to connect their paths. also tropical storms.
!
!----------------------------------------------------------------------
!   input files
!   -----------
!     unit               content
!     ----------------------------------------------------------
!     chzp   22          geopot at e.g.1000hpa, zp, or pressure
!     chvort 23          vorticity at e.g.850hpa, vort
!
!     chtrop 24          tropical basic unit
!            24  itrop=1 wind 10m                  code 171
!            +1  itrop=2 upper level vorticity     code 138
!            +2  itrop=3 temperature               code 130
!
!     choro  29          orography, oro (only of onoro=1)
!
!   output files
!   ------------
!     unit               content
!     -----------------------------------------------------------
!     chpath 30          paths, complete, all vars detected
!     chpgrd 31          grid    (plot with routine of
!     chpdgr 32          degr     f. sielmann)
!     chpar2 33          r2 vs age, log/log, and slope
!     chpand 35          nd (distrib in diica)
!     chppca 36          pca, density p(r,t) r=arc, t=age
!     chxy   37          abs paths seperated by 0 0
!     chserv 38          first data field in service format, if other
!     chcc   39          cat in character form (global)
!     chpalt 40          distribution of life times
!     chgen  41,2        density genesis, in one column
!     chlys  43,4          "       lysis
!     chtot  45,6          "       total
!     chptro 50          tropical data for paths
!     chpos  60          cyclone positions
!------------------------------------------------------------------
!  variables (in all subroutines, in alphabetical order)
!
!     age             age of cyclone in steps
!     ageic(ic)       actual age of cyclone number ic
!     agemh           min age in hours
!     agemin          min age for real cyclones in steps
!     agemx           max age to be stored (in steps, parameter!)
!     anatyp          method to determine cyclones 0=geop,1=vort, 8=trop
!     cat(i,j)        category of point (i,j), cyclone in state 1,2,..
!     cat0,cat1       category of point (i,j), cyclone in state 1,2,..
!     catica(ic,age)  category of cyclone ic at age age
!     ch...           unit numbers, see subroutine inivar
!     czpica(ic,age)  central height in a cyclone
!     datica(ic,age)  date
!     timeica(ic,age) time
!     dattyp          type of data file, =0 service, =1 other
!     delola(i,j)     cyclone detected in i,j
!     diica(ic,age)   last distance travelled by a cyclone
!     dist_fit        range of fit are in fitgauss (radius)
!     dt              time interval in hours between fields
!     dx, dy          grid distance (in 1000km), dx depends on phi
!     gzp             mean gradient of ze around a minimum
!     gzpica(ic,age)  grad z of cyc ic at age
!     gzpmax          max geop/press grad to detect real (intense) vortices
!     gzpmin          min geop/press grad to detect vortices
!     i,j,i1,j1       longitude/latitude  loop variables
!     i0,j0           old positon of cyclone
!     iclola(i,j)     lagragian coord of cyclone at i,j at time t
!     id, idmax       integer number distance, max
!     ih              headers in service format input files
!     itrop           loop variable for tropical fields
!     j               latitiude loop variable
!     jump            0 or 1, =1 if jump in input file date occured
!     lam1, lam2      geographical longitude if search region
!     latica(ic,age)  latitude  of cyclone ic at age age
!     lonica(ic,age)  longitude of cyclone ic at age age
!     lows            +1 or -1, lows or highs are determined
!     misv_fit        fit missing value, for case of no result
!     nc              number of cyclones found
!     ncmx            max. number of cyclones (parameter!)
!     N_corr_xxx      number of corrections in fit
!     newdat          date of new input file
!     nla1,nla2       latitude  grid points of search region
!     nlo,nla         longitude and latitude grid width as found in data
!     nlo1,nlo2       longitude grid points of search region
!     nlomx,nlamx     max longitude and latitude grid points
!     nrealc          number of real (intense and old enough) cyclones
!     nsteps          number of input data time steps
!     nt
!     ntrop           number of variables used for tropical storms
!     olddat          date of last input file
!     onoro           0 or 1, orography used
!     oro(i,j)        orography field, used to exclude e.g.>1000m
!     p(i)            periodic longitude (e.g. i=0 -> max)
!     pca(id,age)     number of cycs with distance id at age
!     phi             latitude angle
!     phi1, phi2      geographical latitudes of search region
!     pp, pp_fit      fit parameters in powell (Numer Recip)
!     r0              earth radius
!     r2(age)         mean square displacement
!     r2lon(age), r2lat(age) in long and latitude
!     realic(ic)      tells if ic is an intense cyclone (life time etc)
!     sr,srmax        search range, max
!     step            number of input data time steps, max
!     t0ic(ic)        time step of first appearance of cyclone ic
!     test            0 or 1, if=1 some additional tests are performed
!     time            time in service format header (in hours)
!     tropica         individual (track) data (ntrop,ic,age)
!     tropf           (ntrop,lon,nla) additional fields for tropical storms
!     vort            vorticity field, input data
!     vortmax         max vorticity to detect real (intense) vortices
!     vortmin         min vorticity to detect vortices
!     xxx_fit         copies of the variables 'xxx' used in fit
!     zrad            radius (std dev of Gaussian)
!     zdep            depth (difference zcenter and zenv)
!     zenv            geopot height of far off the cyclone
!     zenv_max_fit    max value for fit of zenv
!     zmean           zp (geopot) in area around cyclone
!     zp              geopotential or pressure field, zp(i,j)
!     zp_max_fit      Gaussian fit only if zp lower than this
!     zrad_max_fit    max value for fit of zrad
!
!  subroutines and functions
!
!     ancat0   analyzes first category field in a sequence
!     ancat1   analyzes subsequent cat fields and connects tracks
!     chacat   output of cats in character format, easy to read
!     chkcat   check category field (mainly double neighbours)
!     datdmy   yields day,month,year from date in service format
!     datjmp   determines whether there is a jump in the dates
!     ddica    distribution of distances
!     detcat   classifies lows/high using numbers cat=1..9 (categories)
!     detrcy   determines whether the tracks found are intense and long
!     dsmall   distance for small angles
!     fitgauss fit of zrad,zdep (new 2007), only for geopot/pressure
!     findpars used in fitgauss, calls powell (numer recip) (new 2007)
!     func     distance for optimisation in fitgauss/powell/findpars
!     genlys   output of genesis, lysis and total density fields
!     gradzp   gradient of field around a point i,j
!     ijmin1   checks whether field has minimum in i,j
!     imin     integer minimum
!     indat    contrals data input (using insvf or inothf)
!     inidep   init variables which depend on input
!     inivar   init variables, e.g. units (units)
!     inothf   input of other format
!     inpvar   read initial values of variables
!     insr     determines whether point i,j is in the search region
!     insvf    input of service format files
!     outlt    output of lifetimes
!     outsvf   output of data in service format
!     outtr    output of tracks, complete info and for F. Sielmanns  plots
!     prhead   print header of program, output of parameters
!     search   search cyclone in preceding cat field
!     shiftc   shift cat field
!     statr2   output of mean square displacements
!
!----------------------------------------------------------------------

      module xatrm

      implicit none
      
!     change nlomx, e.g.64 for T21

!     service format resolution:
!     T21 64x32, T42 128x64, T63 96x96, T106 320x160

!     ntrop number of variables used for tropical storms
      
!s2.9.15      integer,parameter       :: nlomx=96, nlamx=nlomx/2
      integer,parameter       :: nlomx=144, nlamx=96
!s2.9.15      integer,parameter       :: ncmx=50000
      integer,parameter       :: ncmx=200000
!s2.9.15      integer,parameter       :: agemx=40
      integer,parameter       :: agemx=180
      integer,parameter       :: ntrop=1 ! 1 means not used, max ntrop=3

      real(kind=4),parameter  :: pi=3.1415926
      real(kind=4),parameter  :: r0=6.370 ! x 1000km earth radius

      integer                 :: test=0      
      integer                 :: p(-nlomx:2*nlomx)      ! periodic in i
      
      ! ------------------------
      
      ! NEW 2007 fit radius and depth, variables
      
      logical,parameter       :: do_fit    =.true. ! fit radius and depth
      
      logical,parameter       :: do_fit_out=.false. ! output addit. fit data
      
      integer                 :: i_fit,j_fit
      integer                 :: nlo_fit,nla_fit
      
      integer                 :: N_corr_zenv=0  ! corrected in fitgauss
      integer                 :: N_corr_zdep=0
      integer                 :: N_corr_zrad=0
      integer                 :: N_call_fitgauss=0
      
      real(kind=4)            :: zp_fit(nlomx,nlamx)
      
      real(kind=4)            :: dist_fit=1.0     ! x 1000km fit range
      
      real(kind=4),parameter  :: misv_fit=0.999   ! fit not successful
      real(kind=4),parameter  :: zp_max_fit=4000  ! max zp for fit  
      
      real(kind=4),parameter  :: zenv_max_fit=400  ! max zenv for fit 
      real(kind=4),parameter  :: zrad_max_fit=1.0  ! max radius for fit
      
      real(kind=4),parameter  :: f_rad_out=5000    ! factor zrad output
      
      ! ----------------------------------------------------------------
      
      real(kind=4)            :: phi1,phi2,lam1,lam2  ! region
      
      integer                 :: nlo,nla  ! nlo,nla defined by insvf
      integer                 :: nlo1,nlo2,nla1,nla2
      
      ! input and output units
      
      integer,parameter       :: chzp    = 22
      integer,parameter       :: chvort  = 23
      integer,parameter       :: chtrop  = 24
      integer,parameter       :: choro   = 29

      integer,parameter       :: chpath  = 30
      integer,parameter       :: chpgrd  = 31
      integer,parameter       :: chpdgr  = 32
      integer,parameter       :: chpar2  = 33
      integer,parameter       :: chpand  = 35
      integer,parameter       :: chppca  = 36
      integer,parameter       :: chxy    = 37
      integer,parameter       :: chserv  = 38
      integer,parameter       :: chcc    = 39
      integer,parameter       :: chpalt  = 40
      integer,parameter       :: chgen   = 41
      integer,parameter       :: chlys   = 43
      integer,parameter       :: chtot   = 45
      integer,parameter       :: chptro  = 50
      integer,parameter       :: chpos   = 60
      
      end module xatrm

!----------------------------------------------------------------------

      program xatrs

      use xatrm
      
      implicit none

      real(kind=4) ::    oro  (nlomx,nlamx)
      real(kind=4) ::    zp   (nlomx,nlamx)
      real(kind=4) ::    vort (nlomx,nlamx)
      real(kind=4) ::    tropf(ntrop,nlomx,nlamx)

      real(kind=4) ::    czpica(ncmx,0:agemx)
      real(kind=4) ::    gzpica(ncmx,0:agemx)
      real(kind=4) ::    zradica(ncmx,0:agemx)      
      real(kind=4) ::    zdepica(ncmx,0:agemx)      
      real(kind=4) ::    tropica(ntrop,ncmx,0:agemx)
      real(kind=4) ::    diica (ncmx,0:agemx)

      real(kind=4) ::    vortmin,vortmax,gzpmin,gzpmax
      real(kind=4) ::    wmtrop,vmtrop

      integer   :: cat0(nlomx,nlamx)
      integer   :: cat1(nlomx,nlamx)
      integer   :: lonica(ncmx,0:agemx)
      integer   :: latica(ncmx,0:agemx)
      integer   :: datica(ncmx,0:agemx)
      integer   :: timeica(ncmx,0:agemx)
      integer   :: iclola(nlomx,nlamx)
      integer   :: delola(nlomx,nlamx)
      integer   :: catica(ncmx,0:agemx)
      integer   :: realic(ncmx)
      integer   :: t0ic(ncmx)
      integer   :: ageic(ncmx)
      integer   :: dt,srmx,nc,nrealc,lows
      integer   :: nsteps,step
      integer   :: onoro,anatyp,dattyp,eofile,date,time
      integer   :: agemin,agemh
      integer   :: jump,olddat,newdat

      common /c2/   step

      namelist/namel1/onoro,anatyp,dattyp,nsteps,dt,lam1,lam2,
     &               phi1,phi2,agemh,vortmin,vortmax,gzpmin,gzpmax,
     &               srmx,lows,wmtrop,vmtrop
!----------------------------------------------------------------------

      if(do_fit_out) then
         open(781,file="fit_distzp.sav")
         open(782,file="fit_zpij.sav")
         open(784,file="fit_pp.sav")
         open(785,file="fit_gauss.sav")

         open(791,file="fit_iter.sav")    
         open(792,file="fit_fret.sav")          

         open(801,file="fit_zenv.sav")            
         open(802,file="fit_zrad.sav")            
         open(803,file="fit_zdep.sav")            
         open(804,file="fit_zmean.sav")            
         open(805,file="fit_zcen.sav")            
         open(806,file="fit_gzp.sav")                  

         open(819,file="fit_dep_vs_rad.sav")
         open(820,file="fit_dep_vs_gzp.sav")
         open(821,file="fit_rad_vs_gzp.sav")      
         open(822,file="fit_env_vs_mea.sav")      
         open(823,file="fit_dep_vs_mea.sav")      
         open(824,file="fit_gzp_vs_cen.sav")      
         open(825,file="fit_rad_vs_cen.sav")            
         open(826,file="fit_dep_vs_cen.sav")            

         open(840,file="fit_mea-cen.sav")            
         open(841,file="fit_env-mea.sav")            
      endif
      
      call prhead(nlomx,nlamx,ncmx,agemx)

!      open(22,file="/alphadata04/blumer/Data/Geopotential/Z1000/DJF/Z100
!     &0gp_DJF_${yy}_`expr $yy + 9`.srv",
!                   /alphadata04/blumer/no_Backup/Data/Postprocessed/BPRD_trans_1/Z1000/DJF/

      open(22,file="tmpz1000.srv",
     &     FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL',status='OLD')

!      open(22,file="/alphadata04/blumer/no_Backup/Data/Postprocessed/BPR
!     &D_trans_1/Z1000/DJF/Z1000gp_DJF_${yy}_`expr $yy + 9`.srv",
!     &     FORM='UNFORMATTED',
!     &     ACCESS='SEQUENTIAL',status='OLD')

!      open(22,file="tmp.srv",
!     &     FORM='UNFORMATTED',
!     &     ACCESS='SEQUENTIAL',status='OLD')

      open(29,file="oro.srv",
     &     FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL',status='OLD')
      print*,'>  ok input files '

!     namelist vars

      onoro   = 1
      anatyp  = 1
      dattyp  = 1
      nsteps  = 14600
      dt      = 6
      lam1    = 0.
      lam2    = 360.
      phi1    = 85.
      phi2    = 20.
      agemh   = 24
      vortmin = 1E-5
      vortmax = 4E-5
      gzpmin  = 20
      gzpmax  = 30
      srmx    = 0
      lows    = 1
      wmtrop  = 15.0
      vmtrop  = 3.5E-5

!s2.9.15      print*,'>  read namelist >>namel1<< '
!      read(*,namel1)
!      print*,'>  namelist found '

      call       inivar(step,nsteps,oro,nc,dt,
     &                  ageic,t0ic,lonica,latica,catica,datica,
     &                  timeica,
     &                  diica,czpica,gzpica,iclola,tropica,
     &                  zradica,zdepica,
     &                  agemh,agemin,
     &                  vortmin,gzpmin,vortmax,gzpmax,wmtrop,vmtrop,
     &                  onoro,anatyp,dattyp,srmx,lows)

      print*,'=> step= ',step

      call       indat(zp,vort,tropf,
     &                 dt,time,anatyp,dattyp,date,eofile,step)

      call outsvf(chserv,vort)

      call inidep(dt,srmx,anatyp)

      if(eofile.eq.1 .or. step.gt.nsteps) print*,'**  data/step? **'

      call       detcat(zp,vort,tropf,oro,cat0,lows,
     &                  onoro,anatyp,
     &                  vortmin,gzpmin,wmtrop,vmtrop)

      call       catpos(cat0,date,time)

      call       chacat(cat0)


      call       ancat0(cat0,zp,vort,tropf,date,time,step,
     &                  ageic,t0ic,lonica,latica,catica,datica,
     &                  timeica,
     &                  iclola,czpica,gzpica,
     &                  zradica,zdepica,     
     &                  tropica,anatyp,nc)

      do step=2,nsteps

         print*,'=> step= ',step

         olddat=date
         call       indat(zp,vort,tropf,
     &                    dt,time,anatyp,dattyp,date,eofile,step)

         if(eofile.eq.1) goto 999

         newdat=date
         call datjmp(step,olddat,newdat,time,jump)

         call       detcat(zp,vort,tropf,oro,cat1,lows,
     &                     onoro,anatyp,
     &                     vortmin,gzpmin,wmtrop,vmtrop)

         call       catpos(cat1,date,time)

         if(step.le.10) call       chacat(cat1) 

         if(jump.eq.0) then
            call    ancat1(cat0,cat1,zp,vort,tropf,
     &                     date,time,step,ageic,t0ic,lonica,latica,
     &                     catica,datica,timeica,iclola,delola,
     &                     czpica,gzpica,diica,
     &                     zradica,zdepica,     
     &                     tropica,anatyp,srmx,nc,
     &                     dt)
         else
            print*,'** jump: date not continuous in data file'

            call ancat0(cat1,zp,vort,tropf,date,time,step,
     &                  ageic,t0ic,lonica,latica,catica,datica,
     &                  timeica,
     &                  iclola,czpica,gzpica,
     &                  zradica,zdepica,     
     &                  tropica,anatyp,nc)
         endif

         call shiftc(cat0,cat1)

!        next field (time step)

      enddo !  2000 continue

!     end of file reached:
 999  continue


      call       detrcy(ageic,realic,gzpica,czpica,
!trop&                  tropica, not yet ready
     &                  agemin,anatyp,vortmax,gzpmax,nc,nrealc)

      call       statr2(nc,lonica,latica,ageic,realic)

      call       outxy(lonica,latica,ageic,realic,
     &                 nrealc,nc)

      call       outtr(nc,lonica,latica,czpica,catica,ageic,
     &                 t0ic,datica,timeica,gzpica,diica,realic,
     &                 zradica,zdepica,     
     &                 tropica,nrealc,agemin,anatyp)

      call       genlys(nc,lonica,latica,ageic,realic)

      call       ddica(nc,diica,ageic) ! ,chpand)

      call       outlt(ageic) ! ,chpalt)

!od include deallocate statement
!      deallocate(czpica)
!      deallocate(gzpica)
!      deallocate(zradica)
!      deallocate(zdepica)
!      deallocate(tropica)
!      deallocate(diica)

      write(*,*) 'xatr done'
      write(*,*) ' number corrected in fit N_corr_zenv',N_corr_zenv,
     &      N_corr_zenv/real(N_call_fitgauss)
      
      write(*,*) ' number corrected in fit N_corr_zdep',N_corr_zdep,
     &      N_corr_zdep/real(N_call_fitgauss)
     
      write(*,*) ' number corrected in fit N_corr_zrad',N_corr_zrad,
     &      N_corr_zrad/real(N_call_fitgauss)
     
      write(*,*) ' compare      number N_call_fitgauss',N_call_fitgauss
     
      end  ! main

!----------------------------------------------------------------------
      subroutine prhead(nlomx,nlamx,ncmx,agemx)
!----------------------------------------------------------------------
      integer   :: nlomx,nlamx,ncmx,agemx

      print*,' '
      print*,'=================================================='
      print*,'==   xatr                                       =='
      print*,'==   Analysis and Tracking of (Anti)Cyclones    =='
      print*,'==   Richard Blender May 6 1996/18 9 2006       =='
      print*,'==   Meteorologisches Institut                  =='
      print*,'==   Universitaet Hamburg                       =='
      print*,'==   email: Richard.Blender at zmaw.de          =='
      print*,'==                                              =='
      print*,'==   Main Features                              =='
      print*,'==   + cyclones and anticyclones                =='
      print*,'==   + different data file types                =='
      print*,'==   + use geopotential or vorticity            =='
      print*,'==   + different resolutions (space and time)   =='
      print*,'==   + globally and in regions                  =='
      print*,'==   + specify intensity and minimal lifetime   =='
      print*,'==   + use combined files                       =='
!     print*,'==   New                                        =='
!     print*,'==   + tropical storms                          =='
      print*,'=================================================='
      print*,' '
      print*,'   parameters:'
      print*,'   max longitude           nlomx   =  ',nlomx
      print*,'   max latitude            nlamx   =  ',nlamx
      print*,'   max number of cyclones  ncmx    =  ',ncmx
      print*,'   max age (in time steps) agemx   =  ',agemx
      print*,' '
      print*,'=================================================='

      return
      end

!----------------------------------------------------------------------
      subroutine inivar(step,nsteps,oro,nc,dt,
     &                  ageic,t0ic,lonica,latica,catica,datica,
     &                  timeica,
     &                  diica,czpica,gzpica,iclola,tropica,
     &                  zradica,zdepica,
     &                  agemh,agemin,
     &                  vortmin,gzpmin,vortmax,gzpmax,wmtrop,vmtrop,
     &                  onoro,anatyp,dattyp,srmx,lows)
!----------------------------------------------------------------------
      
      use xatrm
      
      implicit none

      integer   :: step,nsteps,onoro,anatyp,dattyp,ihoro(8)
      integer   :: type,i,j,ic,age,agemh,agemin,nc,itrop
      integer   :: lows,dt,srmx
      integer   :: eofile
      integer   :: ageic(ncmx)
      integer   :: t0ic(ncmx)
      integer   :: lonica(ncmx,0:agemx)
      integer   :: latica(ncmx,0:agemx)
      integer   :: catica(ncmx,0:agemx)
      integer   :: datica(ncmx,0:agemx)
      integer   :: timeica(ncmx,0:agemx)
      integer   :: iclola(nlomx,nlamx)

      real(kind=4) ::    gzpica(ncmx,0:agemx)
      real(kind=4) ::    czpica(ncmx,0:agemx)
      real(kind=4) ::    zradica(ncmx,0:agemx)      
      real(kind=4) ::    zdepica(ncmx,0:agemx)      

      
      real(kind=4) ::    diica(ncmx,0:agemx)
      real(kind=4) ::    tropica(ntrop,ncmx,0:agemx)

      real(kind=4) ::    oro(nlomx,nlamx)
      real(kind=4) ::    vortmin,gzpmin,vortmax,gzpmax
      real(kind=4) ::    wmtrop,vmtrop

      print*,'>> inivar init variables '

      print*,'   input units '
      print*,'   chzp   ',chzp
      print*,'   chvort ',chvort
      print*,'   chtrop ',chtrop
      print*,'   choro  ',choro
      print*,' '
      print*,'   output units '
      print*,'   chpath ',chpath
      print*,'   chpgrd ',chpgrd
      print*,'   chpdgr ',chpdgr
      print*,'   chpar2 ',chpar2
      print*,'   chpand ',chpand
      print*,'   chppca ',chppca
      print*,'   chxy   ',chxy
      print*,'   chserv  ',chserv
      print*,'   chcc   ',chcc
      print*,'   chpalt ',chpalt
      print*,'   chgen  ',chgen
      print*,'   chlys  ',chlys
      print*,'   chtot  ',chtot
      print*,'   chptro ',chptro
      print*,' '

      step=1
      nc=0

      do ic=1,ncmx ! 20
         ageic(ic)=-1
         diica(ic,0)=0
         t0ic(ic)=0
         do age=0,agemx
            lonica(ic,age)=0
            latica(ic,age)=0
            catica(ic,age)=0
            datica(ic,age)=0
            timeica(ic,age)=0
            gzpica(ic,age)=0
            czpica(ic,age)=0
            zradica(nc,age)=-0.1            
            zdepica(nc,age)=-1
            
            do itrop=1,ntrop
               tropica(itrop,ic,age)=0
            enddo
         enddo
      enddo !  20   continue

!     lagrangian coordinate (identification number) ic of cyclone
!     in i,j in the first cat field
!     set  0 first and actualize at each step

      do j=1,nlamx
      do i=1,nlomx
         iclola(i,j)=0
      enddo
      enddo

      print*,'   variables in namelist: '
      print*,' '
      print*,'   onoro   ',onoro
      print*,'   anatyp  ',anatyp
      print*,'   dattyp  ',dattyp
      print*,'   nsteps  ',nsteps
      print*,'   dt      ',dt
      print*,'   lam1    ',lam1  ,'   geographical coordinates'
      print*,'   lam2    ',lam2
      print*,'   phi1    ',phi1
      print*,'   phi2    ',phi2
      print*,'   agemh   ',agemh
      print*,'   vortmin ',vortmin
      print*,'   vortmax ',vortmax
      print*,'   gzpmin  ',gzpmin
      print*,'   gzpmax  ',gzpmax
      print*,'   srmx    ',srmx
      print*,'   lows    ',lows
      if(anatyp.eq.8) then
         print*,'   vmtrop  ',vmtrop
         print*,'   wmtrop  ',wmtrop
      endif

      if(anatyp.eq.8) then
         if(agemh.eq.0) then
            print*,'   min age for tropical storms default '
            agemh=32
            print*,'   agemh=',agemh,' hours'
         endif
      endif

      agemin=agemh/dt
      print*,'   agemin (in dt-steps) ',agemin

!     check input
      print*,' '
      if(onoro.lt.0 .or. onoro.gt.1) then
         print*,'*** onoro wrong STOP ***, onoro=0,1'
         stop
      endif
      if(onoro .eq.1) print*,'   orogrophy > 1000m excluded '

      if(anatyp.lt.1 .or. anatyp.gt.8) then
         print*,'*** anatyp wrong STOP ***, anatyp=1,2,8'
         stop
      endif
      if(anatyp.eq.1) print*,'   geopotential z/pressure used '
      if(anatyp.eq.2) print*,'   vorticity used '
      if(anatyp.eq.8) print*,'   search tropical cyclones '

      if(dattyp.lt.1 .or. dattyp.gt.2) then
         print*,'*** dattyp wrong STOP ***, dattyp =1,2'
         stop
      endif
      if(dattyp.eq.1) print*,'   service format '
      if(dattyp.eq.2) then
         print*,'      non service format'
         print*,'      read with double loop '
         print*,'        do 1 j=1,nlo '
         print*,'        do 1 i=1,nla '
         print*,'      1   read(unit,*) field(i,j) '
         print*,'      correct ?'
      endif

      if(nsteps.lt.1 ) then
         print*,'*** nsteps wrong STOP ***, nsteps>=1'
         stop
      endif
      if(dt.le.0 .or. dt.gt.24) then
         print*,'*** dt wrong STOP ***, 0<=dt<=24 '
         stop
      endif
      if(dt.gt.12) print*,'*** dt large, results not meaningful ***'
      if(phi1.gt. 90) then
         print*,'*** phi1 too large STOP ***, <=90.'
         stop
      endif
      if(phi1.lt.-90) then
         print*,'*** phi1 too small STOP ***, >= -90. '
         stop
      endif
      if(phi2.gt. 90) then
         print*,'*** phi2 too large STOP ***, <= 90. '
         stop
      endif
      if(phi2.lt.-90) then
         print*,'*** phi2 too small STOP ***, >= -90. '
         stop
      endif

      if(srmx.lt.0 .or. srmx.gt.20) print*,'** srmx possibly wrong **'
      if(srmx.eq.0) print*,'   srmx determined automatically '

      if((lows.ne.-1) .and. (lows.ne.1)) then
         print*,'*** lows  wrong **, lows=1,-1'
         stop
      endif
      if(lows.eq.1)   print*,'   cyclones are detected'
      if(lows.eq.-1)  print*,'   anticyclones are detected'

!      if((catmin.lt.0).or.(catmin.gt.9)) then
!         print*,'*** catmin  wrong STOP ***'
!         stop
!      endif

      if(vortmin.eq.0) then
         print*,'   vortmin determined automatically '
         print*,'   depends on anatyp and lows'
         if((anatyp.eq.2).and.(lows.eq. 1)) vortmin=1E-5
         if((anatyp.eq.8).and.(lows.eq. 1)) vortmin=3.5E-5
         if((anatyp.eq.2).and.(lows.eq.-1)) vortmin=1E-5
         print*,'   vortmin =',vortmin
      endif
      if(vortmax.eq.0) then
         print*,'   vortmax determined automatically '
         print*,'   depends on anatyp and lows'
         if((anatyp.eq.2).and.(lows.eq. 1)) vortmax=4E-5
         if((anatyp.eq.8).and.(lows.eq. 1)) vortmin=3.5E-5
         if((anatyp.eq.2).and.(lows.eq.-1)) vortmax=4E-5
         print*,'   vortmax =',vortmax
      endif

      if(gzpmin.eq.0) then
         print*,'   gzpmin determined automatically '
         print*,'   depends on anatyp and lows'
         if((anatyp.eq.1).and.(lows.eq. 1)) gzpmin =20
         if((anatyp.eq.1).and.(lows.eq.-1)) gzpmin =20
         print*,'   gzpmin =',gzpmin
      endif

      if(gzpmax.eq.0) then
         print*,'   gzpmax determined automatically '
         print*,'   depends on anatyp and lows'
         if((anatyp.eq.1).and.(lows.eq. 1)) gzpmax=20
         if((anatyp.eq.1).and.(lows.eq.-1)) gzpmax=20
         print*,'   gzpmax=',gzpmax
      endif

      if(vmtrop.eq.0) then
         if(anatyp.eq.8) then
            print*,'   min. lower level vorticity: '
            print*,'   vmtrop determined automatically '
            vmtrop=3.5E-5
            print*,'   vmtrop=',vmtrop
         endif
      endif

      if(wmtrop.eq.0) then
         if(anatyp.eq.8) then
            print*,'   min. wind 10 m '
            print*,'   wmtrop determined automatically '
            wmtrop=15
            print*,'   wmtrop=',wmtrop
         endif
      endif


!     read orography

      if(onoro.eq.1)  then
         print*,'   try to read orography unit 22'
         if(dattyp.ne.1) then
            print*,'   service format only for orography stop'
            stop
         endif
         call insvf(choro,ihoro,oro,type,eofile)
         if(eofile.eq.1) print*,'** orography file eof !'
         print*,'   orography found '
         print*,'   size',nlo,nla,' type',type
      endif

      print*,' '
      return
      end

!----------------------------------------------------------------------
      subroutine inidep(dt,srmx,anatyp)
!----------------------------------------------------------------------
!     init variables after reading the first data file

      use xatrm

      implicit none
      
      integer   :: anatyp
      integer   :: dt 
      integer   :: srmx,i,inlo,idt

      print*,'>> inidep  init data dependent variables'

      if(srmx.eq.0) then
         print*,'   srmx automatic (diffusive)'
         print*,'      formula for srmx: (for srmx=0 input)'
         print*,'      srmx=int(1.0+nlo/64*(dt/6.)**0.75)'
         print*,'      will be set to 1 if this formula yields 0'
         srmx=int(1.0+nlo/64*(dt/6.)**0.75)
         if(srmx.eq.0) then
            srmx=1
            print*,'   set srmx=1 after formula yielded 0'
         endif

         print*,'      for various nlo and dt: srmx'
         do inlo=1,5
         do idt=1,4
            write(*,'(I10,2I6)') inlo*64,idt*6,
     &              int(0.5+inlo*(idt*1.)**0.75)
         enddo
         enddo
      endif

      print*,'   search region in grid points at the equator'
      print*,'      longitude srmx        ',srmx
      print*,'      latitude  (srmx+1)/2  ',(srmx+1)/2
      print*,' '

!     determine analyzed region and check input

      nla1=int( 0.5 + 1+(0.5-phi1/180.)*(nla-1) )
      nla2=int( 0.5 + 1+(0.5-phi2/180.)*(nla-1) )

      if(nla1.ge.nla2) then
         print*,'*** wrong meridional range *** phi1,2 ',phi1,phi2
         print*,'*** nla1, nla2 ',nla1,nla2
         print*,'*** STOP ***'
         stop
      endif

      nlo1=1 + int( 0.5 + lam1/360.*nlo )
      nlo2=1 + int( 0.5 + lam2/360.*nlo )

      if(nlo2.gt.nlo) then
         print*,'   analyzed range is global '
         nlo2=nlo
      endif

      if(nlo1.lt.nlo2) then
         print*,'   analyzed range between longs ',lam1,lam2
      else
         print*,'   analyzed range between longs ',lam2,lam1
      endif

      print*,' '
      print*,'   analyzed region in grid points '
      print*,'   long: nlo1 nlo2 ',nlo1,nlo2
      print*,'   lat:  nla1 nla2 ',nla1,nla2

      if(anatyp.eq.1) then
         print*,'   geopot cat formula: int(gzp/gzpmin)'
      endif

      if(anatyp.eq.2) then
         print*,'   vorticity cat formula:  '
         print*,'   cat(i,j)=int( abs(vort(i,j)/vortmin) ) '
         print*,' '
      endif

      if(anatyp.eq.8) then
         print*,'   trop cat formula:  '
         print*,'      cat(i,j)=   int(max/wmtrop)'
         print*,' '
      endif

C     longitude periodic map

      do i=-nlomx,2*nlomx
         p(i)=i
         if(p(i).gt.nlo) p(i)=p(i)-nlo
         if(p(i).lt.1  ) p(i)=p(i)+nlo
      enddo ! 10   continue

      print*,' '

      return
      end

!----------------------------------------------------------------------
      subroutine insvf(chann,ih,f,type,eofile)
!----------------------------------------------------------------------
!     read a service format field
!     ih     header, f      field

      use xatrm
      implicit none

      integer   :: chann,ih(8),i,j,eofile 
      integer   :: date,time,type
      real  (kind=4) ::      f(nlomx,nlamx)

      eofile=0
      print*,'>> insvf '

      print*,'   read service format field, unit ',chann
      read(chann,end=9) ih

      type=ih(1)
      date=ih(3)
      time=ih(4)
      nlo=ih(5)
      nla=ih(6)

      print*,type
      print*,date
      print*,time
      print*,nlo
      print*,nla
      print*,ih
      STOP
      if((nlo.gt.nlomx).or.(nla.gt.nlamx)) then
         print*,'   input file too large STOP '
         print*,'   ',ih(5),' x ',ih(6)
         STOP
      endif

      print*,'   field found, type, size ',type,nlo,nla
      print*,'                date, time ',date,time

      read(chann) f
   
      goto 90

!     eofile=1 if end of file
  9   continue
         eofile=1
         print*,' '
         print*,'*  insvf: end of file'
         print*,' '
  90  continue


      print*,' '

      return
      end

!----------------------------------------------------------------------
      subroutine outsvf(chann,f)
!----------------------------------------------------------------------
!     output of f in servive format

      use xatrm
      implicit none

      integer        :: chann
      integer        :: ih(8),i,j
      real  (kind=4) :: f(nlomx,nlamx)

C     header not completely known, only size

      do i=1,8
         ih(i)=0
      enddo
      ih(5)=nlo
      ih(6)=nla

      write(chann) ih

      do j=1,nla
      do i=1,nlo
         write(chann) f(i,j)
      enddo
      enddo

      return
      end

!----------------------------------------------------------------------
      subroutine inothf(chann,ih,f,type,eofile)
!----------------------------------------------------------------------
!     read non-service format field
      
      use xatrm
      
      implicit none
      integer        :: chann,type,eofile
      integer        :: ih(8),i,j,dt,step
      real(kind=4)   :: f(nlomx,nlamx)

      common /c2/ step

      eofile=0

      print*,'>> inothf:  read non-service format file, unit ',chann

      do i=1,8
         ih(i)=0
      enddo

      nlo=144
      nla=96
      type=0
      dt=12

      ih(1)=type
      ih(2)=0
      ih(3)=900101+(step-1)*12/dt
      ih(4)=mod(step*dt,12)
      ih(5)=nlo
      ih(6)=nla
      ih(7)=0
      ih(8)=0

      print*,'   artific. date/time ',ih(3),ih(4),',size ',nlo,nla
      print*,'   type=',type,'  dt=',dt

!     particular form of anna ghelli's file, in rows

      do j=1,nla
      do i=1,nlo
         read(chann,*,end=9) f(i,j)
      enddo
      enddo

! 8378 format(1e10.4)

!     eofile=1 if end of file
      goto 90
  9   continue
         eofile=1
         print*,' '
         print*,'*  inothf: end of file'
         print*,' '
  90  continue

      return
      end

!----------------------------------------------------------------------
      subroutine indat(zp,vort,tropf,
     &                 dt,time,anatyp,dattyp,date,eofile,step)
!----------------------------------------------------------------------

!     read geopot/vorticity/etc  fields until end of file
!     files depend on anatyp and on dattyp
      
      use xatrm
      
      implicit none

      integer        :: anatyp,dattyp,date,time,eofile,dt
      integer        :: step,type
      real(kind=4)   :: zp(nlomx,nlamx)
      real(kind=4)   :: vort(nlomx,nlamx)
      real(kind=4)   :: tropin(nlomx,nlamx)
      real(kind=4)   :: tropf(ntrop,nlomx,nlamx)
      integer        :: ih(8),day,month,year,i,j,itrop

      print*,'>> indat '

!     read zp if anatyp=1
!     ----------------------
      if(anatyp.eq.1) then
         print*,'   read geopotential/pressure zp '
         if(dattyp.eq.1) then
            call insvf(chzp,ih,zp,type,eofile)
            if((type.ne.156).and.(type.ne.134)) then
               print*,'*** warning: type not geopot/surf press **'
            endif
            do i=1,nlo
            do j=1,nla
               zp(i,j) = zp(i,j)
            end do
            end do   
         endif
         if(dattyp.eq.2) then
            call inothf(chzp,ih,zp,type,eofile)
         endif
         if(eofile.eq.1) then
            print*,'   eof of file reached '
            goto 90
         endif
      endif

!     read vorticity if anatyp=2
!     --------------------------
      if(anatyp.eq.2) then
         print*,'   read vorticity'
         if(dattyp.eq.1) then
            call insvf(chvort,ih,vort,type,eofile)
            if(type.ne.138) then
               print*,'*** warning: type not vorticity **',type
            endif
         endif
         if(dattyp.eq.2) then
            call inothf(chvort,ih,vort,type,eofile)
         endif
         if(eofile.eq.1) then
            print*,'   eof of file reached '
            goto 90
         endif
      endif

!     read zp, vort and tropf (wind10,vort,temp) if anatyp=8
!     ------------------------------------------------------
      if(anatyp.eq.8) then
         print*,'   read zp,vort,tropf (=wind10,vort,temp) '
         if(dattyp.eq.1) then
!           surface pressure
            call insvf(chzp,ih,zp,type,eofile)
            if(type.ne.156.or.type.ne.134) then
               print*,'*** warning: type not geopot or surf press **'
            endif

!           lower level vorticity
            call insvf(chvort,ih,vort,type,eofile)
            if(type.ne.138) then
               print*,'*** warning: type not vorticity **'
            endif

!           read additional tropical data
            do itrop=1,ntrop

               call insvf(chtrop+itrop-1,ih,tropin,type,eofile)
               if(itrop.eq.1 .and. type.ne.171) then
                  print*,'*** warning: type not wind 10m **'
               endif
               if(itrop.eq.2 .and. type.ne.138) then
                  print*,'*** warning: type not vorticity **'
               endif
               if(itrop.eq.3 .and. type.ne.130) then
                  print*,'*** warning: type not temperature **'
               endif

               do i=1,nlo
               do j=1,nla
                  tropf(itrop,i,j)=tropin(i,j)
               enddo
               enddo

            enddo

         endif

         if(dattyp.eq.2) then
!           type not checked
            call inothf(chzp  ,ih,zp  ,type,eofile)
            call inothf(chvort,ih,vort,type,eofile)

            do itrop=1,ntrop
               call inothf(chtrop+itrop-1,ih,tropin,type,eofile)
               do i=1,nlo
               do j=1,nla
                  tropf(itrop,i,j)=tropin(i,j)
               enddo
               enddo
            enddo 

         endif

         if(eofile.eq.1) then
            print*,'   eof of file reached '
            goto 90
         endif
      endif

      time=ih(4)

      if(step.eq.0 .and. time.ne.0) then
         print*,'*** '
         print*,'*** first data field at time ',time
         print*,'*** correct ?'
         print*,'*** '
      endif

      if((step.eq.1) .and. ((time.ne.dt).and.(time.ne.0))) then
         print*,'*** '
         print*,'*** possible error ***'
         print*,'*** second data field at time=',time
         print*,'*** but timestep:          dt=',dt
         print*,'*** STOP '
         stop
      endif

      date=ih(3)

      call datdmy(date,day,month,year)

 90   continue

      return
      end

!----------------------------------------------------------------------
      subroutine detcat(zp,vort,tropf,oro,cat,lows,
     &                  onoro,anatyp,
     &                  vortmin,gzpmin,wmtrop,vmtrop)
!----------------------------------------------------------------------

!     this subroutine searches the centers of vortices and marks
!     them in the field cat(i,j) by the integers 1..9, according
!     to the intensity
!     for zp, vort ot tropical depressions
!     ncat   number of cat=1,2,3,...
      
      use xatrm
      
      implicit none

      integer      :: onoro,anatyp,lows
      real(kind=4) :: zp(nlomx,nlamx)
      real(kind=4) :: vort(nlomx,nlamx)
      real(kind=4) :: tropf(ntrop,nlomx,nlamx)
      real(kind=4) :: tropf1(nlomx,nlamx)
      real(kind=4) :: oro(nlomx,nlamx)
      real(kind=4) :: gzp,vmtrop,wmtrop,max,zrad,zdep
      real(kind=4) :: vortmin,gzpmin

      integer   :: cat(nlomx,nlamx),i,j,i1,j1,k,answer,ncat(0:9),sregion

      print*,'>> detcat: determine cat, anatyp= ',anatyp

      do i=0,9
         ncat(i)=0
      enddo

      do j=1,nla
      do i=1,nlo
         cat(i,j)=0
      enddo
      enddo

      if(lows.eq.-1) then
!        highs are detected, change signs
         do j=1,nla
         do i=1,nlo
            if(anatyp.eq.1) zp(i,j)=-zp(i,j)
            if(anatyp.eq.2) vort(i,j)=-vort(i,j)
         enddo
         enddo
      endif

!     change sign of vorticity in nh, minima over the globe!
!     is retransformed at the end

      do j=1,nla/2
      do i=1,nlo
         vort(i,j)=-vort(i,j)
      enddo
      enddo

      open (unit=80,file='oro.txt',STATUS='REPLACE')

      do j=1,nla
	write(80,'(*(F20.7,:,","))') oro(:,j)
      do i=1,nlo

!        drop if orography is too high, 1000m
         if((onoro.eq.1).and.(oro(i,j).ge.9810.)) go to 30 !9810 3.9. sandro
!        --------------------
         if(anatyp.eq.1) then
!           zp minima, min. if answer=1, =0 otherwise
            call ijmin1(zp,i,j,answer)
            if(answer.eq.1) then
               call gradzp(zp,i,j,gzp)
               cat(i,j)=int(gzp/gzpmin)

	       if (cat(i,j).gt.0) then
		  call fitgauss(zp,i,j,zrad,zdep,gzp)
		  if (zrad.lt.0.1) cat(i,j)=0
	       endif
            else
               cat(i,j)=0
            endif
         endif
!        --------------------
         if(anatyp.eq.2) then
!           vorticity min. in i,j? (changed sign in nh!)
            call ijmin1(vort,i,j,answer)
            if(answer.eq.1) then
!              vorticity has changed sign  in NH!
               cat(i,j)=int( abs(vort(i,j)/vortmin) )
            else
               cat(i,j)=0
            endif
         endif
!        --------------------
         if(anatyp.eq.8) then
!           tropical cyclones
!           itrop=1 wind 10m                  code 171
!           itrop=2 upper level vorticity     code 138
!           itrop=3 temperature               code 130
!           position at zp-minimum (surface pressure)
            call ijmin1(zp,i,j,answer)
            if(answer.eq.1) then
!              surface position found
!              search max of vorticity in neighbourhood
               sregion=nlo/45
               if(ntrop.ge.2) then
                  do i1=1,nlo
                  do j1=1,nla
                     tropf1(i1,j1)=tropf(2,i1,j1)
                  enddo
                  enddo
               endif
               call maxreg(tropf1,i,j,sregion,max)
               if(max.gt.vmtrop) then
                  do i1=1,nlo
                  do j1=1,nla
                     tropf1(i1,j1)=tropf(1,i1,j1)
                  enddo
                  enddo
                  call maxreg(tropf1,i,j,sregion,max)
                  if(max.gt.wmtrop) then
                     cat(i,j)= int(max/wmtrop)
                  else
                    cat(i,j)=0
                  endif
               else
                  cat(i,j)=0
               endif

            else
               cat(i,j)=0
            endif
         endif

 8123    format(i5,3g16.8)

!        correct cat if outside 0...9
         if(cat(i,j).gt.9) cat(i,j)=9
         if(cat(i,j).lt.0) cat(i,j)=0
         
  30  continue ! from a goto 30
  
      enddo ! j 30   continue
      enddo ! i
      
      call chkcat(cat) 
      close(80)

! new at August 13,  1997

!c    call modcat(cat) 

      do i=1,nlo
      do j=1,nla
         ncat(cat(i,j))=ncat(cat(i,j)) + 1
      enddo
      enddo
 
      print*,'   no of cycs with cat=0,1,2,3,...'
      write(*,8000) (ncat(k),k=0,9)
 8000 format(i9,10i6)

!     change sign of vorticity back '

      if(anatyp.eq.2) then
         do j=1,nla/2
         do i=1,nlo
            vort(i,j)=-vort(i,j)
         enddo
         enddo
      endif


      return
      end

!----------------------------------------------------------------------
      subroutine chkcat(cat)
!----------------------------------------------------------------------
!     check if two neighbor points have cat>0
!     repair field, set cat=0 for the neighbor detected
      
      use xatrm
      implicit none
      integer   :: cat(nlomx,nlamx)
      integer   :: i,j,di,dj

      do j=1+1,nla-1
      do i=1,nlo
         if(cat(i,j).gt.0) then
            do di=-1,1
            do dj=-1,1
!              scan 1-neighbohood (without i j)
               if((di.ne.0).or.(dj.ne.0)) then
               if(cat(p(i+di),j+dj) .gt. 0) then
                  if(test.gt.0) then
                     print*,'>> chkcat: neighbor cat=cat, repair ',i,j
                  endif
                  cat(p(i+di),j+dj)=0
               endif
               endif
            enddo
            enddo
         endif
      enddo
      enddo

      return
      end

!----------------------------------------------------------------------
      subroutine catpos(cat,date,time)
!----------------------------------------------------------------------
!     output of cyclones found in cat
!     in long/lat coords
      
      use xatrm
      
      implicit none
      integer   :: date,time
      integer   :: cat(nlomx,nlamx)
      integer   :: i,j
      real(kind=4) ::    xlon,xlat

      write(chpos,*) date, time

      do j=1+1,nla-1
      do i=1,nlo
         if(cat(i,j).gt.0) then
            xlat=90.-180*j/(nla+1.)
            xlon=(i*360.)/nlo
            write(chpos,8000) xlon,xlat,cat(i,j)
         endif
      enddo
      enddo

 8000 format(2f12.3,i3)

      return
      end

!----------------------------------------------------------------------
      subroutine modcat(cat) 
!----------------------------------------------------------------------

!     modify cat field to eliminate nearby cyclones
!     in the range sr (in gridpoints)
!     a weaker one next to a strong one is eliminated
      
      use xatrm
      implicit none
      integer   :: cat(nlomx,nlamx)
      integer   :: i,j,di,dj,mi,mj,i0,j0,sr,nelim
      real(kind=4) ::    dphi,phi

!     the range is 5 degrees (1000km in total)
      sr=nlo/64

      nelim=0

      do i=1,nlo
      do j=1,nla

!        skip if not occupied
         if(cat(i,j).eq.0) goto 10

         dphi=pi/(nla+1)
         phi =pi/2-j*dphi
         mi  =int(0.5+ sr/abs(cos(phi)) )

         mj  =sr
         if((nlo.le.64).and.(sr.eq.1)) mj=1

         do 1 di=-mi,mi
         do 1 dj=-mj,mj
!           exclude centre
            if((di.eq.0).and.(dj.eq.0)) goto 1
!           exclude pole boundary crossing
            if( ((j+dj).lt.1).or.((j+dj).gt.nla) ) goto 1
            i0=p(i+di)
            j0=j+dj
            if(cat(i0,j0).gt.0) then
!              another cyc found, keep the stronger one
               if(cat(i,j).gt.cat(i0,j0)) then
                  cat(i0,j0)=0
               else
                  cat(i,j)=0
               endif
               nelim=nelim+1
            endif
 1       continue
 10      continue

      enddo !  i,j  10  continue
      enddo !       

      print*, '>> modcat: ',nelim,' features eliminated in cat'

      return
      end

!----------------------------------------------------------------------
      subroutine chacat(cat)
!----------------------------------------------------------------------
!     output of cat in character form
      
      use xatrm
      
      implicit none
      integer     :: cat(nlomx,nlamx),i,j
      character*1 :: catc(nlomx,nlamx)

      print*,'>> chacat:  character output of cat, unit ',chcc
      print*,'            longitude maximum 64 grid points'


      do j=1,nla
      do i=1,nlo
         if(cat(i,j).eq.0) catc(i,j)='.'
         if(cat(i,j).eq.1) catc(i,j)='1'
         if(cat(i,j).eq.2) catc(i,j)='2'
         if(cat(i,j).eq.3) catc(i,j)='3'
         if(cat(i,j).eq.4) catc(i,j)='4'
         if(cat(i,j).eq.5) catc(i,j)='5'
         if(cat(i,j).eq.6) catc(i,j)='6'
         if(cat(i,j).eq.7) catc(i,j)='7'
         if(cat(i,j).eq.8) catc(i,j)='8'
         if(cat(i,j).eq.9) catc(i,j)='9'
      enddo
      enddo
 
      write(chcc,2013)
 2013 format('  j  12345678901234567890123456789012345678901234567890')
      do j=1,nla
         write(chcc,'(i4,1x,120a1)') j,(catc(i,j),i=1,110)
      enddo

      return
      end

!----------------------------------------------------------------------
      subroutine ijmin1(zp,i,j,answer)
!----------------------------------------------------------------------
!  check whether zp(i,j) is smaller than in the neighborhood
!  then answer=1, =0 otherwise
      
      use xatrm
      implicit none 
      real(kind=4) :: zp(nlomx,nlamx)
      integer      :: i,j,di,dj,mi,mj,answer

!     set answer=0 if (i,j) is near the poles

      if((j.eq.1) .or. (j.eq.nla)) then
         answer=0
         return
      endif

      answer=1

!     neighborhood range
      mi=1
      mj=1

!     set answer=0 if z is smaller somewhere in the neighborhood

      do di=-mi,mi
      do dj=-mj,mj
         if((di.ne.0).or.(dj.ne.0)) then
            if((j+dj.ge.1) .and. (j+dj.le.nla)) then
               if( zp(p(i+di),j+dj) .lt. zp(i,j) ) then
                  answer=0
                  goto 9
               endif
           endif
         endif
      enddo
      enddo

 9    return
      end

!----------------------------------------------------------------------
      subroutine gradzp(zp,i,j,gzp)
!----------------------------------------------------------------------
!  mean gradient in z around i,j in m/1000km
      
      use xatrm
      implicit none
      real(kind=4) :: zp(nlomx,nlamx)
      real(kind=4) :: sumgzp,gzp ,j1,j2
      real(kind=4) :: phi,dphi,dist,dx,dy,nt
      integer      :: di,dj,mi,mj,i,j,nnb,k


!  r0 radius of earth, in 1000km
!  dx, dy distance along long/lat between i,j-points in 1000km
!  in latitude simply dphi=const
!  the range is determined by nt

      nt    =nlo/64
      dphi  =pi/(nla+1)
      phi   =pi/2-(j)*dphi        !change on 27.07.2006 pi/2-(j+0.5)*dphi
      dx    =abs(cos(phi))*2*pi*r0/nlo
      dy    =r0*dphi
      nnb   =0
      sumgzp=0

!  range

      mi=int( nt/abs(cos(phi)) )
      mj=nt
      if(mi.lt.1) mi=1
      if(mj.lt.1) mj=1

      if(j-mj.lt.1 .or. j+mj.gt.nla) then
!        exclude boundaries
         gzp=0
      else
         do di=-mi,mi
         do dj=-mj,mj

!             j1 = 90.-1.125 *j
!             j2 = j1+dj*1.125
             j1 = 0.+180./(nlamx+1) *j
             j2 = j1+dj*180./(nlamx+1)

            if((di.ne.0).or.(dj.ne.0)) then
               dist=sqrt( (di*dx)**2 + (dj*dy)**2 )
               sumgzp=sumgzp +  (zp(p(i+di),j+dj) - zp(i,j))/dist
!            sumgzp=sumgzp +  (zp(p(i+di),j+dj)* (sin(pi*50./180.))
!     &           / sin((pi*(j2))/180.) - zp(i,j)* (sin(pi*50./180.))
!     &            / sin((pi*(j1))/180.) )/dist
               nnb=nnb+1
            endif
         enddo
         enddo

         gzp=sumgzp/nnb
      endif

      return
      end

!----------------------------------------------------------------------
      subroutine ancat0(cat0,zp,vort,tropf,date,time,step,
     &                  ageic,t0ic,lonica,latica,catica,datica,
     &                  timeica,
     &                  iclola,czpica,gzpica,
     &                  zradica,zdepica,     
     &                  tropica,anatyp,nc)
!----------------------------------------------------------------------
!     analyze first cat field, cat0
      
      use xatrm
      implicit none
      
      integer   :: step,anatyp
      integer   :: date,time,age,nc,i,j,itrop
      integer   :: cat0(nlomx,nlamx)
      integer   :: ageic(ncmx),t0ic(ncmx)
      integer   :: lonica(ncmx,0:agemx)
      integer   :: latica(ncmx,0:agemx)
      integer   :: catica(ncmx,0:agemx)
      integer   :: datica(ncmx,0:agemx)
      integer   :: timeica(ncmx,0:agemx)
      integer   :: iclola(nlomx,nlamx)
      real(kind=4) ::    zp(nlomx,nlamx)
      real(kind=4) ::    vort(nlomx,nlamx)
      real(kind=4) ::    tropf(ntrop,nlomx,nlamx)
      real(kind=4) ::    czpica(ncmx,0:agemx)
      real(kind=4) ::    gzpica(ncmx,0:agemx)
      real(kind=4) ::    zradica(ncmx,0:agemx)      
      real(kind=4) ::    zdepica(ncmx,0:agemx)      
      real(kind=4) ::    tropica(ntrop,ncmx,0:agemx)
      real(kind=4) ::    gzp,zdep,zrad
      logical   insr

      print*,'>> ancat0  analyze first field '

!  analyze cat0 and store cyclones

      age=0
!.    nc =0 no more here, data jumps possible

      do j=1,nla
      do i=1,nlo
      
      if(insr(i,j)) then
         if( cat0(i,j).gt.0 ) then
!           new cyclone found, identified
            nc       =nc+1
            ageic(nc)=age
!           store position, start time, lagrangian coord, and categ.
            lonica(nc,age)=i
            latica(nc,age)=j

            if(anatyp.eq.1) then
               czpica(nc,age)=zp(i,j)
               call gradzp(zp,i,j,gzp)
               gzpica(nc,age)=gzp
               if(do_fit) then
                  call fitgauss(zp,i,j,zrad,zdep,gzp)
                  if(do_fit_out) 
     &            write(*,88)' ancat0 fitg',i,j,zrad,zdep,zp(i,j),gzp
  88              format(A,2i4,4f12.3)
                  write(77,'(A,2i5,3e13.4)') 
     &                  '0',nc,age,zrad,zp(i,j),gzp
                  zradica(nc,age)=zrad
                  zdepica(nc,age)=zdep
               endif
            endif

            if(anatyp.eq.2) then
               czpica(nc,age)=vort(i,j)
               call gradzp(vort ,i,j,gzp)
               gzpica(nc,age)=gzp
            endif

            if(anatyp.eq.8) then
!              tropical storms, czp and gzp vorticity and gradient
               czpica(nc,age)=vort(i,j)
               call gradzp(vort ,i,j,gzp)
               gzpica(nc,age)=gzp
!              other parameters
               do itrop=1,ntrop
                  tropica(itrop,nc,age)=tropf(itrop,i,j)
               enddo
            endif

            t0ic(nc)      =step
            datica(nc,age)=date
            timeica(nc,age)=time
            iclola(i,j)   =nc
            catica(nc,age)=cat0(i,j)
         endif
      endif

      enddo
      enddo

      print*,'   number of cyclones found, nc ',nc

      return
      end

!----------------------------------------------------------------------
      subroutine ancat1(cat0,cat1,zp,vort,tropf,
     &                  date,time,step,ageic,t0ic,lonica,latica,
     &                  catica,datica,timeica,iclola,delola,
     &                  czpica,gzpica,diica,
     &                  zradica,zdepica,     
     &                  tropica,anatyp,srmx,nc,
     &                  dt)
!----------------------------------------------------------------------
!  search procedure:
!  look in increasing ranges around i,j for a previous cyclone
!  width sr, from 0...srmax
!  if no cyclone found in the largest neighborhood store as new cyclone
!  the search region is restricted by nla1,2 and nlo1,2 (see insr)
      use xatrm
      implicit none

      integer      :: date,time,step,anatyp

      real(kind=4) :: zp(nlomx,nlamx)
      real(kind=4) :: vort(nlomx,nlamx)
      real(kind=4) :: tropf(ntrop,nlomx,nlamx)
      real(kind=4) :: czpica(ncmx,0:agemx)
      real(kind=4) :: diica(ncmx,0:agemx)
      real(kind=4) :: gzpica(ncmx,0:agemx)
      real(kind=4) :: zradica(ncmx,0:agemx)      
      real(kind=4) :: zdepica(ncmx,0:agemx)      
      real(kind=4) :: tropica(ntrop,ncmx,0:agemx)
      real(kind=4) :: dsmall,gzp,zdep,zrad

      integer      :: cat0(nlomx,nlamx)
      integer      :: cat1(nlomx,nlamx)
      integer      :: iclola(nlomx,nlamx)
      integer      :: delola(nlomx,nlamx)
      integer      :: lonica(ncmx,0:agemx)
      integer      :: latica(ncmx,0:agemx)
      integer      :: datica(ncmx,0:agemx)
      integer      :: timeica(ncmx,0:agemx)
      integer      :: catica(ncmx,0:agemx)
      integer      :: t0ic(ncmx)
      integer      :: ageic(ncmx)
      integer      :: i,j,i0,j0,answer,itrop,dt
      integer      :: icold,age
      integer      :: sr,nc,ncf,ncnew,srmx
      
      logical   insr

      print*,'>> ancat1:  analyze fields '

!  analyze cat1,in a smaller j-range
!  check whether a cyclone is in regions of increasing size sr
!  search from east to west
!  in neighborhoods, depends also on dt
!  sr: search region size, depends on resolution and time step
!  actual region depends on latitute in search, increases at poles
!  and is sr/2 in meridional direction

      do j=1,nla
      do i=1,nlo
         delola(i,j)=0
      enddo
      enddo

      do 204 sr=0,srmx

      ncf  =0

      do 205 j=1,nla
      do 205 i=1,nlo
      if(insr(i,j)) then
!        look for cyclones that are not already identified
         if((cat1(i,j).eq.0).or.(delola(i,j).eq.1)) goto 205
         call search(cat0,i,j,i0,j0,answer,sr,dt)
!        no cyc found:
         if(answer.eq.0) goto 205
!        i0,j0 is outside of search region:
         if(.not.insr(i0,j0)) goto 205
!.          print*,':: ',i,j

!           cyclone found in i0,j0, considered as the same
!           identify, determine lagrangian coordinate, is detected
            delola(i,j)=1
            ncf        =ncf+1
            icold      =iclola(i0,j0)
            if(icold.le.0) then
               print*,'** error, wrong detect of cyc **'
               print*,'**  icold=0!,i,j,i0,j0 ',i,j,i0,j0
               print*,'**  cat1(i,j) cat0(i0,j0)',cat1(i,j),cat0(i0,j0)
               print*,'**  STOP'
               stop
            endif
            age=ageic(icold)+1
            if(age.gt.agemx)  then
               print*,'*  age > agemx, ic=',icold,' in ',i,j
            endif
            ageic(icold)=age
            if(age.le.agemx) then
!              if storage available store date, position, geopot, cat
               datica(icold,age)=date
               timeica(icold,age)=time
               lonica(icold,age)=i
               latica(icold,age)=j
! old          if(anatyp.eq.1) then
! old             czpica(icold,age)=zp(i,j)
! old             call gradzp(zp,i,j,gzp)
! old             gzpica(icold,age)=gzp
! old          endif
               if(anatyp.eq.1) then
                  czpica(icold,age)=zp(i,j)
                  call gradzp(zp,i,j,gzp)
                  gzpica(icold,age)=gzp
                  if(do_fit) then
                     call fitgauss(zp,i,j,zrad,zdep,gzp)
                     if(do_fit_out) 
     &               write(*,88)' ancat1 fitg',i,j,zrad,zdep,zp(i,j),gzp
                     write(77,'(A,2i5,3e13.4)') 
     &                     '1',icold,age,zrad,zp(i,j),gzp
                     zradica(icold,age)=zrad
                     zdepica(icold,age)=zdep
                  endif
               endif
               if(anatyp.eq.2) then
                  czpica(icold,age)=vort(i,j)
                  call gradzp(vort ,i,j,gzp)
                  gzpica(icold,age)=gzp
               endif
               if(anatyp.eq.8) then
!                 tropical storms, czp and gzp vorticity and gradient
                  czpica(nc,age)=vort(i,j)
                  call gradzp(vort ,i,j,gzp)
                  gzpica(nc,age)=gzp
!                 other parameters
                  do itrop=1,ntrop
                     tropica(itrop,nc,age)=tropf(itrop,i,j)
                  enddo
               endif


               catica(icold,age)=cat1(i,j)
!              distance from i,j to j0,j0 (small only)
               diica(icold,age)=dsmall(i,i0,j,j0)
            endif
!           renew lagrangian coordinate
            iclola(i0,j0)=0
            iclola(i,j)  =icold
!           remove this cyclone in cat0 to prevent further detection
            cat0(i0,j0)  =0
      endif
 205  continue

      print*,'   sr nc ncf  ',sr,nc,ncf
 204  continue

!     if not identified up to now store as new cyclone

      ncnew=0

      do j=1,nla
      do i=1,nlo
      if(insr(i,j)) then
         if((cat1(i,j).gt.0).and.(delola(i,j).ne.1)) then
            ncnew=ncnew+1
            nc   =nc+1
            if(nc.gt.ncmx) then
                  print*,' '
                  print*, '** number nc of cycs too large, STOP'
                  print*, '** date=',date,' nc=',nc,' ncmx=',ncmx
                  print*, '** increase ncmx'
                  stop
            endif
            age           =0
            ageic(nc)     =age
            lonica(nc,age)=i
            latica(nc,age)=j
            if(anatyp.eq.1) then
               czpica(nc,age)=zp(i,j)
               call gradzp(zp,i,j,gzp)
               gzpica(nc,age) =gzp
               if(do_fit) then
                  call fitgauss(zp,i,j,zrad,zdep,gzp) 
                  if(do_fit_out) 
     &               write(*,88)' ancat1 fitg',i,j,zrad,zdep,zp(i,j),gzp
                  write(77,'(A,2i5,3e13.4)') '1',nc,age,zrad,zp(i,j),gzp
                  zradica(nc,age)=zrad
                  zdepica(nc,age)=zdep
               endif
  88           format(A,2i4,4f12.3)
            endif
            if(anatyp.eq.2) then
               czpica(nc,age)=vort(i,j)
               call gradzp( vort,i,j,gzp)
               gzpica(nc,age)=gzp
            endif
            t0ic(nc)      =step
            datica(nc,age)=date
            timeica(nc,age)=time
            iclola(i,j)   =nc
            catica(nc,age)=cat1(i,j)
         endif
      endif
      enddo
      enddo

      print*,'   nc new ',ncnew,' number of new cycs'
      print*,' '

      return
      end

!----------------------------------------------------------------------
      subroutine shiftc(cat0,cat1)
!----------------------------------------------------------------------
!     shift cat1 to cat0
      
      use xatrm
      implicit  none
      integer   :: cat0(nlomx,nlamx)
      integer   :: cat1(nlomx,nlamx)
      integer   :: i,j

      do j=1,nla
      do i=1,nlo
         cat0(i,j)=cat1(i,j)
      enddo
      enddo

      return
      end

!----------------------------------------------------------------------
      subroutine detrcy(ageic,realic,gzpica,czpica,
!trop&                  tropica,
     &                  agemin,anatyp,vortmax,gzpmax,nc,nrealc)
!----------------------------------------------------------------------
!  consider only those cyclones ic with age gt agemin and
!  vort > vortmax or gzp > gzpmax at least once
!  these are called real cyclones, realic=1
      
      use xatrm
      implicit none

      integer      :: realic(ncmx)
      integer      :: ageic(ncmx)
      integer      :: age,agemin,nc,nrealc,ic,imin,anatyp
      real(kind=4) :: czpica(ncmx,0:agemx)
      real(kind=4) :: gzpica(ncmx,0:agemx)
      real(kind=4) :: vortmax,gzpmax

      print*,'>> detrcy determine intense (anti) cyclones '

      if(anatyp.eq.8) then
         print*,'>> tropical max condit not clear, set all realc=1'
      endif

      print*,'   agemin ',agemin,' minimum age (in time steps)'
      if(anatyp.eq.2) then
         print*,'   condition for intense vortex > abs(vortmax)',
     &          vortmax
      endif

      do 501 ic=1,nc
         realic(ic)=0
         if( ageic(ic).ge.agemin ) then
            do 502 age=0,imin(ageic(ic),agemx)
               if(anatyp.eq.1) then
                  if(gzpica(ic,age).ge.gzpmax) realic(ic)=1
               endif
               if(anatyp.eq.2) then
                  if(abs(czpica(ic,age)).ge.vortmax) realic(ic)=1
               endif
               if(anatyp.eq.8) then
!                 tropical conditions not yet clear
                  realic(ic)=1
               endif
 502        continue
         endif
 501  continue

      nrealc=0
      do 503 ic=1,nc
 503     if(realic(ic).eq.1) nrealc=nrealc+1

      print*,'   nrealc ',nrealc,'  number of intense (anti) cyclones'
      print*,' '

      return
      end

!----------------------------------------------------------------------
      subroutine ddica(nc,diica,ageic) ! ,chpand)
!----------------------------------------------------------------------
!   determine distribution of distances
!   nd(0..id)    number of distances, id=number of slots
      
      use xatrm
      
      integer   ::  nc,ic,age,idmax
      integer   ::  ageic(ncmx)
      real(kind=4) ::     diica(ncmx,0:agemx),deltar

!     idmax*deltar should be the maximum distance in one time step

      parameter  (idmax=40,deltar=0.1)
      integer   ::  nd(0:idmax),imin

      print*,'>> ddica, determine frequency of distances in diica'

      do 50 id=0,idmax
 50      nd(id)=0

      do ic=1,nc
      do age=1,imin(ageic(ic),agemx)
         id=int(diica(ic,age)/deltar)
!        use upper limit to detect large distances
         if(id.gt.idmax) id=idmax
         nd(id)=nd(id)+1
      enddo
      enddo

      print*,'   write result to file pand, unit=',chpand
      do id=0,idmax
         write(chpand,8000) id*deltar,nd(id)
      enddo
 8000 format(f15.5,i10)

      print*,' '

      return
      end

!----------------------------------------------------------------------
      subroutine outxy(lonica,latica,ageic,realic,nrealc,nc)
!----------------------------------------------------------------------
!   output of trajectories seperated by 2 x seper, (=0.)
!   format: geographical coordinates (angles)
!   longitude in  -180 ... 180
!   can be smoothed bw w>0

      use xatrm
      implicit none
      
      integer       :: ic,age,nrealc,nc,maxage
      integer       :: lonica(ncmx,0:agemx)
      integer       :: latica(ncmx,0:agemx)
      integer       :: ageic(ncmx)
      integer       :: realic(ncmx)
      real(kind=4)  :: seper,xlon(0:agemx),xlat(0:agemx)
      real(kind=4)  :: xlatsm(0:agemx),xlonsm(0:agemx),w

!     weight for smoothing of coordinates, makes nicer plots
!     otherwise all cycs are on the grid points, w should by 0 .. 0.2

      w=0.15
      seper=0

      print*,'>> outxy, paths seperated by (0,0) to unit ',chxy
      print*,'   ',nrealc,'  intense (anti) cyclones'

      write(chxy,8000) seper,seper

      do ic=1,nc
         if(realic(ic).eq.1) then
!        this is the maximum age to be outputted (and stored)
         maxage=min(ageic(ic),agemx)
         do age=0,maxage
            xlat(age)=90.-180*latica(ic,age)/(nla+1.)
            xlon(age)=lonica(ic,age)*360./nlo
            if(xlon(age) .gt. 180.) xlon(age)=xlon(age)-360.
         enddo
  

!        smooth xlon and xlat, do not change end points
         do age=0,maxage
            xlatsm(age)=xlat(age)
            xlonsm(age)=xlon(age)
         enddo
  
         do 40 age=1,maxage-1
!           prevent smoothing ad dateline
            if(xlon(age  )*xlon(age-1).lt.0) goto 40
            if(xlon(age  )*xlon(age+1).lt.0) goto 40
            if(xlon(age-1)*xlon(age+1).lt.0) goto 40
            xlatsm(age)=w*xlat(age-1)+w*xlat(age+1)+(1-2*w)*xlat(age)
            xlonsm(age)=w*xlon(age-1)+w*xlon(age+1)+(1-2*w)*xlon(age)
  40     continue
!           write(chxy,8000) seper,seper

         do age=0,min(ageic(ic),agemx)
            write(chxy,8000) xlonsm(age),xlatsm(age)
         enddo

         write(chxy,8000) seper,seper
         endif
      enddo ! ic

 8000 format(2f8.2)

      print*,' '

      return
      end

!----------------------------------------------------------------------
      subroutine outtr(nc,lonica,latica,czpica,catica,ageic,
     &                 t0ic,datica,timeica,gzpica,diica,realic,
     &                 zradica,zdepica,     
     &                 tropica,nrealc,agemin,anatyp)
!----------------------------------------------------------------------

!   output of trajectories
!   in three different files, two for f. sielmanns plot routine

      use xatrm
      implicit none

      real(kind=4) ::    x,y
      real(kind=4) ::    czpica(ncmx,0:agemx)
      real(kind=4) ::    gzpica(ncmx,0:agemx)
      real(kind=4) ::    diica(ncmx,0:agemx)
      real(kind=4) ::    zradica(ncmx,0:agemx)      
      real(kind=4) ::    zdepica(ncmx,0:agemx)      
      real(kind=4) ::    tropica(ntrop,ncmx,0:agemx)

      integer
     &          nc,ic, age,agemin
     &         ,lonica(ncmx,0:agemx),latica(ncmx,0:agemx)
     &         ,t0ic(ncmx), ageic(ncmx),datica(ncmx,0:agemx)
     &         ,timeica(ncmx,0:agemx)
     &         ,catica(ncmx,0:agemx)
     &         ,realic(ncmx)
     &         ,nrealc,inum,imin,itrop,anatyp

!  output of tracks and categories

      print*,'>> outtr  '
      print*,'   output of path information to unit',chpath
      print*,'      nrealc'
      print*,'      agemin'
      print*,'      ic T0 age Lo La Ca Da Cz Di Re Gz TAg Tim Rad Dep'
      print*,' '

      write(chpath,8012) nrealc
      write(chpath,8012) agemin
 8012 format(i8)

      do ic=1,nc
      do age=0,imin(ageic(ic),agemx)
         if(realic(ic).gt.0) write(chpath,8018)
     &      ic+$zeiger,t0ic(ic),age,
     &      lonica(ic,age),latica(ic,age),
     &      catica(ic,age),datica(ic,age),
     &      czpica(ic,age), diica(ic,age), realic(ic),
     &      gzpica(ic,age), ageic(ic), timeica(ic,age),
     &      zradica(ic,age),zdepica(ic,age)
      enddo ! 2011    continue
      enddo ! 2010 continue

 8018 format(2i12,i4,2i4,i2,i9,e13.4,f7.3,i2,e13.4,i4,i7,f7.3,e13.4)


!     tropical path information
      if(anatyp.eq.8) then
         write(chptro,8012) nrealc
         write(chptro,8012) agemin

         do ic=1,nc
         do age=0,imin(ageic(ic),agemx)
            write(chptro,8019)
     &         ic,t0ic(ic),age,
     &         lonica(ic,age),latica(ic,age),
     &         catica(ic,age),datica(ic,age),
     &         czpica(ic,age), diica(ic,age), realic(ic),
     &         gzpica(ic,age), ageic(ic),
     &         (tropica(itrop,ic,age), itrop=1,ntrop),
     &         timeica(ic,age)
      enddo !  3011       continue
      enddo !  3010    continue
      endif
 8019 format(2i5,i3,2i4,i2,i7,e12.4,f6.2,i2,e12.4,i3,5e12.4,i3)


!  output in the format to plot traj
!  two files, one with grid positions, unit chpgrd
!  the other degrees, chpdgr, for f. sielmann's routine

      do ic=1,nc
         if(realic(ic).eq.1) then

            inum=1+imin(ageic(ic),agemx)
            write(chpgrd,'(i6)') inum
            write(chpdgr,'(i6)') inum

            do age=0,imin(ageic(ic),agemx)
               write(chpgrd,'(2i3,i2,i10)') lonica(ic,age),
     &            latica(ic,age),catica(ic,age),datica(ic,age)
               x=(360.0*(lonica(ic,age)-1))/nlo
               y=90.0-(180.0*latica(ic,age))/(nla+1)
               write(chpdgr,'(2f7.2)') x,y
            enddo !  2013       continue

         endif
      enddo !  2012 continue

      return
      end

!----------------------------------------------------------------------
      subroutine outlt(ageic) ! ,chpalt)
!----------------------------------------------------------------------
!  determine freque ncy distribution of life times  nage
!  to unit chpalt, file palt
      
      use xatrm
      implicit none
      
      integer   :: ic
      integer   :: ageic(ncmx),a
      integer   :: nage(0:agemx)

      print*,'>> outlt: determine life times '

      do a=0,agemx
         nage(a)=0
      enddo

      do ic=1,ncmx
         a=ageic(ic)
         if((a.ge.0).and.(a.le.agemx))  nage(a)=nage(a)+1
      enddo !  100  continue

      print*,'   result to unit ',chpalt
      print*,'   age nage(age)=number of cyclones with lifetime age '

      do a=0,agemx
         write(chpalt,8000) a,nage(a)
      enddo !  200  continue

 8000 format(2I10)

      print*,' '

      return
      end

!----------------------------------------------------------------------
      subroutine genlys(nc,lonica,latica,ageic,realic)
!----------------------------------------------------------------------
!  distribution of genesis, lysis and total density
!  genesis to chgen, lysis to chlys, total to chtot
!  variables
!     genden(i,j) frequency of starting cyclones at i,j
!     lysden(i,j) "            ending   "
!     totden(i,j) "            total    "
!  could be weighted with categories catica
      
      use xatrm
      
      implicit none

      integer      :: nc,ic, age
      integer      :: lonica(ncmx,0:agemx)
      integer      :: latica(ncmx,0:agemx)
      integer      :: ageic(ncmx)
      integer      :: realic(ncmx)
      integer      :: i,j,imin,ih(8)
      real(kind=4) :: genden(nlomx,nlamx)
      real(kind=4) :: lysden(nlomx,nlamx)
      real(kind=4) :: totden(nlomx,nlamx)

      print*,'>> genlys: genesis, lysis and total densities '
      print*,'           units:',chgen,chlys,chtot

      do 5 i=1,8
 5       ih(i)=0
      ih(3)=900001
      ih(5)=nlo
      ih(6)=nla

      do i=1,nlo
      do j=1,nla
         genden(i,j)=0
         lysden(i,j)=0
         totden(i,j)=0
      enddo
      enddo
      
!     totden(1,nla/2)=10

      do ic=1,nc
         if(realic(ic).eq.1) then
!           starts at age=0
            i=lonica(ic,0)
            j=latica(ic,0)
            genden(i,j)=genden(i,j)+1
!           ends
            age=ageic(ic)
            if(age.le.agemx) then
               i=lonica(ic,age)
               j=latica(ic,age)
               lysden(i,j)=lysden(i,j)+1
            endif
            do age=0,imin(ageic(ic),agemx)
               i=lonica(ic,age)
               j=latica(ic,age)
               totden(i,j)=totden(i,j)+1
            enddo !      33        continue
         endif
      enddo !  20   continue

      write(chgen) ih
      write(chgen) genden
      write(chlys) ih
      write(chlys) lysden
      write(chtot) ih
      write(chtot) totden

      do 41 j=1,nla
      do 41 i=1,nlo
         write(chgen+1,*) genden(i,j)
 41   continue

      do 42 j=1,nla
      do 42 i=1,nlo
         write(chlys+1,*) lysden(i,j)
 42   continue

      do 43 j=1,nla
      do 43 i=1,nlo
         write(chtot+1,*) totden(i,j)
 43   continue

      return
      end

!----------------------------------------------------------------------
      subroutine statr2(nc,lonica,latica,ageic,realic)
!----------------------------------------------------------------------
      
      use xatrm
      
      implicit none

      integer,parameter   :: cimax=40

      real(kind=4) ::    a,b,c,sum 
      real(kind=4) ::    r2(agemx),r2lat(agemx),r2lon(agemx)
     &         ,suma(agemx)
     &         ,slope,slopes,slopem
     &         ,pca(0:cimax,0:agemx)

      integer
     &          imin
     &         ,nc,ic, age,ci
     &         ,lonica(ncmx,0:agemx)
     &         ,latica(ncmx,0:agemx)
     &         ,ageic(ncmx)
     &         ,realic(ncmx)
     &         ,lonr(ncmx,0:agemx)
     &         ,latr(ncmx,0:agemx)
     &         ,dlon(ncmx,1:agemx)
     &         ,nslope

      
!----------------------------------------------------------------------
!  variables
!      dlon, dlat  shift of cyclone ic at age in grid distances
!                  jumps removed
!      lonr, latr  relative movement of cyclone ic at age
!      suma        number of cyclones with age
!      r2          mean square displacement vbs age
!      a,b,c       distance in spherical coords (bronstein)
!      slope..     exponent in r2(t)
!----------------------------------------------------------------------
!  determine distances dlon
!  remove jumps over the border at i=1
!  map  into   -nlo/2 .. nlo/2, physical distances, not periodic


      print*,'>> statr2 '

      do ic=1,nc
         if(realic(ic).eq.1) then
         do age=1,imin(ageic(ic),agemx)
            dlon(ic,age)=lonica(ic,age)-lonica(ic,age-1)
            if(dlon(ic,age).gt. nlo/2) dlon(ic,age)=dlon(ic,age)-nlo
            if(dlon(ic,age).lt.-nlo/2) dlon(ic,age)=dlon(ic,age)+nlo
         enddo !  58      continue
         endif
      enddo !  57   continue

!  determine relative paths lonr, latr in real distances
!  (not in periodic) beginning in 0,0

      do ic=1,nc
         if(realic(ic).eq.1) then
            lonr(ic,0)=0
            do age=1,imin(ageic(ic),agemx)
               lonr(ic,age)=dlon(ic,age) + lonr(ic,age-1)
               latr(ic,age)=latica(ic,age) - latica(ic,0)
            enddo !  78         continue
         endif
      enddo !  77   continue

!----------------------------------------------------------------------
!  mean square displacement r2(age)
!
!  spherical trigonometry: cos c=cos a cos b
!  a = arc in longs, b = arc in lats
!  c = arc between initial and endpoint
!  a = deltai*2*pi/nlo
!  approximately is deltaphi=5.5 degrees, for t21, roughly 180/32
!  b = deltaj*2*pi/(2*nla)
!
      print*,'   determine mean square displacements (in arc units)'

      do age=1,agemx
         r2lon(age)=0
         r2lat(age)=0
         r2   (age)=0
      enddo
!  radial
!  density p(r,t) = pca(ci,age) = sum delta(ci-c(age,ic))
!                                  ic

      do ci=0,cimax
      do age=0,agemx
         pca(ci,age)=0
      enddo
      enddo
      
      do ic=1,nc
         if(realic(ic).eq.1) then
         do age=1,imin(ageic(ic),agemx)
!           distance in spherical coords, bronstein, r2 in arc units
            a=(lonr(ic,age)*2*pi)/nlo
            b=(latr(ic,age)*pi)/(nla+1)
            c=acos(cos(a)*cos(b))
            r2   (age)=r2   (age) + c**2
            r2lat(age)=r2lat(age) + b**2
            r2lon(age)=r2lon(age) + a**2
!           determ distrib p(r,age), ci=discrete of r=c (in arcs)
            ci=int(c/0.05)
            if(ci.gt.cimax) ci=cimax
            pca(ci,age)=pca(ci,age)+1
         enddo !  93      continue
         endif
      enddo !  92   continue

!  normalize r2 by the number of cyclones with age age

      do age=1,agemx
         suma(age)=0
      enddo

      do ic=1,nc
         if(realic(ic).eq.1) then
            do age=1,imin(ageic(ic),agemx)
               suma(age)=suma(age) + 1
            enddo
         endif
      enddo 

!od      print*,'age   suma(age)'
      do age=1,agemx
         if(suma(age).ne.0)  then
            r2   (age)=r2   (age)/suma(age)
            r2lat(age)=r2lat(age)/suma(age)
            r2lon(age)=r2lon(age)/suma(age)
!            do 410 ci=0,cimax
! 410           pca(ci,age)=pca(ci,age)/suma(age)
!od            print*,age,suma(age)
         endif
      enddo !    continue

      print*,'   output of r2 r2lat r2lon and slopes to unit 33'

      do 103 age=3,agemx
         if(r2   (age)  .le.0) goto 103
         if(r2   (age-1).le.0) goto 103
         if(r2lat(age)  .le.0) goto 103
         if(r2lat(age-1).le.0) goto 103
         if(r2lon(age)  .le.0) goto 103
         if(r2lon(age-1).le.0) goto 103
         slope=(log(r2(age))-log(r2(age-1)))
     &          /(log(float(age))-log(age-1.))
         nslope=nslope+1
         slopes=slopes+slope
         slopem=slopes/nslope
         write(chpar2,1032) log(float(age)),log(r2(age)),
     &                  log(r2lat(age)),log(r2lon(age)),slope,slopem
 103  continue
      
 1032 format(6f11.4)

!   output of pca

      print*,'   write p(r,t) (=pca(ci,age)) to unit ',chppca

!     normalize pca for each age  (int dr p(r,t)=1)
      do 680 age=0,agemx
         sum=0
!        determine sum over r
         do 682 ci=0,cimax
 682        sum=sum+pca(ci,age)
         if(sum.gt.0) then
            do 684 ci=0,cimax
 684           pca(ci,age)=pca(ci,age)/sum
         endif
 680  continue

      do 700 ci=cimax,0,-1
 700     write(chppca,8700) (pca(ci,age),age=0,agemx)
 8700 format(30f8.2)

      print*,' '

      return
      end

!----------------------------------------------------------------------
      subroutine maxreg(f,i,j,sr,fmax)
!----------------------------------------------------------------------
!  max of field f in region around i,j
!  sr = latitude grid width at equator
     
      use xatrm
      implicit none
      
      integer      :: di,dj,mi,mj,i,j,i0,j0,sr
      integer      :: f(nlomx,nlamx) 
      real(kind=4) :: dphi,phi,fmax

!     search range depends on resolution
!     longit range is larger in higher latitudes
!     latit. range is half of longit.

      dphi=pi/(nla+1)
      phi =pi/2-j*dphi
      mi  =int(0.5+ sr/abs(cos(phi)) )
      if(mi.le.0) mi=1
      mj  =sr

      fmax=-1e30
      do 1 di=-mi,mi
      do 1 dj=-mj,mj
!        exclude boundary crossing
         if( ((j+dj).lt.1).or.((j+dj).gt.nla) ) goto 1
         i0=p(i+di)
         j0=j+dj
         if(f(i0,j0).gt.fmax) fmax=f(i0,j0)
 1    continue

 9    return
      end

!----------------------------------------------------------------------
      subroutine search(cat,i,j,i0,j0,answer,sr,dt)
!----------------------------------------------------------------------
!  search cyclone in cat around i,j.
!  size of region i=-mi...mi, j=-mj...mj, depends on nt
!  if found, answer=1, store former position in i0,j0
      
      use xatrm
      implicit none
      
      integer   :: di,dj,mi,mj,i,j,i0,j0,answer,sr,dt
      integer   :: cat(nlomx,nlamx) ! ,p(-nlomx:2*nlomx)

      real(kind=4) ::    dphi,phi,srmx_r

!     changed, date 23.2.98
!     this variable was introduced to restrict the search region
!     after problems occured in T42/4h with a too large region

!     srmx  =int(0.5+nlo/64*(dt/6.)**0.75)

!     search range depends on resolution
!     longit range is larger in higher latitudes
!     latit. range is half of longit.

!     this is the maximum for srmx_r, for smaller shells, srmx_r=sr
      srmx_r=nlo/64*(dt/6.)**0.75
      if(srmx_r.gt.sr) srmx_r=sr

      dphi=pi/(nla+1)
      phi =pi/2-j*dphi

      mi  =int(0.0+srmx_r/abs(cos(phi)) )
      if(mi.eq.0) mi=1

!     search region is now quadratic:

      mj  =sr
      if((nlo.le.64).and.(sr.eq.1)) mj=1

      answer=0

      do 1 di=-mi,mi
      do 1 dj=-mj,mj
!        exclude boundary crossing
         if( ((j+dj).lt.1).or.((j+dj).gt.nla) ) goto 1
         i0=p(i+di)
         j0=j+dj
         if(cat(i0,j0).gt.0) then
!           cyclone found in cat, position (i0,j0) stop here
            answer=1
            write(*,80) dt,sr,mi,mj,srmx_r
!           the first cyc found is considered to be  the right one,
!           skip remainder
            goto 9
         endif
 1    continue

  80  format(4i6,f12.3)

 9    return
      end
!----------------------------------------------------------------------
      integer function imin(a,b)
!----------------------------------------------------------------------
      implicit none
      integer   :: a,b

      if(b.lt.a) then
         imin=b
      else
         imin=a
      endif

      return
      end
!----------------------------------------------------------------------
      real function dsmall(i,i0,j,j0)
!----------------------------------------------------------------------
      use xatrm
      
      implicit none
      real(kind=4) ::    phi,dphi,dx,dy
      integer   :: i,i0,j,j0,di,dj

      dphi=pi/(nla+1)
      dy  =r0*dphi

!  approximation for the distance for small i-i0, j-j0

      phi=pi/2-j*dphi
      dx =abs(cos(phi))*2*pi*r0/nlo
      di =i-i0
      if(di.gt. nlo/2) di=di-nlo
      if(di.lt.-nlo/2) di=di+nlo
      dj=j-j0
      dsmall=sqrt((di*dx)**2 + (dj*dy)**2)

      return
      end

!----------------------------------------------------------------------
      logical function insr(i,j)
!----------------------------------------------------------------------
!  determines whether i,j is in the region where cyclones are searched


      use xatrm
      implicit none
      integer   :: i,j 

!  search region northern mid-lat: j=nla1=1+nde,...,nla2=nla/3
!                atlantic/europe:
!
!       i<nlo1  30e                or            280e  i>nlo2  360e
!     !----------!--------------------------------!--------------!
!    i=0        nlo1                             nlo2          i=nlo

!      if( ((i.le.nlo1) .or. (i.ge.nlo2)) .and.
!     &    ((j.gt.nla1).and. (j.le.nla2)) ) then
!         insr=.true.
!      else
!         insr=.false.
!      endif

!     global search, everywhere

C      if( ((i.ge.nlo1).and. (i.le.nlo2)) .and.
C     &    ((j.ge.nla1).and. (j.le.nla2)) ) then
C         insr=.true.
C      else
C         insr=.false.
C      endif

      if(nlo1.lt.nlo2) then
C        inside nlo1..nlo2
         if(  ((j.ge.nla1).and.(j.le.nla2)) .and.
     &        ((i.ge.nlo1).and.(i.le.nlo2)) ) then
            insr=.true.
         else
            insr=.false.
         endif
      else
C        outside nlo1..nlo2 (crosses lambda=0)
         if(  ((j.ge.nla1).and.(j.le.nla2)) .and.
     &        ((i.ge.nlo1).or. (i.le.nlo2)) ) then
            insr=.true.
         else
            insr=.false.
         endif
      endif

!     eliminate mediterranean sea in T106

C      if(i.ge.0 .and. i.lt.60 .and. j.gt.36) insr=.false.

      return
      end
!----------------------------------------------------------------------
      subroutine datdmy(date,day,month,year)
!----------------------------------------------------------------------
!     extracts date, form yymmdd

      implicit  none
      integer   :: date,day,month,year

      year =int(date/10000)
      month=mod(date/100,100)
      day  =mod(date,100)

      return
      end
!----------------------------------------------------------------------
      subroutine datjmp(step,olddat,newdat,time,jump)
!----------------------------------------------------------------------
!     jump=1 if newdat not diercetly after olddat
!     only month is checked, not day

      implicit  none
      integer   :: step,olddat,newdat,jump,oldday,oldmon,oldyr
      integer   :: newday,newmon,newyr,time

      call datdmy(olddat,oldday,oldmon,oldyr)
      call datdmy(newdat,newday,newmon,newyr)

      jump=0

!     go back if at beginning and if time not 0 (otherwise 2xjump!)
      if(step .eq.1) return
      if(time.ne.0) return

      if(newday.ne.1) then
         return
      else
!        for day=1:
         if(newmon.eq.1) then
            if(oldmon.ne.12) jump=1
         else
            if(newmon.ne.(oldmon+1)) jump=1
         endif
      endif

      return
      end
      
!----------------------------------------------------------------------
      subroutine fitgauss(zp,i,j,zrad,zdep,gzp)
!----------------------------------------------------------------------
!  routine in xatr  RB 18.9.2006
!
!  fit a Gaussian to zp in neighborhood of cyc in i,j
!
!  minimise  D2 = 1/2 sum_ij ( z(i,j) - zgauss(r(i,j)) )^2
!
!  zgauss(r)=zenv-(zenv-zcen)*exp(-r^2/2*zrad^2)
!
!  centr:      zgauss(r=0) = zcen this is not fitted, zp(i,j)
!
!  depth:      zenv-zcen
!
!  environm    zgauss(r->inf) = zenv
!
!  radius      zrad
!
!
!  store result in
!      real(kind=4) ::    zradica(ncmx,0:agemx)
!      real(kind=4) ::    zdepica(ncmx,0:agemx)
!
!     call everywhere with
!        call gradzp(zp,i,j,gzp)  is called
!
!  works possibly also with vorticity
!
!     radius   r is by:
!     dist=sqrt( (di*dx)**2 + (dj*dy)**2 )
!
!     periodic function in i:  p(i+di)
!
!  phi latitude (in radians) determined approximately
!
!  r0 radius of earth, in 1000km  (r0=6.370)
!  dx, dy distance along long/lat between i,j-points in 1000km
!  in latitude simply dphi=const
!  the range is determined by nt
!
!
!  for fit
!  =======
!  mi,mj range for fit area,  1 grid point in T21 for latit.
!  could be larger
!
!  deth zdep is known: zp in minimum i,j
!
!  parameters in min search 
!
!     1  zenv
!     2  zrad
!
!  be careful with f77 and f90, real*4 and real*8
!  
!  ===========================================
!  overview to fit radii and depth to cyclones
!  ===========================================
!  call fitgauss
!  
!  in ancat0, 1
!  
!              if(anatyp.eq.1) then
!                 czpica(nc,age)=zp(i,j)
!                 call gradzp(zp,i,j,gzp)
!                 gzpica(nc,age)=gzp
!                 call fitgauss(...)
!                 write(*,88)' ancat0 fitg',i,j,zdep,zrad,zp(i,j),gzp
!    88           format(A,2i4,4f12.3)                              
!  c              zdepica(nc,age)=zdep
!  c              zradica(nc,age)=zrad
!              endif
!          
!  produces:   zdep,zrad
!          
!  fitgauss calls
!          
!  findpars(pp_fit,n_fit,np_fit,zenv_max_fit,zrad_max_fit)
!  
!     pp(1) zenv
!     pp(2) zrad
!          
!  findpars calls Num Rec routine
!     
!  powell(pp,xi,n_fit,np_fit,ftol,iter,fret)
!-----------------------------------------------------------------------
      
      use xatrm
      
      implicit none

      integer      :: izrad,Nzrad
      integer      :: di,dj,mi,mj,i,j,Nzmean
      
      integer,parameter :: n_fit=2,np_fit=2
      
      real(kind=4) :: zp(nlomx,nlamx)
      real(kind=4) :: phi,dphi,dist,dx,dy,nt
      real(kind=4) :: zenv,zdep,zcen,zrad,pp_fit(np_fit),zgauss
      real(kind=4) :: zmean,gzp,dz_env_mean
      real(kind=4) :: a,b,D
           
      ! vars for function func in powell
      
      zp_fit(:,:)=zp(:,:)
      nlo_fit=nlo
      nla_fit=nla
      i_fit=i
      j_fit=j
      
      nt    =nlo/64
      dphi  =pi/(nla+1)
      phi   =pi/2-(j)*dphi      
      dx    =abs(cos(phi))*2*pi*r0/nlo
      dy    =r0*dphi

      ! output to see the behaviour of zp
      ! and determine zmean
      
      if(zp(i,j)<zp_max_fit) then ! fit only intense lows
      
         !if(zp(i,j)<-1500 .and. zp(i,j)>-1540) then ! test1
      
         mi=int(1.2*nt/abs(cos(phi)) )
         mj=1.2*nt
         if(mi.lt.1) mi=1
         if(mj.lt.1) mj=1

         if(j-mj.lt.1 .or. j+mj.gt.nla) then
         ! exclude boundaries
         else
            
            ! zmean is the average of zp in the fit range
            ! zmean is approx 550 deeper than zenv
            
            dz_env_mean=550

            Nzmean=0
            zmean=0
            do di=-mi,mi
            do dj=-mj,mj
               dist=sqrt( (di*dx)**2 + (dj*dy)**2 )
               if(dist<dist_fit) then
               !  output to see the behaviour of zp
                  if(do_fit_out) 
     &              write(781,'(3f10.3)') dist,zp(p(i+di),j+dj)
               endif
               Nzmean=Nzmean+1
               zmean=zmean+zp(p(i+di),j+dj)
            enddo
            enddo
            if(do_fit_out) then
               write(782,*) ' mi mj ',mi,mj
               do dj=-mj,mj
                  write(782,'(100g12.4)') (zp(p(i+di),j+dj),di=-mi,mi)
               enddo
            endif ! do_fit_out
         endif
         if(do_fit_out) write(781,*)'//nc'
      endif 
      
      zmean=zmean/Nzmean


!     fit of a Gaussian in (i,j), determine: zenv,zrad
      
      if(zp(i,j)<zp_max_fit) then

         ! test1 if(zp(i,j)<-1500 .and. zp(i,j)>-1540) then
      
         ! initial values for optimization

         ! zenv:
         
         pp_fit(1)=zmean ! +dz_env_mean ! zp(i,j) ! zp_max_fit
         
         ! init zrad by climate formula using zdep
         
         zcen=zp(i,j)  
         zdep=zmean-zcen
         
         !if(zdep<=0) then
         !   zdep=0
         !   zrad=0
         !endif
         
         if(zdep<=0) then
            zrad=0.5
         else
            a=2628.
            b=1.37
            zrad=(zdep/a)**(1./b)
         endif

         pp_fit(2)=zrad
         
         
         call findpars(pp_fit,n_fit,np_fit)
         
         N_call_fitgauss=N_call_fitgauss+1
         
         zenv=pp_fit(1)
         
         ! correct wrong optimised values
         
         if(zenv<zmean .or.zenv>zenv_max_fit) then
            zenv=zmean 
            N_corr_zenv=N_corr_zenv+1
         endif
         
         zcen=zp(i,j)         
         
         zdep=zenv-zcen

         if(zdep<0) then
            zdep=0
            N_corr_zdep=N_corr_zdep+1
         endif
         
         zrad=pp_fit(2)

         if(zrad<0 .or.zrad>=zrad_max_fit) then
            ! ! power-law
            ! a=2628.
            ! b=1.37
            ! zrad=(zdep/a)**(1./b)  ! zdep is >=0
            
            ! if(zrad>=zrad_max_fit) zrad=0.36+0.0308*gzp**(-0.059)
            
            zrad=0
            
            zdep=0
            
!            linear: zrad=zdep/1703  
!            or:
!            D=zdep
!            a=1102
!            b=4755
!            zrad = (SQRT((27*b*D
!     &      **2+4*a**3)/b)/(SQRT(3.)*b)/6.0+D/b/2.0)**(1.0/3.0)
!     &      -a*(SQRT((27*b*D**2+4*a**3)/b)/
!     &      (SQRT(3.)*b)/6.0+D/b/2.0)**((-1.0)/3.0)/b/3.0

            N_corr_zrad=N_corr_zrad+1
         endif

         ! output of fitted Gaussian curves
         
         if(zrad>0) then
            Nzrad=50
            do izrad=0,Nzrad
               dist=0.05*izrad
               zgauss = zenv-(zenv-zcen)*exp(-dist**2/(2*zrad**2))
               if(do_fit_out) write(785,'(2f10.3)') dist,zgauss
            enddo
         endif
         
         if(do_fit_out) then
            write(785,*)'//nc'

            write(801,*) zenv
            write(802,*) zrad
            write(803,*) zdep
            write(804,*) zmean
            write(805,*) zcen
            write(806,*) gzp

            write(819,*) zrad,zdep         
            write(820,*) gzp,zdep
            write(821,*) gzp,zrad
            write(822,*) zmean,zenv
            write(823,*) zmean,zdep
            write(824,*) zcen,gzp
            write(825,*) zcen,zrad
            write(826,*) zcen,zdep

            write(840,*) zmean-zcen ! should be positive
            write(841,*) zenv-zmean
         endif
      else
         zenv=zmean ! misv_fit
         zcen=zp(i,j)
         zrad=misv_fit
         zdep=zmean-zcen ! misv_fit
      endif

      return
      end      

!-----------------------------------------------------------------------
      subroutine findpars(pp,n_fit,np_fit)
!-----------------------------------------------------------------------

      use xatrm
      
      implicit none
      integer      :: n_fit,np_fit
      real(kind=4) :: pp(np_fit),ppout(np_fit)
      real(kind=4) :: xi(np_fit,np_fit),ftol,fret
      real(kind=4) :: func
      integer      :: iter,l,istep,Nstep
      
      ftol=0.00000001
      Nstep=1 !

      ppout(:)=pp(:)

      do istep=1,Nstep

         xi(:,:)=0
         do l=1,np_fit
            xi(l,l)=1.0
         enddo

         call powell(pp,xi,n_fit,np_fit,ftol,iter,fret)
         
         ! correct iterated values
         
         if(pp(2)<0) pp(2)=abs(pp(2))
         
         ppout(:)=pp(:)
         
         if(do_fit_out) write(784,'(A,i4,5e12.4)') 
     &                ' pp ',istep,(ppout(l),l=1,n_fit)
         if(do_fit_out) write(791,'(i9)') iter
         if(do_fit_out) write(792,'(5e12.4)') fret

      enddo ! istep

      return
      end

!-----------------------------------------------------------------------
      real(kind=4) function func(pp)
!-----------------------------------------------------------------------
!      
!     function D2 to minimise
!
!     D2 = sum_ij ( z(i,j) - zgauss(r(i,j)) )^2
!
!     zgauss(r)=zenv-(zenv-zcen)*exp(-r^2/2*zrad^2)
!
!                  |->  zrad=std dev
!     zenv ---         ----
!              \      /
!               \    /
!                ---   zcen=zp(i,j) fix
!
!     depth:      zenv-zcen
!     use all zp values in fit area
!     pp are fit parameters in powell

      use xatrm
      
      implicit none
      
      integer,parameter :: n_fit=2,np_fit=2
      integer         :: mi,mj,i,j,di,dj
      integer         :: nt

      real(kind=4) :: pp(np_fit),r,zgauss,zenv,zcen,zrad
      real(kind=4) :: D2,dx,dy
      real(kind=4) :: dphi,phi
      real(kind=4) :: zp(nlomx,nlamx)

      ! these are the local names, determined by global vars
      
      nlo=nlo_fit
      nla=nla_fit
      zp(:,:)=zp_fit(:,:)
      i=i_fit
      j=j_fit
           
      
      nt    =nlo/64
      dphi  =pi/(nla+1)
      phi   =pi/2-j*dphi      
      dx    =abs(cos(phi))*2*pi*r0/nlo
      dy    =r0*dphi

      zenv=pp(1) 
      zrad=pp(2) 
      
      D2=0
      
      !  fit area, range 

      mi=int(1.2*nt/abs(cos(phi)) )
      mj=1.2*nt
      if(mi.lt.1) mi=1
      if(mj.lt.1) mj=1
      
      zcen=zp(i,j)
      
      if(j-mj.lt.1 .or. j+mj.gt.nla) then
!        
      else
         do di=-mi,mi
         do dj=-mj,mj
            r=sqrt( (di*dx)**2 + (dj*dy)**2 )
            
            ! fit in distfit-area
            if(r<dist_fit) then 
               zgauss = zenv-(zenv-zcen)*exp(-r**2/(2*zrad**2))
               D2=D2+(zp(p(i+di),j+dj) -zgauss)**2
            endif
    
         enddo
         enddo
      endif 

      func=D2

      return
      end

!-----------------------------------------------------------------------
!     Numerical Recipes routines for optimization
!-----------------------------------------------------------------------

      SUBROUTINE powell(p,xi,n,np,ftol,iter,fret)

      implicit none

      INTEGER iter,n,np,NMAX,ITMAX
      REAL(kind=4) :: fret,ftol,p(np),xi(np,np),func
      EXTERNAL func
      PARAMETER (NMAX=20,ITMAX=50)  ! ITMAX was 200
CU    USES func,linmin
      INTEGER i,ibig,j
      REAL(kind=4) :: del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX)

      !write(77,*)' powell ',p
      !write(77,*)' powell ',xi

      fret=func(p)
      do 11 j=1,n
        pt(j)=p(j)
11    continue
      iter=0
1     iter=iter+1
      fp=fret
      ibig=0
      del=0.
      do 13 i=1,n
        do 12 j=1,n
          xit(j)=xi(j,i)
12      continue
        fptt=fret
        call linmin(p,xit,n,fret)
        if(abs(fptt-fret).gt.del)then
          del=abs(fptt-fret)
          ibig=i
        endif
13    continue
      if(2.*abs(fp-fret).le.ftol*(abs(fp)+abs(fret)))return
      ! if(iter.eq.ITMAX) pause 'powell exceeding maximum iterations'
      do 14 j=1,n
        ptt(j)=2.*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
14    continue
      fptt=func(ptt)
      if(fptt.ge.fp)goto 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0.)goto 1
      call linmin(p,xit,n,fret)
      do 15 j=1,n
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
15    continue
      goto 1
      END

!-----------------------------------------------------------------------

      SUBROUTINE linmin(p,xi,n,fret)

      implicit none

      INTEGER n,NMAX
      REAL (kind=4) :: fret,p(n),xi(n),TOL
      PARAMETER (NMAX=50,TOL=1.e-4)
CU    USES brent,f1dim,mnbrak
      INTEGER j,ncom
      REAL (kind=4) :: ax,bx,fa,fb,fx,xmin,xx
      REAL (kind=4) :: pcom(NMAX),xicom(NMAX),brent
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim

      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.
      xx=1.
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END


!-----------------------------------------------------------------------

      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)

      implicit none

      REAL (kind=4) :: ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
      REAL (kind=4) :: dum,fu,q,r,u,ulim

      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END


!-----------------------------------------------------------------------

      FUNCTION f1dim(x)

      implicit none

      INTEGER NMAX
      REAL (kind=4) :: f1dim,func,x
      PARAMETER (NMAX=50)
CU    USES func
      INTEGER j,ncom
      REAL (kind=4) :: pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      END

!-----------------------------------------------------------------------

      FUNCTION brent(ax,bx,cx,f,tol,xmin)

      implicit none

      INTEGER ITMAX
      REAL (kind=4) :: brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
      EXTERNAL f
      PARAMETER (ITMAX=50,CGOLD=.3819660,ZEPS=1.0e-10) ! ITMAX was 100
      INTEGER iter
      REAL (kind=4) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r
      REAL (kind=4) :: tol1,tol2,u,v,w,x,xm


      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x))
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
       endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      ! pause 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      END
EOF

#pgf90 -o xatrs_CCSM_0.9x1.25.exe xatrs_CCSM_0.9x1.25.f

#gfortran -o xatrs_CCSM_0.9x1.25.exe xatrs_CCSM_0.9x1.25.f
#./xatrs_CCSM_0.9x1.25.exe
#/bin/mv fort.30 test_$yy


#/bin/mv fort.30 fort_30_${arr[1]}"_"${arr[10]}

# later step 


cat > extract_slp_prec.f90 <<EOF
      program extract

      implicit none
      integer  ncmx
      parameter(ncmx=600000)
      integer ntime
	integer count
      parameter(ntime=13140)
      integer i,j,k,nn,month(ncmx,6),year(0:ncmx,6),datloc,jj(6),ii
      integer ih(8),islp(ntime),ihour(ntime),ff
      integer  iic,t0,age,ilon,ilat,cat,date,realc,agetot,idumi,tmpvar
      integer tmpvar2,lon5,lon6,lat5,lat6,lat7,kk
      integer agemin,lon1,lon2,lon3,lon4
      real dumi1,dumi2,sum,mean,lon,lat, zrad, zdep
      real     cz,di,gz 
      real dec1,quart1,med,quart2,dec2,rprec(144,96,ntime)
      real slp(144,96),prec(144,96),rslp(144,96,ntime),precmmean,precmmax
      real ws(144,96),rws(144,96,ntime),wsmmean,wsmmax
      real q(144,96), rq(144,96,ntime),qmmean,qmmax
	real precc(144,96),rprecc(144,96,ntime),preccmmean,preccmmax
	real precl(144,96),rprecl(144,96,ntime),preclmmean,preclmmax

      Open (30,File="fort_30_${arr[1]}_${arr[10]}")
      Open (31,File='PSL.srv',FORM='UNFORMATTED',ACCESS='SEQUENTIAL',status='OLD')
      Open (32,File="PRECT.srv",FORM='UNFORMATTED',ACCESS='SEQUENTIAL',status='OLD')
      Open (33,File="WS.srv",FORM='UNFORMATTED',ACCESS='SEQUENTIAL',status='OLD')
      Open (34,File="Q.srv",FORM='UNFORMATTED',ACCESS='SEQUENTIAL',status='OLD')
      Open (35,File="PRECL.srv",FORM='UNFORMATTED',ACCESS='SEQUENTIAL',status='OLD')

      k=0
	count=0

 900  read(31,end=997) ih
      read(32) ih
      read(33) ih
	read(34) ih
	read(35) ih
!      print*,ih
      read(31) slp
      read(32) prec
      read(33) ws
	read(34) q
	read(35) precl
	  k=k+1
      islp(k)= ih(3)
      ihour(k) = ih(4)
!      print*,k,islp(k)
      do i=1,144
         do j=1,96
            if (prec(i,j).lt.0.) prec(i,j)=0.0
!            if (precc(i,j).lt.0) precc(i,j)=0.0
            if (precl(i,j).lt.0) precl(i,j)=0.0
            if (ws(i,j).lt.0) ws(i,j)=0.0
            if (q(i,j).lt.0) q(i,j)=0.0
            rslp(i,j,k) = slp(i,j)
            rprec(i,j,k) = prec(i,j)
            rprecl(i,j,k) = precl(i,j)
!            rprecc(i,j,k) = precc(i,j)
            rws(i,j,k) = ws(i,j)
            rq(i,j,k) = q(i,j)
            
!           print*,i,j,k,prec(i,j), rprec(i,j,k)
         end do
      enddo
      goto 900
997   continue
!     print*,"k",islp
      k=0
      i=0
      j=0
      
      read(30,*) dumi1
      read(30,*) dumi2
      write(36,*) dumi1
      write(36,*) dumi2
      print*,"OK"
800   read(30,*,end=998) iic, t0, age, ilon, ilat, cat, date,cz, di, realc, gz, agetot , idumi, zrad, zdep
!      print*,iic, t0, age, ilon, ilat, cat, date,cz, di, realc, gz, agetot , idumi, zrad, zdep
!      print*,date,islp(i)
      do i=1,ntime
!         print*,i,date,islp(i),idumi,ihour(i)
         if (date.eq.islp(i).and.idumi.eq.ihour(i)) then
         print*,i,islp(i)
         k=i
!         print*,k

         endif
      enddo
!            print*,k,date,islp(k)
	zrad = zrad*1.5
    do i=1,144
         do j=1,96
            ws(i,j) = rws(i,j,k)
            prec(i,j) = rprec(i,j,k) *86400.*1000.
            precl(i,j) = rprecl(i,j,k) *86400.*1000.
!            precc(i,j) = rprecc(i,j,k) *86400.*1000.
            q(i,j) = rq(i,j,k)

 !           if (prec(i,j).lt.0.) print*,i,j,k,prec(i,j), rprec(i,j,k) 
         end do
      enddo

	call varstats(ws,ilon,ilat,zrad,wsmmean,wsmmax)
      call varstats(prec,ilon, ilat, zrad, precmmean,precmmax)
	call varstats(precl,ilon,ilat,zrad,preclmmean,preclmmax)
      call varstats(q,ilon,ilat,zrad,qmmean,qmmax)

!	print*,precmmean,precmmax,preccmmean,preclmmean,wsmmean,wsmmax
	
      write(36,777) iic, t0, age, ilon, ilat, cat, date,cz, di, realc, gz, &
                    agetot , idumi, zrad, zdep, rslp(ilon, ilat,k)/100., precmmean, &
			  precmmax,qmmean,preclmmean,wsmmean,wsmmax
!      print*,ilon, ilat,k,rslp(ilon, ilat,k)/100.
!      write(35,*)cz
!      write(36,*)rslp(ilon, ilat,k)
!	count = count + 1
!	if (count.eq.50) then
!	stop
!	end if
777   format(2i12,i4,2i4,i2,i9,e13.4,f7.3,i2,e13.4,i4,i7,f8.3,8e13.5) 
      goto 800
     
998   continue

      end program extract

      subroutine varstats(svar,ilon, ilat, zrad, varmean, varmax)
      implicit none
      integer i,j,k,ilon, ilat,l,m, icount
      integer ni, nj,nloout, nlaout
      real rlon,rlat,pi
      real svar(144,96), zrad, varmean

      real  dphi,dlam,lonxy , latxy, radxy 
      real varmax
	nloout =144
      nlaout =96

      dlam= 360./real(nloout)
      dphi=180./real(nlaout)

      pi         = 4 * atan(1.0)
      rlon = ilon/144. *360. 
      rlat = 90. - (ilat-1) *180./nlaout
      call radiuscount(ilon, ilat, rlon,rlat, zrad, nloout, nlaout, ni, nj)
!      print*,ilon, ilat, ni,nj,  rlon, rlat 

      if ((ilat-nj .le. 0) .or. (ilat+nj .ge. nlaout)) then
         nj = 0
      endif
      varmean = 0.0
	varmax = 0.0
      icount =0
      

 !     print*,"HA", ilon-ni
      if ((ilon-ni .le. 0)) then
!         lon = longic(ic,age)
!         lat = latiic(ic,age)
!         rcrit = zradic(ic,age)

!        print*,"test"
      
         do l = ilon-ni+nloout, nloout
             do m = ilat-nj, ilat+nj
                 lonxy = rlon+abs(l-ilon)*dlam
                 latxy = rlat+abs(m-ilat)*dphi
!                 print*,"TESTLATLON"
!                 print*,latxy,lonxy
                 call sphdis(rlon,lonxy,rlat,latxy,radxy)
                 if (radxy .lt. zrad) then
                    varmean =varmean + svar(l, m)
                    if (varmax .lt. svar(l,m)) then
				varmax = svar(l,m)
			  end if
			   icount = icount +1
                 end if
             end do 
         end do
         do l = 1, ilon+ni
             do m = ilat-nj, ilat+nj
                 lonxy = rlon+abs(l-ilon)*dlam
                 latxy = rlat+abs(m-ilat)*dphi
!                 print*,"TESTLATLON"
!                 print*,latxy,lonxy
                 call sphdis(rlon,lonxy,rlat,latxy,radxy)
                 if (radxy .lt. zrad) then
                    varmean =varmean + svar(l, m)
                    if (varmax .lt. svar(l,m)) then
                        varmax = svar(l,m)
      		  end if
			   icount = icount +1
			  
!		    print*,"hello world" 
                 end if
             end do 
         end do



      else
         if (ilon+ni .ge. nloout) then 
!               print*,"test2"

               do l = ilon-ni, nloout
                  do m = ilat-nj, ilat+nj
                     lonxy = rlon+abs(l-ilon)*dlam
                     latxy = rlat+abs(m-ilat)*dphi
!                 print*,"TESTLATLON"
!                 print*,latxy,lonxy
                     call sphdis(rlon,lonxy,rlat,latxy,radxy)
                 if (radxy .lt. zrad) then
                    varmean =varmean + svar(l, m)
                    if (varmax .lt. svar(l,m)) then
                        varmax = svar(l,m)
			  end if
			   icount = icount +1
			   
                 end if
                  end do 
               end do
               do l = 1, ilon+ni-nloout
                  do m = ilat-nj, ilat+nj
                     lonxy = rlon+abs(l-ilon)*dlam
                     latxy = rlat+abs(m-ilat)*dphi
!                 print*,"TESTLATLON"
!                 print*,latxy,lonxy
                     call sphdis(rlon,lonxy,rlat,latxy,radxy)
                 if (radxy .lt. zrad) then
                   varmean =varmean + svar(l, m)
                    if (varmax .lt. svar(l,m)) then
                        varmax = svar(l,m)
			  end if
			  icount = icount +1
			  
                 end if
                  end do 
               end do
            else
!               print*,"test3"
!        the above lines were just to make sure not to overshoot bdrs
!        boundaries, now for the real stuff ...
!         lon = longic(ic,age)
!         lat = latiic(ic,age)
!         rcrit = zradic(ic,age)
         do l = ilon-ni, ilon+ni
             do m = ilat-nj, ilat+nj
                 lonxy = rlon+abs(l-ilon)*dlam
                 latxy = rlat+abs(m-ilat)*dphi
!                 print*,"TESTLATLON"
!                 print*,latxy,lonxy
                 call sphdis(rlon,lonxy,rlat,latxy,radxy)
                if (radxy .lt. zrad) then
!			  print*,varmean
!              print*,varmax
!   			  print*,svar(l, m) * sin(pi*(m)/nlaout)
                    varmean =varmean + svar(l, m)
                    if (varmax .lt. svar(l,m)) then
                        varmax = svar(l,m)

!                    print*,l,m,precm
			  end if
                    icount = icount +1
!			  print*,icount
                 end if     
!                 if (radxy .lt. zrad) print*,"Precm",ilon, ilat, precm
             end do 
         end do
! ccr edits       
         endif
         endif
         varmean = varmean/icount
!         print*,"icount",icount,precm,sprec(ilon, ilat),ilon, ilat
      return
    stop
  end
	
      subroutine radiuscount(icenter, jcenter, lon, lat, radius, nlo, nla, ni, nj)
!     New subroutine in 2016: convert gaussian (anti-)cyclone radius to 
!     square of grid points and increase the total counts by one within this square.

      implicit none
    
     
      integer      nlo, nla
      integer      i,j
      integer      icenter,jcenter,ni,nj
      real         lon, lat
      real         radius
      real         pi, rearth, circ, dlam, dphi
      real         latrd, dislat, rlat, dislon
      
!      print*,nlo,nla
!      print*,"RADIUS"
!      print*,radius

      pi         = 4 * atan(1.0)
      rearth     = 6.370
      circ       = 2 * pi * rearth
      latrd      = lat/360. * 2 * pi               !latitude in radian  
      dlam       = 360./real(nlo)
      dphi       = 180./real(nla)
!       dislat     = circ/360. * dphi                !distance 1 gridpoint lat
!       rlat       = rearth * cos(latrd)      !radius at given latitude
!       dislon     = rlat/360. * dlam                !distance one gridpoint lon at given latitude
!     new version with a call to subroutine sphdis, more accurate?
!      print*,dlam,dphi
      call sphdis(lon,lon+dlam,lat,lat,dislon)
      call sphdis(lon,lon,lat,lat+dphi,dislat)

!      print*,dislon,dislat
!     now calculate number of gridpoints with above parameter
      ni = ceiling(radius/dislon)
      nj = ceiling(radius/dislat)
      

      end subroutine


!----------------------------------------------------------------------
      subroutine sphdis(long1,long2,lati1,lati2,distan)
!----------------------------------------------------------------------

!     spherical distance
!     input:  2 points
!             longitude and latitude given in geographical coordinates
!     output: distance in 1000km              ------------
!     use: (Bronstein p. 259)
!           cos c = cos a cos b  + sin a sin b cos gamma    (radian)
!                                                            ------
!     a=distance of point 1 to North Pole, b of 2 to NP.
!     gamma=longitude between 2 and 1
!
!                       . x  lati2
!                     .   |
!                  (c)    |
!                .        |
!               x---------x  lati1
!             long1     long2

      real     long1,long2,lati1,lati2,distan
      real     pi,rearth,a,b,gamma,argum,small

      integer  chin,chin2,chlog,chden,easypl,choro
      common/ch/chin,chin2,chlog,chden,easypl,choro


      small=1d-10
      pi=4*atan(1.0)
      rearth=6.370

      a=pi*(90-lati1)/180.
      b=pi*(90-lati2)/180.
      gamma=2*pi/360.*(long2-long1)

      argum=cos(a)*cos(b) + sin(a)*sin(b)*cos(gamma)

      if(abs(argum).gt.1D0) print*,'argum>1:',argum
      distan=rearth*acos(argum)

      if(distan.lt.-small) then
 !        write(chlog,*) 'xx distan<',real(small),real(distan)
 !        write(*    ,*) 'xx distan<',real(small),real(distan)
      endif
      if(distan.lt.0) distan=0

      return
      end

EOF


gfortran -mcmodel=medium -o  extract_psl_prec.exe extract_slp_prec.f90
./extract_psl_prec.exe




/bin/mv fort.36 fort_36_${arr[1]}"_"${arr[10]}

#/bin/mv fort.36 fort_36_test
#/bin/mv ws.txt ws_${arr[1]}"_"${arr[10]}.txt


#/bin/rm tmpz1000.srv
#/bin/rm PSL.srv
#/bin/rm PRECT.srv
# /bin/rm oro.srvy

#set zeiger = `/usr/bin/expr $zeiger + 1000000`
set y = `/usr/bin/expr $y + 10` 
echo $y $zeiger
end

