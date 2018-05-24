function grid_core,filename,pix_size,T,cls_dist=cls_dist,dp=dp,h=h,r_pix_lim=r_pix_lim

; PURPOSE:
;    Use the gravitational potential from surface density to find cores
;
; CALLING SEQUENCE:
;    output=grid_core(filename,pix_size,T,cls_dist=cls_dist,dp=dp,h=h,r_pix_lim=r_pix_lim)
;
;    e.g. output=grid_core('column.fits', 0.005, 10, dp=0.001)
;
; INPUTS:
;    filename     : the name of the FITS file containing column
;                   density N_H, in units of cm^-2
;
;    pix_size     : the resolution of each pixel in the map, in units of pc
;
;    T            : temperature of the cloud (assumed isothermal) in units of K
;
;    dp           : the interval of potential used for core-finding, 
;                   in units of c_s^2; default is 0.01
;
;    cls_dist     : the closest distance permitted between two local potential,
;                   minima for separate cores, in units of pix_size; default is 2
;                   (closer regions are merged)
;
;    h            : the length at which to adjust the 2D potential, in
;                   units of pix_size; default is 1
;
;    r_pix_lim    : the radius of the smallest core that is considered
;                   resolved, in units of pix_size; default is 2
;
; OUTPUTS:
;    A .fits file containing the identified core areas 
;    
;    A .fits file containing the identified bound core areas
;
;    A .ps file showing the core and bound core regions as well as the
;               calculated gravitational potential contours,
;               overlaid on the surface density map
;
;    A .dat file containing the properties of each core: 
;             core number
;             location i (pixel #)
;             location j (pixel #)
;             total mass (Msun)
;             total mass with background subtraction (Msun)
;             bound mass (Msun)
;             area of whole core region (pixels)
;             area of bound region (pixels)
;             potential well |Phi| at center [(km/s)^2]
;             potential well Phi_max-Phi_min [(km/s)^2]
;
; AUTHORS:
;    Hao Gong & Eve C. Ostriker  2011
;    Department of Astronomy, University of Maryland
;    hgong@astro.umd.edu, ostriker@astro.umd.edu
;
;    Modified by Che-Yu Chen, 2016
;    Department of Astronomy, University of Virginia
;    cheyu.c@gmail.com
;
; NOTES:  Some algorithms are adapted from the CLUMPFIND code of 
;         J. Williams, E. de Geus, and L. Blitz (1994, ApJ 428, 693); see 
;         http://www.ifa.hawaii.edu/~jpw/page16/page4/page4.html


COMMON SHARE100,surfd,data,assign,assign_b,ncore
COMMON SHARE200,G,mp,dx,cs2,msun

on_error,2

if N_params() LT 3 THEN BEGIN
  print,'Wrong number of input arguments ...'
  print,'Syntax - output=grid_core(filename,pix_size,T,cls_dist=cls_dist,dp=dp,h=h,r_pix_lim=r_pix_lim)  '
  print,'Example: output=grid_core("column.fits",0.01,10.,cls_dist=6,dp=0.1,h=1.0,r_pix_lim=3)          '
  print,' filename          : FITS file containing column density N_H, in units of cm^-2 '
  print,' pix_size          : the resolution of each pixel, in units of pc'  
  print,' T                 : temperature of the cloud, in units of K '
  print,' cls_dist          : the minimum distance between core centers, in units of pix_size  '
  print,' dp                : the potential interval for core-finding, in units of c_s^2 '
  print,' h                 : the effective thickness of the layer (smoothing scale for the 2D potential), in units of pix_size '
  print,' r_pix_lim         : the limit of the radius of resolved cores, in units of pix_size '
  return,0
endif

;; start to count time
t0=systime(1)

;; Physical constants and parameters
G    = double(6.6743  *10.^(-8))                ; Gravitational constant
k    = double(1.3807  *10.^(-16))               ; Boltzmann's constant
mp   = double(1.6726  *10.^(-24))               ; mass of proton
msun = double(1.989   * 10.^33)                 ; solar mass
pc   = double(3.086   * 10.^18 )                ; pc/cm conversion

;*************************** setting up parameters *****************************
dx        = double(pix_size*pc)                 
 
cs2       = double(k*T/(2.3*mp))                ; sound speed square, assuming molecular gas
cs2       = cs2/10.^10                          ; convert cs to km/s

IF KEYWORD_SET(cls_dist)  THEN cls_pix = cls_dist           ELSE cls_pix   = double(2.)
IF KEYWORD_SET(dp)        THEN dphi    = double(dp*cs2)     ELSE dphi      = double(0.01*cs2)
IF KEYWORD_SET(h)         THEN hdx     = double(h*dx)       ELSE hdx       = double(dx)
IF KEYWORD_SET(r_pix_lim) THEN r_pix   = r_pix_lim          ELSE r_pix     = 2

;*************************** print out parameters *****************************
print,'The set up is ...'
print, pix_size,  format='("Resolution =                                 ", f6.4, 2x, "pc")'
print, T,         format='("temperature of the cloud =                   ", f6.2, 2x, "K")'
print, cls_pix,   format='("minimum distance between core centers =      ", f6.2, 2x, "pixel")'
print, dphi,      format='("the potential interval for core-finding =  ",   e8.2, 2x, "km^2/s^2")'
print, hdx/dx,    format='("the effective thickness of the layer =       ", f6.2, 2x, "pixel")'
print, r_pix,     format='("the limit of the radius of resolved cores =  ", f6.2, 2x, "pixel")'

;***************************readin data***************************************
;; readin data
surfd=readfits(filename,sdh)

;; test
;i_nan = where(~finite(surfd), /null)
;surfd(i_nan)=0.0
;; end of test

surfd=surfd*1.4*mp                         ; convert number density to surface density

fsize=size(surfd)
nx = fsize(1)
ny = fsize(2)

;; set up coordinates for plotting
x=dindgen(nx)*dx                           ; coordinates for later plot
y=dindgen(ny)*dx

;******************************************************************************

;***********Calculating gravitational potential of the column density**********

k0x=2.*!pi/(2.*nx*dx)
k0y=2.*!pi/(2.*ny*dx)

;; Faking periodic surface density: original map at center + mean N outside
meanN = mean(surfd)
pad0_surfd=dblarr(2*nx,2*ny) + meanN
mnx=fix(0.5*nx)
mny=fix(0.5*ny)
pad0_surfd(mnx:nx+mnx-1, mny:ny+mny-1)=surfd(0:nx-1,0:ny-1)


;; wave number for gravitational potential calculation
kw=dblarr(2*nx,2*ny)

kx = dindgen(2*nx)-nx+1.
kx = shift(kx,-nx+1)

kyt = dindgen(2*ny)-ny+1.
kyt = shift(kyt,-ny+1)
ky = transpose(kyt)

kxx = k0x*rebin(kx,2*nx,2*ny)
kyy = k0y*rebin(ky,2*nx,2*ny)

kw = sqrt(kxx^2+kyy^2)

kw(0,0) = (min(kw)+max(kw))/2.             ; set kw(0,0) as non-zero value

;; FFT of zero-padded surface density
fsurfd=fft(pad0_surfd)

;; calculate gravitational potential
fphi_surfd=-2.*!pi*G*fsurfd/kw/(1.+kw*hdx)

fphi_surfd(0,0) = 0.0

phi_surfd=fft(fphi_surfd,/inverse)

;; gravitational potential
data=real_part(phi_surfd)

;; extract the original region
data = data(mnx:nx+mnx-1, mny:ny+mny-1)

;; make the potential positive for later core-finding
data=-(data-max(data))

;; convert potential to km/s
data=data/10.^10

;******************************************************************************

;*************merge local potential minima if Delta x < cls_pix***************

indx_phimin=extremes_2d(data,1) ; index array of local potential minima
                                ; (these are maxima of "data", since we have
                                ; changed the sign)

n_extremes=n_elements(indx_phimin)      ; number of local potential minima

;; if only one minimum, stop
if (n_extremes eq 1) then return,0 ; one core, no largest closed contour

;; 2d minima coordinates
jpeak=indx_phimin/nx
ipeak=indx_phimin-jpeak*nx

;; generate index for distance between minima
if (n_extremes gt 2) then begin
   ;measure distances between minima
   coord=[transpose(ipeak),transpose(jpeak)]
   dists=distance_measure(coord)

   i1=intarr(n_extremes-1)
   i1(*)=1
   for i=n_extremes-2,1,-1 do begin
       a=intarr(i)
       a(*)=n_extremes-i
       i1=[i1,a]
   endfor

;; get rid of the less-extreme minimum when two are close
   id_dis=where(dists le cls_pix,ct)
   if (ct gt 0) then begin
      indx_phimin(i1(id_dis))=-1
      n_extremes=n_extremes-ct

      indx_phimin=indx_phimin(where(indx_phimin ge 0))
      jpeak=indx_phimin/nx
      ipeak=indx_phimin-jpeak*nx
      print,"*********************************"
      print,ct,format='(5x,i3.3,1x,"minima excluded")'
      print,"*********************************"
   endif
endif

;****************************************************************************

;; if only one minimum, stop
if (n_extremes eq 1) then return,0

;**************************GRID core-finding*********************************
up_limit = 1.5 * max(data)         ; set up the upper threshold for search_2d

assign=intarr(nx,ny)

;; first mark all the extrema with unique labels
for i=0,n_extremes-1 do assign(indx_phimin(i))=i+1

for i=0,n_extremes-1 do begin
    flag_stop = 0
    
    phi_start=data(indx_phimin(i))
    print,i+1,n_extremes,format='("Working on No.",5x,i6.6,1x,"extreme out of",5x,i6.6)'
    while (flag_stop eq 0) do begin
       pix=search2d(data,ipeak(i),jpeak(i),phi_start,up_limit,/diagonal)
       new_pix=assign(pix)
       nc=new_pix(uniq(new_pix,sort(new_pix))) 
       if (n_elements(nc) le 2) then begin
          assign(pix) = i + 1
          phi_start  = phi_start - dphi
       endif else begin
          flag_stop = 1
       endelse
    endwhile
endfor
;***************************************************************************

;***********************Reject unresolved cores******************************
;; reject clumps with fewer than npixmin pixels
ncore=n_extremes
npixmin=fix(!pi*r_pix^2)
destroy_bad,npixmin,nbad
print,nbad,format='(5x,i6.6,1x,"bad cores excluded")'
;***************************************************************************

;*****resort the core ID sequence using the depth of potential wells********
pot_well = fltarr(ncore)
for i=1,ncore do begin
    clp_id = where(assign eq i)
    pot_well(i-1) = abs(max(data(clp_id)) - min(data(clp_id)))
endfor

new_seq = reverse(sort(pot_well))+1
new_assign= intarr(nx,ny)
for i=1,ncore do begin
    clp_id = where(assign eq new_seq(i-1))
    new_assign(clp_id)= i
endfor

assign = new_assign
new_assign = 0
;***************************************************************************

;***********************background subtraction******************************
;; if the original column density is background subtracted, no need to do this step
;; subtract the surface density background (the mean of bottom 10% of the non-core region)

atemp=where(assign eq 0)
btemp=where(assign ne 0)

tt=surfd(atemp)
tt=tt(sort(tt))
n_tt=n_elements(tt)
bg=mean(tt(0:n_tt/10))

;; surface density with background subtracted
surfd_rb = surfd-bg
;***************************************************************************

;***********************core parameters*************************************
core_para=dblarr(ncore,9)
assign_b=intarr(nx,ny)

boundcore2d,npixmin,surfd_rb,core_para,x,y
;***************************************************************************

;; time
delta_t=(systime(1)-t0)/60.0
print,format='(f9.1," minutes elapsed")',delta_t

;***********************output identified core areas***********************
lccname='lcc.fits'
lccname_b='lcc_b.fits'
outfn="core_parameters.dat"
coreplt="cores_phi_on_surfd.ps"

writefits,lccname,assign
writefits,lccname_b,assign_b
;***************************************************************************

;*********************output identified core parameters*********************
openw,100,outfn,/append
printf,100,"   #    i    j      mass   mass_bs  mass_bound n_tot n_bound  phi_max    D_phi"
printf,100,"                   M_sun    M_sun      M_sun    pix     pix  [km/s]^2   [km/s]^2"
for i = 0, ncore - 1 do begin
    coret=reform(core_para(i,*))
    printf,100,i+1,coret,format='(3(1x,i4.4),3(1x,g9.5),2(1x,i6.6),2(1x,g9.5))'
endfor
close,100
;***************************************************************************

;*********************plot identified core regions**************************

lev = [0.9,1.1]
indx=where(assign gt 0)
assign(indx) = 1
indx=where(assign_b gt 0,ct)
if (ct gt 0) then assign_b(indx) = 1

;; export the core+phi+surfd plot to a postscript file
x=x/pc
y=y/pc
set_plot,'ps'
device,filename=coreplt,/color,bits_per_pixel=8,/landscape
loadct,13
imdisp,alog10(surfd),/axis,xr=[min(x),max(x)],yr=[min(y),max(y)],xs=1,ys=1,range=[min(alog10(surfd)),max(alog10(surfd))],/noresize,xtitle="x [pc]",ytitle="y [pc]"
philev = findgen(20)*(max(data)-min(data))*0.05 + min(data)
contour,data,x,y,xstyle=1,ystyle=1,/overplot,levels=philev,c_colors=0
contour,assign,x,y,xstyle=1,ystyle=1,/overplot,levels=lev,c_colors=210,thick=8
contour,assign_b,x,y,xstyle=1,ystyle=1,/overplot,levels=lev,c_colors=255,thick=8
device,/close
set_plot,'x'

;***************************************************************************

print,"********************************************"
print,ncore,format='(5x,i6.6,1x,"cores found")'
print,"********************************************"
print,"corefind_surfd exited sucessfully"
;return core mass to main program
return,core_para

END
