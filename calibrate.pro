pro calibrate, field, chip, band


t_exp = 3.32
dmax_sirius = 3.  ; max HAWK-I pixel distance between SIRIUS and HAWK-I for photometric calibration star selection
                  ; 1 pixel = 0.053" 
dmax_vvv = 2.  
srej = 2.0 ; sigma rejection for RESISTANT_MEAN
ngrid = 60 ; Defines size of star grid for alignment 
mag_lim_VVV = 16 ; only use stars brighter than this from VVV (Ks-magnitude)

;dmax_near_stars = 10 ; eliminate nearby stars within this amount of HAWK-I pixels (10 pix = 0.53")
; Any secondary star within within delta_r of a photometric reference
; star needs to be at least delta_mag magnitudes fainter
;delta_r = 10. ; ; with 0.053"  pixel scale 20 = 1"
;delta_mag = 5
delta_r = 10. ; ; with 0.106"  pixel scale 10 = 1" 
delta_mag = 2.5
SIGMA_CUT = 5.0 ; sigma cut to be applied in resistant_mean

magerr_si = 0.03 ; max acceptable sirius magnitude uncertainty for photometric calibration stars
magerr_hawki = 0.03 ; max acceptable HAWK-I magnitude uncertainty for photometric calibration stars

; max acceptable astrometric uncertainty for astrometric reference stars (pixels)
sastro_max = 0.05

astro_dist = 'none' ; astrometric distortion for SOLVE_ASTRO (can be left empty)
mag_max = 14 ; Bright magnitude cut off for calibration stars in given band
mag_min = 18 ; Faint magnitude cut off for calibration stars in given band
rebfac = 1
scale_holo = 0.106/rebfac ; pixel scale of holography image (" /pixel)

; Radius for creating median ZP
rmed = 1800. ; with 0.053"  pixel scale 1200 ~1 arcmin
;rmed = 900. ; with 0.053"  pixel scale 1200 ~1 arcmin
minexp = 0.2 ; valid stars should be covered by at leasy 20% of all exposures
minexp_cal = 0.2 ; valid calibration stars should be covered by at leasy 40% of all exposures

; Define input and output paths
; and file names
; -----------------------------

data_path = '/home/data/GNS/2021/' + band + '/' + strn(field) + 'HB/data/'
indir = '/home/data/GNS/2021/' + band + '/' + strn(field) + 'HB/cubeims/chip' + strn(chip) + '/'
indir_list = '/home/data/GNS/2021/' + band + '/' + strn(field) + 'HB/photo/chip' + strn(chip) + '/lists/'
star_list = 'astrophot.txt' ; from merge_sublists.pro

; output
res_dir = '/home/data/GNS/2021/' + band + '/' + strn(field) + 'HB/photo/chip' + strn(chip) + '/lists/'
res_dir_im = '/home/data/GNS/2021/' + band + '/' + strn(field) + 'HB/cubeims/chip' + strn(chip) + '/'

; These data will be used to establish astrometry (via VVV - J is used in GNS1 - therefore the variable name)
indir_list_J = '/home/data/GNS/2021/H' + '/' + strn(field) + 'HB/photo/chip' + strn(chip) + '/lists/' 
star_list_J = 'astrophot.txt'
res_dir_imJ = '/home/data/GNS/2021/H' + strn(field) + 'HB/photo/chip' + strn(chip) + '/'

tmpdir =  '/home/data/GNS/2021/' + band + '/' + strn(field) + 'HB/tmp/' 

; Calibration files and images are taken from VVV J-band
vvv_dir = '/home/data/VVV/Fields/J/' ; new VVV
;vvv_dir = '/home/data/working/GNS_2/VVV/Fields/J/' ; old VVV
ref = readfits(vvv_dir + 'Field' + strn(field)  + '.fits.gz',vvvheader)
ref_file = vvv_dir + 'Field' + strn(field)  + '_stars.txt'
;ref_file = vvv_dir + 'Field' + strn(field)  + '_stars_fine.txt'

; Holography image is similarly aligned as lnx_aligned (see alignquadrants), although there may be an offset
; To transfrom HAWK-I positions from holography image
; to VVV we merely need to scale and shift (x_off, y_off from alignquadrants)
scale = (scale_holo/0.34) ; from HAWK-I to VVV


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (1) Calibrate astrometry
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Read images
im = readfits(indir +  strn(field) + '_chip' + strn(chip) + '_holo.fits.gz', im_head)
noise = readfits(indir +  strn(field) + '_chip' + strn(chip) + '_noise.fits.gz', noise_head)
wt = readfits(indir +  strn(field) + '_chip' + strn(chip) + '_nsubims.fits.gz', wt_head)
exp = readfits(indir +  strn(field) + '_chip' + strn(chip) + '_exp.fits.gz', exp_head)

sz = size(im)
xsize_hawki = sz[1]
ysize_hawki = sz[2]
sz_vvv = size(ref)
xsize_vvv = sz_vvv[1]
ysize_vvv = sz_vvv[2]
sub_size = 1200 ; size of sub-fields for check of ZP variations across the field

; Read star lists and sort by decreasing brightness 
; (important when applying compare_lists)
readcol, ref_file, x_vvv, y_vvv, J_vvv, Ks_vvv, a_vvv, d_vvv, a_pm_vvv, d_pm_vvv, FORMAT='F,F,F,F,F,F,F,F' 
vvv_good = where(Ks_vvv lt mag_lim_VVV)
x_vvv = x_vvv[vvv_good]
y_vvv = y_vvv[vvv_good]
a_vvv = a_vvv[vvv_good]
d_vvv = d_vvv[vvv_good]
m_vvv = Ks_vvv[vvv_good]
ord = sort(m_vvv)
x_vvv = x_vvv[ord]
y_vvv = y_vvv[ord]
a_vvv = a_vvv[ord]
d_vvv = d_vvv[ord]
m_vvv = m_vvv[ord]

n_stars = n_elements(m_vvv)
openw, out, vvv_dir + 'VVV_all.reg', /get_lun
printf, out, 'global color=magenta dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
printf, out, 'image'
for s = 0, n_stars-1 do begin
   printf, out, 'circle('+ strn(x_vvv[s])+','+strn(y_vvv[s])+',5)'
endfor
free_lun, out

; Read HAWK-I list for field to be calibrated
; correct for offset inside LXP as created by findstars_holo.pro
; --------------------------------------------------------------
x_large = 4800  ; xaxis length of large cube
y_large = 4800  ; yaxis length of large cube
xw = x_large/2
yw = y_large/2
if (chip eq 1) then begin
  x_s = 0
  y_s = 0
endif 
if (chip eq 2) then begin
  x_s = xw
  y_s = 0
endif 
if (chip eq 4) then begin
  x_s = xw
  y_s = yw
endif 
if (chip eq 3) then begin
  x_s = 0
  y_s = yw
endif 

readcol, indir_list + star_list, x, y, f, sx, sy, sf, corr, ndet   
; input lists are in frame of large long exposure
; but positions from align_VVV.py are in frame of chip
; this must be corrected to identify the calibration stars
x = x - x_s
y = y - y_s

; sort by decreasing flux and select only stars with small astrometric uncertainty
n_stars = n_elements(f)
ord = REVERSE(sort(f))
x = x[ord]
y = y[ord]
f = f[ord]
sx = sx[ord]
sy = sy[ord]
sf = sf[ord]
good = where((sx lt sastro_max) and (sy lt sastro_max))
x = x[good]
y = y[good]
f = f[good]
sx = sx[good]
sy = sy[good]
sf = sf[good]

n_stars = n_elements(f)
print, 'Selected ' + strn(n_stars) + ' stars from HAWK-I image for astrometric calibration. '

; make a map of the sources
dat = ptr_new({X_size: rebfac*20, Y_size: rebfac*20, Sigma_x: 1.5, Sigma_y:1.5, Angle: 0.0})
map = image_model(x,y,f,xsize_hawki,ysize_hawki,'gaussian', dat)
writefits, indir + 'map.fits', map

; find common stars VVV J and HAWK-I
; We can use the common stars found by ../scripts/align_VVV.py
; to initiate the search
; =================================================================
readcol, data_path + 'aa_stars_hawki_' + strn(chip) + '.txt', x_ref, y_ref       ; positions in HAWK-I image
readcol, data_path + 'aa_stars_vvv_' + strn(chip) + '.txt', x_ref_vvv, y_ref_vvv ; positions of HAWK-I transformed to VVV

n_astro = n_elements(x_ref)
print, 'Using ' + strn(n_astro) +  ' stars to compute preliminary astrometric solution.'
polywarp, x_ref_vvv, y_ref_vvv, x_ref, y_ref, 1, Kx, Ky
print, Kx
print, Ky
xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y
yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y
compare_lists, x_vvv, y_vvv, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax_vvv, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
nc = n_elements(subc1)
print, 'Initial alignment.'
print, 'Found ' + strn(nc) + ' common stars.'

openw, out, indir + 'HAWKI_ini.reg', /get_lun
printf, out, 'global color=cyan dashlist=8 3 width=2 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
printf, out, 'image'
for s = 0, nc-1 do begin
   printf, out, 'circle('+ strn(x[subc2[s]])+','+strn(y[subc2[s]])+',2)'
endfor
free_lun, out
;STOP

; ;Now compute transformation with HAWK-I positions, not the
; transformed HAWK-I positions
n_astro = n_elements(subc1)
print, 'Using ' + strn(n_astro) +  ' stars to compute preliminary astrometric solution.'
polywarp, x_vvv[subc1], y_vvv[subc1], x[subc2], y[subc2], 1, Kx, Ky
print, Kx
print, Ky
xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y
yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y
compare_lists, x_vvv, y_vvv, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax_vvv, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
nc = n_elements(subc1)
print, 'Iteration 0'
print, 'Found ' + strn(nc) + ' common stars.'
;STOP

; iterative degree search for alignment stars
; deselect outliers
; ------------------------------------------

n_stars = n_elements(f)
degree = 2

dx_tmp =  (x1c-x2c)*0.34
dy_tmp = (y1c-y2c)*0.34
dpos_tmp = sqrt(dx_tmp^2 + dy_tmp^2)
RESISTANT_Mean, dpos_tmp, 3.0, dpos_mean, dpos_sigma, Num_RejECTED, GOODVEC=goodvec

for it = 1, 3 do begin
  grid_idx = gridselect(x[subc2[goodvec]], y[subc2[goodvec]],f[subc2[goodvec]],gridsize=ngrid)
  polywarp, x_vvv[subc1[goodvec[grid_idx]]], y_vvv[subc1[goodvec[grid_idx]]], x[subc2[goodvec[grid_idx]]], y[subc2[goodvec[grid_idx]]], degree, Kx, Ky
  print, Kx
  print, Ky
  xi = replicate(0.0,n_stars)
  yi = replicate(0.0,n_stars)
  for k = 0, degree do begin
    for m = 0, degree do begin
      xi = xi + Kx[m,k]*x^k*y^m
      yi = yi + Ky[m,k]*x^k*y^m
    endfor
  endfor

  compare_lists, x_vvv, y_vvv, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax_vvv, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
  nc = n_elements(subc1)
  print, 'Iteration ' + strn(it)
  print, 'Found ' + strn(nc) + ' common stars.'

  dx_tmp =  (x1c-x2c)*0.34
  dy_tmp = (y1c-y2c)*0.34
  dpos = sqrt(dx_tmp^2 + dy_tmp^2)
  RESISTANT_Mean, dpos, 3.0, dpos_mean, dpos_sigma, Num_RejECTED, GOODVEC=goodvec
  print, 'Rejected ' + strn(Num_RejECTED) + ' stars.'

endfor


;compare_lists, x_vvv, y_vvv, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax_vvv, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
;dx_tmp =  (x1c-x2c)*0.34
;dy_tmp = (y1c-y2c)*0.34
;dpos = sqrt(dx_tmp^2 + dy_tmp^2)
;RESISTANT_Mean, dpos, 3.0, dpos_mean, dpos_sigma, Num_RejECTED, GOODVEC=goodvec
;print, 'Rejected ' + strn(Num_RejECTED) + ' stars.'

;STOP

; Plot distributions of differences of
; astrometric  positions between reference stars in HAWK-I and VVV
; -----------------------------------------------------------------

dx_tmp =  (x1c[goodvec]-x2c[goodvec])*0.34
dy_tmp = (y1c[goodvec]-y2c[goodvec])*0.34
set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME =  res_dir + 'position_uncertainty_' + strn(chip) + '.eps', XSIZE=20., YSIZE=12., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[0,2,1]
!P.CHARSIZE=1.2
!P.THICK=4.0
!P.CHARTHICK=4.0

binsize = 0.01
cgHistoplot, dx_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Distance in X-axis (")', YTITLE = 'Number of Stars',xrange = [-0.5,0.5], xticks = 4

binCenters1 = loc + (binsize / 2.0)
yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.17, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.17, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

cgHistoplot, dy_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Distance in Y-axis (")' , YTITLE = '' ,xrange = [-0.5,0.5], xticks = 4
binCenters1 = loc + (binsize / 2.0)

yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.65, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.65, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

device, /close

; Now compute and apply the astrometric solution 
; ----------------------------------------------
; for the correct astrometric solution we must return into the 
; frame of the long exposure, so add the offsets
x = x + x_s
y = y + y_s
; +1 MUST BE ADDED TO THE PIXEL COORDINATES IN SOLVE_ASTRO, otherwise there is an offset.
astro = SOLVE_ASTRO(a_vvv[subc1[goodvec]], d_vvv[subc1[goodvec]], x[subc2[goodvec]]+1, y[subc2[goodvec]]+1)
print, astro
;STOP

putast, im_head, astro;, crpix, crval, ctype
putast, noise_head, astro;, crpix, crval, ctype
putast, wt_head, astro;, crpix, crval, ctype
putast, exp_head, astro;, crpix, crval, ctype
sxaddpar, im_head, 'Band',band,'HAWK-I VLT'
sxaddpar, noise_head, 'Band',band,'HAWK-I VLT'
sxaddpar, wt_head, 'Band',band,'HAWK-I VLT'
sxaddpar, exp_head, 'Band',band,'HAWK-I VLT'

; Write astrometrically calibrated images to disc
writefits, res_dir_im + strn(field) + '_chip' + strn(chip) + '_holo_cal.fits', im, im_head, /COMPRESS
writefits, res_dir_im + strn(field) + '_chip' + strn(chip) + '_noise_cal.fits', noise, noise_head, /COMPRESS
writefits, res_dir_im + strn(field) + '_chip' + strn(chip) + '_nsubims_cal.fits', wt, wt_head, /COMPRESS
writefits, res_dir_im + strn(field) + '_chip' + strn(chip) + '_exp_cal.fits', exp, exp_head, /COMPRESS


; Plot distributions of differences of CALIBRATED
; astrometric  positions between reference stars in HAWK-I and VVV
; -----------------------------------------------------------------

XY2AD, x, y, astro, a, d 

a_hawki = a[subc2[goodvec]]
d_hawki = d[subc2[goodvec]]
dra_tmp =  cos(d_vvv[subc1[goodvec]] * !PI/180.) * (a_vvv[subc1[goodvec]] - a_hawki) * 3600.
ddec_tmp = (d_vvv[subc1[goodvec]] - d_hawki) * 3600.

set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME =  res_dir + 'dpos_calibrated_' + strn(chip) + '.eps', XSIZE=20., YSIZE=12., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[0,2,1]
!P.CHARSIZE=1.2
!P.THICK=4.0
!P.CHARTHICK=4.0

binsize = 0.01
cgHistoplot, dra_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Difference in R.A. (")', YTITLE = 'Number of Stars',xrange = [-0.5,0.5], xticks = 4

binCenters1 = loc + (binsize / 2.0)
yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.17, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.17, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

cgHistoplot, ddec_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Difference in Dec. (")' , YTITLE = '' ,xrange = [-0.5,0.5], xticks = 4
binCenters1 = loc + (binsize / 2.0)

yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.65, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.65, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

device, /close


; Save list of positions and uncertainties for later analysis
f_hawki = f[subc2[goodvec]]
forprint, TEXTOUT =  indir_list + 'astrometric_uncertainty_chip' + strn(chip) + '.txt', a_hawki, d_hawki, f_hawki, dra_tmp, ddec_tmp, format='(7F18.8)',/NOCOMMENT


; Read full list of stars and  Write astrometrically calibrated list to disc
; ============================================================================
readcol, indir_list + star_list, x, y, f, sx, sy, sf, corr, ndet   
XY2AD, x, y, astro, a, d 
star_list = 'calibrated_' + band + '_' + strn(field) + '_' + strn(chip) + '.txt'
forprint, TEXTOUT =  indir_list + star_list, a, d, x, y, f, sx, sy, sf, format='(8F18.8)',/NOCOMMENT

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; (2) Photometric calibration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Read astrometrically calibrated star list
readcol, indir_list + star_list,  a, d, x, y, f, sx, sy, sf, FORMAT='D,D,D,D,D,D,D,D'

; Read in SIRIUS photometric reference catalogue 
; -------------------------------------------------------

readcol, '/home/data/working/GNS_2/SIRIUS/list_all_Sirius.txt',$
  a_si, dec_si, J_si, dJ_si, H_si, dH_si, K_si, dK_si, Format = 'D,D,D,D,D,D,D,D'   

; check the uncerinties of SIRIUS to get an idea
;test = where(H_si gt 0 and H_si lt 15)
;print, strn(mean(dH_si[test])) + ', ' + strn(stddev(dH_si[test]))
;test = where(K_si gt 0 and K_si lt 13)
;print, strn(mean(dK_si[test])) + ', ' + strn(stddev(dK_si[test]))
;test = where(J_si gt 0 and J_si lt 17)
;print, strn(mean(dJ_si[test])) + ', ' + strn(stddev(dj_si[test]))
;STOP
;  Select corresponding band 
; -------------------------------

if band eq 'J' then begin
  m_si = J_si
  dm_si = dJ_si
endif

if band eq 'H' then begin
 m_si = H_si
 dm_si = dH_si
endif

if band eq 'Ks' then begin
 m_si = K_si
 dm_si = dK_si
endif

; use only valid measurements
valid = where(m_si gt 0)
a_si = a_si[valid]
dec_si = dec_si[valid]
m_si = m_si[valid]
dm_si = dm_si[valid]

; Select calibration stars by brightness
if band eq 'J' then begin
  good = where(m_si gt 10.0 and dm_si lt magerr_si) 
endif

if band eq 'H' then begin
;  good = where(m_si gt 11.0  and m_si lt 14 and dm_si lt 0.04)
  good = where(m_si gt 11.0  and dm_si lt magerr_si)
endif

if band eq 'Ks' then begin
  good = where(m_si gt 11.0 and dm_si lt magerr_si)    
endif

; appendix _si means all valid stars in SIRIUS catalogue
; appendic _ref means potential photmetric reference stars
a_ref = a_si[good]
dec_ref = dec_si[good]
m_ref = m_si[good]
dm_ref = dm_si[good]


; Compute pixel positions of reference stars
AD2XY, a, d, astro, x, y                        ; stars in HAWK-I image
AD2XY, a_ref, dec_ref, astro, x_ref, y_ref      ; selected SIRIUS reference stars
;x_ref = x_ref - x_s
;y_ref = y_ref - y_s

; Select  stars within HAWK-I image
; ----------------------------------------------
cutx_min = min(x)
cutx_max = max(x)
cuty_min = min(y)
cuty_max = max(y)
good = where(x_ref gt cutx_min and x_ref lt cutx_max and y_ref gt cuty_min and y_ref lt cuty_max, n_ref)
x_ref = x_ref[good]
y_ref = y_ref[good]
m_ref = m_ref[good]
dm_ref = dm_ref[good]
print, strn(n_ref) + ' potential photometric reference stars in the field.'

; Now avoid calibration stars in low S/N areas, perticularly near the
; borders
; Use only areas with > 30% of max exposure number
good = []
exp = exp/max(exp)
for i = 0, n_ref -1 do begin
  xx = round(x_ref[i]) 
  yy = round(y_ref[i])
  if exp[xx,yy] gt minexp_cal then good = [good,i]
endfor
n_good = n_elements(good)
x_ref = x_ref[good]
y_ref = y_ref[good]
m_ref = m_ref[good]
dm_ref = dm_ref[good]
print, strn(n_good) + ' high SNR photometric reference stars in the field.'
;STOP
; Create map of all potential reference stars for this field
dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 2, Sigma_y: 2, Angle: 0.0})
im = image_model(x_ref,y_ref,10^((-m_ref)/2.5),xsize_hawki,ysize_hawki,'gaussian', dat)
writefits, tmpdir +'sirius_' + band + '_0.fits', im, /COMPRESS

; Find all potential reference stars common to HAWK-I and  SIRIUS
compare_lists, x, y, x_ref, y_ref, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax_sirius,$
   SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2 
print, 'Found ' + strn(n_elements(subc1)) + ' common stars for SIRUS and HAWK-I.'

; Create map of all potential reference stars for this field
dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 2, Sigma_y: 2, Angle: 0.0})
im = image_model(x_ref[subc2],y_ref[subc2],10^((-m_ref[subc2])/2.5),xsize_hawki,ysize_hawki,'gaussian', dat)
writefits, tmpdir +'sirius_' + band + '_common.fits', im, /COMPRESS
;STOP

x_common_hawki = x[subc1]
y_common_hawki = y[subc1]
f_common_hawki = f[subc1]
df_common_hawki = sf[subc1]
mag_common_hawki = - 2.5*alog10(f_common_hawki) ; instrumental magnitude of common stars
mag = -2.5*alog10(f)                            ; instrumental magnitude of all HAWKI stars
;STOP
; Find all stars common stars that are isolated
ind_iso = []
ISOLATED_STARS, x, y, mag, x_common_hawki, y_common_hawki, mag_common_hawki, delta_mag, delta_r, ind_iso
print, 'Found ' + strn(n_elements(ind_iso)) + ' isolated common stars in HAWK-I.'
x_ref_hawki = x_common_hawki[ind_iso]
y_ref_hawki = y_common_hawki[ind_iso]
f_ref_hawki = f_common_hawki[ind_iso]
df_ref_hawki = df_common_hawki[ind_iso]
;STOP

; Now make sure that there is no invalid pixel within +- mask pixels
; of the potential calibration stars
;sz = size(wt)
;valid = []
;mask = 50
;for i=0, n_elements(x_ref_hawki)-1 do begin
;  xlo = 0 > round(x_ref_hawki[i]-mask)
;  xhi = (xsize_hawki-1) < round(x_ref_hawki[i]+mask)
;  ylo = 0 > round(y_ref_hawki[i]-mask)
;  yhi = (ysize_hawki-1) < round(y_ref_hawki[i]+mask)
;  value = total(where(exp[xlo:xhi,ylo:yhi] eq 0))
;  if value eq -1 then begin
;    valid = [valid, i]  
;  endif
;endfor
;x_ref_hawki = x_ref_hawki[valid]
;y_ref_hawki = y_ref_hawki[valid]
;f_ref_hawki = f_ref_hawki[valid]
;df_ref_hawki = df_ref_hawki[valid]


; Find common stars with potential SIRIUS reference sources
compare_lists, x_ref_hawki, y_ref_hawki, x_ref, y_ref, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax_sirius,$
   SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2 
print, 'Found ' + strn(n_elements(subc1)) + ' stars for photometric calibration.'
x_ref_hawki = x_ref_hawki[subc1]
y_ref_hawki = y_ref_hawki[subc1]
f_ref_hawki = f_ref_hawki[subc1]
df_ref_hawki = df_ref_hawki[subc1]
mag_ref_sirius = m_ref[subc2]

; Select sources with photometric uncertainties magerr_hawki,
; fainter than mag_max, and brighter than mag_min
; need to get first ZP estimate for this purpose
zp0 = median(mag_ref_sirius + 2.5*alog10(f_ref_hawki/t_exp))
m0 = zp0 - 2.5*alog10(f_ref_hawki/t_exp)
dm0 = 2.5/alog(10) * df_ref_hawki/f_ref_hawki
good = where(dm0 lt magerr_hawki and m0 gt mag_max and m0 lt mag_min, count)
x_ref_hawki = x_ref_hawki[good]
y_ref_hawki = y_ref_hawki[good]
f_ref_hawki = f_ref_hawki[good]
df_ref_hawki = df_ref_hawki[good]
print, 'Found ' + strn(count) + ' with right magnitude and sufficiently low uncertainty.'
;STOP

; Find common stars with potential SIRIUS reference sources
compare_lists, x_ref_hawki, y_ref_hawki, x_ref, y_ref, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax_sirius,$
   SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2 
print, 'Found ' + strn(n_elements(subc1)) + ' stars for photometric caibration.'
x_ref_hawki = x_ref_hawki[subc1]
y_ref_hawki = y_ref_hawki[subc1]
f_ref_hawki = f_ref_hawki[subc1]
df_ref_hawki = df_ref_hawki[subc1]
x_ref_sirius = x_ref[subc2]
y_ref_sirius = y_ref[subc2]
mag_ref_sirius = m_ref[subc2]
dmag_ref_sirius = dm_ref[subc2]

; Make map of calibration stars
dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 2, Sigma_y: 2, Angle: 0.0})
im = image_model(x_ref_sirius,y_ref_sirius,10^(-0.4 * mag_ref_sirius),xsize_hawki,ysize_hawki,'gaussian', dat)
;writefits, res_dir + 'phot_ref_' + band + '.fits', im, /COMPRESS
writefits, tmpdir + 'phot_ref_' + band + '_chip' + strn(chip) +  '.fits', im, /COMPRESS

; DETERMINE ZP
; ----------------------------------------------------------------
zp = mag_ref_sirius + 2.5*alog10(f_ref_hawki/t_exp)
resistant_mean, zp, SIGMA_CUT, zp_mean, zp_sigma, rej, goodvec=good_cal
print, '*****'
print, 'ZP from resistant mean: ' + strn(zp_mean) + ' +- ' +strn(zp_sigma)
;zp_fin = zp_mean
 
n_ref = n_elements(good_cal)
x_ref_hawki = x_ref_hawki[good_cal]
y_ref_hawki = y_ref_hawki[good_cal]
mag_ref_sirius = mag_ref_sirius[good_cal]
f_ref_hawki = f_ref_hawki[good_cal]
 

; Map of zero points
; ---------------------------
zpoints = mag_ref_sirius + 2.5*alog10(f_ref_hawki/t_exp)
redfac = 20 ; reduction factor for faster computation
xs_map = xsize_hawki/redfac
ys_map = ysize_hawki/redfac
zpmap = fltarr(xs_map,ys_map)
dzpmap = fltarr(xs_map, ys_map)
for ix = 0L, xs_map-1 do begin
for iy = 0L, ys_map-1 do begin
  r = sqrt((redfac*ix-x_ref_hawki)^2 + (redfac*iy-y_ref_hawki)^2)
  nearby = where(r le rmed, n_near)
  if (n_near gt 3) then begin
    resistant_mean, zpoints[nearby], SIGMA_CUT, zp_mean, zp_sigma, rej, goodvec=good_cal
    zpmap[ix,iy] = zp_mean
    dzpmap[ix,iy] = zp_sigma
  endif
endfor
;print, ix
endfor
zpmap = filter_image(zpmap, median=10,/ALL)
zpmap = CREBIN(zpmap,xsize_hawki,ysize_hawki)
dzpmap = CREBIN(dzpmap,xsize_hawki,ysize_hawki)
writefits, tmpdir + 'zpmap_' + band + '_chip' + strn(chip) +  '.fits', zpmap, /COMPRESS
writefits, tmpdir + 'dzpmap_' + band + '_chip' + strn(chip) +  '.fits', dzpmap, /COMPRESS


; Histogram of zero points
; ---------------------------
zpoints = mag_ref_sirius + 2.5*alog10(f_ref_hawki/t_exp)
set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
    FILENAME= res_dir + 'ZP_hist_' + band + '.eps', XSIZE=20, YSIZE=15., $
    /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
  !P.MULTI=0
  !P.CHARSIZE=1.3
  !P.THICK=4.0
  !X.THICK=4.0
  !Y.THICK=4.0
  !P.CHARTHICK=4.0
  darkblue = cgColor("Blue")
  blue = cgColor("Dodger Blue")
  red = cgColor("Dark Red")
  green = cgColor("Green")
  black =  cgColor("Black")

  xaxrange = [25.8,26.8]
  if band eq 'Ks' then  xaxrange = [25.0,26.0] else  xaxrange = [25.8,26.8]
  binsize = 0.03
  cgHistoplot, zpoints, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'ZP', xrange = xaxrange, xticks = 4, YTITLE='N'
  binCenters1 = loc + (binsize / 2.0)

  yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
  cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  print, '*****'
  print, 'ZP from Gauss fit to histogram: ' + strn(coeff[1]) + ' +- ' +strn(coeff[2]) + ' ' + strn(coeff[2]/sqrt(n_elements(zpoints)))
  print, '*****'
 
  zp_fin = coeff[1]

device, /close

print, 'ZERO POINT: ' + strn(zp_fin) + ' +- ' + strn(zp_sigma)
;print, 'used ' + strn(n_elements(good_cal)) + ' stars for photometric calibration.'


;Comparison of calibrated HAWK-I reference star magnitudes with SIRIUS magnitudes to estimate the uncertainty
; -------------------------------------------------------------------------------------------------------------

; Scatter plot
; -----------

set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
    FILENAME= res_dir + 'ZP_calib_' + band + '.eps', XSIZE=30, YSIZE=15., $
    /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
;  !P.MULTI=[0,2,1]
  !P.MULTI=[0,1,1]
  !P.CHARSIZE=1.3
  !P.THICK=4.0
  !X.THICK=4.0
  !Y.THICK=4.0
  !P.CHARTHICK=4.0
  darkblue = cgColor("Blue")
  blue = cgColor("Dodger Blue")
  red = cgColor("Dark Red")
  green = cgColor("Green")
  black =  cgColor("Black")

  ; Calibration with single ZP
  m_hawki = zp_fin - 2.5*alog10(f_ref_hawki/t_exp)
  dm = m_hawki - mag_ref_sirius
  plot, mag_ref_sirius, dm,  XRAN = [10,16], XTITLE='[' + band + ']', YTITLE='[d' + band + ']', XSTYLE=1, PSYM=1, YSTYLE=1, YRAN=[-0.25,0.25]
  oplot, [10.0, 16.0], [0.0,0.0], PSYM=-3, color = darkblue

  ; plot uncertainty in bins
  min_m = round(min(mag_ref_sirius))
  max_m = round(max(mag_ref_sirius))
  nbin = 4
  mags = 11.5 + findgen(nbin)
  mmag = fltarr(nbin)
  sigmag = fltarr(nbin)
  for j = 0, nbin - 1 do begin
    thisbin = where(abs(mag_ref_sirius - mags[j]) le 0.5)
    vals = dm[thisbin]
    if n_elements(vals) gt 1 then begin
    mmag[j] = mean(vals)
    sigmag[j] = stddev(vals);/sqrt(n_elements(thisbin))
    xyouts, mags[j]-0.5, -0.3, strn(sigmag[j],FORMAT='(f5.3)'), color = blue, charsize = 0.8
    xyouts, mags[j]-0.5, -0.35, strn(n_elements(vals),FORMAT='(i)'), color = black, charsize = 1
    endif
  endfor
  oploterror, mags, mmag, sigmag, errcolor = blue, errthick = 4, psym = 3
  xyouts, 11., 0.45, 'Number of detected stars: ' + strn(n_elements(dm),format='(I)')
  xyouts, 11, 0.35, 'ZP =  ' + strn(zp_fin,FORMAT='(f6.3)') + ' +- ' + strn(zp_sigma,FORMAT='(f6.3)')

; Calibration with variable ZP
m_hawki = fltarr(n_ref)
for s = 0, n_ref-1 do begin
  xx = round(x_ref_hawki[s])
  yy = round(y_ref_hawki[s])
  zp_loc = zpmap[xx,yy]
  m_hawki[s] = zp_loc - 2.5*alog10(f_ref_hawki[s]/t_exp)
endfor
dm = m_hawki - mag_ref_sirius

;  plot, mag_ref_sirius, dm,  XRAN = [10,16], XTITLE='[' + band + ']', YTITLE='[d' + band + ']', XSTYLE=1, PSYM=1, YSTYLE=1, YRAN=[-0.25,0.25]
  oplot, mag_ref_sirius, dm,  PSYM = 1, color = red

device, /close


; FINAL PHOTOMETRIC CALIBRATION
; ==============================
; Read astrometrically calibrated star list
readcol, indir_list + star_list,  a, d, x, y, f, sx, sy, sf, FORMAT='D,D,D,D,D,D,D,D'
n_stars = n_elements(f)
sm = 2.5/alog(10) * (sf/f)

; Calibration with variable ZP
; ---------------------------
m = fltarr(n_stars)
for s = 0, n_stars-1 do begin
   xx = round(x[s])
   yy = round(y[s])
   zp_loc = zpmap[xx,yy]
   m[s] = zp_loc - 2.5*alog10(f[s]/t_exp)
   dzp_loc = dzpmap[xx,yy]
   sm[s] = sqrt((sm[s])^2 +dzp_loc^2) 
 endfor

; Calibration with single ZP
; COMMENT IF YOU WISH TO USE A VARIABLE ZP
; ---------------------------
;m = zp_fin - 2.5 * alog10(f/t_exp)



  ; map with only calibrated stars in HAWK-I
dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 1.7, Sigma_y: 1.7, Angle: 0.0})
im = image_model(x,y,f,xsize_hawki,ysize_hawki,'gaussian', dat)
writefits, res_dir +'calibrated_stars_' + strn(chip) + '.fits', im, im_head, /COMPRESS
   
; Determine number of exposures for each star
exp = readfits(indir +  strn(field) + '_chip' + strn(chip) + '_exp.fits.gz', exp_head)
expos = lonarr(n_stars)
for s = 0, n_stars -1 do begin
  xx = round(x[s])
  yy = round(y[s])
  expos[s] = exp[xx,yy]
endfor

;Write calibrated list of stars to file
forprint, TEXTOUT= res_dir + 'stars_calibrated_' + band + '_chip' + strn(chip) + '.txt', a ,d, x, y, f, m, sx, sy, sf, sm, expos, format='(11F18.8)', /NOCOMMENT

; Map of photometric uncertainty
; ---------------------------
redfac = 20 ; reduction factor for faster computation
rmed = 100 
bright = where(m lt 18)
x = x[bright]
y = y[bright]
sm = sm[bright]
xs_map = xsize_hawki/redfac
ys_map = ysize_hawki/redfac
smmap = fltarr(xs_map,ys_map)
for ix = 0L, xs_map-1 do begin
for iy = 0L, ys_map-1 do begin
  r = sqrt((redfac*ix-x)^2 + (redfac*iy-y)^2)
  nearby = where(r le rmed, n_near)
  if (n_near gt 3) then begin
    resistant_mean, sm[nearby], 3.0, sm_mean, sm_sigma, rej, goodvec=good_cal
    smmap[ix,iy] = sm_mean
  endif
endfor
;print, ix
endfor
smmap = CREBIN(smmap,xsize_hawki,ysize_hawki)
writefits, tmpdir + 'smmap_' + band + '_chip' + strn(chip) +  '.fits', smmap, /COMPRESS

; PLots of relative astrometric and photometric uncertainties
; ============================================================== 

readcol,  res_dir + 'stars_calibrated_' + band + '_chip' + strn(chip) + '.txt', a ,d, x, y, f, m, sx, sy, sf, sm, FORMAT='D,D,D,D,D,D,D,D,D,D'

set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
    FILENAME = res_dir + 'uncertainty_phot_' + strn(chip) + '.eps', XSIZE=20., YSIZE=15., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
    !P.MULTI=0
    !P.CHARSIZE=1
    !P.THICK=4.0
    !P.CHARTHICK=1.0
    cgplot, m, sm,  PSYM=4, Color = 'black', XTITLE='[' + band + strn(chip) +  ']', YTITLE='[d' + band + ']',$
      THICK = 0.5, XRANGE = [12,23], YRANGE = [0,0.2], SYMSIZE=0.2
device, /close 


sx = sx * scale_holo * 1000
sy = sy * scale_holo * 1000
set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
    FILENAME = res_dir + 'uncertainty_xpos_' + strn(chip) + '.eps', XSIZE=20., YSIZE=15., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
    !P.MULTI=0
    !P.CHARSIZE=1
    !P.THICK=4.0
    !P.CHARTHICK=1.0
    cgplot, m, sx,  PSYM=4, Color = 'black', XTITLE='[' + band + ', Chip ' + strn(chip) +  ']', YTITLE='dpos, x-axis [mas]',$
            THICK = 0.5, XRANGE = [11,23], YRANGE = [0,30], SYMSIZE=0.2
;    nsteps = 18
;    dm_bins = fltarr(nsteps)
;    dm_sigma_bins = fltarr(nsteps)
;    m_bins = 12 + 0.5 * findgen(nsteps)
;    for im = 0, nsteps-1 do begin
;       thisbin = m_bins[im]
;       vals = sx[where(m ge (thisbin-0.5) and m lt (thisbin+0.5),nvals)]
;       resistant_mean, vals, 3.0, dm_mean, dm_sigma, rej, goodvec=goodvals
;       dm_bins[im] = dm_mean
;       dm_sigma_bins[im] = dm_sigma * sqrt(n_elements(goodvals))
;    endfor
;    oploterror, m_bins, dm_bins, dm_sigma_bins, color = cgColor('Dodger Blue')
    device, /close

set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
    FILENAME = res_dir + 'uncertainty_ypos_' + strn(chip) + '.eps', XSIZE=20., YSIZE=15., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
    !P.MULTI=0
    !P.CHARSIZE=1
    !P.THICK=4.0
    !P.CHARTHICK=1.0
   cgplot, m, sy,  PSYM=4, Color = 'black', XTITLE='[' + band + strn(chip) +  ']', YTITLE='dpos, y-axis [mas]',$
            THICK = 0.5, XRANGE = [11,23], YRANGE = [0,30], SYMSIZE=0.2
;    nsteps = 18
;    dm_bins = fltarr(nsteps)
;    dm_sigma_bins = fltarr(nsteps)
;    m_bins = 12 + 0.5 * findgen(nsteps)
;    for im = 0, nsteps-1 do begin
;       thisbin = m_bins[im]
;       vals = sy[where(m ge (thisbin-0.5) and m lt (thisbin+0.5),nvals)]
;       resistant_mean, vals, 3.0, dm_mean, dm_sigma, rej, goodvec=goodvals
;       dm_bins[im] = dm_mean
;       dm_sigma_bins[im] = dm_sigma * sqrt(n_elements(goodvals))
;    endfor
;    oploterror, m_bins, dm_bins, dm_sigma_bins, color = cgColor('Dodger Blue')
device, /close 


end
