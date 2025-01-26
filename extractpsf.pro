FUNCTION EXTRACTPSF, im, noise, mask, maskrad, psf_fwhm, tmpdir, UNWEIGHTED=unweighted

DIT = 3.32
ZP = 26.3  ; rough ZP

n_ref_min = 5
n_ref_max = 40

; Parameters to select isolated stars
delta_r = psf_fwhm * 2.
delta_mag = 5.

psf_size = 2 * (maskrad + 5) + 1
n_pix = float(psf_size)^2

sz = size(im)
n1 = sz[1]
n2 = sz[2]

; Parameters for PSF estimation and StarFinder
; ------------------------------------------------
min_correlation = 0.8
correl_mag = 4.0
deblend = 0
deblost = 0
niter = 2
rel_thresh = 1
guide_x = ""
guide_y = ""
weigh_med = 0
norm_max = 1


; 1) First estimate of PSF
; ----------------------------

psf = psf_gaussian(NPIXEL=psf_size,FWHM=psf_fwhm,/NORMALIZE,/DOUBLE)

threshold = 10.*noise
background = estimate_background(im,2*maskrad)
search_objects, im, LOW_SURFACE = background, threshold, $
	                PRE_SMOOTH = 1, MINIF = 2, $
	                n, x, y, f
good = where(f gt 0,n)
x = x[good]
y = y[good]
f = f[good]
m = ZP - 2.5 * alog10(f/DIT)


; -----------------------------
; Select reference stars
; -----------------------------

 ; stars must be sorted by decreasing brightness
 ; so that the brightest ones will be used to extract the PSF
 ord = sort(m)
 x_psf = x[ord]
 y_psf = y[ord]
 m_psf = m[ord]
 n_ref = n_elements(m_psf)
 print, 'Found '+ strn(n_ref) + ' stars.'

; Reference stars must lie in support region
; ----------------------------------------------------------------
 accept = replicate(1,n_ref)
 for s = 0, n_ref-1 do begin
   xx = round(x_psf[s])
   yy = round(y_psf[s])
   xx = 0 > xx & xx = (n1 - 1) < xx
   yy = 0 > yy & yy = (n2 - 1) < yy
   if (mask[xx,yy] lt 1) then accept[s] = 0
 endfor
good = where(accept gt 0, n_ref)
x_psf = x_psf[good]
y_psf = y_psf[good]
m_psf = m_psf[good] 
print, 'Found '+ strn(n_ref) + ' supported reference stars.'


if (n_ref lt n_ref_min) then begin
  print, 'Could not find enough reference stars with sufficient support.'
  STOP
endif

; Select isolated stars
isolated_stars, x, y, m, x_psf, y_psf, m_psf, delta_mag, delta_r, ind_iso
x_psf = x_psf[ind_iso]
y_psf = y_psf[ind_iso]
m_psf = m_psf[ind_iso]
n_ref = n_elements(m_psf)
print, 'Found '+ strn(n_ref) +  ' supported, unsaturated and isolated reference stars.'
;STOP
; Use the brightest stars
; At least n_ref_min stars
if n_ref gt n_ref_max then begin
   n_ref = n_ref_max
endif else begin
  if n_ref gt n_ref_min then n_ref = n_ref
  if n_ref lt n_ref_min then begin
    print, 'Cannot find enough PSF reference stars.'
    STOP
  endif
 endelse
 x_psf = x_psf[0:n_ref-1]
 y_psf = y_psf[0:n_ref-1]
 m_psf = m_psf[0:n_ref-1]  
 print, 'using '+ strn(n_ref) + ' PSF reference stars.'

; print
 print, 'Magnitudes of reference stars: '
 print, m_psf

; -----------------------------

debug = 0
nrad = 2 * psf_fwhm
iter = 0                        ; more iterations do not necessarily make it better
mindist = 2* psf_fwhm
MAKEPSF, x_psf, y_psf, x, y, f, im, nrad, FOVMASK = fov_mask, PSF=psf,  BACKGROUND=background, DEBUG = debug, ITER = iter, MINDIST = mindist, NOISE_PSF = psf_noise, MASKRAD = maskrad, UNWEIGHTED=unweighted, TMPDIR = tmpdir, LOCAL_SKY=1

 psf = centroider(psf)
 MASK_PSF, psf, maskrad, PSF_MASKED=psf_masked
 psf = psf_masked  
 psf = psf/total(psf)  ; normalization of PSF
 writefits, tmpdir + 'tmppsf0.fits', psf
;STOP

 ; Run StarFinder and iterate search for PSF reference stars and PSF extraction
 ; --------------------------------------------------------------------------

  Threshold = [10.]
  estim_bg = 1
  back_box = maskrad
  starfinder, im, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
        ESTIMATE_BG = estim_bg, DEBLEND = deblend, N_ITER = niter, SILENT=1, $
        GUIDE_X = guide_x, GUIDE_Y = guide_y, $
        SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
        x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC;, /NO_SLANT
    writefits, tmpdir + 'stars.fits', stars
    writefits, tmpdir + 'bg.fits', background
    
 ; 2) Repeat PSF extraction
 ; ---------------------

 
; -----------------------------
; Select reference stars
; -----------------------------

 ; stars must be sorted by decreasing brightness
 ; so that the brightest ones will be used to extract the PSF
 m = ZP - 2.5 * alog10(f)
 ord = sort(m)
 x_psf = x[ord]
 y_psf = y[ord]
 m_psf = m[ord]
 n_ref = n_elements(m_psf)

; Reference stars must lie in support region
; ----------------------------------------------------------------
 accept = replicate(1,n_ref)
 for s = 0, n_ref-1 do begin
   xx = round(x_psf[s])
   yy = round(y_psf[s])
   xx = 0 > xx & xx = (n1 - 1) < xx
   yy = 0 > yy & yy = (n2 - 1) < yy
   if (mask[xx,yy] lt 1) then accept[s] = 0
 endfor
good = where(accept gt 0, n_ref)
x_psf = x_psf[good]
y_psf = y_psf[good]
m_psf = m_psf[good] 
print, 'Found '+ strn(n_ref) + '  reference stars in support region.'
;STOP

if (n_ref lt n_ref_min) then begin
  print, 'Could not find enough reference stars with sufficient support.'
  STOP
endif
 
isolated_stars, x, y, m, x_psf, y_psf, m_psf, delta_mag, delta_r, ind_iso
x_psf = x_psf[ind_iso]
y_psf = y_psf[ind_iso]
m_psf = m_psf[ind_iso]
n_ref = n_elements(m_psf)
print, 'Found  '+ strn(n_ref) + ' supported, unsaturated and isolated reference stars.'

; Use the brightest stars
; At least n_ref_min stars
if n_ref gt n_ref_max then begin
  n_ref = n_ref_max
endif else begin
  if n_ref gt n_ref_min then n_ref = n_ref
  if n_ref lt n_ref_min then begin
    print, 'Cannot find enough PSF reference stars.'
    STOP
   endif
endelse
       
 x_psf = x_psf[0:n_ref-1]
 y_psf = y_psf[0:n_ref-1]
 m_psf = m_psf[0:n_ref-1]  
 print, 'Found '+ strn(n_ref) + ' PSF reference stars.'
 dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 2., Sigma_y: 2., Angle: 0.0})
 f_psf = 1.e6 * 10^(-0.4*m_psf)
 refstars = image_model(x_psf,y_psf,f_psf,n1,n2,'gaussian', dat)
 writefits, tmpdir + 'refstars.fits', refstars
; forprint, TEXTOUT= dir + outnam + 'psfstars.txt', x_psf, y_psf, m_psf, /NOCOMMENT

 print
 print, 'Magnitudes of reference stars: '
 print, m_psf

; -----------------------------

  debug = 0
  nrad =  2* psf_fwhm
  mindist = 2 * psf_fwhm
  iter = 1  ; more iterations do not necessarily make it better
  MAKEPSF, x_psf, y_psf, x, y, f, im, nrad, FOVMASK = fov_mask, PSF=psf,  BACKGROUND=background, DEBUG = debug, ITER = iter, MINDIST = mindist, NOISE_PSF = psf_noise, MASKRAD = maskrad, UNWEIGHTED=unweighted, TMPDIR = tmpdir, LOCAL_SKY=1


  psf = centroider(psf)
;  writefits, tmpdir + 'tmppsf_nomask.fits', psf
  MASK_PSF, psf, maskrad, PSF_MASKED=psf_masked
  psf = psf_masked  
;  writefits, tmpdir + 'tmppsf_masked.fits', psf  
  psf = psf/total(psf)  ; normalization of PSF
  RETURN, psf


END
