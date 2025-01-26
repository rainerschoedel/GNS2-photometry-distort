pro merge_subims, field_nr, chip_nr, Band


; PURPOSE: Piece the holographically reduced sub-images together
;          use input from align_lists.pro/merge_sublists.pro
;


field = strn(field_nr)
chip = 'chip' + strn(chip_nr)

; input and output directories
indir = '/home/data/GNS/2021/'+ Band + '/' + field + 'HB'
outdir = indir
photdir = indir + '/photo/chip' + strn(chip_nr) + '/lists/'
imdir = indir + '/cubeims/chip' + strn(chip_nr) + '/'

border = 20. ; border of images that is to be masked


; Basic parameters
; size of final image
xaxis = 4800L
yaxis = 4800L  
sub_size_x0 = 600
sub_size_y0 = 600
edge = 20

rebfac = 1   ; rebin factor
border = border * rebfac
xaxis = xaxis * rebfac
yaxis = yaxis * rebfac
sub_size_x0 = sub_size_x0 * rebfac
sub_size_y0 = sub_size_y0 * rebfac
; Here, the assumption is that each half of a sub-image
; overlaps with one half of its neighbour (see subcubes.pro)
x_sub_shift = sub_size_x0/2
y_sub_shift = sub_size_y0/2
nx = xaxis/x_sub_shift - 1 
ny = yaxis/y_sub_shift - 1
       
; Read image offsets
; obtained by running 
; merge_lists.pro on H-band
RESTORE, 'x_offsets.SAV'
RESTORE, 'y_offsets.SAV'

; initialise mosaics
im =  fltarr(xaxis,yaxis)
nsubims =  fltarr(xaxis,yaxis)
exp =  fltarr(xaxis,yaxis)
sigma =  fltarr(xaxis,yaxis)

mask_borders = fltarr(sub_size_x0, sub_size_y0)
mask_borders[*,*] = 0
mask_borders[border:sub_size_x0-border-1, border:sub_size_y0-border-1] = 1     


for i_x = 0, nx-1 do begin
  for i_y = 0, ny-1 do begin 

;     name = 'mosaic_' + strn(i_x) + '_' +  strn(i_y)
     name = 'holo_' + strn(i_x) + '_' +  strn(i_y)

   if FILE_TEST(imdir + name + '.fits.gz') then begin
     ; for J and Ks use lists pre-aligned with H
     ; that are produced by align_fields.pro
     if (Band ne 'H') then name = 'offH_' + name    

     tmp_subfield = fltarr(xaxis,yaxis)
     tmp_subfield_exp = fltarr(xaxis,yaxis)
     tmp_subfield_nsubims = fltarr(xaxis,yaxis)
     tmp_subfield_noise = fltarr(xaxis,yaxis)
     
     subim = readfits(imdir + name + '.fits.gz') * mask_borders
     subnoise = readfits(imdir + name + '_sigma.fits.gz') * mask_borders
     subexp = readfits(imdir + name + '_expmap.fits.gz')  * mask_borders


    ; Apply offsets for this sub-field
    ; ----------------------------------------
     xoff_0 =  i_x * x_sub_shift
     yoff_0 =  i_y * y_sub_shift 
     x_off = xoff_0 + x_offsets[i_x,i_y]
     y_off = yoff_0 + y_offsets[i_x,i_y]

     tmp_subfield[0:sub_size_x0-1,0:sub_size_y0-1] = subim            
     tmp_subfield_noise[0:sub_size_x0-1,0:sub_size_y0-1] = subnoise
     tmp_subfield_exp[0:sub_size_x0-1,0:sub_size_y0-1] = subexp

     accept = where(tmp_subfield_exp gt 0, complement=reject)
     tmp_subfield_nsubims[accept] = 1
     tmp_subfield_nsubims[reject] = 0
      
     tmp_subfield = image_shift(tmp_subfield, x_off, y_off)  
     tmp_subfield_noise = image_shift(tmp_subfield_noise, x_off, y_off) 
     tmp_subfield_exp = image_shift(tmp_subfield_exp, x_off, y_off) 
     tmp_subfield_nsubims = image_shift(tmp_subfield_nsubims, x_off, y_off) 
      
      ; Sub-pixel shifts can lead to strange weights
      ; because of interpolation near edges.
      ; Discard those pixels!
      accept = where(tmp_subfield_nsubims gt 0.99, complement=reject)
      tmp_subfield_nsubims[reject] = 0
      tmp_subfield[reject] = 0
      tmp_subfield_noise[reject] = 0
      tmp_subfield_exp[reject] = 0

      im =  im + tmp_subfield 
      sigma =  sigma + (tmp_subfield_noise)^2
      exp =  exp + tmp_subfield_exp
      nsubims =  nsubims + tmp_subfield_nsubims
      
   endif ; FILE_TEST  
      endfor
   endfor
      
good = where(nsubims gt 0, complement=bad)
im[good] = im[good]/nsubims[good]      
exp[good] = exp[good]/nsubims[good]
sigma[good] = sigma[good]/nsubims[good]
im[bad] = 0  
sigma[bad] = 0
exp[bad] = 0      
  
writefits, imdir + Field + '_' + chip + '_holo' + '.fits.gz', im, /COMPRESS         
writefits, imdir + Field + '_' + chip + '_noise' + '.fits.gz', sqrt(sigma), /COMPRESS 
writefits, imdir + Field + '_' + chip + '_exp' +'.fits.gz', exp, /COMPRESS 
writefits, imdir + Field + '_' + chip + '_nsubims' +'.fits.gz', nsubims, /COMPRESS 

END
