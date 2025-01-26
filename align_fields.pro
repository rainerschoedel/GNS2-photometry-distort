pro align_fields, field_nr, chip_nr


; PURPOSE: Take the input lists from compute_uncertainties.pro
;      
;          Align all images and lists with the ones from the H-band 
;

field = strn(field_nr)
chip = 'chip' + strn(chip_nr)
max_distance = 2 ; Correponds to dmax in search for 
              ; common stars between a subimmages
              ;  for fine-alignment of subimages
              ; with rebfac = 2, tolerance = 2
              ; corresponds to half the resolution limit
 
; Numbers of first and last sub-images (usually not to be changed, see
; runholo.pro)
i_x0 = 0
i_x1 = 2

i_y0 = 0
i_y1 = 2

; input and output directories
basedir = '/home/data/GNS/2021/'

photdirh = basedir + 'H/' + field + '/photo/chip' + strn(chip_nr) + '/lists/'
imdirh = basedir + 'H/' + field + '/cubeims/chip' + strn(chip_nr) + '/'
photdirj = basedir + 'J/' + field + '/photo/chip' + strn(chip_nr) + '/lists/'
imdirj = basedir + 'J/' + field + '/cubeims/chip' + strn(chip_nr) + '/'
photdirk = basedir + 'Ks/' + field + '/photo/chip' + strn(chip_nr) + '/lists/'
imdirk = basedir + 'Ks/' + field + '/cubeims/chip' + strn(chip_nr) + '/'

for i_x = i_x0, i_x1 do begin
  for i_y = i_y0, i_y1 do begin 

     print, 'Now working on sub-image : ' + strn(i_x) + ', ' + strn(i_y)

     name = strn(i_x) + '_' +  strn(i_y)
     readcol, photdirh + 'list_' + name + '.txt', xh, yh, fh, sxh, syh, sfh, ch

     ; Determine offsets for J
     ; Apply offsets to lists and images
     ; and write them to file
     ; --------------------------
;     photdir = photdirj
;     imdir = imdirj
;     readcol, photdir + 'list_' + name + '.txt', x, y, f, sx, sy, sf, c
;     compare_lists, xh, yh, x, y, x1c, y1c, x2c, y2c, $
;           MAX_DISTANCE = max_distance, $
;           SUBSCRIPTS_1 = subc1, SUBSCRIPTS_2 = subc2
;     x_off = median(x1c - x2c)
;     y_off = median(y1c - y2c)
;     x = x + x_off
;     y = y + y_off
;     print, 'Offsets in x and y of J relative to H: ' + strn(x_off) + ', ' + strn(y_off)
;     forprint, TEXTOUT=photdir + 'offH_list_' + name + '.txt', x, y, f, sx, sy, sf, c, FORMAT='(7F13.3)', $
;              COMMENT='#x[pix]  y[pix]  f[ADU]  dx[pix]  dy[pix]  c0[corr]  output from align_fields.pro'

;     im = readfits(imdir + 'holo_' + name + '.fits.gz')
;     noise = readfits(imdir + 'holo_' + name + '_sigma.fits.gz')
;     wt = readfits(imdir + 'holo_' + name + '_wt.fits.gz')
;     writefits, imdir +  'offH_holo_' + name + '.fits.gz', image_shift(im,x_off,y_off)
;     writefits, imdir +  'offH_holo_' + name + '_sigma.fits.gz', image_shift(noise,x_off,y_off)
;     writefits, imdir +  'offH_holo_' + name + '_wt.fits.gz', image_shift(wt,x_off,y_off)

     ; Determine offsets for Ks
     ; Apply offsets to lists and images
     ; and write them to file
     ; --------------------------
;     photdir = photdirk
;     imdir = imdirk
;     readcol, photdir + 'list_' + name + '.txt', x, y, f, sx, sy, sf, c
;     compare_lists, xh, yh, x, y, x1c, y1c, x2c, y2c, $
;           MAX_DISTANCE = max_distance, $
;           SUBSCRIPTS_1 = subc1, SUBSCRIPTS_2 = subc2
;     x_off = median(x1c - x2c)
;     y_off = median(y1c - y2c)
;     x = x + x_off
;     y = y + y_off
;     print, 'Offsets in x and y of Ks relaative to H: ' + strn(x_off) + ', ' + strn(y_off)
;     forprint, TEXTOUT=photdir + 'offH_list_' + name + '.txt', x, y, f, sx, sy, sf, c, FORMAT='(7F13.3)', $
;              COMMENT='#x[pix]  y[pix]  f[ADU]  dx[pix]  dy[pix]  c0[corr]  output from align_fields.pro'

;     im = readfits(imdir + 'holo_' + name + '.fits.gz')
;     noise = readfits(imdir + 'holo_' + name + '_sigma.fits.gz')
;     wt = readfits(imdir + 'holo_' + name + '_wt.fits.gz')
;     writefits, imdir +  'offH_holo_' + name + '.fits.gz', image_shift(im,x_off,y_off)
;     writefits, imdir +  'offH_holo_' + name + '_sigma.fits.gz', image_shift(noise,x_off,y_off)
;     writefits, imdir +  'offH_holo_' + name + '_wt.fits.gz', image_shift(wt,x_off,y_off)

  endfor
endfor

END
