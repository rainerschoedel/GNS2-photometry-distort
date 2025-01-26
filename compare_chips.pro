PRO COMPARE_CHIPS, band, field

; PURPOSE: Compare chips pairwise to check 
;          the magnitude of systematic offsets

field = strn(field)
n_fields = 4 ; number of fields to compare

maxmag = 14.0

minmag = 18 ; comparison stars must be brighter than minmag
if band eq 'Ks' then minmag = 16


SIGMA_CUT = 3.0 ; parameter for RESISTANT_MEAN

; read astrometry structure
;im = readfits('/home/data/working/HAWK-I/2018/J/' + field + '/photo/chip1/Field' + field + '_chip1holo_2_cal.fits', new_header)
;EXTAST, new_header, astr   

dmax = 0.02 ; max distance for list comparison in arcseconds
magerr = 0.03 ; valid comparison stars must have smaller magnitude uncertainties

resdir = '/home/data/GNS/2021/' + band + '/' + field + '/photo/'
tmpdir = '/home/data/GNS/2021/' + band + '/' + field + '/tmp/'

openw, out, resdir + 'offsets_' + band + '.txt', /get_lun


row_nr = 0L  ; refers to row in N x N matrix of Dong et al. (2011), equ. (A2)

for ix = 1, n_fields do begin
   
    readcol, resdir + 'chip' + strn(ix) + '/lists/stars_calibrated_' + band + '_chip' + strn(ix) + '.txt',  a_row ,d_row, x, y, f, m_row, sx, sy, sf, sm_row, FORMAT='D,D,D,D,D,D,D,D,D,D';, /SILENT
    ord = reverse(sort(f))
    a_row = a_row[ord]
    d_row = d_row[ord]
    m_row = m_row[ord]
    sm_row = sm_row[ord]

    good = where(sm_row lt magerr and m_row lt minmag and m_row gt maxmag, nlist1) 
    a_row = a_row[good]
    d_row = d_row[good]
    m_row= m_row[good]
    sm_row = sm_row[good]

    ; Compare list with lists from the other fields
    ;----------------------------------------------------------


    for j = ix+1, n_fields do begin

      if (not(j eq ix)) then begin
       readcol, resdir + 'chip' + strn(j) + '/lists/stars_calibrated_' + band + '_chip' + strn(j) + '.txt',  a_col ,d_col, x, y, f, m_col, sx, sy, sf, sm_col, FORMAT='D,D,D,D,D,D,D,D';, /SILENT
       ord = reverse(sort(f))
       a_col = a_col[ord]
       d_col = d_col[ord]
       m_col = m_col[ord]
       sm_col = sm_col[ord]

        good = where(sm_col lt magerr and m_col lt minmag and m_col gt maxmag,nlist2) 
        a_col = a_col[good]
        d_col = d_col[good]
        m_col = m_col[good]
        sm_col = sm_col[good]

        ncommon = 0L
        subc1 = [] & subc2 = []
        for l = 0, nlist1-1 do begin
          d = sphdist(a_row[l],d_row[l],a_col,d_col, /DEGREES) * 3600.0
          good = where(d lt dmax,ngood)
          if (ngood eq 1) then begin
            ncommon = ncommon + 1
            if ncommon eq 1 then begin
              subc1 = [l]
              subc2 = [good[0]]
            endif else begin
              subc1 = [subc1,l]
              subc2 = [subc2,good[0]]
            endelse
           endif else if (ngood gt 1) then begin
              print, 'Several matches...'
              print, ngood
              print, a_col[good]
              print, d_col[good]
              print, f[good]
              print, 'Using the one closest in flux.'
              dm = abs(m_col[good] - m_row[l])
              m_min = min(dm,min_sub)
              print, 'Minimum magnitude difference is: ' + strn(m_min)
              subc1 = [subc1,l]
              subc2 = [subc2,good[min_sub]]
           endif
        endfor
        print, 'Fields ' + strn(ix) + ' and ' + strn(j)

        if ncommon gt 1 then begin


         RESISTANT_Mean,(a_row[subc1] - a_col[subc2]),SIGMA_CUT,Mean,Sigma,Num_Rej
          a_off = Mean * 3600.
          sa_off = Sigma * 3600.
          
          RESISTANT_Mean,(d_row[subc1] - d_col[subc2]),SIGMA_CUT,Mean,Sigma,Num_Rej
          d_off = Mean * 3600.
          sd_off = Sigma * 3600.
          
          RESISTANT_Mean, (m_row[subc1] - m_col[subc2]),SIGMA_CUT,Mean,Sigma,Num_Rej, goodvec= goodvec
          m_off = mean
          sm_off = sigma

          print, 'Number of stars used: ' + strn(n_elements(goodvec))
          print, 'Offsets: ' + strn(a_off,FORMAT='(F9.4)') + '+-' + strn(sa_off,FORMAT='(F9.4)') + ', ' + strn(d_off,FORMAT='(F9.4)') + '+-' + strn(sd_off,FORMAT='(F9.4)') + ', ' + strn(m_off,FORMAT='(F9.4)') + '+-' + strn(sm_off,FORMAT='(F9.4)')

          printf, out, ix, j, a_off, sa_off, d_off, sd_off, m_off, sm_off, FORMAT='(2I,6F8.3)'

          ; Create ds0 regions
          openw, l1, resdir + 'Fields_' + strn(ix) + '_and_' + strn(j) + '.reg', /get_lun
          printf, l1,'# Region file format: DS9 version 4.1'
          printf, l1, 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
          printf, l1,'fk5'
          for i=0,n_elements(goodvec)-1 do begin
            astrn =  strn(a_row[subc1[goodvec[i]]])
            dstrn = strn(d_row[subc1[goodvec[i]]])
            printf, l1, 'circle('+ astrn + ',' + dstrn + ',0.5")'
          endfor
          free_lun, l1
          ; Ploto histogram of differences
          ; -------------------------------
         set_plot,'PS',/interpolate
         device, XOFFSET=0, YOFFSET=0, $
         FILENAME= tmpdir + 'Comp_' + strn(ix) + '_' + strn(j) + '.eps', XSIZE=20, YSIZE=15., $
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
         dm = m_row[subc1] - m_col[subc2]
         binsize = 0.01
         cgHistoplot, dm, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'ZP', xrange = [-0.1,0.1], xticks = 4, YTITLE='N', TITLE='Fields ' + strn(ix) + ' and '  + strn(j)
         binCenters1 = loc + (binsize / 2.0)
  ;       yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
  ;       cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  ;       print, '*****'
  ;       print, 'Gauss fit to histogram: ' + strn(coeff[1]) + ' +- ' +strn(coeff[2]) + ' ' + strn(coeff[2]/sqrt(n_elements(dm)))
  ;       print, '*****'
        device, /close

        endif
      endif
      ;    endfor
    endfor
    ; endfor
  endfor

  free_lun, out


END

