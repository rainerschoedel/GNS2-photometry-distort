FUNCTION GRIDSELECT, xpos, ypos, flux, gridsize=gridsize

; From an input list, define a grid and select 
; one star inside each grid box


; Define grid
; =============
xlo = round(min(xpos))
xhi = round(max(xpos))
ylo = round(min(ypos))
yhi = round(max(ypos))

xw = xhi - xlo
yw = yhi - ylo
width = xw > yw
step = round(width/float(gridsize))

; Select a star inside of each grid box
; ======================================

n_grid = 0
grid_idx = []
x0 = xlo
x1 = xlo + step
while (x1 lt xhi) do begin
  y0 = ylo
  y1 = ylo + step
  while (y1 lt yhi) do begin
     in_box = where((xpos ge x0) and (xpos lt x1) and (ypos ge y0) and (ypos lt y1),count)
     if (count gt 0) then begin
       fbright = max(flux[in_box],idx_bright)
       grid_idx = [grid_idx, in_box[idx_bright]]
     endif
     y0 = y1
     y1 = y1 + step
     n_grid = n_grid + 1
  endwhile
  x0 = x1
  x1 = x1 + step
endwhile

print, 'Wodth of x- and y-axes: ' + strn(xw) + ', ' + strn(yw)
print, 'Number of grid points: ' + strn(n_grid)
print, 'Number of grid points with stars: ' + strn(n_elements(grid_idx))

return, grid_idx

END
