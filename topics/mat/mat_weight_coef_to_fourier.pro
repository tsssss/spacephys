; show that at higher order, the weight resembles the
; frequency response of a delta function in frequency space.

ofn = shomedir()+'/weight_coef.pdf'
;ofn = 0

w = 5d          ; width of smoothing.
wp = (w-1)/2    ; half width.
wc = 1          ; center width.

p = 1d/w
norder = 4      ; order.
cptrs = ptrarr(norder,/allocate_heap)   ; the weight coef of each order.
mptrs = ptrarr(norder,/allocate_heap)   ; the matrix of each order.

mx0 = [dblarr(wp),dblarr(wc)+1,dblarr(wp)]-p
for i = 0, norder-1 do begin
    nr = 1+i*2*wp
    nc = w+i*2*wp
    tmx = dblarr(nc,nr)
    for j = 0, nr-1 do tmx[j,j] = mx0
    *mptrs[i] = tmx
endfor

*cptrs[0] = mx0
for i = 1, norder-1 do *cptrs[i] = (*cptrs[i-1]) ## (*mptrs[i])


sgopen, ofn, xsize = 3, ysize = 4, /inch


tpos = [0.2,0.15,0.95,0.90]
poss = sgcalcpos(norder, position = tpos)

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

xr = [-1,1]*9
xtickv = smkarthm(xr[0],xr[1],4,'n')
xticks = n_elements(xtickv)-1
xminor = 6
xticklen = -0.02

yr = [-0.6,1.2]
ytickv = [0,1]
yticks = n_elements(ytickv)-1
yminor = 5
yticklen = -0.01
ytitle = ''

gridcolor = sgcolor('grey')
gridstyle = 2

letters = ['a','b','c','d','e','f','g']

for i = 0, norder-1 do begin
    nrec = wc+(i+1)*2*wp
    tx = findgen(nrec)-(i+1)*wp
    ty = *cptrs[i]
    tpos = poss[*,i]
    
    ; set-up coord and box.
    plot, tx, ty, position = tpos, /noerase, /nodata, $
        xrange = xr, xstyle = 1, xtickformat='(A1)', $
        xtickv=xtickv, xticks=xticks, xminor=xminor, xticklen=xticklen, $
        yrange = yr, ystyle = 1, ytitle = ytitle, $
        ytickv=ytickv, yticks=yticks, yminor=yminor, yticklen=yticklen
    ; plot the weight.
    oplot, tx, ty, psym = -1, symsize=0.5, color=sgcolor('red')
    
    xyouts, tpos[0]-xchsz*5, tpos[3]-ychsz*1, /normal, letters[i]+'.'
    xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*1.2, /normal, 'order='+sgnum2str(i+1)
    
    ; plot the theoretical sinc.
    u = 8/1.5/!dpi
    u = 8/!dpi
    u = 7.5/!dpi
    u = 2.5
    ; interpolate to avoid aliasing.
    ttx = smkarthm(min(tx),max(tx),1000,'n')
    ttx = smkarthm(xr[0], xr[1] ,1000,'n')
    tty = sin(ttx*u)/(ttx*u)
    oplot, ttx, tty, color=gridcolor, thick=0.5
    oplot, xr, [0,0], linestyle = 2, color=gridcolor


endfor

plots, (poss[0,0]+poss[2,0])*0.5+[0,0], [poss[1,norder-1],poss[3,0]], /normal, $
    linestyle=gridstyle, color=gridcolor
tpos = poss[*,norder-1]
for i=0, xticks do xyouts, /normal, alignment=0.5, $
    tpos[0]+(tpos[2]-tpos[0])/xticks*i, tpos[1]-ychsz*1, sgnum2str(xtickv[i])
xyouts, /normal, alignment=0.5, (tpos[0]+tpos[2])*0.5, tpos[1]-ychsz*2.5, $
    'Number of record'
tpos = poss[*,0]
xyouts, /normal, alignment=0.5, (tpos[0]+tpos[2])*0.5, tpos[3]+ychsz*0.5, $
    'Weight of '+sgnum2str(w)+'-point MAT'


sgclose




for i = 0, norder-1 do begin
    ptr_free, cptrs[i]
    ptr_free, mptrs[i]
endfor

end
