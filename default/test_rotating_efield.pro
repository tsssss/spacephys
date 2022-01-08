
ofn = 0
sgopen, ofn, xsize = 5, ysize = 5, /inch

sgtruecolor
white = sgcolor('white')
red = sgcolor('red')
blk = sgcolor('black')

e0 = 1
r0 = 1
v0 = 2*!dpi*r0*0.05

; where to draw E field.
nrec = 30
xs = smkarthm(-2,2,nrec,'n')
ys = smkarthm(-2,2,nrec,'n')

; times.
nt = 50
ts = smkarthm(0,1,nt,'n')

; electron orbit.
rxs = r0*cos(2*!dpi*ts)
rys = r0*sin(2*!dpi*ts)
vxs =-v0*sin(2*!dpi*ts)
vys = v0*cos(2*!dpi*ts)

; user symbol.
tmp = findgen(11)*2*!dpi/10
txs = cos(tmp)
tys = sin(tmp)
usersym, txs, tys, color = red, /fill


for i = 0, nt-1 do begin
    exs = fltarr(nrec,nrec)+e0*cos(2*!dpi*ts[i])
    eys = fltarr(nrec,nrec)+e0*sin(2*!dpi*ts[i])
    
    velovect, exs, eys, xs, ys, /isotropic, color = blk, background = white
    oplot, rxs, rys, color = red, linestyle = 2
    plots, rxs[i], rys[i], psym = 8, symsize = 0.8
    arrow, rxs[i], rys[i], rxs[i]+vxs[i], rys[i]+vys[i], /data, $
        color = red, hsize = 5
    wait, 0.1
endfor

sgclose

end
