
np = 20         ; # of cycles.
nrec = 1e4      ; # of rec.
p0 = nrec/np    ; period in rec.
c0 = 109          ; phase difference in rec.
r0 = findgen(nrec)/(nrec-1)
xs = nrec/p0*2*!dpi*r0
fs = sin(xs + c0/p0*2*!dpi)    ; f(x) = sin(2*pi*np*t)

; down sample.
ps = p0*1.08    ; sample period in rec, slower than real wave.
idx = smkarthm(0,nrec, ps, 'dx')

sgopen, 0

plot, xs, fs, xstyle = 1


f1s = interpol(fs[idx],xs[idx],xs,/spline)
;oplot, xs, f1s, color = sgcolor('red')
oplot, xs[idx], fs[idx], psym = 1, color = sgcolor('red')

pa = abs(1d/(1d/ps-1d/p0))
f2s = sin(nrec/pa*2*!dpi*r0 + c0/p0*2*!dpi)
oplot, xs, f2s, color = sgcolor('blue')

ofn = shomedir()+'/linear_vs_log_scale.pdf'
;ofn = 0
sgopen, ofn, xsize = 5, ysize = 7, /inch

txs = smkarthm(0,1, nrec, 'n')
tys = smkarthm(16384,1e-1, nrec, 'n')

poss = sgcalcpos(2, ypad = 8)
plot, txs, tys, ylog = 1, xrange = [-0.2,1.2], yrange = [100, 2e4], $
    position = poss[*,0], /noerase, $
    xstyle = 1, ystyle = 1, title = 'Log scale', $
    xtitle = 'Time', ytitle = 'Freq (Hz)'
plot, txs, tys, ylog = 0, xrange = [-0.2,1.2], yrange = [100, 2e4], $
    position = poss[*,1], /noerase, $
    xstyle = 1, ystyle = 1, title = 'Linear scale', $
    xtitle = 'Time', ytitle = 'Freq (Hz)'
    
sgclose

end