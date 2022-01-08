; leap frog solver for lorentz equation.


pro test_leap_frog, rs, vs, q = q, m = m

; initial E/B field.
e0 = [1d,0,0]
b0 = [0,0,1d]

; charge and mass.
if n_elements(q) eq 0 then q = 1d
if n_elements(m) eq 0 then m = 1d

; step info.
t1 = 50d
dt = abs(2d*!dpi/(q*snorm(b0)/m)/(50))
nt = t1/dt
r0 = [0d,0,0]   ; initial position.
v0 = [0d,0,0]   ; initial velocity.

rs = dblarr(nt,3)
vs = dblarr(nt,3)

rs[0,*] = r0
vs[0,*] = v0

for i = 1, nt-1 do begin
    tb = b0
    te = e0
    vo = reform(vs[i-1,*])  ; old velocity.
    ro = reform(rs[i-1,*])  ; old position.
    tomega = q*tb/m
    tsigma = q*te/m
    ta = tomega*(0.5*dt)
    tc = vo+dt*(tsigma+scross(vo,tomega*0.5))
    vs[i,*] = (tc+sdot(ta,tc)*ta-scross(ta,tc))/(sdot(ta,ta)+1d)
    rs[i,*] = ro+vs[i,*]*dt
endfor

;p = plot3d(rs[*,0],rs[*,1],rs[*,2], $
;    axis_style=2, /perspective, xy_shadow=1, yz_shadow=1, xz_shadow=1, $
;    xtitle = 'x', ytitle = 'y', ztitle = 'z', shadow_color='deep sky blue')
;
;ax = p.axes
;ax[2].hide = 1
;ax[6].hide = 1
;ax[7].hide = 1

end


poss = sgcalcpos(1,2)

test_leap_frog, ers, evs, q =-1, m = 1
test_leap_frog, irs, ivs, q = 1, m = 4

ofn = 0
sgopen, ofn, xsize = 4, ysize = 3, /inch

xr = [-12,0]
yr = [-28,0]
plot, ers[*,0], ers[*,1], /noerase, position = poss[*,0], $
    title = 'Electron', /isotropic, $
    xstyle = 1, xtitle = 'x', xrange = xr, xticks = 2, xminor = 6, $
    ystyle = 1, ytitle = 'y';, yrange = yr

xr = [0,12]
yr = [-28,0]
plot, irs[*,0], irs[*,1], /noerase, position = poss[*,1], $
    title = 'Ion', /isotropic, $
    xstyle = 1, xtitle = 'x', xrange = xr, xticks = 2, xminor = 6, $
    ystyle = 1, ytitle = 'y';, yrange = yr

sgclose
    
end