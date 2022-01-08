
    device, decompose=1
    red = sgcolor('red')

    nrec = 1000
    dr0 = 0.1
    dur = nrec*dr0
    uts = findgen(nrec)*dr0
    f0 = sin(uts)
    
    s0 = dr0
    s1 = dur
    dj = 1d/8
    j1 = floor(alog(s1/s0)/alog(2)/dj)  ; # of powers-of-two with dj
    s1 = s0*2d^(dj*j1)            
    ns = j1+1

    mor = wavelet(f0, dr0, /pad, s0=s0, dj=dj, j=j1, $
        mother='Morlet', param=w0, $
        period=ps, scale=ss, coi=coi)

    cdelta = 0.776d
    psi0 = !dpi^(-0.25)
    print, cdelta, psi0
    f1 = dblarr(nrec)
    for i = 0, ns-1 do f1 += real_part(mor[*,i]/sqrt(ss[i]))
    f1 *= (dj*sqrt(dr0))/cdelta/psi0
    
    mor = wavelet(f0, dr0, /pad, s0=s0, dj=dj, j=j1, $
        mother='Morlet', param=w0, $
        period=ps, scale=ss, coi=coi)

stop
    cdelta = 0.776d
    psi0 = !dpi^(0.25)
    print, cdelta, psi0
    f1 = dblarr(nrec)
    for i = 0, ns-1 do f1 += real_part(mor[*,i]/sqrt(ss[i]))
    f1 *= (dj*sqrt(dr0))/cdelta/psi0

    plot, f0
    oplot, f1, color=red

stop



; **** make up data for test.
T = 80d         ; 60 sec long.
dr0 = 0.001d     ; data rate 1 msec, or 1 kHz sample/sec.
nrec = long(T/dr0)+1    ; # of record.

ns = findgen(nrec)
ts = ns*dr0

; wave 1: 1.8 mV/m, 3 sec period, initial phase 0 rad.
a0 = 1
p0 = 2d
t0 = 0
x0s = a0*sin(ts*2*!dpi/p0+t0)

; wave 2: 0.5 mV/m, 11 sec period, initial phase 10 rad.
a1 = 1
p1 = 8d
t1 = 10
x1s = a1*sin(ts*2*!dpi/p1+t1)

; wave 2: 0.5 mV/m, 11 sec period, initial phase 10 rad.
a1 = 1
p1 = 4d
t1 = 10
x2s = a1*sin(ts*2*!dpi/p1+t1)


f0 = x0s+x1s+x2s


; set scales.
s0 = dr0
s1 = t
dj = 1d/8
j1 = floor(alog(s1/s0)/alog(2)/dj)  ; # of powers-of-two with dj
s1 = s0*2d^(dj*j1)
ns = j1+1


mor = wavelet(f0, dr0, /pad, s0=s0, dj=dj, j=j1, $
    mother='Morlet', param=w0, $
    period=ps, scale=ss, coi=coi)

cdelta = 0.776d
psi0 = !dpi^(-0.25)
f1 = dblarr(nrec)
for i = 0, ns-1 do f1 += real_part(mor[*,i]/sqrt(ss[i]))
f1 *= (dj*sqrt(dr0))/cdelta/psi0

plot, f0
oplot, f1, color=red




end
