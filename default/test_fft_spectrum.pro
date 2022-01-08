
; test fft spectrum.


dr = 1d/800
uts = smkarthm(0,8,dr, 'dx')
nrec = n_elements(uts)

fs = findgen(nrec)/(nrec*dr)

freq = 200
dat = sin(2*!dpi*uts*freq)+sin(2*!dpi*uts*0.25*freq)+randomn(1,nrec)*0.5

fp = abs(fft(dat))

plot, uts, dat
stop

plot, fs[0:nrec/2], fp[0:nrec/2]

; peak at 50 and 200 Hz.

end
