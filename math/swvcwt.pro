;+
; Type: procedure.
; Purpose: Perform wavelet transform, and/or inverse transform and filtering.
; Parameters:
;   f0, in, dblarr[n], required. Input signal.
;   cwt, in/out, dblarr[n,m], optional if f0 is set, otherwise required.
;       Input/Output continuous wavelet transform.
;   waveform, out, dblarr[n,k], optional. Output filtered waveforms.
; Keywords:
;   family, in, string, optional. Can be ['gaussian','paul','morlet'].
;   order, in, int, optional. Order for certain wavelet.
;
;   scale = scl, in/out, type = dblarr[m], optional.
;       Set scales for transform or return the default scale.
;       Scale refers to number of record in Fourier period.
;
;   filter = ftr, in/out, type = dblarr[m], optional.
;       Set filter for inverse tranaform.
;       Filter refers to number of record.
;
;   inverse = inverse, in, type = boolean, optional.
;       Set inverse to do inverse transform or filter.
;
;   waveform = waveform, in, type = boolean, optional.
;       Set waveform to return cwt converted to wave form color-coded.
;
; Return:
;   return, out, dblarr[n] or dblarr[n,m].
; 
; Example:
;   t1 = swvcwt(f0).
;   f1 = swvcwt(f0, filter = recftr, /inverse).
;   f1 = swvcwt(t0, filter = recftr, /inverse).
; 
; Notes: The cwt is time-period spectrogram, NOT scalogram.
;   none.
; 
; Dependence:
;   none.
; 
; Author:
;   Sheng Tian.
; 
; History:
;   2013-03-26, Sheng Tian, create.
;-

function swvcwt, f0, order, family = wvtype, scale = recscl, $
    filter = ftr, inverse = inverse, waveform = waveform, _extra = extra
    compile_opt idl2

    ; check f0.
    sz = size(f0) & ndim = sz[0] & nrec = sz[1]

    ; check wavelet type.
    ; lamb is the "Fourier wavelength lambda" in table 1.
    ; ecoi is the "e-folding time tau_s' in table 1.
    ; cdel is C_delta * Psi_0(0) in table 2.
    if n_elements(wvtype) eq 0 then wvtype = 'gaussian'
    case wvtype of
        'morlet': begin
            info = wv_fn_morlet(order)      ; auto set undefined order.
            ecoi = info.efolding
            lamb = info.fourier_period
            cdel = 0.776D*!dpi^(-0.25D)
            dj = 0.6
        end
        'paul': begin
            info = wv_fn_paul(order)        ; auto set undefined order.
            ecoi = info.efolding
            lamb = info.fourier_period
            cdel = 1.132D*1.079
            dj = 1.5
        end
        'gaussian': begin 
            info = wv_fn_gaussian(order)    ; auto set undefined order.
            ecoi = info.efolding
            lamb = info.fourier_period
            case order of
                2: begin
                    cdel = 3.541D*0.867
                    psi0 = 9.867
                    dj = 1.4 & end
                6: begin
                    cdel = 1.966D*0.884
                    psi0 = 0.884
                    dj = 0.97 & end
            endcase
        end
    endcase

    ; check scale. s = s0*2^{ds*[0,1,...,N]}.
    if n_elements(recscl) ne 0 then begin
        nscl = n_elements(scl)
        s0 = recscl[0]/lamb
        ds = alog(recscl[1]/recscl[0])/alog(2)
    endif

    if ndim eq 1 then begin
        ; do transformation.
        cwt = wv_cwt(f0, wvtype, order, $
            start_scale = s0, dscale = ds, nscale = nscl, $
            scale = recscl, /pad, _extra = extra)
        if n_elements(scl0) eq 0 then scl = recscl*lamb $
        else scl = scl0
        if keyword_set(waveform) then begin
            if wvtype eq 'morlet' then $
                cwt = abs(cwt)*cos(atan(cwt, /phase))
            cwt *= (dj*sqrt(lamb)/cdel)
            for i = 0, nscl-1 do cwt[*,i] *= (1D/sqrt(scl[i]))
        endif
        if ~keyword_set(inverse) then return, cwt
    endif

    ; check ftr.
    if n_elements(ftr) eq 0 then ftr = scl[[0,nscl-1]]

    ; inverse transform or filter.
    if ndim eq 2 then cwt = f0
    idx = where(scl ge ftr[0] and scl le ftr[1])
    if idx[0] eq -1 then message, 'wrong filter ...'
    return, (dj*sqrt(lamb)/cdel)*(cwt[*,idx] # (1D/sqrt(scl[idx])))

end

; f(t) = sin(x)+sin(6*x).
; period: 10 sec and 60 sec.
nrec = 12001
dr = 0.01    ;sec
x = 4*!dpi/(nrec-1)*findgen(nrec)
f0 = sin(x)+sin(6*x)
t0 = time_double('1997-05-27/18:00')+findgen(nrec)*dr   ; 2 min.

nscl = 50
s0 = 4
dj = alog((0.125D*nrec)^(1D/(nscl-1)))/alog(2)

c = 2*!dpi/sqrt(2.5)    ; for dog m = 2.

; idl wv_cwt.
cwt = wv_cwt(f0, 'gaussian', $
;    start_scale = s0, dscale = dj, nscale = nscl, $
    scale = recscl, /pad)
timescl = recscl*c*dr
store_data, 'idlcwt', data = {x:t0, y:cwt, v:timescl}

; sheng swvcwt.
recscl = recscl*c
cwt = swvcwt(f0, scale = recscl)
timescl = recscl*dr
store_data, 'swvcwt', data = {x:t0, y:cwt, v:timescl}

f1 = swvcwt(cwt, scale = recscl, filter = [0,15]/dr, /inverse)
store_data, 'fcwt', data = {x:t0, y:f1}

; lynn wavelet.
cwt = wavelet(f0, mother = 'dog', dr, dj = 0.25, $
;    s0 = s0, dj = dj, j = nscl-1, $
    scale = recscl, /pad)
timescl = recscl*c
store_data, 'lyncwt', data = {x:t0, y:real_part(cwt), v:timescl}

vars = ['idlcwt', 'swvcwt', 'lyncwt']
options, vars, 'spec', 1
options, vars, 'x_no_interp', 1
options, vars, 'y_no_interp', 1
ylim, vars, 0.2, 200, 1
zlim, vars, -60, 60, 0

device, decomposed = 0
loadct, 39

vars = ['idlcwt', 'swvcwt', 'lyncwt', 'fcwt']
tplot, vars

end
