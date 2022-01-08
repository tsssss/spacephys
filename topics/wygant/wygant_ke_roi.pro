
; test time.
utr = ['2013-05-01/07:00','2013-05-01/08:00']
prb = 'b'

utr = ['2012-11-14/04:00','2012-11-14/05:00']
prb = 'b'

utr = time_double(utr)
pre0 = 'rbsp'+prb+'_'

rootdir = srootdir()
fn = rootdir+'/hope_l3_ion.dat.tplot'

if file_test(fn) eq 0 then begin
    vars = []
    ; read hope l3.
    hopel3 = sread_rbsp_hope_l3(utr, probe = prb)
    
    uts = sfmepoch(hopel3.epoch_ele,'unix')
    val = hopel3.hope_energy_ele
    
    lim = {spec:1, no_interp:1, ylog:1, zlog:1, $
        ystyle:1, yrange:[20,5e4], zrange:[1e4,1e10], $
        ytitle:'RBSP-'+strupcase(prb)+'!CHope!Ce (eV)'}
    dat = total(hopel3.fedu, 3, /nan)
    tvar = pre0+'e_spec'
    store_data, tvar, uts, dat, val, limits = lim
    vars = [vars,tvar]
    
    uts = sfmepoch(hopel3.epoch_ion,'unix')
    val = hopel3.hope_energy_ion
    
    lim = {spec:1, no_interp:1, ylog:1, zlog:1, $
        ystyle:1, yrange:[20,5e4], zrange:[1e4,1e7], $
        ytitle:'RBSP-'+strupcase(prb)+'!CHope!CH (eV)'}
    dat = total(hopel3.fpdu, 3, /nan)
    tvar = pre0+'h_spec'
    store_data, tvar, uts, dat, val, limits = lim
    vars = [vars,tvar]
    
    lim = {spec:1, no_interp:1, ylog:1, zlog:1, $
        ystyle:1, yrange:[20,5e4], zrange:[1e4,1e7], $
        ytitle:'RBSP-'+strupcase(prb)+'!CHope!CO (eV)'}
    dat = total(hopel3.fodu, 3, /nan)
    tvar = pre0+'o_spec'
    store_data, tvar, uts, dat, val, limits = lim
    vars = [vars,tvar]

;    ; plot each angle.
;    pas = [4.5,18,36,175.5,162,144]
;    pas = pas[sort(pas)]
;    npa = n_elements(pas)
;    vars = []
;    for i = 0, npa-1 do begin
;        ; electrons.
;        uts = sfmepoch(hopel3.epoch_ele,'unix')
;        val = hopel3.hope_energy_ele
;        idx = where(hopel3.pitch_angle eq pas[i])
;        
;        lim = {spec:1, no_interp:1, ylog:1, zlog:1, $
;            ystyle:1, yrange:[10,5e4], $
;            ytitle:'RBSP-'+strupcase(prb)+'!CHope!Ce (eV)!CPA='+$
;            sgnum2str(pas[i])+' deg'}
;        dat = reform(hopel3.fedu[*,*,idx])
;        tvar = pre0+'e_spec_'+string(idx,format='(I0)')
;        store_data, tvar, uts, dat, val, limits = lim
;        vars = [vars,tvar]
;
;        ; ions.
;        uts = sfmepoch(hopel3.epoch_ion,'unix')
;        val = hopel3.hope_energy_ion
;        idx = where(hopel3.pitch_angle eq pas[i])
;        
;        lim = {spec:1, no_interp:1, ylog:1, zlog:1, $
;            ystyle:1, yrange:[10,5e4], zrange:[1e4,1e7], $
;            ytitle:'RBSP-'+strupcase(prb)+'!CHope!CH (eV)!CPA='+$
;            sgnum2str(pas[i])+' deg'}
;        dat = reform(hopel3.fpdu[*,*,idx])
;        tvar = pre0+'p_spec_'+string(idx,format='(I0)')
;        store_data, tvar, uts, dat, val, limits = lim
;        vars = [vars,tvar]
;        
;        lim = {spec:1, no_interp:1, ylog:1, zlog:1, $
;            ystyle:1, yrange:[10,5e4], zrange:[1e3,1e6], $
;            ytitle:'RBSP-'+strupcase(prb)+'!CHope!CO (eV)!CPA='+$
;            sgnum2str(pas[i])+' deg'}
;        dat = reform(hopel3.fodu[*,*,idx])
;        tvar = pre0+'o_spec_'+string(idx,format='(I0)')
;        store_data, tvar, uts, dat, val, limits = lim
;        vars = [vars,tvar]
;
;        lim = {spec:1, no_interp:1, ylog:1, zlog:1, $
;            ystyle:1, yrange:[10,5e4], zrange:[1e3,1e6], $
;            ytitle:'RBSP-'+strupcase(prb)+'!CHope!CHe (eV)!CPA='+$
;            sgnum2str(pas[i])+' deg'}
;        dat = reform(hopel3.fhedu[*,*,idx])
;        tvar = pre0+'he_spec_'+string(idx,format='(I0)')
;        store_data, tvar, uts, dat, val, limits = lim
;        vars = [vars,tvar]
;    endfor

    tplot_save, vars, filename = fn
endif else begin
    tplot_restore, filename = fn
endelse


;get_data, 'rbspb_o_spec_8', t0, dat, val, limits = lim
;tdat = alog10(dat)
;zrng = lim.zrange
;tdat = bytscl(tdat, min = alog10(zrng[0]), max = alog10(zrng[1]))

;xsz = 500*2
;ysz = 300*2
;tdat = congrid(tdat,xsz,ysz)

device, decompose = 0
loadct2, 43
tvlct, rr, gg, bb, /get

tplot, 'rbspb_'+['h_spec','o_spec','e_spec']
stplot_add_line, vxs, vys



stop

end
