;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

pro rbsp_hope_decompose_sector, vname, id, newname = newname
    
    get_data, vname, t0, v0, e0, limit = lim
    
    nrec = n_elements(t0)

    nen = n_elements(e0)
    if nen/nrec ne 0 then nen = nen/nrec
    
    ; pixel id: 0,1,2,3,4.
    ; # of sec: 4,8,16,8,4.
    nsec = 16/2^abs(id-2)
    nunt = nsec+2
    
    t1 = dblarr(nrec*nunt)
    v1 = dblarr(nrec*nunt,nen)
    e1 = dblarr(nrec*nunt,nen)
    
    secidx = fix(smkarthm(0,16/nsec,nsec,'x0'))
    dtsec = 12d/nsec
    
    for i = 0ull, nrec-1 do begin
        idx1 = i*nunt
        idx2 = (i+1)*nunt-1
        t1[idx1:idx2] = smkarthm(t0[i]-dtsec,dtsec,nunt,'x0')
        v1[idx1+1:idx2-1,*] = transpose(reform(v0[i,*,secidx,id]))
        v1[[idx1,idx2]] = !values.d_nan
        e1[idx1+1:idx2-1,*] = e0[i,*] ## (intarr(nsec)+1)
        e1[[idx1,idx2]] = !values.d_nan
    endfor

    if n_elements(newname) eq 0 then newname = vname+'_px'+string(id,format='(I0)')
    store_data, newname, t1, v1, e1, limits = lim
end

utr = time_double(['2013-03-14/00:00','2013-03-14/00:10'])
type = 'level2'
prob = 'a'

hope = sread_rbsp_hope(utr)

ion_uts = sfmepoch(hope.epoch_ion, 'unix')
nrec = n_elements(ion_uts)

uts = sfmepoch(hope.epoch, 'unix')
ion_idx = intarr(nrec)
for i = 0, nrec-1 do ion_idx[i] = where(uts eq ion_uts[i])

fpdu = hope.fpdu
ion_ens = hope.hope_energy_ion
;    detectors = *hope[where(vars eq 'HOPE_DETECTOR')].value
;    sectors = *hope[where(vars eq 'HOPE_SECTOR')].value

tvar = 'fpdu'
store_data, tvar, ion_uts, fpdu, ion_ens
options, tvar, 'spec', 1
options, tvar, 'no_interp', 1
options, tvar, 'zlog', 1
ylim, tvar, 1, 1e5, 1
for i = 0, 4 do rbsp_hope_decompose_sector, tvar, i


tvar = 'en'
tmp = total(total(fpdu,3),3)
store_data, tvar, ion_uts, tmp, ion_ens
options, tvar, 'spec', 1
options, tvar, 'no_interp', 1
options, tvar, 'zlog', 1
ylim, tvar, 1, 1e5, 1

ion_collaps = (hope.sector_collapse_cntr)[ion_idx,*]
en_collaps = (hope.energy_collapsed)[ion_idx]

get_data, 'fpdu', t0, fpdu
get_data, 'h_en', t0, en

; load L, MLT
tvar = 'lshell'
ion_lshell = hope.l_ion
store_data, tvar, ion_uts, ion_lshell
options, tvar, 'ytitle', 'L'

tvar = 'mlt'
ion_mlt = hope.mlt_ion
store_data, tvar, ion_uts, ion_mlt
options, tvar, 'ytitle', 'MLT'


; load B field.
emfisis = sread_rbsp_emfisis(utr)
vars = emfisis.name

b_uts = sfmepoch(*emfisis[where(vars eq 'Epoch')].value,'unix',/tt2000)
b_gse = *emfisis[where(vars eq 'Mag')].value

store_data, 'b_gse', b_uts, b_gse
tvar = 'b_gse'
options, tvar, 'labels', ['x','y','z']

rbsp_uvw2gse, 'b_gse', newname = 'b_uvw', probe = prob, /inverse
tvar = 'b_uvw'
options, tvar, 'labels', ['u','v','w']
get_data, tvar, t0, dat
store_data, tvar+'_mag', t0, snorm(dat), limits = {ynozero:1}



; load s/c pos.
probe = 'a'
coord = 'gse'
rbsp_load_spice_state, probe = probe, coord = coord, times = b_uts, /no_spice_load
tvar = 'r_gse'
stplot_renew, 'rbsp'+probe+'_state_pos_'+strlowcase(coord), newname = tvar
get_data, tvar, t0, tmp
tmp*= (1d/6378)
store_data, tvar, t0, tmp
options, tvar, 'labels', 'GSE '+['x','y','z']+' (Re)'



vars = ['b_gse','b_uvw','r_gse']
options, vars, 'colors', [6,4,2]

vars = ['fpdu_px?','b_gse','b_uvw','b_uvw_mag']
labs = ['lshell','mlt']
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'xmargin', [20,10]
tplot, vars, var_label = labs, title = 'RBSP-'+strupcase(probe)+' HOPE Proton', trange = utr

end
