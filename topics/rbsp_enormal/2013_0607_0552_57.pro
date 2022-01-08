


; settings.
utr = time_double(['2013-06-07/04:52:57.00','2013-06-07/04:52:57.90'])
utr0 = time_double(['2013-06-07/04:50','2013-06-07/04:55'])   ; pad time.

tprobe = 'a'
pre0 = 'rbsp'+tprobe+'_'
b1dr = 1d/4096
period0 = 20e-3 ; sec.
boomlens = [100d,100,12]    ; tip-to-tip, m.
freq0 = 1d/period0
omega0 = 2*!dpi*freq0
load = 0
reload = 1


; **** local proton cyclotron frequency is eB/mp = 22.86 Hz, |B| = 238.6 nT.
; wave frequency is 50 Hz, omega = 314 rad/s, k = 0.077*(-0.91,-0.22,0.32) in GSE.
; phase velocity is 422.2 km/s.
; use n = 1 cc, Ti = 10 keV, Te = 1 keV,
; Va = 5249 km/s, Vthi = 979 km/s, Vthe = 41931 km/s.


; one settings.
vars = ['b_gse','pos_gse','mlt','lshell','dis','vsc','vsvy','vb1', $
    'vsvy_e_uvw','vb1_e_uvw','m_uvw2gse']
datvars = pre0+vars    ; basic variables.


device, decomposed = 0
loadct2, 43
tplot_options, 'constant', 0
tplot_options, 'labflag', -1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1



datfn = shomedir()+'/psbl_de_32hz/2013_0607_0552_57_data.tplot'
if tnames(datvars[0]) eq '' then if file_test(datfn) ne 0 then begin
    store_data, '*', /delete
    tplot_restore, filename = datfn
endif
foreach tvar, datvars do if tnames(tvar) eq '' then reload = 1
if load eq 1 then reload = 1

if reload eq 1 then begin
    if file_test(datfn) ne 0 then begin
        store_data, '*', /delete
        tplot_restore, filename = datfn
    endif
    defsysv,'!rbsp_spice', exists=flag
    if flag eq 0 then rbsp_load_spice_kernels, trange = utr0, probes = tprobe
endif



if reload ne 0 then begin
    store_data, '*', /delete
    
; B in 32Hz.
    emfisis = sread_rbsp_emfisis_l3(utr0, type = 'hires', probe = tprobe, coord = 'gse')
    tvar = pre0+'b_gse'
    store_data, tvar, sfmepoch(emfisis.epoch,'unix'), emfisis.mag, $
        limit = {ytitle:'B GSE!C(nT)', labels:['x','y','z'], colors:[6,4,2]}
    get_data, tvar, uts, dat
    idx = where(uts ge utr0[0] and uts le utr0[1])
    store_data, tvar, uts[idx], dat[idx,*]
    
    
    
; position in 32Hz.
    tvar = pre0+'pos_gse'
    tuts = uts[idx]
    rbsp_load_spice_state, probe = tprobe, coord = 'gse', times = tuts, /no_spice_load
    get_data, pre0+'state_pos_gse', tmp, posgse
    posgse = posgse*re1
    tdat = posgse
    store_data, pre0+'state_*', /delete
    
    store_data, tvar, tuts, tdat, $
        {ytitle:'R GSE (Re)', colors:[6,4,2], labels:'GSE '+['x','y','z']}
        
        
    tvar = pre0+'mlt'
    tdat = atan(posgse[*,1],posgse[*,0])*deg
    tdat = (tdat+360) mod 360   ; convert to 0-360.
    tdat = (tdat/15 + 12) mod 24
    store_data, tvar, tuts, tdat, limits = {ytitle:'MLT (hr)'}
    
    
    tvar = pre0+'lshell'
    possm = sgse2sm(posgse, stoepoch(tuts,'unix'))
    mlat = atan(possm[*,2],sqrt(possm[*,0]^2+possm[*,1]^2)) ; in rad.
    dis = sqrt(possm[*,0]^2+possm[*,1]^2+possm[*,2]^2)
    tdat = dis/(cos(mlat)^2)
    store_data, tvar, tuts, tdat, limits = {ytitle:'L-shell'}
    
    store_data, pre0+'dis', tuts, dis
    
    
    
; vsc in 32Hz.

    tvar = pre0+'vsc'
    dat = sread_rbsp_efw_l2(utr0, probes = tprobe, type = 'vsvy')
    uts = sfmepoch(dat.epoch, 'unix')
    nrec = n_elements(uts)
    vsvy = dat.vsvy
    ; calibrate vsvy, tried rbsp_efw_get_cal_params, but really no effect
    ; since offset is 0, gain is not used.
    ; end up using V12, may consider more complex ways.
    vsc = mean(vsvy[*,0:1], dimension = 2)
    store_data, tvar, uts, vsc, limits = {ytitle:'Vsc!C(V)', labels:'Vsc'}
    store_data, pre0+'vsvy', uts, vsvy, limits = $
        {ytitle:'Vsvy!C(V)', colors:[1,2,3,4,5,6], labels:'V'+['1','2','3','4','5','6']}



; vb1 in 4096Hz.
    psbl_de_trim_b1_data, utr, tprobe
    
    
    vars = pre0+['vsvy','vb1']
    options, vars, 'yrange', [-21,-17]
    foreach tvar, vars do begin
        get_data, tvar, uts, vsvy
        dr0 = sdatarate(uts)
        ; V -> V/m -> mV/m.
        eu = (vsvy[*,0]-vsvy[*,1])/boomlens[0]*1e3
        ev = (vsvy[*,2]-vsvy[*,3])/boomlens[1]*1e3
        ew = (vsvy[*,4]-vsvy[*,5])/boomlens[2]*1e3
        store_data, tvar+'_e_uvw', uts, [[eu],[ev],[ew]], limits = $
            {ytitle:'E UVW!C(mV/m)', labels:['Eu','Ev','Ew'], colors:[6,4,2], yrange:[-200,200]}
    endforeach
    

    
    
    
    ; convert all data to uniform time.
    tvar = pre0+'vb1_e_uvw'
    get_data, tvar, tuts, euvw
    
    dt0 = (utr[1]-utr[0])
    nrec = floor(dt0/b1dr)
    uts = utr[0]+findgen(nrec)*b1dr
    
    foreach tvar, datvars do begin
        if tvar eq pre0+'m_uvw2gse' then continue
        get_data, tvar, tuts, dat
        store_data, tvar, uts, sinterpol(dat, tuts, uts)
    endforeach
    
    
    ; **** rotate E uvw to gse.
    get_data, pre0+'vb1_e_uvw', uts, dat
    
    tmp = time_string(uts[0],tformat='YYYY-MM-DDThh:mm:ss.ffffff')
    cspice_str2et, tmp, tet0
    tets = tet0+uts-uts[0]
    
    scid = strupcase(pre0+'science')
    cspice_pxform, scid, 'GSE', tets, muvw2gse
    
    tvar = pre0+'m_uvw2gse'
    store_data, tvar, uts, muvw2gse

endif


; update data file.
if reload eq 1 then tplot_save, datvars, filename = datfn




; remove background for Euvw and VB1.
if tnames(pre0+'vb1_e_uvw_bg') eq '' then begin
    tvar = pre0+'vb1_e_uvw'
    get_data, tvar, uts, dat
    nrec = n_elements(uts)
    bg = dat
    for l = 0, 2 do begin
        bg[*,l] = smooth(dat[*,l],period0/b1dr,/edge_truncate)
        dat[*,l] = dat[*,l]-bg[*,l]
    endfor
    store_data, tvar, uts, dat
    store_data, tvar+'_bg', uts, bg, limits = {ytitle:'B1 E UVW Bg!C(mV/m)', $
        labels:'E'+['u','v','w'], colors:[6,4,2]}
    
    options, tvar, 'yrange', [-20,20]
    options, pre0+['vsvy_e_uvw','vb1_e_uvw_bg'], 'yrange', [-15,15]
    options, pre0+'vsvy_e_uvw', 'ytitle', 'Svy E UVW!C(mV/m)'
    options, pre0+'vb1_e_uvw', 'ytitle', 'B1 E UVW!C(mV/m)'
endif

tvar = pre0+'vb1'
get_data, tvar, uts, dat
nrec = n_elements(uts)
bg = dat
for l = 0, 5 do begin
    bg[*,l] = smooth(dat[*,l],period0/b1dr,/edge_truncate)
    dat[*,l] = dat[*,l]-bg[*,l]
endfor
store_data, tvar+'_bg', uts, bg, limits = {ytitle:'VB1 Bg!C(V)', $
    labels:'V'+['1','2','3','4','5','6'], colors:[1,2,3,4,5,6], yrange:[-21,-17]}

tdat = dblarr(nrec,6)
for l = 0, 5 do dat[*,l] = dat[*,l]/boomlens[l/2]*1e3*2
store_data, pre0+'eb1', uts, dat, limits = {ytitle:'EB1!C(mV/m)', $
    labels:'E'+['1','2','3','4','5','6'], colors:[1,2,3,4,5,6]}





; **** rotate Euvw into FAC.
get_data, pre0+'vb1_e_uvw', uts, euvw
get_data, pre0+'m_uvw2gse', uts, muvw2gse

egse = euvw
egse[*,0] = euvw[*,0]*muvw2gse[0,0,*] + euvw[*,1]*muvw2gse[1,0,*] + euvw[*,2]*muvw2gse[2,0,*]
egse[*,1] = euvw[*,0]*muvw2gse[0,1,*] + euvw[*,1]*muvw2gse[1,1,*] + euvw[*,2]*muvw2gse[2,1,*]
egse[*,2] = euvw[*,0]*muvw2gse[0,2,*] + euvw[*,1]*muvw2gse[1,2,*] + euvw[*,2]*muvw2gse[2,2,*]

store_data, pre0+'e_gse', uts, egse, limits = $
    {ytitle:'E GSE!C(mV/m)', colors:[6,4,2], labels:'GSE '+['x','y','z']}


get_data, pre0+'b_gse', tuts, bgse
bgse = sinterpol(bgse,tuts,uts,/quadratic)
get_data, pre0+'pos_gse', tuts, rgse
rgse = sinterpol(rgse,tuts,uts,/quadratic)

rhat = sunitvec(rgse)
bhat = sunitvec(bgse)
phat = sunitvec(scross(rhat,bhat))
vhat = scross(bhat,phat)

efac = [[sdot(egse,bhat)],[sdot(egse,phat)],[sdot(egse,vhat)]]
store_data, pre0+'e_fac', uts, efac, limits = $
    {ytitle:'E FAC!C(mV/m)', colors:[6,4,2], labels:'FAC '+['b','p','v']}





; **** cross correlation to determine the time shift b/w V1 and -V2, etc.
    
    drec0 = fix(period0*0.5/b1dr)    ; max rec shift.
    tnrec = 2*drec0+1       ; # time shift to be tested.
    corr = dblarr(tnrec)
    drecs = -drec0+findgen(tnrec)

    corrdts = dblarr(3)     ; in msec.

    get_data, pre0+'vb1', uts, vb1
    get_data, pre0+'vb1_bg', uts, vb1bg
    vb1 = vb1-vb1bg
    
    poss = sgcalcpos(3, ypad = 8)
    ofn = shomedir()+'/fig_cross_correlation.pdf'
; ofn = 0
    sgopen, ofn, xsize = 11, ysize = 8.5, /inch
    device, decompose = 0
    loadct2, 43
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
            
    for k = 0, 2 do begin
        v1 = vb1[*,k*2]
        v2 = -vb1[*,k*2+1]
        for l = 0, tnrec-1 do corr[l] = c_correlate(v1,v2,drecs[l])
        tmp = max(corr, idx)
        drec = drecs[idx]
        corrdts[k] = drec*b1dr*1e3
        
        
        tpos = poss[*,k]
        tpos[2] = tpos[0]+(tpos[2]-tpos[0])*0.6
        tdat = [[v1],[shift(v2,-drec)]]
        tdat[0:abs(drec),*] = !values.d_nan
        tdat[nrec-1-abs(drec):nrec-1,*] = !values.d_nan
        store_data, tvar+'_tmp', uts, tdat, limits = $
            {labels:'V'+string(k*2+[1,2],format='(I0)'),colors:k*2+[1,2],ytitle:'V'+sgnum2str(k*2+1)+sgnum2str(k*2+2)+' shfted!C(V)'}
        tplot, tvar+'_tmp', trange = tutr, position = tpos, /noerase, /novtitle
        store_data, tvar+'_tmp', /delete
        
        tpos = poss[*,k]
        tpos[0] = tpos[0]+(tpos[2]-tpos[0])*0.75
        plot, drecs*b1dr*1e3, corr, /noerase, position = tpos, $
            yrange = [0,1], ystyle = 1, $
            xtitle = 'Time shift (msec)', ytitle = 'Correlation', $
            title = 'dT = '+sgnum2str(corrdts[k],ndec=3)+' msec'
        plots, corrdts[k]+[0,0], !y.crange, linestyle = 1
    endfor
    sgclose
        
    k0uvw = omega0*corrdts*1e-3/boomlens ; in /m.
    k0mag = snorm(k0uvw)                  ; in /m.
    vphuvw = boomlens/corrdts  ; in km/s.
    vphmag = snorm(vphuvw)     ; in km/s.
    
    ; to gse.    
    muvw2gse = reform(muvw2gse[*,*,nrec/2])
    k0gse = reform(k0uvw # muvw2gse)
    vphgse = reform(vphuvw # muvw2gse)
    
    ; to fac.
    k0fac = [sdot(k0gse,bhat[nrec/2,*]),sdot(k0gse,phat[nrec/2,*]),sdot(k0gse,vhat[nrec/2,*])]
    vphfac = [sdot(vphgse,bhat[nrec/2,*]),sdot(vphgse,phat[nrec/2,*]),sdot(vphgse,vhat[nrec/2,*])]

stop







; **** frequency is 50 Hz, 20 msec period, |E| 15-20 mV/m.
; plot the FFT power spectrum for VB1 and EB1.

freqs = findgen(nrec)/(nrec*b1dr)

ofn = 0
ofn = shomedir()+'/fig_fft.pdf'
sgopen, ofn, xsize = 11, ysize = 8.5, /inch
device, decomposed = 0
loadct2, 43

poss = sgcalcpos(3,3, ypad = 5, xpad = 5, region = [0,0,1,0.9])
xr = [10,200]

get_data, pre0+'vb1', uts, vb1
get_data, pre0+'vb1_bg', uts, vb1bg
vb1 = vb1-vb1bg
for i = 0, 2 do begin
    for j = 0, 1 do begin
        tdat = vb1[*,i*2+j]
        fftpow = (abs(fft(tdat)))^2*2*!dpi*(utr[1]-utr[0])
        plot, freqs[0:nrec/2], fftpow[0:nrec/2], /noerase, position = poss[*,i,j], $
            /xlog, xrange = xr, xstyle = 1, $
            ytitle = '(V!U2!N/Hz)', $
            xtitle = '', title = 'V'+sgnum2str(i*2+j+1)
    endfor
endfor

get_data, pre0+'vb1_e_uvw', uts, euvw
uvwlabs = 'E'+['u','v','w']
for i = 0, 2 do begin
    tdat = euvw[*,i]
    fftpow = (abs(fft(tdat)))^2*2*!dpi*(utr[1]-utr[0])
    plot, freqs[0:nrec/2], fftpow[0:nrec/2], /noerase, position = poss[*,i,2], $
        /xlog, xrange = xr, xstyle = 1, $
        ytitle = '(mV!U2!N/m!U2!N-Hz)', $
        xtitle = 'Frequency (Hz)', title = uvwlabs[i]
endfor

xyouts, 0.5, 0.9, /normal, alignment = 0.5, charsize = 1.25, $
    'FFT power spectrum of Burst V[123456] and E[uvw]'
    
sgclose

stop



end
