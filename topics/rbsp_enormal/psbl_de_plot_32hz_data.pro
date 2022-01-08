;+
; download rbsp 'raw' data, save plot to disk.
;-
pro psbl_de_plot_32hz_data, id, trange = utr, rootdir = rootdir

    if n_elements(id) eq 0 then message, 'no id ...'
    if n_elements(tprobe) eq 0 then tprobe = strmid(id,strlen(id)-1)
    if n_elements(utr) eq 0 then $
        utr = time_double(strmid(id,0,14),tformat='YYYY_MMDD_hhmm')+[-2,8]*60
    if n_elements(rootdir) eq 0 then rootdir = shomedir()+'/psbl_de_32hz'

    plotvar = []
    lablvar = []
    
    ; expand utr to at least 10 min.
    theutr = utr
    utr = mean(utr) & utr = utr-(utr mod 60) & utr+= [-1,1]*0.5*600

    
    pre0 = 'rbsp'+tprobe+'_'
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    minratio = 0.1  ; Bx/B ratio in MGSE>
    rgb = [6,4,2]
    dt0 = 86400d
    utr0 = utr-(utr mod dt0)+[0,dt0]
    timespan, utr0[0], dt0, /second
    
    tplot_options, 'constant', 0
    tplot_options, 'num_lab_min', 10
    tplot_options, 'labflag', -1
    tplot_options, 'ynozero', 1
    tplot_options, 'ymargin', [5,5]
    tplot_options, 'xmargin', [25,15]
    tplot_options, 'no_interp', 1
    
    
; **** load E field. 32Hz dE MGSE.
; rbspx_de_mgse, rbspx_de_mgse_survey.

    tvar = pre0+'de_mgse'
    dat = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'esvy')
    if size(dat,/type) ne 8 then message, 'no E data ...'
    tmp = sfmepoch(dat.epoch, 'unix')
    dr = sdatarate(tmp[sort(tmp)])
    uts = smkarthm(min(tmp), max(tmp), dr, 'dx')    ; get uniform time.
    dat = sinterpol(dat.efield_mgse, tmp, uts)
    dat[*,0] = 0.   ; set x to 0.
    store_data, tvar, uts, dat[*,[1,2]], limits = $
        {ytitle:'dE MGSE!C(mV/m)', colors:rgb[[1,2]], $
        labels:'MGSE '+['y','z'], constant:[-25,0,25]}    
    plotvar = [plotvar,tvar]
    
    
    tvar = pre0+'de_mgse_survey'
    dat = sread_rbsp_efw_l3(utr, probes = tprobe)
    if size(dat,/type) ne 8 then message, 'no E data ...'
    tmp = sfmepoch(dat.epoch, 'unix')
    dr = sdatarate(tmp[sort(tmp)])
    uts = smkarthm(min(tmp), max(tmp), dr, 'dx')    ; get uniform time.
    tdat = sinterpol(dat.efield_inertial_frame_mgse, tmp, uts)
    tdat[*,0] = 0.   ; set x to 0.
    store_data, tvar, uts, tdat[*,[1,2]], limits = $
        {ytitle:'dE MGSE!C(mV/m)', colors:rgb[[1,2]], $
        labels:'MGSE '+['y','z'], constant:[-25,0,25]}    
    plotvar = [plotvar,tvar]
    

; **** load position.
; rbspx_[pos_gsm,dis,mlt,lshell].

    tvar = pre0+'pos_gse'
    tdat = sinterpol(dat.pos_gse, tmp, uts)*re1
    store_data, tvar, uts, tdat, $
        {ytitle:'R GSE (Re)', colors:rgb, labels:'GSE '+['x','y','z']}
    store_data, tvar+'x', uts, tdat[*,0], limits = {ytitle:'X GSE (Re)'}
    store_data, tvar+'y', uts, tdat[*,1], limits = {ytitle:'Y GSE (Re)'}
    store_data, tvar+'z', uts, tdat[*,2], limits = {ytitle:'Z GSE (Re)'}
    lablvar = [lablvar,tvar+['x','y','z']]
    
    tvar = pre0+'mlt'
    tdat = sinterpol(dat.mlt_lshell_mlat, tmp, uts)
    store_data, tvar, uts, tdat[*,0], limits = {ytitle:'MLT (hr)'}
    lablvar = [lablvar,tvar]
    
    tvar = pre0+'lshell'
    store_data, tvar, uts, tdat[*,1], limits = {ytitle:'L-shell'}
    lablvar = [lablvar,tvar]



; **** load B field. l3, gse, emfisis, hires 171 MB for 1 day.
; rbspx_b_gse, rbspx_b.

    emfisis = sread_rbsp_emfisis_l3(utr, probes = tprobe, type = 'hires')
    if size(emfisis,/type) ne 8 then message, 'no B data ...'
    uts = sfmepoch(emfisis.epoch,'unix')
    store_data, pre0+'b_gse', uts, emfisis.mag, $
        limits = {ytitle:'B!C(nT)',colors:rgb, labels:'GSE '+['x','y','z']}
    store_data, pre0+'b', uts, sqrt(total(emfisis.mag^2,2)), $
        limits = {ytitle:'B magnitude!C(nT)'}    
    plotvar = [plotvar,pre0+'b_gse']

    
; **** load HOPE spectrograms.
; rbspx_[e,p]_enspec.

    dat = sread_rbsp_hope_l3(utr0, probes = tprobe, type = 'pa')
    if size(dat,/type) eq 8 then begin
        npa = n_elements(dat.pitch_angle)
        uts = sfmepoch(dat.epoch_ele, 'unix')
        store_data, pre0+'e_enspec', uts, reform(total(dat.fedu,3,/nan)), $
            dat.hope_energy_ele, limits = {ytitle:'Electron!CEnergy (eV)',spec:1, $
            no_interpol:1, ylog:1, zlog:1, yrange:[10,4e4], zrange:[1e7,1e9], $
            ztitle:'(s!E-1!Ncm!E-2!Nster!E-1!NkeV!E-1!N)', ystyle:1}
        uts = sfmepoch(dat.epoch_ion, 'unix')
        store_data, pre0+'p_enspec', uts, reform(total(dat.fpdu,3,/nan)), $
            dat.hope_energy_ion, limits = {ytitle:'Proton!CEnergy (eV)',spec:1, $
            no_interpol:1, ylog:1, zlog:1, yrange:[10,4e4], zrange:[1e5,1e7], $
            ztitle:'(s!E-1!Ncm!E-2!Nster!E-1!NkeV!E-1!N)', ystyle:1}
    endif else begin
        store_data, pre0+'e_enspec', 0, 0
        store_data, pre0+'p_enspec', 0, 0
    endelse    
    plotvar = [plotvar,pre0+['e_enspec','p_enspec']]
    
    ; use the z-value to set the zrange.
    tvar = pre0+'e_enspec'
    get_data, tvar, uts, dat, val
    idx = where(uts ge utr[0] and uts le utr[1], cnt)
    if cnt gt 0 then begin
        dat = dat[idx,*,*]
        tmp = where(dat eq 0, cnt)
        if cnt ne 0 then dat[tmp] = !values.d_nan
        val = reform(val[idx[0],*])
        minen = 100     ; eV.
        tmp = min(val-minen, tidx, /absolute)
        vmax = median(dat[*,tidx])
        if finite(vmax) then vmax = 10d^(ceil(alog10(vmax))) else begin
            vmax = max(dat,/nan)
            if finite(vmax) then vmax = 10d^(floor(alog10(vmax))) else vmax = 1e8
        endelse
        maxen = 10000   ; eV.
        tmp = min(val-maxen, tidx, /absolute)
        vmin = median(dat[*,tidx])
        if finite(vmin) then vmin = 10d^(floor(alog10(vmin))) else begin
            vmin = min(dat,/nan)
            if finite(vmin) then vmin = 10d^(ceil(alog10(vmin))) else vmin = 1e4
        endelse
        options, tvar, 'zrange', [vmin,vmax]
;        scl = alog10(vmax/vmin)
;        if scl le 2 or scl le 4 then options, tvar, 'zrange', [1e4,1e8]
        
        tvar = pre0+'p_enspec'
        get_data, tvar, uts, dat, val
        dat = dat[idx,*,*]
        tmp = where(dat eq 0, cnt)
        if cnt ne 0 then dat[tmp] = !values.d_nan
        vmax = max(dat, /nan)
        if finite(vmax) then vmax = 10d^(floor(alog10(vmax))) else vmax = 1e7
        vmin = min(dat, /nan)
        if finite(vmin) then vmin = 10d^(ceil(alog10(vmin))) else vmin = 1e5
        options, tvar, 'zrange', [vmin,vmax]
;        scl = alog10(vmax/vmin)
;        if scl le 2 or scl le 4 then options, tvar, 'zrange', [1e5,1e7]
    endif
    


; **** load HOPE moments.
; rbspx_[n,t]
    hopemom = sread_rbsp_hope_l3(utr0, probes = tprobe, type = 'mom')
    if size(hopemom,/type) eq 8 then begin
        uts = sfmepoch(hopemom.epoch_ele,'unix')
        store_data, pre0+'n', uts, hopemom.dens_e_200, $
            limits = {ytitle:'N!Ie!N!C(cm!E-3!N)', ylog:1, constant:[0.1,1]}

        store_data, pre0+'t', uts, [[hopemom.tperp_e_200],[hopemom.tpar_e_200]], $
            limits = {ytitle:'T!Ie!I!N!C(eV)', ylog:1, colors:[6,0], $
            labels:['Tperp','Tpara'], constant:1000}
        tvar = pre0+'t'
        get_data, tvar, uts, dat
        idx = where(dat eq 1e20, cnt)
        if cnt ne 0 then dat[idx] = !values.d_nan
        store_data, tvar, uts, dat
    endif else begin
        store_data, pre0+'n', 0, 0
        store_data, pre0+'t', 0, 0
    endelse        
    plotvar = [plotvar,pre0+['n','t']]
    
    ; use the dynamic range to set yrange.
    get_data, pre0+'n', uts
    idx = where(uts ge utr[0] and uts le utr[1], cnt)
    if cnt ne 0 then begin
        foreach tvar, pre0+['n','t'] do begin
            get_data, tvar, uts, dat
            dat = dat[idx]
            uts = uts[idx]
            tmin = min(dat,/nan)
            tmax = max(dat,/nan)
            scl = 10^ceil(alog10(tmax/tmin))
            ymax = double(ceil(sqrt(tmin*tmax*scl)))
            ymin = ymax/scl
            options, tvar, 'yrange', [ymin,ymax]
            options, tvar, 'ystyle', 1
        endforeach
    endif


; **** omni Dst and AE.
; [dst,ae].
    omni = sread_omni(utr)
    uts = sfmepoch(omni.epoch,'unix')
    store_data, 'dst', uts, omni.sym_h, limits = {ytitle:'Dst (nT)'}
    store_data, 'ae', uts, omni.ae_index, limits = {ytitle:'AE (nT)'}    
    lablvar = [lablvar,['dst','ae']]


; **** beta.
; rbspx_beta.
    ; plasma beta. 2mu0*nkT/B^2
    tmp = 2*4*!dpi*1e-7*1e6*1.6e-19*1e18
    get_data, pre0+'n', t0, dat & dat*= tmp
    get_data, pre0+'t', t0, tmp & dat*= tmp[*,0]
    get_data, pre0+'b', t1, tmp & tmp = interpol(tmp, t1, t0) & dat/= tmp^2
    tvar = pre0+'beta'
    store_data, tvar, t0, dat, limits = {ytitle:'Beta', $
        ylog:1, constant:0.01}
    plotvar = [plotvar,tvar]

    
    
; **** save data and plot to disk.
    ofn = rootdir+'/psbl_de_32hz_survey_plot_'+id+'.pdf'
;    ofn = 0
    sgopen, ofn, xsize = 8.5, ysize = 11, /inch
    device, decomposed = 0
    loadct2, 43
    tplot, plotvar, trange = utr, var_label = lablvar
    timebar, theutr, color = 6
    sgclose

end


rootdir = shomedir()+'/psbl_de_32hz/'

runlog = rootdir+'/run_process.log'
dir = file_dirname(runlog)
if file_test(dir) eq 0 then file_mkdir, dir
if file_test(runlog) eq 0 then stouch, runlog

idinfos = psbl_de_id('detect_step')
nidinfo = n_elements(idinfos)

;rbsp_load_spice_kernels, probe = 'a', /all
;rbsp_load_spice_kernels, probe = 'b', /all


for i = 0, nidinfo-1 do begin
    
    id = idinfos[i].id
    utr = idinfos[i].utr
    
    ; write process in log file.
    openw, lun, runlog, /get_lun, /append
    printf, lun, 'processing '+id+' ...'
    printf, lun, 'time range: '+time_string(utr[0])+' to '+time_string(utr[1])
    printf, lun, ''
    free_lun, lun
    
    store_data, '*', /delete
    psbl_de_plot_32hz_data, id, trange = utr, rootdir = rootdir
endfor

end
