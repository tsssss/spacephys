;+
; download rbsp 'raw' data, save plot to disk.
;-
pro psbl_de_gen_vsc_survey_plot, id, trange = utr, rootdir = rootdir

    if n_elements(id) eq 0 then message, 'no id ...'
    tprobe = strmid(id,strlen(id)-1)
    if n_elements(utr) eq 0 then message, 'no time range ...'
    if n_elements(rootdir) eq 0 then rootdir = shomedir()+'/psbl_de_32hz'
    

    plotvar = []
    lablvar = []

    
    ; expand utr to at least 10 min.
    theutr = utr
    dt = 15*60    ; sec
    utr = mean(utr) & utr = utr-(utr mod 60) & utr+= [-1,1]*dt
    
    
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
    
    
    
    
    
    ; **** load Vsc.
    ; rbspx_[vsc].
    
    tvar = pre0+'vsc'
    dat = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'vsvy')
    if size(dat,/type) ne 8 then begin
        message, 'no V data ...', /continue
        return
    endif else begin
        tmp = sfmepoch(dat.epoch, 'unix')
        dr = sdatarate(tmp[sort(tmp)])
        uts = smkarthm(min(tmp), max(tmp), dr, 'dx')    ; get uniform time.
        dat = sinterpol(dat.vsvy, tmp, uts)
        
        cp0 = rbsp_efw_get_cal_params(utr[0])
        case strlowcase(tprobe) of
            'a': cp = cp0.a
            'b': cp = cp0.b
            else: dprint, 'Invalid probe name. Calibration aborted.'
        endcase
        gain = cp.adc_gain_vdc
        offset = cp.adc_offset_vdc
        for i = 0, 5 do dat[*,i] =  (dat[*,i]-offset[i])*gain[i]
        vsc = mean(dat[*,0:3], dimension = 2)
        store_data, tvar, uts, vsc, limits = {ytitle:'Vsc!C(V)', labels:'Vsc'}
        plot, vsc, /nodata & yr = !y.crange & dyr = abs(yr[1]-yr[0])
        if dyr le 0.1 then yr = mean(yr)+[-1,1]*dyr ; 0.2 V range.
        options, tvar, 'yrange', yr
    endelse
    plotvar = [plotvar,tvar]


; **** load E field. 32Hz dE MGSE.
; rbspx_de_mgse, rbspx_de_mgse_survey.
    
    tvar = pre0+'de_mgse'
    dat = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'esvy')
    if size(dat,/type) ne 8 then begin
        message, 'no E data ...', /continue
        return
    endif else begin
        tmp = sfmepoch(dat.epoch, 'unix')
        dr = sdatarate(tmp[sort(tmp)])
        uts = smkarthm(min(tmp), max(tmp), dr, 'dx')    ; get uniform time.
        dat = sinterpol(dat.efield_mgse, tmp, uts)
        dat[*,0] = 0.   ; set x to 0.
        store_data, tvar, uts, dat[*,[1,2]], limits = $
            {ytitle:'dE MGSE!C(mV/m)', colors:rgb[[1,2]], $
            labels:'MGSE '+['y','z'], constant:[-25,0,25,-100,100], $
            yrange:[-200,200], ystyle:1, yticks:4, yminor:5}
    endelse
    plotvar = [plotvar,tvar]
    
    
    tvar = pre0+'de_mgse_survey'
    dat = sread_rbsp_efw_l3(utr, probes = tprobe)
    if size(dat,/type) ne 8 then begin
        message, 'no E survey data ...', /continue
        store_data, tvar
    endif else begin
        tmp = sfmepoch(dat.epoch, 'unix')
        dr = sdatarate(tmp[sort(tmp)])
        uts = smkarthm(min(tmp), max(tmp), dr, 'dx')    ; get uniform time.
        tdat = sinterpol(dat.efield_inertial_frame_mgse, tmp, uts)
        tdat[*,0] = 0.   ; set x to 0.
        store_data, tvar, uts, tdat[*,[1,2]], limits = $
            {ytitle:'dE MGSE!C(mV/m)', colors:rgb[[1,2]], $
            labels:'MGSE '+['y','z'], constant:[-10,0,10], $
            yrange:[-50,50], ystyle:1, yticks:4, yminor:5}
    endelse
    plotvar = [plotvar,tvar]

    
    
; **** load position.
; rbspx_[pos_gse,dis,mlt,lshell].

    if size(dat,/type) ne 8 then begin
        rbsp_load_spice_kernels, trange = utr, probes = tprobe
        tvar = pre0+'pos_gse'
        rbsp_load_spice_state, probe = tprobe, coord = 'gse', times = uts, /no_spice_load
        get_data, pre0+'state_pos_gse', tmp, posgse
        posgse = posgse*re1
        tdat = posgse
        store_data, pre0+'state_*', /delete
        rbsp_load_spice_kernels, trange = utr, probes = tprobe, /unload
        
        store_data, tvar, uts, tdat, $
            {ytitle:'R GSE (Re)', colors:rgb, labels:'GSE '+['x','y','z']}
        store_data, tvar+'x', uts, tdat[*,0], limits = {ytitle:'X GSE (Re)'}
        store_data, tvar+'y', uts, tdat[*,1], limits = {ytitle:'Y GSE (Re)'}
        store_data, tvar+'z', uts, tdat[*,2], limits = {ytitle:'Z GSE (Re)'}
        lablvar = [lablvar,tvar+['x','y','z']]


        tvar = pre0+'mlt'
        ets = stoepoch(uts,'unix')
        tdat = atan(posgse[*,1],posgse[*,0])*deg
        tdat = (tdat+360) mod 360   ; convert to 0-360.
        tdat = (tdat/15 + 12) mod 24
        store_data, tvar, uts, tdat, limits = {ytitle:'MLT (hr)'}
        lablvar = [lablvar,tvar]
        
        
        tvar = pre0+'lshell'
        possm = sgse2sm(posgse, ets)
        mlat = atan(possm[*,2],sqrt(possm[*,0]^2+possm[*,1]^2)) ; in rad.
        dis = sqrt(possm[*,0]^2+possm[*,1]^2+possm[*,2]^2)
        tdat = dis/(cos(mlat)^2)
        store_data, tvar, uts, tdat, limits = {ytitle:'L-shell'}
        lablvar = [lablvar,tvar]        
    endif else begin
        tdat = sinterpol(dat.pos_gse, tmp, uts)*re1
        tvar = pre0+'pos_gse'
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
    endelse



; **** save data and plot to disk.
    ofn = rootdir+'/psbl_de_32hz_vsc_survey_'+id+'.pdf'
;    ofn = 0
    
    titl = id+': Vsc, dE MGSE 32 Hz and 0.1 Hz'
    sgopen, ofn, xsize = 8.5, ysize = 11, /inch
    device, decomposed = 0
    loadct2, 43
    tplot, plotvar, trange = utr, var_label = lablvar, title = titl
    timebar, theutr, color = 6
    sgclose
    
end

idinfos = psbl_de_id('all_e')
nidinfo = n_elements(idinfos)

;idx = where(idinfos.id eq '2015_0131_2153_a')
;idinfos = idinfos[idx[0]:*]
;nidinfo = n_elements(idinfos)


foreach tinfo, idinfos do psbl_de_gen_vsc_survey_plot, tinfo.id, trange = tinfo.utr

end
