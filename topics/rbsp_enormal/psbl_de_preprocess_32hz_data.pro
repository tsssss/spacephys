;+
; read the tplot vars created by psbl_de_save_32hz_data.
; check E direction, and register the event.
; calc Poynting flux.
; write event info array.
;-

pro psbl_de_preprocess_32hz_data, id, trange = utr, rootdir = rootdir

    if n_elements(id) eq 0 then message, 'no id ...'
    if n_elements(utr) eq 0 then message, 'no time range ...'
    
    tprobe = strmid(id,strlen(id)-1)
    if n_elements(rootdir) eq 0 then rootdir = shomedir()+'/psbl_de_32hz'

    ; load the tplot vars.
    datadir = shomedir()+'/Google Drive/works/data/rbsp_de/32hz_data'
    datafn = datadir+'/psbl_de_32hz_data_'+id+'.tplot'
    tplot_restore, filename = datafn
    
    
    datavar = []
    
testplot = 0

    figdir = rootdir+'/backup_plot/'+id
    datdir = rootdir+'/data'
    if file_test(figdir,/directory) eq 0 then file_mkdir, figdir
    
    pre0 = 'rbsp'+tprobe+'_'    
    rgb = [6,4,2]
    deg = 180/!dpi
    theta = '!9'+string(113b)+'!X'
    sdeg = '!9'+string(176b)+'!X'
    
    tplot_options, 'constant', 0
    tplot_options, 'num_lab_min', 10
    tplot_options, 'labflag', -1
    tplot_options, 'ynozero', 1
    tplot_options, 'ymargin', [5,5]
    tplot_options, 'xmargin', [25,15]
    tplot_options, 'no_interp', 1
    
    
    
; **** separate B0 and dB.
; rbspx_[b0,db]_gse
    get_data, pre0+'b_gse', uts, bgse
    b0gse = bgse
    for i = 0, 2 do b0gse[*,i] = scalcbg(bgse[*,i])
    store_data, pre0+'b0_gse', uts, b0gse, limits = $
        {ytitle:'B0!C(nT)', colors:rgb, labels:'GSE '+['x','y','z']}
    store_data, pre0+'db_gse', uts, bgse-b0gse, limits = $
        {ytitle:'dB!C(nT)', colors:rgb, labels:'GSE '+['x','y','z']}
    vars = pre0+['b0_gse','db_gse']
    datavar = [datavar,vars]


; **** reconstruct 3-D field. use B0.
; rbspx_b0_angle,rbspx_[de,de_survey,dedot0,dedot0_survey]_mgse.
    bang0 = 20  ; deg.
    bvar = pre0+'b0_'
    get_data, bvar+'gse', uts, dat
    get_data, pre0+'w_gse', tmp, wgse
    wgse = sinterpol(wgse, tmp, uts)
    dat = sgse2mgse(dat, wsc = wgse)
    store_data, bvar+'mgse', uts, dat, limits = $
        {ytitle:'B0!C(nT)', colors:rgb, labels:'MGSE '+['x','y','z']}
    bang = atan(dat[*,0],sqrt(dat[*,1]^2+dat[*,2]^2))*deg
    store_data, pre0+'b0_angle', uts, bang, limits = $
        {ytitle:'Angle!C(deg)', labels:'<B!Dw!N,|B!Du,v!N|>', $
        constant:bang0*[-1,1]}
    datavar = [datavar,pre0+'b0_angle']
    

    vars = pre0+['de_mgse','de_survey_mgse']
    foreach tvar, vars do begin
        get_data, tvar, uts, dat
        nrec = n_elements(uts)
        dat = [[dblarr(nrec)], [dat]]
        store_data, tvar, uts, dat, limits = $
            {ytitle:'dE!C(mV/m)', colors:rgb, labels:'MGSE '+['x','y','z']}
            
        get_data, bvar+'mgse', tmp, bmgse
        bmgse = sinterpol(bmgse, tmp, uts)
        dat[*,0] = -(dat[*,1]*bmgse[*,1]+dat[*,2]*bmgse[*,2])/bmgse[*,0]
        bang = atan(bmgse[*,0],sqrt(bmgse[*,1]^2+bmgse[*,2]^2))*deg
        idx = strpos(tvar, 'de')
        tvar = strmid(tvar,0,idx)+'dedot0'+strmid(tvar,idx+2)
        store_data, tvar, uts, dat, limits = $
            {ytitle:'dE!C(mV/m)', colors:rgb, labels:'MGSE '+['x','y','z']}
    endforeach
    datavar = [datavar,pre0+['dedot0_mgse','dedot0_survey_mgse']]
    options, pre0+'de_mgse', 'ytitle', 'dE!C(mV/m)'
    options, pre0+'de_survey_mgse', 'ytitle', 'dE survey!C(mV/m)'
    options, pre0+'dedot0_mgse', 'ytitle', 'dE dot0!C(mV/m)'
    options, pre0+'dedot0_survey_mgse', 'ytitle', 'dE dot0 survey!C(mV/m)'
    

; **** use wgse to rotate E mgse to gse.
    vars = pre0+['de','de_survey','dedot0','dedot0_survey']
    foreach tvar, vars do begin
        get_data, tvar+'_mgse', uts, dat, limits = lims
        get_data, pre0+'w_gse', tmp, wgse
        wgse = sinterpol(wgse, tmp, uts)
        dat = smgse2gse(dat, wsc = wgse)
        store_data, tvar+'_gse', uts, dat, limits = $
            {ytitle:lims.ytitle, colors:rgb, labels:'GSE '+['x','y','z']}
    endforeach

    get_data, pre0+'w_gse', uts, wgse
    dat = sgse2mgse(wgse, wsc = wgse)
    store_data, pre0+'w_mgse', uts, dat, limits = $
        {ytitle:'Spin Axis!Cdirection', colors:rgb, labels:'MGSE '+['x','y','z']}
    
    
    ; ** plot: shows the B0 and dB decomposition, rotation from MGSE to GSE.
    posl = [0.10,0.1,0.3,0.9]
    posm = [0.43,0.1,0.63,0.9]
    posr = [0.75,0.1,0.95,0.9]
    
    
    ofn = figdir+'/'+id+'_bsep_erot.pdf'
    if testplot then ofn = 0
    sgopen, ofn, xsize = 11.5, ysize = 8, /inch
    device, decomposed = 0
    loadct2, 43
    tutr = utr+[-1,1]*60
    
    get_data, pre0+'b_gse', uts, dat
    bmag = snorm(dat)
    get_data, pre0+'b0_gse', uts, dat
    bmag = [[bmag],[snorm(dat)]]
    store_data, pre0+'bmag', uts, bmag, limits = $
        {ytitle:'|B|!C(nT)', labels:['B!DTotal','B!DBackground'], colors:[0,2]}
    
    vars = pre0+['b_gse','db_gse','b0_gse','bmag','b0_angle']
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, position = posl)
    tplot, vars, trange = tutr, position = poss, /noerase, title = 'B Decomposition', vlab_margin = 12
    timebar, utr, linestyle = 1
    
    vars = pre0+['dedot0','dedot0_survey','de','de_survey','w']
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, position = posm)
    tplot, vars+'_mgse', trange = tutr, position = poss, /noerase, title = 'E&B fields in MGSE', /novtitle
    timebar, utr, linestyle = 1
    
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, position = posr)
    tplot, vars+'_gse', trange = tutr, position = poss, /noerase, title = 'E&B fields in GSE', /novtitle
    timebar, utr, linestyle = 1
    
    xyouts, 0.5, 0.95, /normal, alignment = 0.5, charsize = 1.25, $
        'B decomposition and rotation between MGSE and GSE'
    
    sgclose
    


; **** rotate E/B fields into FAC coord.
; no need to save the rotation angles, They can be calculated from b0gse.
    vars = pre0+['de_gse','de_survey_gse','dedot0_gse','dedot0_survey_gse',$
        'db_gse','b0_gse']
    foreach tvar, vars do begin
        get_data, tvar, tuts, dat, limits = lims

        get_data, pre0+'b0_gse', uts, b0gse
        b0gse = sinterpol(b0gse, uts, tuts)
        bhat = sunitvec(b0gse)
        get_data, pre0+'pos_gse', uts, r0gse
        r0gse = sinterpol(r0gse, uts, tuts)
        
        rhat = sunitvec(r0gse)
        phat = scross(rhat,bhat)
        vhat = scross(bhat,phat)        
        dat = [[sdot(dat,bhat)],[sdot(dat,phat)],[sdot(dat,vhat)]]
        
        
;        p = atan(bhat[*,1],bhat[*,0])   ; angle in GSE x-y plane.
;        t = atan(bhat[*,2],sqrt(bhat[*,1]^2+bhat[*,0]^2))   ; angle out of GSE x-y plane.
;        srotate, dat,-p, 2
;        srotate, dat, t, 1

        tmp = strmid(tvar,0,strpos(tvar,'gse'))+'fac'
        store_data, tmp, tuts, dat, limits = $
            {colors:rgb, labels:'FAC '+['b','p','v'], ytitle:lims.ytitle, constant:0}
    endforeach
    vars = pre0+['de_fac','de_survey_fac','dedot0_fac','dedot0_survey_fac','db_fac']
    datavar = [datavar,vars]



    ; ** plot: shows the rotation from GSE to FAC.
    posl = [0.1,0.1,0.4,0.9]
    posr = [0.6,0.1,0.9,0.9]
    
    ofn = figdir+'/'+id+'_gse_to_fac.pdf'
    if testplot then ofn = 0
    sgopen, ofn, xsize = 11.5, ysize = 8, /inch
    device, decomposed = 0
    loadct2, 43
    tutr = utr+[-1,1]*60
    
    vars = pre0+['dedot0','dedot0_survey','de','de_survey','db','b0']
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, position = posl)
    tplot, vars+'_gse', trange = tutr, position = poss, /noerase, title = 'E&B fields in GSE', vlab_margin = 12
    timebar, utr, linestyle = 1

    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, position = posr)
    tplot, vars+'_fac', trange = tutr, position = poss, /noerase, title = 'E&B fields in FAC', vlab_margin = 12
    xyouts, 0.5, 0.95, /normal, alignment = 0.5, 'Rotation from GSE to FAC', charsize = 1.25
    timebar, utr, linestyle = 1
    sgclose



; **** use e_survey to determine E direction, bipolar or unipolar.

    @tplot_com.pro
    
    vars = pre0+['de','dedot0']
    nvar = n_elements(vars)
    
    ; ** plot: shows the polarization.
    ofn = figdir+'/'+id+'_e_polarization.pdf'
    if testplot then ofn = 0
    sgopen, ofn, xsize = 6, ysize = 8, /inch
    device, decomposed = 0
    loadct2, 43

    poss = sgcalcpos(nvar, ypad = 8, tmargin = 8)
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    for i = 0, nvar-1 do begin
        tvar = vars[i]+'_survey_fac'
        get_data, tvar, uts, dat
        tutr = utr+[-1,1]*60
        idx = where(uts gt tutr[0] and uts le tutr[1])
        uts = uts[idx]
        dat = dat[idx,*]
        maxde = max(dat[*,2], tmp) & maxdeut = uts[tmp]
        minde = min(dat[*,2], tmp) & mindeut = uts[tmp]
        
        polratio = abs(maxde/minde)
        if polratio lt 1 then polratio = 1d/polratio
        
        
        color = 2
        tpos = poss[*,i]
        
        tplot, tvar, trange = tutr, position = tpos, vlab_margin = 12, /noerase, title = vars[i]
        timebar, utr, linestyle = 1

        
        tx = maxdeut
        ty = maxde
        xr = tutr
        yr = tplot_vars.settings.y.crange
        tx = (tx-xr[0])/(xr[1]-xr[0])
        tx = tpos[0]+(tpos[2]-tpos[0])*tx
        ty = (ty-yr[0])/(yr[1]-yr[0])
        ty = tpos[1]+(tpos[3]-tpos[1])*ty
        plots, tx, ty, /normal, psym = 2, symsize = 0.5, color = color
        xyouts, tx+xchsz, ty-ychsz, /normal, 'max Ez '+sgnum2str(maxde,nsgn=3)+' mV/m', color = color
        
        tx = mindeut
        ty = minde
        xr = tutr
        yr = tplot_vars.settings.y.crange
        tx = (tx-xr[0])/(xr[1]-xr[0])
        tx = tpos[0]+(tpos[2]-tpos[0])*tx
        ty = (ty-yr[0])/(yr[1]-yr[0])
        ty = tpos[1]+(tpos[3]-tpos[1])*ty
        plots, tx, ty, /normal, psym = 2, symsize = 0.5, color = color
        xyouts, tx+xchsz, ty+0.5*ychsz, /normal, 'min Ez '+sgnum2str(minde,nsgn=3)+' mV/m', color = color
        
        tdat = snorm(dat)
        maxde = max(tdat, tmp) & maxdeut = uts[tmp]
        polang = atan(dat[tmp,1],dat[tmp,2])*deg
        timebar, maxdeut
        tx = maxdeut
        xr = tutr
        tx = (tx-xr[0])/(xr[1]-xr[0])
        tx = tpos[0]+(tpos[2]-tpos[0])*tx
        ty = (tpos[1]+tpos[3])*0.5
        xyouts, tx+xchsz, ty+2*ychsz, /normal, 'max |E| '+ $
            sgnum2str(maxde,nsgn=3)+' mV/m!C<E,normal>!N '+sgnum2str(polang,nsgn=3)+sdeg
    endfor
    
    xyouts, 0.5, 0.95, /normal, alignment = 0.5, $
        'E FAC polarization and asymmetry', charsize = 1.25
        
    sgclose
        
    
; **** save derived data.
    datfn = datdir+'/psbl_de_32hz_deriv_data_'+id+'.tplot'
    tplot_save, datavar, filename = datfn

end

rootdir = shomedir()+'/psbl_de_32hz/'

idinfos = psbl_de_id('detect_step')
nidinfo = n_elements(idinfos)


for i = 0, nidinfo-1 do begin
    id = idinfos[i].id
    utr = idinfos[i].utr
;    if id ne '2013_0501_0736_b' then continue
;    if id ne '2013_0126_2125_a' then continue
;    if id ne '2013_1002_0459_b' then continue
    
    store_data, '*', /delete
    psbl_de_preprocess_32hz_data, id, trange = utr, rootdir = rootdir
endfor

end
