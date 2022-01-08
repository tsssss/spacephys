;+
; on daily bases, preprocess efield for certain probe.
; subtract background E field near earth (<3 Re)
; remove bad E field when Eu and Ev do not match and
; when they are asymetric around 0.
; 
; calc E GSM from V, save the quaternion of rotation matrix uvw2gsm, 
; from which wsc_gsm can be extracted.
;-


pro rbsp_preprocess_efield, ut0, probes = tprobe, no_spice_load = no_spice_load

;---Constant.
    secofday = 86400d
    dr0 = 1d/16     ; data rate, 16 Hz.
    spinrate = 12   ; sec, 10.9378 sec.
    blen0s = [100d,100,12]  ; m, twice of boom lengths.
    fshrts = [1d,1,1]       ; use 1 temperarily.
    blen1s = blen0s*fshrts
    maxerng = 400d
    minerng = 5d

    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re


;---Settings.
    tbuff = 3600d   ; 1 hour of data buffer.
    if n_elements(tprobe) eq 0 then tprobe = 'a'
    pre0 = 'rbsp'+tprobe+'_'

    rootdir = shomedir()+'/rbsp_efield_preprocess'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
    
    utr = ut0-(ut0 mod secofday)+[0,secofday]
    utr1 = utr+[-1,1]*tbuff
    ut2016 = time_double('2015-10-01')
    date = time_string(utr[0],tformat='YYYY_MMDD')        
    
    
    rgb = [6,4,2]
    xyz = ['x','y','z']
    uvw = ['u','v','w']


    device, decomposed = 0
    loadct2, 43

    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 10
    tplot_options, 'labflag', -1
    tplot_options, 'xticklen', -0.03
    tplot_options, 'yticklen', -0.005
    tplot_options, 'constant', 0
    tplot_options, 'xmargin', [18,10]
    
;goto, startplot

    
;---load spice product.
; rbspx_[pos_gsm,mlt,lshell,dis]

    spice = sread_rbsp_spice_product(utr1, probe=tprobe)
    tuts = spice.ut_pos
    store_data, pre0+'pos_gsm', tuts, spice.pos_gsm, limits={ytitle:'(Re)', colors:rgb, labels:'GSM '+xyz}
    store_data, pre0+'mlt', tuts, spice.mlt, limits={ytitle:'(hr)',labels:'MLT'}
    store_data, pre0+'lshell', tuts, spice.lshell, limits={ytitle:'(Re)',labels:'L-shell'}
    store_data, pre0+'mlat', tuts, spice.mlat, limits = {ytitle:'MLat (deg)'}
    store_data, pre0+'dis', tuts, spice.dis, limits = {ytitle:'Dist (Re)'}
    tuts = spice.ut_cotran
    store_data, pre0+'quvw2gsm', tuts, spice.q_uvw2gsm
    
    
;---load vsvy and calc esvy (instead of load esvy).
; rbspx_[vsc,vsvy]
    
    tvar = pre0+'vsc'    
    dat = sread_rbsp_efw_l2(utr1, probes=tprobe, type='vsvy')
    if size(dat,/type) ne 8 then begin
        message, 'no V data ...', /continue
        return
    endif
    
    uts = sfmepoch(dat.epoch, 'unix')
    nrec = n_elements(uts)
    vsvy = dat.vsvy
    ; calibrate vsvy, tried rbsp_efw_get_cal_params, but really no effect
    ; since offset is 0, gain is not used.  
    ; end up using V12, may consider more complex ways.
    vsc = mean(vsvy[*,0:1], dimension = 2)
    if utr[0] gt ut2016 then vsc = mean(vsvy[*,2:3], dimension = 2)
    store_data, tvar, uts, vsc, limits = {ytitle:'(V)', labels:'Vsc'}
    store_data, pre0+'vsvy', uts, vsvy[*,0:3], limits = $
        {ytitle:'(V)', colors:[1,2,3,4], labels:'V'+['1','2','3','4']}
  
    
    ; calc E uvw.
    ; rbspx_euvw

    tvar = pre0+'euvw'    
    ; the following is based on linear regression on v1,v2 and eu, v3,v4 and ev.
    eu = (vsvy[*,0]-vsvy[*,1])/blen1s[0]*1e3   ; V -> V/m -> mV/m.
    ev = (vsvy[*,2]-vsvy[*,3])/blen1s[1]*1e3
    ;ew = (vsvy[*,4]-vsvy[*,5])/blen1s[2]*1e3
    ew = dblarr(nrec)
    ; remove dc-offset.
    nspin = 1
    tnrec = nspin*spinrate/dr0
    eu = eu-smooth(eu, tnrec, /edge_truncate, /nan)
    ev = ev-smooth(ev, tnrec, /edge_truncate, /nan)
    store_data, tvar, uts, [[eu],[ev],[ew]], limits = $
        {ytitle:'(mV/m)', labels:'E'+uvw, colors:rgb}
    tvar = pre0+'esvy'  ; this will be further modified.
    store_data, tvar, uts, [[eu],[ev],[ew]], limits = $
        {ytitle:'(mV/m)', labels:'E'+uvw, colors:rgb}


    ; coefficient to remove E field near earth (< 3 Re).
    ; rbspx_[e0,coef_e0].

;    ; make a decaying window around 3 Re at 0.5 Re width.
    mindis0 = 3.50  ; Re.
    deldis0 = 0.25  ; Re.
    get_data, pre0+'dis', tuts, dis
    dis = interpol(dis, tuts, uts)
    drec = atan((dis-mindis0)/deldis0)
    ; norm to 50:1, 1-> 1 spin.
    drec = ((drec-min(drec))/(max(drec)-min(drec))*600+1)*spinrate/dr0*0.2
    
    ; remove a min background.
    emag = sqrt(eu^2+ev^2)    
    emagbg = [] & emagut = []
    i1 = 0
    while i1 lt nrec-1 do begin
        i2 = i1+drec[i1]<(nrec-1)
        emagbg = [emagbg,min(emag[i1:i2-1],/nan)]
        emagut = [emagut,uts[(i1+i2)*0.5]]
        i1 = i2
    endwhile
    emagbg = interpol(emagbg, emagut, uts, /nan)
    
    ec0 = 0>(1-emagbg/emag)<1
    ec1 = [] & eut = []
    drec = spinrate*4/dr0
    i1 = 0
    while i1 lt nrec-1 do begin
        i2 = i1+drec<(nrec-1)
        ec1 = [ec1,max(ec0[i1:i2-1],/nan)]
        eut = [eut,uts[(i1+i2)*0.5]]
        i1 = i2
    endwhile
    ec = interpol(ec1, eut, uts, /nan)
    
    store_data, pre0+'coef_e0', uts, ec
    store_data, pre0+'e0', uts, [[emag],[emagbg]], limits=$
        {ytitle:'(mV/m)',colors:[6,2],labels:['|E|','|E| BG']}

;options, pre0+'e0', 'yrange', [0,2]
;tplot, pre0+['euvw','coef_e0','e0','dis']
;stop


    ; apply the coefficient to E field.
    tvar = pre0+'esvy'
    get_data, tvar, tuts, esvy
    for i=0, 1 do esvy[*,i] *= ec
    store_data, tvar, tuts, esvy
    eu = esvy[*,0]
    ev = esvy[*,1]
    
;    tplot, pre0+['e0','esvy','dis']


;---Make flags to tell bad E field.
    maxde0 = 2d ; mV/m.
    padt = 300d ; sec.
    
        
    ; get the upper and lower envelope over several spins.
    ampu = sqrt(eu^2+shift(eu,spinrate*0.25/dr0)^2)
    ampv = sqrt(ev^2+shift(ev,spinrate*0.25/dr0)^2)
    ampu = smooth(ampu, spinrate/dr0*0.2, /nan, /edge_truncate) ; remove small spikes.
    ampv = smooth(ampv, spinrate/dr0*0.2, /nan, /edge_truncate)
    ampu = smooth(ampu, spinrate/dr0, /nan, /edge_truncate)     ; remove spin tone.
    ampv = smooth(ampv, spinrate/dr0, /nan, /edge_truncate)
    store_data, pre0+'eamp', uts, [[ampu],[ampv]], limits={ytitle:'(mV/m)',labels:['u','v'],colors:[6,2]}

    ; check the amplitude of high-freq waves.    
    tf0s = abs(eu-smooth(eu, spinrate/dr0*0.2, /nan, /edge_truncate))
    drec = spinrate/dr0
    tf1s = [] & tuts = []
    i1 = 0
    while i1 lt nrec-1 do begin
        i2 = i1+drec<(nrec-1)
        if i2-i1 lt 1 then continue
        tdat = tf0s[i1:i2]
        tdat = tdat[sort(tdat)]
        tf1s = [tf1s,max(tdat[0:(i2-i1)*0.8])]
        tuts = [tuts,uts[(i1+i2)*0.5]]
        i1 = i2
    endwhile
    henv = interpol(tf1s, tuts, uts, /nan)
    henv = smooth(henv, spinrate/dr0*10, /nan, /edge_truncate)+1    ; make it always loose the criteria.  

    damp = ampu-ampv
    ;damp = smooth(damp, spinrate/dr0, /nan, /edge_truncate)
    flags = (abs(damp) ge maxde0*henv)
    store_data, pre0+'damp', uts, damp
    store_data, pre0+'flag_euv', uts, flags, limits={yrange:[-0.5,1.5]}
    ;tplot, [pre0+['euvw','damp','henv','flag_euv']]


    ; check eclipse time.
    tmp = sread_rbsp_eclipse_time(utr, probes = tprobe)
    flags = interpol(double(tmp.flags),tmp.uts,uts)
    idx = where(flags ne 0, cnt)
    if cnt ne 0 then flags[idx] = 1
    store_data, pre0+'flag_eclipse', uts, flags, limits = $
        {ytitle:'', yrange:[-0.5,1.5], yticks:1, yminor:0, ystyle:1, $
        labels:['1:eclipse']}
        
    ; check sdt time.
    tmp = sread_rbsp_efw_sdt_time(utr[0], probes = tprobe)
    flags = interpol(double(tmp.flags),tmp.uts,uts)
    idx = where(flags ne 0, cnt)
    if cnt ne 0 then flags[idx] = 1
    store_data, pre0+'flag_sdt', uts, flags, limits = $
        {ytitle:'', yrange:[-0.5,1.5], yticks:1, yminor:0, ystyle:1, $
        labels:['1:sdt']}
    

    ; combine flags.
    vars = pre0+'flag_'+['euv','eclipse','sdt']
    flags = bytarr(nrec)
    foreach tvar, vars do begin
        get_data, tvar, uts, dat
        flags = flags or dat
    endforeach

    idx = where(flags eq 1, cnt)
    drec = padt/dr0
    if cnt ne 0 then begin
        for i=0, cnt-1 do flags[(idx[i]-drec)>0:(idx[i]+drec)<(nrec-1)] = 1
        idx = where(flags eq 1)
        tvar = pre0+'esvy'
        get_data, tvar, tuts, esvy
        esvy[idx,*] = !values.d_nan
        store_data, tvar, uts, esvy
    endif
    store_data, pre0+'bade_flag', uts, flags, limits=$
        {labels:'1: bad E', yrange:[-0.5,1.5], ystyle:1, yticks:1, ytickv:[0,1], panel_size:0.4, ytitle:''}
    
    idx = where(flags eq 1, cnt)
    if cnt ne 0 then begin
        eu[idx] = !values.d_nan
        ev[idx] = !values.d_nan
    endif
    

;    ;---test!!! correct for boom shorting and dc-offset.   
;    ; balance amplitude.    
;    tamp = (ampu+ampv)*0.5
;    bsu = tamp/ampu
;    bsv = tamp/ampv
;    eu = eu*bsu
;    ev = ev*bsv
;
;    store_data, pre0+'bs', uts, [[bsu],[bsv]], limits=$
;        {labels:['u','v'], colors:[6,4]}
;    get_data, pre0+'esvy', limits=lim
;    store_data, pre0+'esvy', uts, [[eu],[ev],[ew]], limits=lim
;    store_data, pre0+'emag', uts, [[sqrt(eu^2+ev^2)],[emag*ec]], $
;        limits={ytitle:'(mV/m)',labels:['|E|','|E|2'],colors:[6,2]}

    ;tplot, pre0+['bs','eamp','esvy','emag','eamp']
    ;stop


;tplot, pre0+['esvy','euvw']
;stop

    ; clean up and trim data to utr.
    vars = pre0+['vsvy','euv_mag','euv_mag_diff_abs','euv_avg','e0']
    store_data, vars, /delete
    
    vars = pre0+['euvw','esvy']
    foreach tvar, vars do begin
        get_data, tvar, uts, dat
        idx = where(uts ge utr[0] and uts le utr[1])
        store_data, tvar, uts[idx], dat[idx,*]
    endforeach
    
    

        
    ;---rotate uvw into gsm.
        get_data, pre0+'esvy', uts, dat
        eu = dat[*,0]
        ev = dat[*,1]
        ew = dat[*,2]
        
        get_data, pre0+'quvw2gsm', tuts, quvw2gsm
        quvw2gsm = qslerp(quvw2gsm, tuts, uts)
        muvw2gsm = transpose(qtom(quvw2gsm))
        
        ex = eu*muvw2gsm[0,0,*] + ev*muvw2gsm[1,0,*] + ew*muvw2gsm[2,0,*]
        ey = eu*muvw2gsm[0,1,*] + ev*muvw2gsm[1,1,*] + ew*muvw2gsm[2,1,*]
        ez = eu*muvw2gsm[0,2,*] + ev*muvw2gsm[1,2,*] + ew*muvw2gsm[2,2,*]
        
        tvar = pre0+'egsm'
        store_data, tvar, uts, [[ex],[ey],[ez]], limits = $
            {ytitle:'(mV/m)', colors:[6,4,2], labels:'GSM E'+['x','y','z']}
        
        vars = pre0+['flag_'+['euv','eclipse','sdt'],$
            'bade_flag','eamp','quvw2gsm','vsc','esvy','coef_e0']
        store_data, vars, /delete


;---Plot to show: remove bad E, co-rotation E.
;        xsz = 11
;        ysz = 8.5
;        
;
;        tdir = rootdir+'/rbsp'+tprobe+'/'+time_string(ut0,tformat='YYYY')
;        if file_test(tdir,/directory) eq 0 then file_mkdir, tdir
;        ofn = tdir+'/'+date+'_rbsp'+tprobe+'_preprocess_field.pdf'
;ofn = 0
;        
;        sgopen, ofn, xsize = xsz, ysize = ysz, /inch
;        
;        device, decomposed = 0
;        loadct2, 43
;        
;        xchsz = double(!d.x_ch_size)/!d.x_size
;        ychsz = double(!d.y_ch_size)/!d.y_size
;
;
;        titl = 'RBSP-'+strupcase(tprobe)+' '+time_string(utr[0],tformat='YYYY-MM-DD')+' E field preprocess: remove co-rotation E, bad E'        
;        
;        vars = pre0+['eamp','bade_flag','euvw'+uvw[0:1],'egsm'+xyz]
;        flabs = ['a','b','c','d','e','f','g']+'. '+$
;            ['|Euv|','Flag','Eu','Ev','GSM E'+xyz]
;        
;        
;        tvar = pre0+'euvw'
;        stplot_split, tvar, newnames=tvar+uvw, colors=rgb
;        tvar = pre0+'egsm'
;        stplot_split, tvar, newnames=tvar+xyz, colors=rgb
;        
;        get_data, pre0+'euvw', uts, dat
;        erng = max(abs(dat),/nan)
;        erng = ceil(erng/10)*10
;        erng = erng<maxerng
;        tvar = pre0+['eamp','euvw'+uvw] 
;        options, tvar, 'yrange', [-1,erng]
;        options, tvar, 'ystyle', 1
;        options, tvar, 'yticks', 4
;        options, tvar, 'yminor', 5
;        
;        tvar = pre0+'egsm'
;        get_data, tvar, uts, dat
;        erng = max(abs(dat),/nan)
;        erng = ceil(erng/10)*10
;        erng = erng>minerng
;        
;        tvar = tvar+xyz
;        options, tvar, 'yrange', [-1,1]*erng
;        options, tvar, 'ystyle', 1
;        options, tvar, 'yticks', 4
;        options, tvar, 'yminor', 5
;        
;         
;        tplot, vars, trange=utr, figlabel=flabs, get_plot_position=poss
;        xyouts, (poss[0,0]+poss[2,0])*0.5, poss[3,0]+ychsz*0.8, /normal, alignment=0.5, titl, charsize=1.2
;
;        idx = where(vars eq pre0+'bade_flag')
;        tpos = reform(poss[*,idx])
;        plot, utr, [-0.5,1.5], /nodata, /noerase, position=tpos, $
;            xstyle=5, ystyle=5
;        vars = pre0+'flag_'+['euv','eclipse','sdt']
;        foreach tvar, vars, i do begin
;            get_data, tvar, uts, dat
;            oplot, uts, (dat<(1-(i+1)*0.1)), color=rgb[i]
;        endforeach
;stop
;        tvar = pre0+['euvw'+uvw,'egsm'+xyz]
;        store_data, tvar, /delete
;
;
;        sgclose
    
end

probes = ['a']
utr0 = time_double(['2012-09-25','2016-12-31'])
utr0 = time_double(['2013-01-12','2013-01-13']) ; steps in ec.
utr0 = time_double(['2013-11-27','2013-11-28']) ; should reduce 2fsp.
;utr0 = time_double(['2013-05-01','2013-05-02']) & probes = 'b'
;utr0 = time_double(['2012-11-14','2012-11-15']) 
;utr0 = time_double(['2012-09-25','2012-09-26'])
utr0 = time_double(['2013-11-27','2013-11-28']) ; trouble day.
;utr0 = time_double(['2015-03-16','2015-03-17']) ; trouble day, 180 deg phase diff.
;utr0 = time_double(['2013-03-17','2013-03-18']) ; trouble day.
;utr0 = time_double(['2012-10-09','2012-10-10']) ; trouble day.
;utr0 = time_double(['2012-09-25','2012-09-26'])
;utr0 = time_double(['2015-01-10','2015-01-11']) ; contain sdt.

;utr0 = time_double(['2013-06-07','2013-06-07']) ; good wave & bad data.
;probes = ['b']

utr0 = time_double(['2013-06-07','2013-06-07'])

    

;---Loop through each day.
    uts = smkarthm(utr0[0], utr0[1], 86400, 'dx')
    foreach tut, uts do foreach tprobe, probes do $
        rbsp_preprocess_efield, tut, no_spice_load=no_spice_load, probes=tprobe


end
