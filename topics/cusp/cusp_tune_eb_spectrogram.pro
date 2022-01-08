; Tune up the scale and zrange for E and B field.
; Need to run polar_sdt_prep_poynting_flux or fast_sdt_prep_poynting_flux.

pro cusp_tune_eb_spectrogram, dename, dbname, tr = tr, ctr = ctr

    ; plot settings.
    sgopen, 0, xsize = 9, ysize = 7, /inch
    device, decomposed = 0
    loadct2, 43
    posde = [0.1,0.15,0.45,0.95]
    posdb = [0.6,0.15,0.95,0.95]
    faclabs = ['v','p','b']
    rgb = [6,4,2]
    !p.font = 1 & !p.charsize = 1.5
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', 1
    time_stamp, /off
    tlimit, /full

    ; read data.
    pre = strmid(dename,0,2)+'_'
    get_data, dename, t0
    dr = sdatarate(t0)
    nscale = 60
    timsc0 = dr*50 & timsc1 = (t0[-1]-t0[0])*0.25
    if n_elements(ctr) ne 0 then timsc1 = 2*(ctr[1]-ctr[0])
    if n_elements(tr) eq 0 then tr = [t0[0],t0[-1]]
    
    ; split dE to plot.
    stplot_split, dename, labels = faclabs
    devname = dename+'_comp1'
    get_data, dename, t0, tmp
    tmp = smooth(tmp, [30/dr,1], /edge_mirror, /nan)
    desmname = dename+'1'
    store_data, desmname, t0, tmp
    stplot_split, desmname, labels = faclabs
    tmp = tmp[where(t0 ge tr[0] and t0 le tr[1]),*]
    ylim, tnames(desmname+'_comp?'), min(tmp), max(tmp), 0
    erange = [-1,1]*max(tmp)/nscale*5
    for i = 0, 2 do begin
        tmp = desmname+'_comp'+string(i+1,format='(I0)')
        options, tmp, 'ytitle', 'dE'+faclabs[i]+'!C(mV/m)'
        options, tmp, 'labels', faclabs[i]
        options, tmp, 'colors', 0
    endfor
    
    ; split dB to plot.
    stplot_split, dbname, labels = faclabs
    dbpname = dbname+'_comp2'
    get_data, dbname, t0, tmp
    tmp = smooth(tmp, [30/dr,1], /edge_mirror, /nan)
    dbsmname = dbname+'1'
    store_data, dbsmname, t0, tmp
    stplot_split, dbsmname, labels = faclabs
    tmp = tmp[where(t0 ge tr[0] and t0 le tr[1]),*]
    ylim, tnames(dbsmname+'_comp?'), min(tmp), max(tmp), 0
    brange = [-1,1]*max(tmp)/nscale*5
    for i = 0, 2 do begin
        tmp = dbsmname+'_comp'+string(i+1,format='(I0)')
        options, tmp, 'ytitle', 'dB'+faclabs[i]+'!C(nT)'
        options, tmp, 'labels', faclabs[i]
        options, tmp, 'colors', 0
    endfor
    
    ; calculate the spectrogram.
    timscs = smkgmtrc(timsc0, timsc1, nscale, 'n')
    devmatname = devname+'_mat'
    dbpmatname = dbpname+'_mat'
    stplot_mat, devname, scale = timscs, newname = devmatname
    stplot_mat, dbpname, scale = timscs, newname = dbpmatname

    vars = [devmatname,dbpmatname]
    options, vars, 'spec', 1
    options, vars, 'ytitle', 'period (s)'
    options, devmatname, 'ztitle', 'dEv (mV/m)'
    options, devmatname, 'zrange', erange
    options, dbpmatname, 'ztitle', 'dBp (nT)'
    options, dbpmatname, 'zrange', brange

    ; plot the spectrogram.
    erase
    vars = [tnames(desmname+'_comp?'),devmatname]
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position=posde)
    tplot, vars, position = pos, tr = tr, /noerase
    if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
    
    vars = [tnames(dbsmname+'_comp?'),dbpmatname]
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position=posdb)
    tplot, vars, position = pos, tr = tr, /noerase
    if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
    
    ; start to tune the spectrogram.
    nudge = 'y'
    while nudge eq 'y' do begin
        read, nudge, prompt = 'continue to tune up scale? (y/n): '
        if nudge eq 'n' then break
        read, timsc0, timsc1, prompt = 'scale in time (min, max): '
        read, erange, prompt = 'de zrange (min, max): '
        read, brange, prompt = 'db zrange (min, max): '
        
        timscs = smkgmtrc(timsc0, timsc1, nscale, 'n')    
        stplot_mat, devname, scale = timscs, newname = devmatname
        stplot_mat, dbpname, scale = timscs, newname = dbpmatname
    
        vars = [devmatname,dbpmatname]
        options, vars, 'spec', 1
        options, vars, 'ytitle', 'period (s)'
        options, devmatname, 'ztitle', 'dEv (mV/m)'
        options, devmatname, 'zrange', erange
        options, dbpmatname, 'ztitle', 'dBp (nT)'
        options, dbpmatname, 'zrange', brange
        
        erase
        vars = [tnames(desmname+'_comp?'),devmatname]
        nvar = n_elements(vars)
        pos = sgcalcpos(nvar, position=posde)
        tplot, vars, position = pos, /noerase
        if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
        
        vars = [tnames(dbsmname+'_comp?'),dbpmatname]
        nvar = n_elements(vars)
        pos = sgcalcpos(nvar, position=posdb)
        tplot, vars, position = pos, /noerase
        if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
    endwhile
    
    ; filter.
    nudge = 'y' & timfts = [0]
    while nudge eq 'y' do begin
        print, timfts
        read, nudge, prompt = 'continue to tune up filter? (y/n): '
        if nudge eq 'n' then break
        ctime, tmp, timfts, /exact
        ; do the filtering.
        filters = stofilter(timfts) & nfilter = n_elements(filters)/2
        ifilters = string(indgen(nfilter)+1, format='(I0)')     ; filter id string.
        for i = 0, nfilter-1 do begin
            stplot_filter, devmatname, 'mat', filter = filters[i,*], ifilter = ifilters[i]
            stplot_filter, dbpmatname, 'mat', filter = filters[i,*], ifilter = ifilters[i]
        endfor
        erase
        vars = tnames(devmatname+['','f?'])
        nvar = n_elements(vars)
        pos = sgcalcpos(nvar, position=posde)
        tplot, vars, position = pos, /noerase, title = ''
        plot, /noerase, /nodata, ystyle = 5, xstyle = 5, /ylog, $
            position = pos[*,0], tr, [timsc0,timsc1]
        for j = 0, nfilter do plots, tr, timfts[[j,j]]
        if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
        vars = tnames(dbpmatname+['','f?'])
        nvar = n_elements(vars)
        pos = sgcalcpos(nvar, position=posdb)
        tplot, vars, position = pos, /noerase, title = ''
        plot, /noerase, /nodata, ystyle = 5, xstyle = 5, /ylog, $
            position = pos[*,0], tr, [timsc0,timsc1]
        for j = 0, nfilter do plots, tr, timfts[[j,j]]
        if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
    endwhile
    
    print, 'scale in time (min, max): ', timsc0, timsc1
    print, 'de zrange (min, max): ', erange
    print, 'db zrange (min, max): ', brange
    print, 'filter in time: ', timfts
    
    ; final plot.
    erase
    case pre of
        'po_': pre = 'polar_'
        'fa_': pre = 'fast_'
    endcase
    for i = 0, 2 do begin
        tmp = string(i+1,format='(I0)')
        get_data, dename+'_comp'+tmp, t0, f0
        get_data, dename+'1_comp'+tmp, t0, f1
        store_data, dename+'0_comp'+tmp, t0, [[f0],[f1]], $
            limits = {colors:[-1,6], ytitle:'dE'+faclabs[i]+'!C(mV/m)', labels:faclabs[i]}
        get_data, dbname+'_comp'+tmp, t0, f0
        get_data, dbname+'1_comp'+tmp, t0, f1
        store_data, dbname+'0_comp'+tmp, t0, [[f0],[f1]], $
            limits = {colors:[-1,6], ytitle:'dB'+faclabs[i]+'!C(nT)', labels:faclabs[i]}
    endfor
    
    vars = [tnames(dename+'0_comp?'),tnames(dename+'1_comp?'),devmatname]
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position=posde)
    tplot, vars, position = pos, /noerase
    if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
    
    vars = [tnames(dbname+'0_comp?'),tnames(dbname+'1_comp?'),dbpmatname]
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position=posdb)
    tplot, vars, position = pos, /noerase
    if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
    
    ofn = '' & read, ofn, prompt = 'output to eps/pdf? (y/n):'
    if ofn eq 'n' then return
    ; **** first plot.
    if ofn eq 'y' then ofn = shomedir()+'/'+$
        time_string(t0[0],tformat='YYYY_MMDD_')+pre+'mat_spectrogram.eps'
    sgopen, ofn
    sgindexcolor, ct = 43
    vars = [tnames(dename+'0_comp?'),tnames(dename+'1_comp?'),devmatname]
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position=posde)
    tplot, vars, position = pos, /noerase
    if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
    
    vars = [tnames(dbname+'0_comp?'),tnames(dbname+'1_comp?'),dbpmatname]
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position=posdb)
    tplot, vars, position = pos, /noerase
    if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
    sgpsclose, /pdf
    ; **** second plot.
    if ofn eq 'y' then ofn = shomedir()+'/'+$
        time_string(t0[0],tformat='YYYY_MMDD_')+pre+'filtered_field.eps'
    sgpsopen, ofn
    sgindexcolor, ct = 43
    vars = tnames(devmatname+['','f?'])
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position=posde)
    tplot, vars, position = pos, /noerase, title = ''
    plot, /noerase, /nodata, ystyle = 5, xstyle = 5, /ylog, $
            position = pos[*,0], tr, [timsc0,timsc1]
    for j = 0, nfilter do plots, tr, timfts[[j,j]]
    if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
    vars = tnames(dbpmatname+['','f?'])
    nvar = n_elements(vars)
    pos = sgcalcpos(nvar, position=posdb)
    tplot, vars, position = pos, /noerase, title = ''
    plot, /noerase, /nodata, ystyle = 5, xstyle = 5, /ylog, $
            position = pos[*,0], tr, [timsc0,timsc1]
    for j = 0, nfilter do plots, tr, timfts[[j,j]]
    if n_elements(ctr) ne 0 then timebar, ctr, color = 6, thick = 2
    sgclose
    
    store_data, '*', /delete
end




;eventid = '1998_0925_05'
;
;rootdir = shomedir()+'/Google Drive/works'
;fn = rootdir+'/works/cusp/cusp_list_of_conjun_9_10.log'
;info = cusp_read_conjun_list(fn, event = eventid)
;
;potr = info.polar.plot_time & poctr = info.polar.cusp_time
;if potr[1] lt potr[0] then potr[1]+= 86400d
;if poctr[1] lt poctr[0] then poctr[1]+= 86400d
;
;fatr = info.fast.plot_time & factr = info.fast.cusp_time
;if factr[1] lt factr[0] then factr[1]+= 86400d
;if fatr[1] lt fatr[0] then fatr[1]+= 86400d
;
;orb = string(info.fast.orbit,format='(I05)')
;pofn = rootdir+'/data/cusp/po_sdt_fld_'+eventid+'.sdt'
;fafn = rootdir+'/data/cusp/fa_sdt_fld_'+eventid+'_'+orb+'.tplot'
;
;pre = 'po_'
;polar_sdt_prep_poynting_flux, pofn, tr = potr, /noplot
;cusp_tune_eb_spectrogram, pre+'de_fac', pre+'db_fac', tr = potr, ctr = poctr
;
;pre = 'fa_'
;fast_sdt_prep_poynting_flux, fafn
;cusp_tune_eb_spectrogram, pre+'de_fac', pre+'db_fac', tr = fatr, ctr = factr




eventid = '2000_1024_10'

rootdir = shomedir()+'/Google Drive/works'
fn = rootdir+'/works/cusp/cusp_list_of_polar_2-4Re.log'
info = cusp_read_conjun_list(fn, event = eventid, /nofast)

potr = info.polar.plot_time & poctr = info.polar.cusp_time
if potr[1] lt potr[0] then potr[1]+= 86400d
if poctr[1] lt poctr[0] then poctr[1]+= 86400d

pofn = rootdir+'/data/cusp/po_sdt_fld_'+eventid+'.sdt'

pre = 'po_'
polar_sdt_prep_poynting_flux, pofn, tr = potr, /noplot
cusp_tune_eb_spectrogram, pre+'de_fac', pre+'db_fac', tr = potr, ctr = poctr

end
