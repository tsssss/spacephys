;+
; on daily bases, survey plot includes both -A and -B
; vars: vsc, |E| high res, B1 and B2 availability (bars).
;-


pro psbl_de_survey_on_bad_e_removal, ut0, no_spice_load = no_spice_load


    probes = ['a','b']
    rootdir = shomedir()+'/psbl_de_32hz/vsc_survey2'
    if file_test(rootdir,/directory) eq 0 then file_mkdir, rootdir
    
    dt0 = 86400d
    utr = ut0-(ut0 mod dt0)+[0,dt0]
    dr0 = 1d/16     ; data rate, 16 Hz.
    ut2016 = time_double('2015-10-01')
    spinrate = 11   ; sec.
    
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    
    rgb = [6,4,2]
    device, decomposed = 0
    loadct2, 43
    tplot_options, 'labflag', -1
    
;goto, startplot

    foreach tprobe, probes do begin
        pre0 = 'rbsp'+tprobe+'_'
        
        ; load vsc, calc euvw, egse.
        rbsp_preprocess_efield, utr[0], probes = tprobe, /no_spice_load
        
        ; from pos_gse to pos_gse_[x,y,z].
        tvar = pre0+'pos_gse'
        stplot_split, tvar, newnames = tvar+['x','y','z']
        
        ; from euvw calc emag.
        tvar = pre0+'emag'
        get_data, pre0+'euvw', uts, dat
        store_data, tvar, uts, snorm(dat), limits = $
            {ytitle:'|E|!C(mV/m)', labels:'RBSP-'+strupcase(tprobe)+'!C  |E| 16 hz'}
            
        ; get flag for burst data.
        ; **** load B1 and B2 times.
        tvar = pre0+'flag_burst'
        remroot = 'http://rbsp.space.umn.edu/data/rbsp/burst_playback'
        locroot = spreproot('rbsp')
        prb = tprobe
        vsn = 'v[0-9.]{2}'
        ext = 'txt'
        types = ['vb1','vb2']
        
        nrec = n_elements(uts)
        flags = bytarr(nrec,2)
        for j = 0, 1 do begin
            type = types[j]
            
            baseptn = 'rbsp'+prb+'_efw_'+type+'_playback_YYYYMMDD_'+vsn+'.'+ext
            rempaths = [remroot,'rbsp'+prb,type+'_playback','YYYY',baseptn]
            locpaths = [locroot,'rbsp'+prb,'efw/burst_playback',type,'YYYY',baseptn]
            
            remfns = sprepfile(utr, paths = rempaths)
            locfns = sprepfile(utr, paths = locpaths)
            nfn = n_elements(locfns)
            for i = 0, nfn-1 do begin
                basefn = file_basename(locfns[i])
                locpath = file_dirname(locfns[i])
                rempath = file_dirname(remfns[i])
                locfns[i] = sgetfile(basefn, locpath, rempath)
            endfor
            idx = where(locfns ne '', nfn)
            if nfn eq 0 then continue
            
            locfns = locfns[0]
            nline = file_lines(locfns)
            lines = strarr(nline)
            openr, lun, locfns, /get_lun
            readf, lun, lines
            free_lun, lun
            lines = lines[4:*]  ; exclude header.
            nline = n_elements(lines)
            
            for i = 0, nline-1 do begin
                tutr = time_double(strsplit(lines[i],' ',/extract))
                idx = where(uts ge tutr[0] and uts le tutr[1], cnt)
                if cnt ne 0 then flags[idx,j] = 1
            endfor
        endfor
        
        store_data, tvar, uts, flags, limits = $
            {ytitle:'', yrange:[-0.5,1.5], yticks:1, yminor:0, $
            colors:[6,4], labels:['B1','B2']}

    endforeach


;startplot:

    ; **** generate plot.
    xsz = 11
    ysz = 8.5
    ofn = 0
    ofn = rootdir+'/rbsp_daily_bad_e_removal_'+time_string(utr[0],tformat='YYYY_MMDD')+'.pdf'

    sgopen, ofn, xsize = xsz, ysize = ysz, /inch

    device, decomposed = 0
    loadct2, 43
    tplot_options, 'version', 0
    tplot_options, 'num_lab_min', 10
    tplot_options, 'labflag', -1


    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    panx0 = 0.15
    panx1 = 0.90


pos1 = [panx0,0.7,panx1,0.95]    
    pandy = pos1[3]-pos1[1]
    pany1 = pos1[1]+0.3*pandy
    pany2 = pos1[1]+0.7*pandy
    txr = utr
    colors = [6,2]
    
    for i = 0, 1 do begin
        tprobe = probes[i]
        
        get_data, 'rbsp'+tprobe+'_vsc', uts, dat

        tpos = [panx0,pany1,panx1,pany2]+[0,1,0,-1]*ychsz*0.2
        tyr = [-10,0]
        plot, txr, tyr, xrange = txr, yrange = tyr, $
            position = tpos, /noerase, /nodata, $
            xstyle = 1, ystyle = 1, xtickformat = '(A1)', $
            ytitle = 'Vsc!C(V)', $
            yticks = 2, yminor = 6
        oplot, uts, dat, color = colors[i]
        
        tpos = [panx0,pany2,panx1,pos1[3]]
        tyr = [0,200]
        plot, txr, tyr, xrange = utr, yrange = tyr, $
            position = tpos, /noerase, /nodata, $
            xstyle = 1, ystyle = 1, xtickformat = '(A1)', $
            min_value = 0, $
            yticks = 2, yminor = 5
        oplot, uts, dat, color = colors[i]
        
        tpos = [panx0,pos1[1],panx1,pany1]
        tyr = [-110,-10]
        plot, txr, tyr, xrange = utr, yrange = tyr, $
            position = tpos, /noerase, /nodata, $
            xstyle = 1, ystyle = 1, xtickformat = '(A1)', $
            yticks = 2, yminor = 5
        oplot, uts, dat, color = colors[i]
        
        xyouts, panx1+2*xchsz, pany2-ychsz*(i+1), /normal, $
            'RBSP-'+strupcase(tprobe), color = colors[i]
    endfor

    xyouts, 0.5, pos1[3]+ychsz, /normal, alignment = 0.5, $
        'Vsc, |E|, '+time_string(utr[0], tformat='YYYY_MMDD'), charsize = 1.25



    emagyr = [0,10]
    barthick = 16


pos2 = [panx0,0.5,panx1,0.65]
    pandy = pos2[3]-pos2[1]
    pany1 = pos2[1]+0.5*pandy

    tpos = [panx0,pos2[1],panx1,pany1]
    tvar = 'rbspa_emag'
    labs = 'rbspa_'+['pos_gse'+['z','y','x'],'mlt','lshell']
    options, tvar, 'labels', 'RBSP-A!C  |E| 32 Hz'
    options, tvar, 'ylog', 0
    options, tvar, 'yrange', emagyr
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'constant', [5]
    tplot, tvar, var_label = labs, trange = utr, position = tpos, /noerase
    
    tpos = [panx0,pany1,panx1,pos2[3]]+[0,1,0,0]*ychsz*0.2
    tvar = 'rbspa_flag_burst'
    get_data, tvar, uts, dat
    colors = [4,6]
    tys = tpos[3]-ychsz*[0,1]
    labs = ['B1','B2']
    plot, utr, [0,1], /nodata, /noerase, xstyle = 5, ystyle = 5, position = tpos
    dat = double(dat)
    for i = 0, 1 do begin
        idx = where(dat[*,i] eq 1, cnt)
        if cnt ne 0 then begin
            dat[idx,i] = (tys[i]-tpos[1])/(tpos[3]-tpos[1])
            idx = where(dat[*,i] eq 0, cnt)
            dat[idx,i] = !values.d_nan
            plots, uts, dat[*,i], thick = barthick, color = colors[i]
        endif
        xyouts, panx1+2*xchsz, tys[i], /normal, labs[i], color = colors[i]
    endfor
    
    tvar = 'rbspa_emag'
    get_data, tvar, uts, dat
    tyr = [10,200]
    plot, txr, tyr, xrange = utr, yrange = tyr, $
        position = tpos, /noerase, /nodata, $
        xstyle = 1, ystyle = 1, xtickformat = '(A1)', $
        ytickv = [100,200], yticks = 1, yminor = 4
    oplot, uts, dat
    oplot, utr, 50+[0,0], linestyle = 1
    oplot, utr, 100+[0,0], linestyle = 1
    
    

pos3 = [panx0,0.17,panx1,0.32]

pandy = pos3[3]-pos3[1]
    pany1 = pos3[1]+0.5*pandy

    tpos = [panx0,pos3[1],panx1,pany1]
    tvar = 'rbspb_emag'
    labs = 'rbspb_'+['pos_gse'+['z','y','x'],'mlt','lshell']
    options, tvar, 'labels', 'RBSP-B!C  |E| 32 Hz'
    options, tvar, 'ylog', 0
    options, tvar, 'yrange', emagyr
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 5
    options, tvar, 'constant', [5]
    tplot, tvar, var_label = labs, trange = utr, position = tpos, /noerase
    
    
    tpos = [panx0,pany1,panx1,pos3[3]]+[0,1,0,0]*ychsz*0.2
    tvar = 'rbspb_flag_burst'
    get_data, tvar, uts, dat
    colors = [4,6]
    tys = tpos[3]-ychsz*[0,1]
    labs = ['B1','B2']
    plot, utr, [0,1], /nodata, /noerase, xstyle = 5, ystyle = 5, position = tpos
    dat = double(dat)
    for i = 0, 1 do begin
        idx = where(dat[*,i] eq 1, cnt)
        if cnt ne 0 then begin
            dat[idx,i] = (tys[i]-tpos[1])/(tpos[3]-tpos[1])
            idx = where(dat[*,i] eq 0, cnt)
            dat[idx,i] = !values.d_nan
            plots, uts, dat[*,i], thick = barthick, color = colors[i]
        endif
        xyouts, panx1+2*xchsz, tys[i], /normal, labs[i], color = colors[i]
    endfor
    
    tvar = 'rbspb_emag'
    get_data, tvar, uts, dat
    tyr = [10,200]
    plot, txr, tyr, xrange = utr, yrange = tyr, $
        position = tpos, /noerase, /nodata, $
        xstyle = 1, ystyle = 1, xtickformat = '(A1)', $
        ytickv = [100,200], yticks = 1, yminor = 4
    oplot, uts, dat
    oplot, utr, 50+[0,0], linestyle = 1
    oplot, utr, 100+[0,0], linestyle = 1
    
    
    
    sgclose 

    
end


utr0 = time_double(['2012-09-25','2016-12-31'])
utr0 = time_double(['2013-08-02','2016-12-31'])
;utr0 = time_double(['2012-11-14','2012-11-15']) 
;utr0 = time_double(['2012-09-25','2012-09-26'])

; load spice kernel for all.
no_spice_load = 1
rbsp_load_spice_kernels, trange = utr0

;utr0 = time_double(['2012-11-14','2012-11-15'])
uts = smkarthm(utr0[0], utr0[1], 86400, 'dx')
foreach tut, uts do psbl_de_survey_on_bad_e_removal, tut, no_spice_load = no_spice_load

rbsp_load_spice_kernels, trange = utr0, /unload

end



;; **** check relation b/w V and euvw.
;tvar = pre0+'euvw'
;
;dat = sread_rbsp_efw_l2(utr, probes = tprobe, type = 'euvw')
;euvw = sinterpol(dat.efield_uvw, sfmepoch(dat.epoch,'unix'), uts)
;store_data, tvar, uts, euvw, limits = $
;    {ytitle:'Euvw!C(mV/m)', labels:['Eu','Ev','Ew'], colors:[6,4,2]}
;    
;get_data, pre0+'vsvy', uts, vsvy
;tutr = time_double(['2012-11-14/02:04:50','2012-11-14/02:05:20'])
;idx = where(uts ge tutr[0] and uts le tutr[1])
;eu = euvw[idx,0]
;vsc0 = vsvy[idx,0]
;vsc1 = vsvy[idx,1]
;
;idx = where(finite(eu) and finite(vsc0) and finite(vsc1))
;x = [transpose(vsc0[idx]),transpose(vsc1[idx])]
;y = eu[idx]
;
;result = regress(x, y, sigma = sigma, const = const, measure_errors = err)
;print, 'Eu = a*V0 + b*V1, a and b equals to: '
;print, result
;print, 'offset: ', const
;
;
;idx = where(uts ge tutr[0] and uts le tutr[1])
;ev = euvw[idx,1]
;vsc2 = vsvy[idx,2]
;vsc3 = vsvy[idx,3]
;
;idx = where(finite(ev) and finite(vsc2) and finite(vsc3))
;x = [transpose(vsc2[idx]),transpose(vsc3[idx])]
;y = ev[idx]
;
;result = regress(x, y, sigma = sigma, const = const, measure_errors = err)
;print, 'Ev = c*V2 + d*V3, c and d equals to: '
;print, result
;print, 'offset: ', const
