
pro cusp_add_fpoint_fig, locroot

    if strmid(locroot,strlen(locroot)-1) eq '/' then $      ; remove trailing '/'.
        locroot = strmid(locroot,0,strlen(locroot)-1)

    fns = file_search(locroot+'/*')
    for i = 0, n_elements(fns)-1 do begin
        tfn = fns[i]
        if ~file_test(tfn, /directory) then continue
        id = strmid(tfn, strlen(locroot)+1)
        ofn = tfn+'/'+strmid(id,0,9)+'_footpoint.pdf'
        if file_test(ofn) then file_delete, ofn
        ofn = tfn+'/'+id+'_footpoint.pdf'
        if file_test(ofn) then continue
        ; need to add footpoint figure.
        rootdir = shomedir()+'/Google Drive/works'
        if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')
        logfile = rootdir+'/works/cusp/cusp_list_of_conjun.log'
        info = cusp_read_conjun_list(logfile, event = id)
        if size(info,/type) ne 8 then continue      ; no enough info.
        orootdir = locroot
        
        ; prepare.
        potr = info.polar.plot_time
        fatr = info.fast.plot_time
        potrcusp = info.polar.cusp_time
        fatrcusp = info.fast.cusp_time
        if n_elements(fatrcusp) eq 1 then continue  ; no fast cusp tr.
        if potr[1] lt potr[0] then potr[1]+= 86400d
        if fatr[1] lt fatr[0] then fatr[1]+= 86400d
        if potrcusp[1] lt potrcusp[0] then potrcusp[1]+= 86400d
        if fatrcusp[1] lt fatrcusp[0] then fatrcusp[1]+= 86400d
        
        ; polar pos.
        tinfo = info.polar
        fn = rootdir+'/data/cusp/po_sdt_fld_'+id+'.sdt'
        if file_test(fn) eq 0 then message, 'file does not exist ...'
        sdt = ssdtread(fn)
        pre = 'po_'
        t1  = sdt.var.polarinvariantlatitude.depend_0
        mlat = sdt.var.polarmaglat.value
        store_data, pre+'mlat', data = {x:t1, y:mlat}, limits = {ytitle:'MLat'}
        ilat = sdt.var.polarinvariantlatitude.value
        tmp = where(mlat lt 0)
        if tmp[0] ne -1 then ilat[tmp] *= -1
        store_data, pre+'ilat', data = {x:t1, y:ilat}, limits = {ytitle:'ILat'}
        mlt  = sdt.var.polarmlt.value
        store_data, pre+'mlt', data = {x:t1, y:mlt}, limits = {ytitle:'MLT'}
        dis  = sdt.var.polarspcraftdist.value
        store_data, pre+'dis', data = {x:t1, y:dis}, limits = {ytitle:'Dist (Re)'}

        ; fast pos.
        tinfo = info.fast
        fn = rootdir+'/data/cusp/fa_sdt_esa_'+strmid(id,0,9)+'_'+ $
            string(tinfo.orbit,format='(I05)')+'.tplot'
        if file_search(fn) eq '' then $
            fn = rootdir+'/data/cusp/fa_sdt_esa_'+id+'_'+$
            string(tinfo.orbit,format='(I05)')+'.tplot'
        tplot_restore, filename = fn
        
        ; plot footpoint.
        satnames = ['polar (*,b)','fast (+,r)']
        satcolors = [sgcolor('blue'),sgcolor('red')]
        satsyms = [2,1]     ; [*,+].
        pos = [0.1,0.1,0.9,0.9]
        thick = 2
        ofn = tfn+'/'+strmid(id,0,9)+'_footpoint.eps'
        sgpsopen, ofn, xsize = 6, ysize = 3, /inch
        sgtruecolor
        ; determine hemisphere.
        get_data, 'po_ilat', t0, poilat
        if interpol(poilat, t0, potrcusp[0]) ge 0 then begin    ; north hem.
            sgset_map, xrange = [90,270], pos = pos, color = sgcolor('black'), $
                ytickv = [50,60,70,80], ytickpos = 225, yticknudge = [-1.4,-0.4], $
                xtickpos = 47
        endif
        ; use cusp time, expand 1 cusp time on both sides.
        ; polar ilat, mlt.
        j = 0
        dt = 600        ; 10 min.
        tmp = potrcusp
        ttr = 0.5*(tmp[1]+tmp[0])+[-1,1]*(tmp[1]-tmp[0])*2.5
        ttr = ttr-(ttr mod dt) & ttr[1]+= dt
        tuts = smkarthm(ttr[0], ttr[1], dt, 'dx')
        get_data, 'po_ilat', t0, tmp
        tilat = interpol(tmp, t0, tuts)
        get_data, 'po_mlt', t0, tmp
        tmlt = interpol(tmp, t0, tuts)*15
        plots, tmlt[0:-2], tilat[0:-2], color = satcolors[j]
        plots, tmlt[0:-2], tilat[0:-2], color = satcolors[j], $
            psym = satsyms[j], symsize = 0.4
        arrow, tmlt[-2], tilat[-2], tmlt[-1], tilat[-1], /data, $
            color = satcolors[j], /solid, thick = thick
        plots, interpol(tmlt, tuts, potrcusp), interpol(tilat, tuts, potrcusp), $
            color = satcolors[j], thick = thick*5
        get_data, 'po_dis', t0, tmp
        tmp = string(interpol(tmp, t0, potrcusp[0]), format='(F3.1)')
        xyouts, 0.6, 0.25, /normal, 'Polar: *, R = '+tmp+' Re', color = satcolors[j]
        ; fast ilat, mlt.
        j = 1
        dt = 60        ; 5 min.
        tmp = fatrcusp
        ttr = 0.5*(tmp[1]+tmp[0])+[-1,1]*(tmp[1]-tmp[0])*2.5
        ttr = ttr-(ttr mod dt) & ttr[1]+= dt
        tuts = smkarthm(ttr[0], ttr[1], dt, 'dx')
        get_data, 'fa_ilat', t0, tmp
        tilat = interpol(tmp, t0, tuts)
        get_data, 'fa_mlt', t0, tmp
        tmlt = interpol(tmp, t0, tuts)*15
        plots, tmlt[0:-2], tilat[0:-2], color = satcolors[j]
        plots, tmlt[0:-2], tilat[0:-2], color = satcolors[j], $
            psym = satsyms[j], symsize = 0.4
        arrow, tmlt[-2], tilat[-2], tmlt[-1], tilat[-1], /data, $
            color = satcolors[j], /solid, thick = thick
        plots, interpol(tmlt, tuts, fatrcusp), interpol(tilat, tuts, fatrcusp), $
            color = satcolors[j], thick = thick*5
        tmp = string((fatrcusp[0]-potrcusp[0])/3600d, format='(F4.1)')
        tmp = strtrim(tmp,2)
        xyouts, 0.6, 0.2, /normal, 'FAST: +, dT = '+tmp+' hr', color = satcolors[j]
        sgpsclose, /pdf
    endfor
end

locroot = shomedir()+'/Google Drive/works/works/cusp/cusp list conjun'
cusp_add_fpoint_fig, locroot
end