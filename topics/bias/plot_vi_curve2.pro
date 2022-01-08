;+
; plot v-i curve for RBSP-A and -B, opposite booms in one panel.
;-

pro plot_vi_curve2, date, irange = irng, vrange = vrng, probes = probes, reload = reload

    if n_elements(date) eq 0 then message, 'no date info ...'
    ut0 = time_double(date)

    if n_elements(probes) eq 0 then probes = ['a','b']
    nprobe = n_elements(probes)
    nboom = 6
    suffix = ['1','2','3','4','5','6']
    pairs = ['12','34','56']
    rootdir = shomedir()+'/bias'

    ; **** load data.
    ; rbsp[ab]_V[123456], rbsp[ab]_ibias[123456].
    ; rbsp[ab]_efw_[esvy].

    load = 0

    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'

        ; check these variables.
        ; load data if the variable doesn't exist, or at a diff time.
        vars = pre0+['V'+suffix,'IEFI_'+['IBIAS','GUARD','USHER'],'efw_esvy']
        for j = 0, n_elements(vars)-1 do begin
            get_data, vars[j], data = dat
            if n_elements(dat) eq 1 and size(dat,/type) ne 8 then begin
                load = 1            ; no data is stored in tplot.
                break
            endif else begin        ; a different date from what's stored.
                if (min(dat.x) gt ut0+86400d) or (max(dat.x) lt ut0) then load = 1
            endelse
        endfor
        vars = pre0+['IV'+suffix]
        for j = 0, n_elements(vars)-1 do begin
            get_data, vars[j], 0, dat
            if size(dat,/type) ne 8 then load = 1
        endfor

        if keyword_set(reload) then load = 1

        if load eq 1 then bias_param_load_data, date, probes = probes[i]
    endfor

    
    ; **** make plots.
    for i = 0, nprobe-1 do begin
        pre0 = 'rbsp'+probes[i]+'_'
        vars = pre0+'IV'+suffix
        if n_elements(usher0) ne 0 then begin
            u0 = usher0.(where(tag_names(usher0) eq strupcase(probes[i])))
            u0 = reform(u0 ## [1,1], nboom)
        endif
        if n_elements(guard0) ne 0 then begin
            g0 = guard0.(where(tag_names(guard0) eq strupcase(probes[i])))
            g0 = reform(g0 ## [1,1], nboom)
        endif

        ; get the usher and guard values.
        get_data, vars[0], tmp, infos
        ushers = round(infos.usher)
        ushers = ushers[uniq(ushers,sort(ushers))]
        if n_elements(u0) ne 0 then ushers = u0[j]
        nusher = n_elements(ushers)
        guards = round(infos.guard)
        guards = guards[uniq(guards,sort(guards))]
        if n_elements(g0) ne 0 then guards = g0[j]
        nguard = n_elements(guards)

        ; how many combinations.
        nband = nusher*nguard
        c1 = sgcolor('red')
        c2 = sgcolor('blue')

        for j = 0, nband-1 do begin
            tmp = array_indices([nusher,nguard],j,/dimensions)
            k = tmp[0] & l = tmp[1]
            tusher = ushers[k]
            tguard = guards[l]
            
            ofn = rootdir+'/bias_'+time_string(ut0,tformat='YYYY_MMDD')+ $
                '_rbsp'+probes[i]+'_iv_g'+string(tguard,format='(I0)')+ $
                '_u'+string(tusher,format='(I0)')+'.pdf'
;    ofn = 0
            

            sgopen, ofn, xsize = 5, ysize = 8, /inch

            poss = sgcalcpos(nboom/2)
            xchsz = double(!d.x_ch_size)/!d.x_size
            ychsz = double(!d.y_ch_size)/!d.y_size

            for m = 0, nboom/2-1 do begin
                tpos = poss[*,m]
                xtickfmt = (m eq nboom/2-1)? '': '(A1)'
                xtitl = (m eq nboom/2-1)? 'IBias (nA)': ''
                titl = (m eq 0)? 'RBSP-'+strupcase(probes[i])+', '+date+', Guard: '+ $
                    string(tguard,format='(I0)')+', Usher: '+ $
                    string(tusher,format='(I0)'): ''
                
                m1 = 2*m
                get_data, vars[m1], tmp, info
                idx = where(round(info.usher) eq tusher and $
                    round(info.guard) eq tguard, cnt)
                i1 = info.i[idx]
                v1 = info.v[idx]

                m2 = 2*m+1
                get_data, vars[m2], tmp, info
                idx = where(round(info.usher) eq tusher and $
                    round(info.guard) eq tguard, cnt)
                i2 = info.i[idx]
                v2 = info.v[idx]

                plot, [i1,i2], [v1,v2], /nodata, /noerase, position = tpos, $
                    title = titl, xtitle = xtitl, xtickformat = xtickfmt, $
                    xticklen = 1, yticklen = 1, xgridstyle = 1, ygridstyle = 1, $
                    ytitle = 'V'+pairs[m]+'!C(V)', $
                    yrange = vrng, ystyle = 1, xrange = irng, xstyle = 1
                oplot, i1, v1, psym = -1, color = c1
                oplot, i2, v2, psym = -1, color = c2
                xyouts, tpos[2]+2*xchsz, tpos[3]-1*ychsz, /normal, $
                    'Boom'+string(m1+1,format='(I0)'), color = c1
                xyouts, tpos[2]+2*xchsz, tpos[3]-2.2*ychsz, /normal, $
                    'Boom'+string(m2+1,format='(I0)'), color = c2
            endfor
            sgclose
        endfor
    endfor
    
end

date = '2015-11-14'
;date = '2013-02-27'
;date = '2014-12-17'
vrng = [-8,-1]
irng = [-120,-10]
probes = ['a','b']
plot_vi_curve2, date, vrange = vrng, irange = irng, probes = probes
end
