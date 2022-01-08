;+
; plot v-i curve for RBSP-A and -B, one boom in one panel.
;-

pro plot_vi_curve, date, irange = irng, vrange = vrng, probes = probes, reload = reload, $
    usher0 = usher0, guard0 = guard0

    if n_elements(date) eq 0 then message, 'no date info ...'
    ut0 = time_double(date)

    if n_elements(probes) eq 0 then probes = ['a','b']
    nprobe = n_elements(probes)
    nboom = 6
    suffix = ['1','2','3','4','5','6']
    pairs = ['12','34','56']
    spinrate = 12   ; sec.
    maxni = 20      ; max # of steps in ibias.
    dibias = 5      ; nA, error in ibias.
    dusher = 1      ; V.
    dguard = 1      ; V.

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
                
        ofn = shomedir()+'/'+time_string(ut0,tformat='YYYY_MMDD_')+pre0+'iv_curve'+'.pdf'
;    ofn = i
        sgopen, ofn, xsize = 8, ysize = 10, /inch
        
        device, decomposed = 0
        loadct2, 43
        
        poss = sgcalcpos(nboom, lmargin=10, rmargin=15)
        ychsz = double(!d.y_ch_size)/!d.y_size
        erase
        
        for j = 0, nboom-1 do begin
            xtickfmt = (j eq nboom-1)? '': '(A1)'
            xtitl = (j eq nboom-1)? 'IBias (nA)': ''
            titl = (j eq 0)? 'RBSP-'+strupcase(probes[i])+' V-I curve, '+date: ''
            get_data, vars[j], tmp, infos
            
            iminmax = [-60,0]
            vminmax = [-10,0]
            
            if size(infos,/type) eq 8 then begin
                iminmax = minmax(infos.i)
                vminmax = minmax(infos.v)
                iminmax = ceil(iminmax/10)*10+[-10,0]
            endif
            
            plot, iminmax, vminmax, position = poss[*,j], /noerase, /nodata, $
                xtickformat = xtickfmt, ytitle = 'V'+suffix[j]+'!C(V)', $
                xtitle = xtitl, xrange = irng, yrange = vrng, title = titl, $
                xticklen = 1, yticklen = 1, xgridstyle = 1, ygridstyle = 1
            
            if size(infos,/type) ne 8 then continue
            
            ushers = round(infos.usher)
            ushers = ushers[uniq(ushers,sort(ushers))]
            if n_elements(u0) ne 0 then ushers = u0[j]
            nusher = n_elements(ushers)
            guards = round(infos.guard)
            guards = guards[uniq(guards,sort(guards))]
            if n_elements(g0) ne 0 then guards = g0[j]
            nguard = n_elements(guards)
            
            nband = nusher*nguard
            is = dblarr(maxni,nusher,nguard)+!values.d_nan
            vs = dblarr(maxni,nusher,nguard)+!values.d_nan
            for k = 0, nusher-1 do begin
                for l = 0, nguard-1 do begin
                    idx = where(round(infos.usher) eq ushers[k] and $
                        round(infos.guard) eq guards[l], cnt)
                    is[0:cnt-1,k,l] = infos.i[idx]
                    vs[0:cnt-1,k,l] = infos.v[idx]
                endfor
            endfor

            for m = 0, nband-1 do begin
                tmp = array_indices([nusher,nguard],m,/dimensions)
                k = tmp[0] & l = tmp[1]
                tx = is[*,k,l]
                ty = vs[*,k,l]
                
                if nband ne 1 then begin
                    tcolor = double(m)*254/nband
                    tnb = ((nband mod 2) eq 0)? nband: nband+1
                    tposx = poss[2,j] & if m ge (tnb/2) then tposx+= 0.1
                    tposy = poss[3,j]-ychsz*0.5-(poss[3,j]-poss[1,j])*(m mod (tnb/2))/(tnb/2)
                endif else begin
                    tcolor = 0
                    tposx = poss[2,j]
                    tposy = (poss[3,j]+poss[1,j])*0.5-ychsz*0.5
                endelse
                
                idx = where(finite(tx),cnt)
                if cnt ne 0 then $
                    oplot, tx[idx], ty[idx], psym = -1, symsize = 0.5, color = tcolor
                tmp = '    U:'+string(ushers[k],format='(I4)')+$
                    '    G:'+string(guards[l],format='(I4)')
                xyouts, tposx, tposy, $
                    tmp, /normal, color = tcolor, charsize = 0.8
            endfor
        endfor
        sgclose
    endfor
    
end

probes = ['a','b']

date = '2015-11-14'
usher0 = {a:[1e31,-15,1e31],b:[-15,-5,1e31]}
guard0 = {a:[1e31,-10,1e31],b:[-15,-10,1e31]}

date = '2015-11-27'
vrng = [-10,5]
irng = [-50,50]

;date = '2013-02-27'

date = '2014-12-17'
vrng = [-10,0]
irng = [-60,-10]

plot_vi_curve, date, vrange = vrng, irange = irng, probes = probes;, usher0 = usher0, guard0 = guard0;, /reload
end
