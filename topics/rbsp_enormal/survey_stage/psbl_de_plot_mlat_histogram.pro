


tvar = 'psbl_de_info'
if tnames(tvar) ne 0 then begin
   get_data, tvar, tmp, infos $
endif else begin    ; load infos from events.

    ilogfn = shomedir()+'/Google Drive/works/works/rbsp_de/'+$
        'dipolarization/list_large_de_round3.log'

    
    tmp = {ut:0d, bipolar:-1d, data:dblarr(3)+!values.d_nan}
    tinfo = {id:'', $           ; event id YYYY_MMDD_hhmm.
        tsta:0d, tend:0d, $     ; star/end time of the event.
        dst:0d, ae:0d, $        ; min Dst, max AE.
        efac:tmp, edot0fac:tmp, $   ; max abs value (with sign), defined when ebipolar is 0.
        jumpinbeta:-1d, jumpindensity:-1d, $    ; 0/1/-1 for no/yes/no value.
        bzincrease:-1d, bincrease:-1, $         ; 0/1/-1 for no/yes/no value.
        fptmlat:0d, fptmlt:0d, mlt:0d}          ; location info.
    
    nheader = 3
    headers = strarr(nheader)
    nline = file_lines(ilogfn)-nheader
    lines = strarr(nline)
    openr, lun, ilogfn, /get_lun
    readf, lun, headers
    readf, lun, lines
    free_lun, lun
    
    ninfo = nline
    infos = replicate(tinfo, ninfo)

    for i = 0, ninfo-1 do begin        
        tline = lines[i]
        id = strmid(tline,0,16)
        tprobe = strmid(id,strlen(id)-1)
        pre0 = 'rbsp'+tprobe+'_'
        
        print, 'reading event: '+id+' ...'
        
        ; set id, tsta, tend.
        infos[i].id = id
        t1 = time_double(strmid(id,0,9)+strmid(tline,20,5),tformat='YYYY_MMDDhh:mm')
        t2 = time_double(strmid(id,0,9)+strmid(tline,29,5),tformat='YYYY_MMDDhh:mm')
        if t2 lt t1 then t2+= 86400d
        infos[i].tsta = t1
        infos[i].tend = t2
        utr = [t1,t2]
        
        store_data, '*', /delete
        psbl_de_load_data, id, trange = utr
        
        ; set dst, ae.
        get_data, 'dst', t0, dat
        idx = where(t0 ge utr[0] and t0 le utr[1], cnt)
        infos[i].dst = (cnt eq 0)? !values.d_nan: min(dat[idx])
        get_data, 'ae', t0, dat
        idx = where(t0 ge utr[0] and t0 le utr[1], cnt)
        infos[i].ae = (cnt eq 0)? !values.d_nan: max(dat[idx])
        
        ; set ebipolar, efac[,ut], edot0fac[,ut].
        ratio0 = 1.5
        vars = pre0+['de_fac','de_dot0_fac']
        ttinfo = replicate({ut:0d, bipolar:1d, data:dblarr(3)+!values.d_nan},2)
        for j = 0, n_elements(vars)-1 do begin
            get_data, vars[j], t0, dat
            idx = where(t0 ge utr[0] and t0 le utr[1], cnt)
            tidx = (max(abs(dat[idx,1])) gt max(abs(dat[idx,2])))? 1: 2
            tidx = 2 ; normal component.
            maxde = max(dat[idx,tidx])
            minde = min(dat[idx,tidx])
            tmp = abs(maxde/minde) & if tmp lt 1 then tmp = 1d/tmp
            if tmp ge ratio0 then begin     ; uni-directional.
                thede = max([maxde,minde],/absolute)
                ttidx = where(dat[*,tidx] eq thede) ; ttidx marks thede.
                ttinfo[j].ut = t0[ttidx]
                ttinfo[j].bipolar = 0
                ttinfo[j].data = dat[ttidx,*]
            endif
        endfor
        infos[i].efac = ttinfo[0]
        infos[i].edot0fac = ttinfo[1]
        
        ; set mlt, fptmlt, fptmlat.
        deut = infos[i].efac.ut
        if deut eq 0 then begin
            infos[i].mlt = !values.d_nan
            infos[i].fptmlt = !values.d_nan
            infos[i].fptmlat = !values.d_nan
        endif else begin
            get_data, pre0+'mlt', data = dat
            infos[i].mlt = sinterpol(dat.y, dat.x, deut)
            get_data, pre0+'fpt_mlt', data = dat
            infos[i].fptmlt = sinterpol(dat.y, dat.x, deut)
            get_data, pre0+'fpt_mlat', data = dat
            infos[i].fptmlat = sinterpol(dat.y, dat.x, deut)
        endelse
        
        print, infos[i].mlt, infos[i].fptmlat
    endfor
    
    store_data, 'psbl_de_info', 0, infos
endif


ninfo = n_elements(infos)

devs = dblarr(ninfo)
deps = dblarr(ninfo)
mlats = dblarr(ninfo)
mlts = dblarr(ninfo)
mlts[where(mlts gt 12)]-= 24

for i = 0, ninfo-1 do begin
    devs[i] = infos[i].efac.data[2]
    deps[i] = infos[i].efac.data[1]
    mlats[i] = infos[i].fptmlat
    mlts[i] = infos[i].mlt
endfor

poss = sgcalcpos(3, ypad = 1)
binsz = 0.2

ofn = shomedir()+'/fig_mlat_histogram.pdf'
;ofn = 0
sgopen, ofn, xsize = 6, ysize = 5, /inch

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

ys = histogram(abs(mlats), binsize = binsz, min = 56, max = 66, locations = xs)
plot, xs, ys, psym = 10, position = poss[*,0], ytitle = 'Count', xtickformat = '(A1)'
xyouts, poss[0,0]+2*xchsz, poss[3,0]-1.5*ychsz, /normal, 'All Events'

ys = histogram(abs(mlats[where(mlats gt 0)]), binsize = binsz, min = 56, max = 66, locations = xs)
plot, xs, ys, psym = 10, position = poss[*,1], /noerase, ytitle = 'Count', xtickformat = '(A1)', yrange = [0,20]
xyouts, poss[0,1]+2*xchsz, poss[3,1]-1.5*ychsz, /normal, 'Northern Hemisphere'

ys = histogram(abs(mlats[where(mlats le 0)]), binsize = binsz, min = 56, max = 66, locations = xs)
plot, xs, ys, psym = 10, position = poss[*,2], /noerase, ytitle = 'Count', xtitle = 'Mapped MLat (deg)', yrange = [0,20]
xyouts, poss[0,2]+2*xchsz, poss[3,2]-1.5*ychsz, /normal, 'Southern Hemisphere'

sgclose

end
