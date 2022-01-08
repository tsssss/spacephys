

ilogfn = shomedir()+'/Google Drive/works/works/rbsp_de/'+$
    'dipolarization/list_large_de_round3.log'

ofn = shomedir()+'/Google Drive/works/data/rbsp_de/'+$
        'event_info.tplot'
ofigdir = shomedir()+'/psbl_de/'


tmp = {ut:0d, bipolar:-1d, data:dblarr(3)+!values.d_nan}
tinfo = {id:'', $           ; event id YYYY_MMDD_hhmm.
    tsta:0d, tend:0d, $     ; star/end time of the event.
    dst:0d, ae:0d, $        ; min Dst, max AE.
    efac:tmp, edot0fac:tmp, $   ; max abs value (with sign), defined when ebipolar is 0.
    jumpinbeta:-1d, jumpindensity:-1d, $    ; 0/1/-1 for no/yes/no value.
    bzincrease:-1d, bincrease:-1, $         ; 0/1/-1 for no/yes/no value.
    maxpfb:-1d, pfbratio:-1d, pfasym:-1d, $ ; pflux down mapped, ratio of pflux b/total, asym b/w down/up.
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

tfilters = [20,60,250]
filterids = ['1','2']
nfilter = 2
coordlabs = ['b','p','v']
rgb = [6,4,2]
dt = 20*60
sclinfo = [10,dt/2,40]

maxpfs = dblarr(ninfo,nfilter+1,3)
maxdot0pfs = dblarr(ninfo,nfilter+1,3)

for i0 = 0, ninfo-1 do begin
    
    tline = lines[i0]
    id = strmid(tline,0,16)
    tprobe = strmid(id,strlen(id)-1)
    pre0 = 'rbsp'+tprobe+'_'
    
    print, 'reading event: '+id+' ...'
    
    ; set id, tsta, tend.
    infos[i0].id = id
    t1 = time_double(strmid(id,0,9)+strmid(tline,20,5),tformat='YYYY_MMDDhh:mm')
    t2 = time_double(strmid(id,0,9)+strmid(tline,29,5),tformat='YYYY_MMDDhh:mm')
    if t2 lt t1 then t2+= 86400d
    infos[i0].tsta = t1
    infos[i0].tend = t2
    utr = [t1,t2]                   ; smaller time range for large field.
    utr = mean(utr)+[-0.5,0.5]*dt   ; larger time range for pflux calc.
    
    
    store_data, '*', /delete
    psbl_de_load_data, id, trange = utr
    
    denames = pre0+['de_fac','de_dot0_fac']
    dbnames = pre0+['db_fac','db_fac']
    pfnames = pre0+['pf_fac','pf_dot0_fac']
    suffixs = '_mat'+filterids
    
    
    ; statistics.
    ; max pflux in 3-d.
    maxpf = dblarr(nfilter+1,3)
    vars = pfnames[0]+['_mat',suffixs]+'_map'
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], t0, dat
        idx = where(t0 ge t1 and t0 le t2)
        maxpf[i,0] = max(dat[idx,0], maxidx, /absolute)
        maxpf[i,1] = dat[idx[maxidx],1]
        maxpf[i,2] = dat[idx[maxidx],2]
    endfor
    maxpfs[i0,*,*] = maxpf
    
    get_data, vars[0], t0, dat
    idx = where(t0 ge t1 and t0 le t2)
    infos[i0].maxpfb = max(dat[idx,0], maxidx, /absolute)
    infos[i0].pfbratio = abs(maxpf[0,0])/snorm(dat[idx[maxidx],*])
    infos[i0].pfasym = total(dat[idx,0])/total(abs(dat[idx,0]))
    
    get_data, pre0+'fpt_mlat', t0, dat
    idx = where(t0 ge t1)
    hem = (dat[idx[0]] gt 0)? 1: -1
    infos[i0].maxpfb *= hem
    infos[i0].pfasym *= hem


    ; max dot0 pflux in 3-d.
    maxdot0pf = dblarr(nfilter+1,3)
    vars = pfnames[1]+['_mat',suffixs]+'_map'
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], t0, dat
        idx = where(t0 ge t1 and t0 le t2)
        maxdot0pf[i,0] = max(dat[idx,0], maxidx, /absolute)
        maxdot0pf[i,1] = (dat[idx,1])[maxidx]
        maxdot0pf[i,2] = (dat[idx,2])[maxidx]
    endfor
    maxdot0pfs[i0,*,*] = maxdot0pf
    
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
    
    
    ; save plots.
    if keyword_set(save_plot) then begin
        ymax = 0.90
        ymin = 0.10
        nvar = nfilter+1
        depos = sgcalcpos(nvar, position = [0.15,ymin,0.3,ymax])
        dbpos = sgcalcpos(nvar, position = [0.45,ymin,0.6,ymax])
        pfpos = sgcalcpos(nvar, position = [0.75,ymin,0.9,ymax])
        labs = pre0+['fpt_mlat','mlt']
        for i = 0, n_elements(denames)-1 do begin
            if i eq 0 then ofigfn = ofigdir+'/'+id+'_field_and_poynt_freq_band.pdf' $
            else ofigfn = ofigdir+'/'+id+'dot0_field_and_poynt_freq_band.pdf'
            
            tplot_options, 'labflag', -1
            tplot_options, 'num_lab_min', 4
            
            sgopen, ofigfn, xsize = 11, ysize = 8, /inch
            device, decomposed = 0
            loadct2, 43
            
            vars = denames[i]+[suffixs,'_mat']
            options, vars, 'ytitle', 'dE FAC!C(mV/m)'
            titl = 'RBSP-'+strupcase(tprobe)+' dE'
            tplot, vars, var_label = labs, trange = utr, position = depos, title = titl, /noerase, vlab_margin = 10
            timebar, [t1,t2], color = 6
            
            vars = dbnames[i]+[suffixs,'_mat']
            options, vars, 'ytitle', 'dB FAC!C(nT)'
            titl = 'RBSP-'+strupcase(tprobe)+' dB'
            tplot, vars, var_label = labs, trange = utr, position = dbpos, title = titl, /noerase, vlab_margin = 10
            timebar, [t1,t2], color = 6
            
            vars = pfnames[i]+[suffixs,'_mat']+'_map'
            options, vars, 'ytitle', 'S FAC!C(mW/m!U2!N)'
            options, vars, 'labels', coordlabs+' map'
            titl = 'RBSP-'+strupcase(tprobe)+' S map'
            tplot, vars, var_label = labs, trange = utr, position = pfpos, title = titl, /noerase, vlab_margin = 10
            timebar, [t1,t2], color = 6
            
            sgclose
        endfor
    endif
    
endfor

store_data, 'psbl_de_info', 0, infos
vars = 'psbl_de_info'
tplot_save, vars, filename = ofn

end