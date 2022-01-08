

datroot = shomedir()+'/Google Drive/works/data/rbsp_de'
ofn = shomedir()+'/Google Drive/works/data/rbsp_de/'+$
        '32hz_event_info.tplot'
ofigdir = shomedir()+'/psbl_de_32hz/'


deinfo = {ut:0d, $
    de_fac:[0d,0,0], $
    db_fac:[0d,0,0], $
    pf_fac:[0d,0,0], $
    fptmlat:0d, fptmlt:0d, mlt:0d, lshell:0d}
tinfo = {id:'', $           ; event id YYYY_MMDD_hhmm.
    utr:[0d,0], $           ; time range of the event.
    dst:0d, ae:0d, $        ; min Dst, max AE.
    maxde: deinfo, $        ; the max de info, de or dedot0 depend on usedot0
    usedot0: -1d, $         ; 1 for dedot0, 0 for de.
    dev:{min:0d, max:0d, minut:0d, maxut:0d, polratio:0d}, $
    fptmlat:0d, fptmlt:0d, mlt:0d, lshell:0d}  ; location info.


tfilters = [20,60,250]
filterids = ['1','2']
nfilter = 2
coordlabs = ['b','p','v']
rgb = [6,4,2]
dt = 20*60
sclinfo = [10,dt/2,40]


idinfos = psbl_de_id('detect_step')
nidinfo = n_elements(idinfos)

infos = replicate(tinfo,nidinfo)

for i = 0, nidinfo-1 do begin
    id = idinfos[i].id
    utr = idinfos[i].utr 
    tprobe = strmid(id,strlen(id)-1)
    pre0 = 'rbsp'+tprobe+'_'
    
;    if id ne '2013_0414_0702_a' then continue
    
    print, 'reading event: '+id+' ...'
    
    ; set id, tsta, tend.
    infos[i].id = id
    infos[i].utr = utr
    
    store_data, '*', /delete
    tplot_restore, filename = $
        datroot+'/derived_data/psbl_de_32hz_deriv_data_'+id+'.tplot'
    tplot_restore, filename = $
        datroot+'/32hz_data/psbl_de_32hz_data_'+id+'.tplot'
    

    ; collect general info the the event.
    ; dst and ae, extreme value in =/+ 5 min.
    tutr = mean(utr)+[-1,1]*300     ; -/+ 5 min.
    get_data, 'dst', uts, dat
    infos[i].dst = min(dat[where(uts ge tutr[0] and uts le tutr[1])])
    get_data, 'ae', uts, dat
    infos[i].ae = max(dat[where(uts ge tutr[0] and uts le tutr[1])])

    ; position at the center time.
    tut = mean(utr)
    get_data, pre0+'fpt_mlat', uts, dat
    infos[i].fptmlat = sinterpol(dat, uts, tut)
    get_data, pre0+'fpt_mlt', uts, dat
    infos[i].fptmlt = sinterpol(dat, uts, tut)
    get_data, pre0+'mlt', uts, dat
    infos[i].mlt = sinterpol(dat, uts, tut)
    get_data, pre0+'lshell', uts, dat
    infos[i].lshell = sinterpol(dat, uts, tut)
    
    
    ; collect info on dE polarization.
    b0ang0 = 20 ; 20 deg.
    get_data, pre0+'b0_angle', uts, dat
    idx = where(uts ge tutr[0] and uts le tutr[1])
    tmp = where(abs(dat[idx]) lt b0ang0, cnt)
    usedot0 = (cnt eq 0)? 1: 0
    infos[i].usedot0 = usedot0
    
    evar = (usedot0 eq 1)? 'dedot0': 'de'
    
    tvar = pre0+evar+'_survey_fac'
    tvar = pre0+evar+'_fac'
    get_data, tvar, uts, dat
    dr = sdatarate(uts)
    for j = 0, 2 do dat[*,j] = smooth(dat[*,j],11/dr)
    tutr = utr+[-1,1]*60
    idx = where(uts gt tutr[0] and uts le tutr[1])
    uts = uts[idx]
    dat = dat[idx,*]
    infos[i].dev.max = max(dat[*,2], tmp)
    infos[i].dev.maxut = uts[tmp]
    infos[i].dev.min= min(dat[*,2], tmp)
    infos[i].dev.minut = uts[tmp]
    polratio = abs(infos[i].dev.max/infos[i].dev.min)
    if polratio lt 1 then polratio = 1d/polratio
    infos[i].dev.polratio = polratio


    tdat = snorm(dat)
    maxde = max(tdat, idx)
    infos[i].maxde.de_fac = dat[idx,*]
    infos[i].maxde.ut = uts[idx]
    tut = uts[idx]
    get_data, pre0+'fpt_mlat', uts, dat
    infos[i].maxde.fptmlat = sinterpol(dat, uts, tut)
    get_data, pre0+'fpt_mlt', uts, dat
    infos[i].maxde.fptmlt = sinterpol(dat, uts, tut)
    get_data, pre0+'mlt', uts, dat
    infos[i].maxde.mlt = sinterpol(dat, uts, tut)
    get_data, pre0+'lshell', uts, dat
    infos[i].maxde.lshell = sinterpol(dat, uts, tut)

    get_data, pre0+'db_fac', uts, dat
    infos[i].maxde.db_fac = sinterpol(dat, uts, tut)
    
endfor

tvar = 'psbl_de_32hz_info'
store_data, tvar, 0, infos
tplot_save, tvar, filename = ofn

end
