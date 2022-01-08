; plot the events along dst for 2012-2015.

utr = ['2012-09-25','2016-01-01']
;utr = ['2013-01-01','2016-01-01']
utr = time_double(utr)
probes = ['a','b']
nprobe = n_elements(probes)
symh0 = -40

load = 0

symhfn = shomedir()+'/Google Drive/works/data/rbsp_de/symh.tplot'
if file_test(symhfn) eq 1 then tplot_restore, filename = symhfn


if load then begin
    ; load dst.
;    dat = sread_omni(utr, vars=['Epoch','SYM_H'])
;    
;    tvar = 'symh'
;    uts = sfmepoch(dat.epoch,'unix')
;    store_data, tvar, uts, dat.sym_h, limits = {constant:[0,symh0], ytitle:'Sym/H!C(nT)'}
;    tplot_save, tvar, filename = symhfn
    
    ; load the times for all events.
    logfn = shomedir()+'/Google Drive/works/works/rbsp_de/dipolarization/list_large_de_round3.log'
    nheader = 3
    headers = strarr(nheader)
    nline = file_lines(logfn)-nheader
    lines = strarr(nline)
    openr, lun, logfn, /get_lun
    readf, lun, headers
    readf, lun, lines
    free_lun, lun
    
    utrs = dblarr(nline,2)
    prbs = strarr(nline)
    for i = 0, nline-1 do begin
    
        tfn = lines[i]
        id = strmid(tfn,0,16)
        tprobe = strmid(id,strlen(id)-1)
        
        t0 = time_double(strmid(id,0,9)+strmid(tfn,20,5),tformat='YYYY_MMDDhh:mm')
        t1 = time_double(strmid(id,0,9)+strmid(tfn,29,5),tformat='YYYY_MMDDhh:mm')
        if t1 lt t0 then t1+= 86400
        
        utrs[i,*] = [t0,t1]
        prbs[i] = tprobe
    endfor
    for i = 0, nprobe-1 do begin
        tprobe = probes[i]
        idx = where(prbs eq tprobe, cnt)
        if cnt eq 0 then continue
        store_data, 'rbsp'+tprobe+'_utrs', utrs[idx,*], bytarr(cnt)+1, $
            limits = {psym:1, symsize:0.2, yrange:[0,2], panel_size:0.5, $
            yticks:1, yminor:1, labels: 'RBSP-'+strupcase(tprobe), $
            ytickformat:'(A1)', ytitle:''}
    endfor
    
    ; get apogee mlt.
    uts = smkarthm(utr[0],utr[1],10*86400,'dx')
    nut = n_elements(uts)
    mlts = dblarr(nut,nprobe)
    for j = 0, nprobe-1 do begin
        tprobe = probes[j]
        pre0 = 'rbsp'+tprobe+'_'
        uts = smkarthm(utr[0],utr[1],nut,'n')
        
        for i = 0, nut-1 do begin
            tmp = sread_rbsp_efw_l3(uts[i], probes = tprobe, $
                vars = ['epoch','pos_gse','mlt_lshell_mlat'])
            if size(tmp,/type) ne 8 then begin
                uts[i] = !values.d_nan
                mlts[i,j] = !values.d_nan
            endif else begin
                dis = sqrt(total(tmp.pos_gse^2,2))
                mlt = tmp.mlt_lshell_mlat[*,0]
                idx = where(dis eq max(dis)) & idx = idx[0]
                uts[i] = sfmepoch(tmp.epoch[idx], 'unix', /epoch16)
                mlts[i,j] = mlt[idx]
            endelse
        endfor
    endfor
    ; [0,24] -> [-12,12]
    idx = where(mlts gt 12)
    mlts[idx] = mlts[idx]-24
    store_data, 'mlt', uts, mlts, limits = $
        {ytitle:'Apogee MLT!C(hr)', yrange:[-12,12], $
        yticks:4,yminor:6,constant:[-6,0,6], labels:'RBSP-'+['A','B'], $
        colors:sgcolor(['red','blue'])}
endif

tplot_options, 'ystyle', 1
tplot_options, 'labflag', -1
tplot_options, 'version', 2


get_data, 'mlt', uts, mlts
mlts = mlts[*,0]
nut = n_elements(uts)
mltdelims = [0,-6,6,12]
mltdelims = [0,-6,6]
utdelims = []
for j = 0, n_elements(mltdelims)-1 do begin
    if mltdelims[j] eq 12 then begin
        tmp1 = mlts[1:*]*mlts[0:-2]
        tmp2 = mlts[1:*]-tmp[0:-2]
        idx = where(tmp1 le 0 and abs(tmp2) gt 12, cnt)
        for i = 0, cnt-1 do utdelims = [utdelims,mean(uts[idx[i]:idx[i]+1])]
        continue
    endif
    tmp = mlts-mltdelims[j]
    tmp1 = tmp[1:*]*tmp[0:-2]
    tmp2 = mlts[1:*]-tmp[0:-2]
    idx = where(tmp1 le 0 and abs(tmp2) le 12, cnt)
    if cnt ne 0 then utdelims = [utdelims,uts[idx]]
endfor

utdelims = [utr, utdelims]
utdelims = utdelims[sort(utdelims)]
nutdelim = n_elements(utdelims)
nsec = nutdelim-1
secuts = dblarr(nsec)
seccnts = dblarr(nsec)
secnormcnts = dblarr(nsec)
secdts = dblarr(nsec)
get_data, 'symh', uts, symh
dr = sdatarate(uts)
for i = 0, nutdelim-2 do begin
    idx = where(uts ge utdelims[i] and uts lt utdelims[i+1], cnt)
    if cnt eq 0 then continue
    tuts = uts[idx]
    tsymh = symh[idx]
    idx = where(tsymh le symh0, cnt)
    tuts = tuts[idx]
    
    symhcnt = 0
    for j = 0, cnt-2 do begin
        if tuts[j+1]-tuts[j] le 86400/2 then continue
        symhcnt+= 1
    endfor
    print, ''
    print, time_string(utdelims[i:i+1])
    print, symhcnt
    print, symhcnt/(utdelims[i+1]-utdelims[i])*1e6
    
    secuts[i] = (utdelims[i+1]+utdelims[i])*0.5
    secdts[i] = utdelims[i+1]-utdelims[i]
    seccnts[i] = symhcnt
    secnormcnts[i] = symhcnt/(utdelims[i+1]-utdelims[i])*30d*86400
endfor
store_data, 'storm_cnt', secuts, seccnts, limits = $
    {psym:1, yrange:[0,40], ytitle:'(#)', labels:'Storm Count', $
    yticks:2, yminor: 4, panel_size:0.5, constant: [20]}
store_data, 'storm_cnt_norm', secuts, secnormcnts, limits = $
    {psym:2, yrange:[0,6], ytitle:'', labels: 'Storm Count!C  per 30 days', $
    yticks:1, yminor:3, panel_size:0.5, constant: [2,4]}

options, 'rbspa_utrs', 'colors', sgcolor('red')
options, 'rbspb_utrs', 'colors', sgcolor('blue')
options, 'rbsp?_utrs', 'symsize', 0.5
options, 'storm_cnt_norm', 'psym', 1

;timebar, tmp, color = sgcolor('red')
;get_data, 'rbspb_utrs', tmp
;timebar, tmp, color = sgcolor('blue')


; **** plot for full time range.
;ofn = 0
;;ofn = shomedir()+'/fig_dst_distr.pdf'
;sgopen, ofn, xsize = 8, ysize = 5, /inch
;
;xchsz = double(!d.x_ch_size)/!d.x_size
;ychsz = double(!d.y_ch_size)/!d.y_size
;
;vars = ['symh','mlt','storm_cnt_norm','rbspa_utrs','rbspb_utrs']
;nvar = n_elements(vars)
;tplot, vars, trange = utr, get_plot_position = poss
;
;timebar, utdelims[where(utdelims ne utr[0] and utdelims ne utr[1])]
;
;tx = poss[0,nvar-1]-1*ychsz
;ty = (poss[1,nvar-1]+poss[3,nvar-2])*0.5
;xyouts, tx, ty, /normal, 'Large E Events', alignment = 0.5, orientation = 90
;
;sgclose


; **** plot for the pre-midnight sector.


; **** plot for the pre-midnight sector.
ofn = 0
ofn = shomedir()+'/fig_e_dst_correlation.pdf'
sgopen, ofn, xsize = 8, ysize = 8, /inch


xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size


tutrs = [[utdelims[1:2]],[utdelims[2:3]],[utdelims[4:5]],[utdelims[5:6]]]
tits = ['Post-midnight I', 'Pre-midnight I', 'Post-midnight II', 'Pre-midnight II']
ntutr = n_elements(tutrs)/2
poss = sgcalcpos(ntutr, ypad = 5)

tvar = 'mlt'
options, tvar, 'yrange', [-6,6]
options, tvar, 'yticks', 2
options, tvar, 'yminor', 3
options, tvar, 'ytitle', 'Apogee!CMLT!C(hr)'
tvar = 'symh'
options, tvar, 'yrange', [-100,50]
options, tvar, 'yticks', 2
options, tvar, 'yminor', 3
options, tvar, 'constant', 0

for j = 0, ntutr-1 do begin
    vars = ['symh','mlt']
    nvar = n_elements(vars)
    tpos = sgcalcpos(nvar, position = poss[*,j])
    tplot, vars, trange = tutrs[*,j], position = tpos, /noerase, title = tits[j]
    
    vars = ['rbspa','rbspb']+'_utrs'
    cs = sgcolor(['red','blue'])
    for i = 0, nvar-1 do begin
        get_data, vars[i], tmp
;        timebar, tmp, color = cs[i]
        idx = where(tmp ge tutrs[0,j] and tmp le tutrs[1,j])
        tmp = tmp[idx]
        tmp = (tmp-tutrs[0,j])/(tutrs[1,j]-tutrs[0,j])*(tpos[2,nvar-1]-tpos[0,nvar-1])+tpos[0,nvar-1]
        for k = 0, n_elements(tmp)-1 do plots, tmp[k]+[0,0], tpos[[1,3],nvar-1], /normal, color = cs[i]
    endfor
endfor


sgclose


end