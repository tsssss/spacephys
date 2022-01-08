; plot the events along dst for 2012-2015.

utr = ['2012-09-25','2016-01-01']
utr = time_double(utr)
probes = ['a','b']
nprobe = n_elements(probes)
symh0 = -40

load = 0

omnifn = shomedir()+'/Google Drive/works/data/rbsp_de/omni.tplot'
apomltfn = shomedir()+'/Google Drive/works/data/rbsp_de/apogee_mlt.tplot'

if file_test(omnifn) eq 1 then tplot_restore, filename = omnifn
if file_test(apomltfn) eq 1 then tplot_restore, filename = apomltfn


if load then begin
    
    ; load dst.
    dat = sread_omni(utr, vars=['Epoch','SYM_H','AE_INDEX'])

    uts = sfmepoch(dat.epoch,'unix')
    store_data, 'symh', uts, dat.sym_h, limits = {constant:[0,symh0], ytitle:'Sym/H!C(nT)'}
    store_data, 'ae', uts, dat.ae_index, limits = {ytitle:'AE!C(nT)'}
    tplot_save, ['symh','ae'], filename = omnifn

        
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
    store_data, 'rbspa_apogee_mlt', uts, mlts[*,0], limits = $
        {ytitle:'Apogee MLT!C(hr)', yrange:[-12,12], $
        yticks:4,yminor:6,constant:[-6,0,6], labels:'RBSP-A'}
    store_data, 'rbspb_apogee_mlt', uts, mlts[*,1], limits = $
        {ytitle:'Apogee MLT!C(hr)', yrange:[-12,12], $
        yticks:4,yminor:6,constant:[-6,0,6], labels:'RBSP-B'}
    
endif

tplot_options, 'ystyle', 1
tplot_options, 'labflag', -1
tplot_options, 'version', 2

nmltbin = 48
mltbins = smkarthm(-12,12,nmltbin+1,'n')
mltvals = 0.5*(mltbins[1:nmltbin]+mltbins[0:nmltbin-1])

dstlims = [-80,-40,-20,-10]
ndstlim = n_elements(dstlims)
aelims = [1000,700,400,100]
naelim = n_elements(aelims)


tinfo = {cmlt:0d, mltcnt:0d, $
    dstcnts:dblarr(ndstlim), dstmed:0d, dstavg:0d, $
    aecnts:dblarr(naelim), aemed:0d, aeavg:0d}
bininfos = replicate(tinfo, nmltbin)


get_data, 'ae', uts, aes
get_data, 'symh', uts, dsts
nrec = n_elements(uts)

get_data, 'rbspa_apogee_mlt', tmp, tmlts
tmlt1s = interpol(tmlts, tmp, uts)
tmlts[where(tmlts le 0)]+= 24
tmlt2s = interpol(tmlts, tmp, uts)
tmlts = tmlt1s
idx = where(tmlt1s-tmlt2s ne 0)
tmlts[idx] = tmlt2s[idx]-24
for i = 0, nrec-2 do begin      ; patch up the connecting part.
    if tmlt2s[i+1] lt tmlt2s[i] then continue
    tmlts[i] = tmlt1s[i]
endfor


for i = 0, nmltbin-1 do begin
    bininfos[i].cmlt = mltvals[i]
    idx = where(tmlts ge mltbins[i] and tmlts lt mltbins[i+1], cnt)
    bininfos[i].mltcnt = cnt
    if cnt ne 0 then begin
        tdat = dsts[idx]
        bininfos[i].dstmed = median(tdat)
        bininfos[i].dstavg = mean(tdat)
        for j = 0, ndstlim-1 do begin
            tmp = where(abs(tdat) ge abs(dstlims[j]), cnt)
            bininfos[i].dstcnts[j] = cnt
        endfor
        
        tdat = aes[idx]
        bininfos[i].aemed = median(tdat)
        bininfos[i].aeavg = mean(tdat)
        for j = 0, naelim-1 do begin
            tmp = where(abs(tdat) ge abs(aelims[j]), cnt)
            bininfos[i].aecnts[j] = cnt
        endfor
    endif
endfor


; prepare statistics.
premidx = where(bininfos.cmlt ge -6 and bininfos.cmlt le 0)
posmidx = where(bininfos.cmlt ge 0 and bininfos.cmlt le 6)

for i = 0, ndstlim-1 do begin
    tcnts = bininfos.dstcnts[i]
    print, dstlims[i], mean(tcnts[premidx])/mean(tcnts[posmidx]), $
        mean(tcnts[premidx]), mean(tcnts[posmidx])
endfor

for i = 0, naelim-1 do begin
    tcnts = bininfos.aecnts[i]
    print, aelims[i], mean(tcnts[premidx])/mean(tcnts[posmidx]), $
        mean(tcnts[premidx]), mean(tcnts[posmidx])
endfor


; **** plot the MLT histogram.
ofn = shomedir()+'/fig_mlt_histogram.pdf'
sgopen, ofn, xsize = 7, ysize = 7, /inch

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

poss = sgcalcpos(4)

tpos = poss[*,0]
plot, bininfos.cmlt, bininfos.mltcnt, psym = 10, xrange = [-6,6], /noerase, $
    yrange = [5e3,5e4], ylog = 0, xtickformat = '(A1)', $
    ytitle = 'Count (min)', position = tpos

labs = ['MLT']
cs = [0]
nlab = n_elements(labs)
txs = tpos[2]+dblarr(nlab)+2*xchsz
tys = tpos[3]-ychsz-findgen(nlab)*ychsz*1.2
for i = 0, nlab-1 do xyouts, /normal, txs[i], tys[i], labs[i], color = cs[i]




tpos = poss[*,1]
cidx = smkarthm(254,10,ndstlim, 'n')

plot, bininfos.cmlt, bininfos.dstcnts[ndstlim-1], /nodata, psym = 10, xrange = [-6,6], /noerase, $
    yrange = [30,3e4], ylog = 1, xtickformat = '(A1)', $
    ytitle = 'Count (min)', position = tpos
for i = 0, ndstlim-1 do oplot, bininfos.cmlt, bininfos.dstcnts[i], psym = 10, color = sgcolor(cidx[i], ct=40)
plot, bininfos.cmlt, bininfos.dstcnts[ndstlim-1], /nodata, psym = 10, xrange = [-6,6], /noerase, $
    yrange = [30,3e4], ylog = 1, xtickformat = '(A1)', $
    ytitle = 'Count (min)', position = tpos

nlab = ndstlim
txs = tpos[2]+dblarr(nlab)+2*xchsz
tys = reverse(tpos[3]-ychsz-findgen(nlab)*ychsz*1.2)
for i = 0, ndstlim-1 do xyouts, /normal, txs[i], tys[i], 'Dst<'+sgnum2str(dstlims[i]), color = sgcolor(cidx[i], ct=40)




tpos = poss[*,2]
cidx = smkarthm(254,10,naelim, 'n')
plot, bininfos.cmlt, bininfos.aecnts[naelim-1], /nodata, psym = 10, xrange = [-6,6], /noerase, $
    yrange = [30,3e4], ylog = 1, xtickformat = '(A1)', $
    ytitle = 'Count (min)', position = tpos
for i = 0, naelim-1 do oplot, bininfos.cmlt, bininfos.aecnts[i], psym = 10, color = sgcolor(cidx[i], ct=40)
plot, bininfos.cmlt, bininfos.aecnts[naelim-1], /nodata, psym = 10, xrange = [-6,6], /noerase, $
    yrange = [30,3e4], ylog = 1, xtickformat = '(A1)', $
    ytitle = 'Count (min)', position = tpos

nlab = naelim
txs = tpos[2]+dblarr(nlab)+2*xchsz
tys = reverse(tpos[3]-ychsz-findgen(nlab)*ychsz*1.2)
for i = 0, naelim-1 do xyouts, /normal, txs[i], tys[i], 'AE>'+sgnum2str(aelims[i]), color = sgcolor(cidx[i], ct=40)




; load event info.
infofn = shomedir()+'/Google Drive/works/data/rbsp_de/round1_dat_info.tplot'
tvar = 'psbl_de_info'
store_data, tvar, /delete
if file_test(infofn) eq 1 then tplot_restore, filename = infofn


if tnames(tvar) ne '' then begin
    get_data, tvar, tmp, infos
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
    tplot_save, tvar, filename = infofn
endelse


mlts = infos.mlt
mlts[where(mlts gt 12)]-= 24
ninfo = n_elements(mlts)
mltcnts = dblarr(nmltbin-1)

for i = 0, nmltbin-2 do begin
    idx = where(mlts gt mltbins[i] and mlts le mltbins[i+1], cnt)
    mltcnts[i] = cnt
endfor


tpos = poss[*,3]
plot, bininfos.cmlt, mltcnts, psym = 10, xrange = [-6,6], /noerase, $
    xtitle = 'MLT (hr)', $
    ytitle = 'Count #', position = tpos
labs = ['|E|>10 mV/m']
cs = [0]
nlab = n_elements(labs)
txs = tpos[2]+dblarr(nlab)+2*xchsz
tys = tpos[3]-ychsz-findgen(nlab)*ychsz*1.2
for i = 0, nlab-1 do xyouts, /normal, txs[i], tys[i], labs[i], color = cs[i]


tcnts = mltcnts
print, '|E|>10 mV/m', mean(tcnts[premidx])/mean(tcnts[posmidx]), $
    mean(tcnts[premidx]), mean(tcnts[posmidx])
    

sgclose



; **** the MLT-Mlat plot to be compared to Keiling's.
mlats = abs(infos.fptmlat)
mlatvs = [80,70,60]
minlat = 60
mltvs = [0,6,12,18]


rad = !dpi/180
n0 = 50     ; # of azimuthal points for the circles.
xr = [-1,1]
yr = [-1,1]
sdeg = '!9'+string(176b)+'!X'


; convert to polar coord, and normalize.
rs = (90-mlats)/(90d - 60)
rvs = (90-mlatvs)/(90d - 60)
ts = mlts*15
ts = ts*rad
tvs = mltvs*15
tvs = tvs*rad

ofn = shomedir()+'/fig_mlat_mlt_distr.pdf'
;ofn = 0
sgopen, ofn, xsize = 4, ysize = 4, /inch
sgtruecolor

tpos = [0.15,0.15,0.85,0.85]

plot, xr, yr, /polar, /nodata, $
    xrange = xr, yrange = yr, /isotropic, xstyle = 5, ystyle = 5, $
    position = tpos, symsize = 0.8

; add axis.
oplot, xr, [0,0]
oplot, [0,0], yr

foreach tmlat, rvs do $
    oplot, smkarthm(tmlat,tmlat,n0,'n'), smkarthm(0,2*!dpi,n0,'n'), /polar

td = 5e-2   ; nudge value.
tr = 1.1
xyouts, tr*cos(tvs), tr*sin(tvs)-td, string(mltvs,format='(I2)'), /data, alignment = 0.5
tt = -120*rad
xyouts, rvs*cos(tt)+1.5*td, rvs*sin(tt), string(mlatvs,format='(I2)')+sdeg, /data, alignment = 0.5

tmp = findgen(11)*2*!dpi/10
txs = cos(tmp)
tys = sin(tmp)
usersym, txs, tys, /fill
plots, rs*cos(ts), rs*sin(ts), psym = 8, /data, symsize = 0.4

sgclose


end
