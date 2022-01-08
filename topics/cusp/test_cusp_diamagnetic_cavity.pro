
pro test_cusp_pmag_ratio

ids = cusp_id('all')
nid = n_elements(ids)


; read info out of the log file.
rootdir = shomedir()+'/Google Drive/works'
if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')

fns = file_search(rootdir+'/data/cusp/po_sdt_fld_????_????_??.sdt', count=nfn)

; read info.
infos = []
nheader = 4
logfile = rootdir+'/works/cusp/cusp_list_of_conjun.log'
nlines = file_lines(logfile)-nheader
headers = strarr(nheader)
lines = strarr(nlines)
openr, lun, logfile, /get_lun
readf, lun, headers
readf, lun, lines
free_lun, lun
infos = [infos,lines]

logfile = rootdir+'/works/cusp/cusp_list_of_polar_2-4Re.log'
nlines = file_lines(logfile)-nheader
headers = strarr(nheader)
lines = strarr(nlines)
openr, lun, logfile, /get_lun
readf, lun, headers
readf, lun, lines
free_lun, lun
infos = [infos,lines]


dt = 30 ; sec.

uts = []
diss = []
mlts = []
ilats = []
btots = []
bmods = []
bt96s = []
dbp = []
dbv = []
dbb = []
xgsm = []
ygsm = []
zgsm = []

for i = 0, nfn-1 do begin
    fn = fns[i]
    id = strmid(file_basename(fn),11,12)
    
;    idx = where(stregex(infos,id) eq 0, cnt)
;    if cnt eq 0 then continue
;    
;    tinfo = infos[idx]
;    tr = strmid(id,0,9)+'/'+strsplit(strmid(tinfo,27,11),',',/extract)
;    tr = time_double(tr, tformat='YYYY_MMDD/hh:mm')

;    loginfo = cusp_read_conjun_list(logfile, event = id)
;    tr = loginfo.polar.cusp_time
;    if tr[1] lt tr[0] then tr[1]+= 86400d
;    tr = (tr[1]-tr[0])*0.2*[-1,1]+tr    ; expand a little.

    ; read data.
    if file_test(fn) eq 0 then message, 'file does not exist ...'
    sdt = ssdtread(fn)
    pre = 'po_'


    ;**** get uniform time.
    maxnrec = 100000ul   ; prevent memory overflow.
    t0 = sdt.var.polar_b_spc_z.depend_0
    dr = sdatarate(t0) & nrec = n_elements(t0)
    if nrec gt maxnrec then dr = (t0[nrec-1]-t0[0])/maxnrec
    tmp = minmax(sdt.var.polar_e_spc_z.depend_0)
    t0 = smkarthm(max([tmp[0],t0[0]]),min([tmp[1],t0[nrec-1]]), dr, 'dx')
    nrec = n_elements(t0)
    tstr = time_string(t0[0], tformat='YYYY_MMDD')
    if n_elements(eventid) ne 0 then tstr = eventid
    print, 'data rate: ', dr
    
    tr = minmax(t0)
    tr = tr+[1d,-1]*(tr[1]-tr[0])/4
    
   
    ;**** original b field and spike removal.
    ft  = sdt.var.polar_b_spc_z.depend_0
    fxy = sdt.var.polar_b_spc_x_y.value
    f56 = sdt.var.polar_b_spc_56.value
    fz  = sdt.var.polar_b_spc_z.value
    b_spc = [[temporary(fxy)],[temporary(fz)],[temporary(f56)]]
    b_spc = sinterpol(b_spc, ft, t0)
    sdespike, t0, b_spc, _extra = extra

    
    ;**** total b field. 'po_b'
    btotal = sqrt(total(b_spc^2,2))
    
    ;**** t96 model.
    ft  = sdt.var.polar_model_b_t96_spc_z.depend_0
    fxy = sdt.var.polar_model_b_t96_spc_x_y.value
    f56 = sdt.var.polar_model_b_t96_spc_56.value
    fz  = sdt.var.polar_model_b_t96_spc_z.value
    bt96_spc = [[temporary(fxy)],[temporary(fz)],[temporary(f56)]]
    bt96_spc = sinterpol(bt96_spc, ft, t0)
    bt96 = sqrt(total(bt96_spc^2,2,/nan))
    
    ;**** model b field. 'po_b0_spc'.
    db_spc = b_spc-bt96_spc
    bmod_spc = scalcbg(db_spc)+bt96_spc
    db_spc = b_spc-bmod_spc
    bmod = sqrt(total(bmod_spc^2,2))

    ;**** dummy de_spc.
    de_spc = dblarr(nrec,3)

    ;**** ilat, mlt, dis.
    t1  = sdt.var.polarinvariantlatitude.depend_0
    mlat = sdt.var.polarmaglat.value
    ilat = sdt.var.polarinvariantlatitude.value
    tmp = where(mlat lt 0)
    if tmp[0] ne -1 then ilat[tmp] *= -1
    mlt  = sdt.var.polarmlt.value
    dis  = sdt.var.polarspcraftdist.value


    ;**** rotate from SPC to FAC.
    polar_sdt_spc2fac, bmod_spc, db_spc, de_spc, db_fac, de_fac, ilat, lon, lat
    hem = (ilat[0] gt 0)? 1: -1
    db_spc*= hem


    ;**** load pos gsm, in Re.
    txgsm = sdt.var.polarpositiongsm_x.value
    tygsm = sdt.var.polarpositiongsm_y.value
    tzgsm = sdt.var.polarpositiongsm_z.value
    

    ; want bmod and btotal, bt96, mlt, ilat, dis.
    tuts = smkarthm(tr[0],tr[1],dt,'dx')
    tbmod = interpol(bmod, t0, tuts)
    tbt96 = interpol(bt96, t0, tuts)
    tbtot = interpol(btotal, t0, tuts)
    tilat = interpol(ilat, t1, tuts)
    tdis = interpol(dis, t1, tuts)
    tmlt = interpol(mlt, t1, tuts)
    tdbp = interpol(db_fac[*,1], t0, tuts)
    tdbv = interpol(db_fac[*,0], t0, tuts)
    tdbb = interpol(db_fac[*,2], t0, tuts)
    txgsm = interpol(txgsm, t1, tuts)
    tygsm = interpol(tygsm, t1, tuts)
    tzgsm = interpol(tzgsm, t1, tuts)

    uts = [uts,tuts]
    diss = [diss,tdis]
    mlts = [mlts,tmlt]
    ilats = [ilats,tilat]
    btots = [btots,tbtot]
    bmods = [bmods,tbmod]
    bt96s = [bt96s,tbt96]
    dbp = [dbp,tdbp]
    dbv = [dbv,tdbv]
    dbb = [dbb,tdbb]
    xgsm = [xgsm,txgsm]
    ygsm = [ygsm,tygsm]
    zgsm = [zgsm,tzgsm]
endfor

dat = {uts:uts,diss:diss,mlts:mlts,ilats:ilats,btots:btots,bmods:bmods,bt96s:bt96s,$
    dbp:dbp, dbv:dbv, dbb:dbb, xgsm:xgsm, ygsm:ygsm, zgsm:zgsm}
store_data, 'dmag', 0, dat

end



test_cusp_pmag_ratio
get_data, 'dmag', 0, dat

ofn = shomedir()+'/dmag.pdf'
;ofn = 0

sgopen, ofn, xsize = 5, ysize = 5, /inch

poss = sgcalcpos(2,1)
xr = [1,8]
symsize = 0.2
ct = 39

plot, dat.diss, (dat.dbp)/sqrt(dat.bmods), position = poss[*,0], /normal, /noerase, $
    ystyle = 1, ytitle = 'dB!Ceast-west!C(nT)', $
    xstyle = 1, xtickformat='(A1)', xrange = xr, $
    psym = 4, symsize = symsize
oplot, xr, [0,0], linestyle = 1

yr = [-50,50]
plot, dat.diss, (dat.dbb), position = poss[*,1], /normal, /noerase, $
    ystyle = 1, ytitle = 'dB!Cparallel!C(nT)', yrange = yr, $
    xstyle = 1, xtitle = 'Dist (Re)', xrange = xr, $
    psym = 4, symsize = symsize
oplot, xr, [0,0], linestyle = 1

sgclose


ofn = shomedir()+'/dmag_2d.pdf'
;ofn = 1

sgopen, ofn, xsize = 5, ysize = 5, /inch

sgtruecolor

tpos = [0.2,0.2,0.55,0.9]
zr = 10*[-1,1]

symx = [1,0,-1,0,1]
symy = [0,1,0,-1,0]

symx = cos(smkarthm(0,2*!dpi,30,'n'))
symy = sin(smkarthm(0,2*!dpi,30,'n'))

tx = dat.xgsm
ty = abs(dat.zgsm)
tz = bytscl(dat.dbb, min=zr[0], max=zr[1])

plot, tx, ty, /nodata, /noerase, /isotropic, position = tpos, $
    xrange = [0,4], xstyle = 1, xticks = 4, xminor = 5, xticklen = 0.01, xtitle = 'GSM x (Re)', $
    yrange = [0,8], ystyle = 1, yticks = 8, yminor = 5, yticklen = 0.04, ytitle = 'GSM y (Re)'
tmp = smkarthm(0,!dpi/2,100,'n')
oplot, cos(tmp), sin(tmp)
for i = 0, n_elements(tz)-1 do begin
    if tx[i] le 0 then continue
    usersym, symx, symy, color = sgcolor(tz[i], ct=ct), /fill
    plots, tx[i], ty[i], psym = 8, symsize = symsize
endfor

tpos = [0.6,0.2,0.62,0.9]
loadct, ct
sgcolorbar, indgen(256), position = tpos, zrange = zr, ztitle = 'dB!I||!N (nT)'

sgclose

end
