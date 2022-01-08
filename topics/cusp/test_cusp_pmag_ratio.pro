
ids = cusp_id('all')
nid = n_elements(ids)


; read info out of the log file.
rootdir = shomedir()+'/Google Drive/works'
if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')

logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'


dt = 30 ; sec.

uts = []
diss = []
mlts = []
ilats = []
btots = []
bmods = []
bt96s = []

for i = 0, nid-1 do begin
    id = ids[i]
    fn = rootdir+'/data/cusp/po_sdt_fld_'+id+'.sdt'

    loginfo = cusp_read_conjun_list(logfile, event = id)
    tr = loginfo.polar.cusp_time
    if tr[1] lt tr[0] then tr[1]+= 86400d
    tr = (tr[1]-tr[0])*0.2*[-1,1]+tr    ; expand a little.

    ; read data.
    if file_test(fn) eq 0 then message, 'file does not exist ...'
    sdt = ssdtread(fn)
    pre = 'po_'
    
    ; some quantities and settings.
    rgb = [6,4,2]   ; r,g,b in colortable 43,45.
    labfac = ['v','p','b']   ; for north sunward meridian cross
    labspc  = ['xy','z','56']       ; v~xy, p(vxb)~56, b~z.
    defactor = 1.3d ; correction to E due to shielding.

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
    bmod = sqrt(total(bmod_spc^2,2))


    ;**** ilat, mlt, dis.
    t1  = sdt.var.polarinvariantlatitude.depend_0
    mlat = sdt.var.polarmaglat.value
    ilat = sdt.var.polarinvariantlatitude.value
    tmp = where(mlat lt 0)
    if tmp[0] ne -1 then ilat[tmp] *= -1
    mlt  = sdt.var.polarmlt.value
    dis  = sdt.var.polarspcraftdist.value



    ; want bmod and btotal, bt96, mlt, ilat, dis.
    tuts = smkarthm(tr[0],tr[1],dt,'dx')
    tbmod = interpol(bmod, t0, tuts)
    tbt96 = interpol(bt96, t0, tuts)
    tbtot = interpol(btotal, t0, tuts)
    tilat = interpol(ilat, t1, tuts)
    tdis = interpol(dis, t1, tuts)
    tmlt = interpol(mlt, t1, tuts)

    uts = [uts,tuts]
    diss = [diss,tdis]
    mlts = [mlts,tmlt]
    ilats = [ilats,tilat]
    btots = [btots,tbtot]
    bmods = [bmods,tbmod]
    bt96s = [bt96s,tbt96]
endfor

dat = {uts:uts,diss:diss,mlts:mlts,ilats:ilats,btots:btots,bmods:bmods,bt96s:bt96s}
store_data, 'pmag_ratio', 0, dat

end
