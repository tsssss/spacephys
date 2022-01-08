
rootdir = shomedir()+'/psbl_de_32hz'
utr0 = time_double(['2012-09-25','2016-12-31'])
utabad = time_double('2015-10-01')
defsysv,'!rbsp_spice', exists=flag
if flag eq 0 then rbsp_load_spice_kernels, trange = utr0
logfn = rootdir+'/asi_conjunction_full_list.log'
listfn = rootdir+'/asi_conjunction_list.log'
stouch, logfn
stouch, listfn

dir = -1
par = 2
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
alt = 110   ; km.
r0 = 1+alt*re1

maxdr = sqrt(r0^2-1)    ; in Re.

ifn = rootdir+'/epacket/rbsp_e_ge_50.log'
nline = file_lines(ifn)
lines = strarr(nline)
openr, lun, ifn, /get_lun
readf, lun, lines
free_lun, lun

tinfo = {id:'', utr:dblarr(2), probe:'', maxde:0d}
infos = replicate(tinfo,nline)

for i = 0, nline-1 do begin
    tline = lines[i]
    date = strmid(tline,0,9)
    tutr = time_double(date+[strmid(tline,13,8),strmid(tline,25,8)], $
        tformat='YYYY_MMDDhh:mm:ss')
    if tutr[1] le tutr[0] then tutr[1]+= 86400d
    tprobe = strlowcase(strmid(tline,44,1))
    infos[i].id = time_string(tutr[0],tformat='YYYY_MMDD_hhmm_')+tprobe
    infos[i].utr = tutr
    infos[i].probe = tprobe
    infos[i].maxde = double(strmid(tline,37,3))
endfor



; load site position, mlon, mlat.
asc = sread_thg_asc(0)
sites = strlowcase(tag_names(asc))
nsite = n_elements(sites)
glats = dblarr(nsite)
glons = dblarr(nsite)
for i = 0, nsite-1 do begin
    glats[i] = asc.(i).glat
    glons[i] = asc.(i).glon
endfor
glats = glats*rad
glons = glons*rad


foreach tinfo, infos do begin
    id = tinfo.id
    utr = tinfo.utr    
    tprobe = strmid(id, 0,1, /reverse_offset)
    pre0 = 'rbsp'+tprobe+'_'
    
    if tprobe eq 'a' and time_double(utr[0]) ge utabad then continue
    
    ; the new commands to be printed.
    cmd1s = []
    
    print, 'processing '+id+' ...'
    
    ; get the basic results.
    sread_thg_check_data, utr, result = cmds
    nline = n_elements(cmds)
    
    
    ; pick out header line.
    header = cmds[0]
    cmd1s = [cmd1s,header]
    
    ; pick out time line.
    timeidx = where(stregex(cmds, '[0-9]{4}-[0-9]{2}-[0-9]{2}') ne -1, cnt)
    if cnt lt 1 then message, 'no time? something is wrong ...'
    timelines = cmds[timeidx]
    ntimeline = n_elements(timelines)
    
    ; check site location for each time.
    lineidx = [timeidx,nline]
    for j = 0, ntimeline-1 do begin
        cmd1s = [cmd1s, timelines[j]]
        sitelines = cmds[lineidx[j]+1:lineidx[j+1]-1]
        ; remove those site with no data.
        idx = where(stregex(sitelines, 'file exists') ne -1, cnt)
        if cnt eq 0 then begin
            cmd1s = [cmd1s,'no site has data ...']
            continue
        endif
        sitelines = sitelines[idx]
        idx = where(stregex(sitelines, 'no moon ...') ne -1, cnt)
        if cnt eq 0 then begin
            cmd1s = [cmd1s,'all sites have data have moon ...']
            continue
        endif
        sitelines = sitelines[idx]
        nsiteline = n_elements(sitelines)
        
        ; load s/c position and map to get tglat,tglon.
        defsysv,'!rbsp_spice', exists=flag
        if flag eq 0 then rbsp_load_spice_kernels, trange = utr, probes = tprobe
        tvar = pre0+'pos_gsm'
        tut = mean(utr)
        rbsp_load_spice_state, probe = tprobe, coord = 'gsm', times = tut, /no_spice_load
        
        get_data, pre0+'state_pos_gsm', tmp, posgsm
        posgsm = posgsm*re1
        tet = stoepoch(tut,'unix')
        geopack_epoch, tet, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
        geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001d, /date
        
        xp = posgsm[0] & yp = posgsm[1] & zp = posgsm[2]        
        geopack_trace, xp, yp, zp, dir, par, xf, yf, zf, r0 = r0, $
            /refine, /ionosphere, /t89
        fptgsm = [xf,yf,zf]
        
        ; convert from gsm to mag.
        geopack_conv_coord, xf,yf,zf, /from_gsm, txf,tyf,tzf, /to_geo
;        tglat = asin(tzf/r0)*deg
;        tglon = (atan(tyf,txf)*deg+360) mod 360
        
        ; check if any site in conjunction with the s/c.
        tflags = bytarr(nsiteline)   ; 0: no conjunction, 1: conjunction.
        tsites = strarr(nsiteline)
        tdrs = dblarr(nsiteline)
        for k = 0, nsiteline-1 do begin
            tsites[k] = strmid(sitelines[k],0,4)
            idx = where(sites eq tsites[k])
            tpos = [cos(glats[idx])*cos(glons[idx]),cos(glats[idx])*sin(glons[idx]), sin(glats[idx])]
            tdrs[k] = snorm(tpos-[txf,tyf,tzf])
            if tdrs[k] le maxdr then tflags[k] = 1
        endfor
        
        cmd1s = [cmd1s,strjoin(tsites,',')+' has data ...']

        idx = where(tflags eq 1, cnt)
        if cnt eq 0 then begin
            cmd1s = [cmd1s,'no site in conjunction with RBSP-'+strupcase(tprobe)+' ...']
            continue
        endif
        cmd1s = [cmd1s,strjoin(tsites[idx],',')+' in conjunction with RBSP-'+strupcase(tprobe)+' ...']
        tdrs = tdrs[idx]/r0*deg
        tmp = strarr(cnt)
        for k = 0, cnt-1 do tmp[k] = sgnum2str(tdrs[k],ndec=1)
        cmd1s = [cmd1s,strjoin(tmp,',')+' deg from site center ...']
    endfor
    
    openw, loglun, logfn, /get_lun, /append
    printf, loglun, ''
    printf, loglun, '**** processing '+id+' ...'
    printf, loglun, ''
    foreach tmp, cmd1s do printf, loglun, tmp
    free_lun, loglun
    
    ; only when conjunction was found.
    idx = where(stregex(cmd1s,' deg from site center') ne -1, cnt)
    if cnt eq 0 then continue
    openw, loglun, listfn, /get_lun, /append
    printf, loglun, ''
    printf, loglun, '**** processing '+id+' ...'
    printf, loglun, ''
    foreach tmp, cmd1s do printf, loglun, tmp
    free_lun, loglun
endforeach


end