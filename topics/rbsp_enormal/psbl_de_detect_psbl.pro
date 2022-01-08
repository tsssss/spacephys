
pro psbl_de_detect_psbl, tvar, probe = tprobe, trange = utr, $
    filename = logfn, id = id, ratio0 = ratio0, nsigma = nsigma, $
    position = pos0, steputr = steputr, infos = infos

    pre0 = 'rbsp'+tprobe+'_'
    
    steputr = []
    infos = []
    
    if n_elements(ratio0) eq 0 then ratio0 = 3
    if n_elements(nsigma) eq 0 then nsigma = 5
    if n_elements(pos0) eq 0 then pos0 = [0,0,1,1]

    get_data, tvar, uts, dat
    idx = where(finite(dat), cnt)
    if cnt eq 0 then begin
        ; use HOPE eletron spectrogram to calc density.
    endif else dat = sinterpol(dat[idx], uts[idx], uts)
    store_data, tvar, uts, dat
    
    idx = where(uts ge utr[0] and uts le utr[1], cnt)
    if cnt eq 0 then begin
        message, 'no data in given time range ...', /continue
        return
    endif
    nrec = cnt
    uts = uts[idx]
    dat = dat[idx]
    f0 = alog10(dat)
    df = [0,f0[1:nrec-1]-f0[0:nrec-2]]
    
    tdf = abs(df) & tdf = tdf[sort(tdf)] & tdf = tdf[0:0.8*(nrec-1)]
    sgm = stddev(abs(df[where(abs(df) le max(tdf))]))   ; stddev of abs(df).
    tsgm = nsigma*sgm
    idx = where(abs(df) gt tsgm, cnt)
    
    store_data, tvar+'_df', uts, abs(df), limits = $
        {ytitle:'Abs(df/dt)', constant:tsgm}

    device, decomposed = 0
    loadct2, 43
    
    tpos = sgcalcpos(2, region = pos0)
    
    if cnt eq 0 then return
    
    tplot, tvar+['','_df'], trange = utr, /noerase, position = tpos

    
    ; check all the points.
    for i = 0, cnt-1 do begin
        j = i
        while j+1 lt cnt-1 do $
            if idx[j+1]-idx[i] eq 1 and df[idx[j+1]]*df[idx[i]] gt 0 then $
                j+= 1 else break
        idx1 = idx[i]   ; start of the step.
        idx2 = idx[j]   ; end of the step.
        if idx1 eq idx2 then idx1-= 1
        i = j   ; i will be update in each for loop.
        ; find the mono-increase or decrease section.
        while idx1 gt 0 do $
            if df[idx1-1]*df[idx1] gt 0 then idx1-= 1 else break
        while idx2 lt nrec-1 do $
            if df[idx2+1]*df[idx2] gt 0 then idx2+= 1 else break
        ratio = dat[idx1]/dat[idx2]
        if ratio lt 1 then ratio = 1d/ratio
        print, 'jump ratio: ', ratio
        if ratio le ratio0 then continue
        tutr = uts[[idx1,idx2]]
        tdat = dat[idx1:idx2]

; this only works on density data.
;        if max(tdat) lt 0.01 then continue   ; all density < 0.01 cc.
;        if min(tdat) gt 100 then continue    ; all density > 100 cc.

        if n_elements(logfn) ne 0 then begin
            openw, loglun, logfn, /append, /get_lun
            printf, loglun, '** step in '+tvar+': '+ $
                time_string(tutr[0])+' - '+time_string(tutr[1])
            printf, loglun, '** jump ratio: ', ratio
            free_lun, loglun
        endif
        tinfo = {utr: tutr, var: tvar, jumpratio: ratio}
        infos = [infos,tinfo]
        timebar, tutr, color = 6
        steputr = [steputr,tutr]
    endfor
    
end



tprobe = 'a'


rootdir = shomedir()+'/psbl_de_detect_step/rbsp'+tprobe
ofn = rootdir+'/psbl_de_detect_step_rbsp'+tprobe+'.log'
if n_elements(ofn) ne 0 then if file_test(ofn) eq 0 then stouch, ofn


; read all the event info.
logfn = shomedir()+'/Google Drive/works/works/rbsp_de/'+$
    'psbl_de_round2/rbsp'+tprobe+'/list_rbsp'+tprobe+'_large_efield_32hz_quick.log'
nhead = 3   ; 3 lines of header.
heads = strarr(nhead)
nline = file_lines(logfn)-nhead
lines = strarr(nline)
openr, lun, logfn, /get_lun
readf, lun, heads
readf, lun, lines
free_lun, lun

infos = []

for i = 0, nline-1 do begin
    cols = strsplit(lines[i],' ',/extract)
    id = cols[0]+'_'+tprobe
    pre0 = 'rbsp'+tprobe+'_'

    
;    if id ne '2012_1114_0406_b' then continue
    
    openw, loglun, ofn, /append, /get_lun
    printf, loglun, ''
    printf, loglun, id
    free_lun, loglun
    
    date = strmid(cols[0], 0, 10)
    utr = time_double(date+cols[1:2], tformat='YYYY_MMDD_hh:mm:ss')
    utr = utr-(utr mod 60) & utr[1] += 60   ; floor and ceil to 1 min.
    if utr[0] gt utr[1] then utr[0] -= 86400d   ; a time in previous day.
    utr0 = utr-(utr mod 86400d)+[0,86400]
    
    ; load hope t, n.
    hopemom = sread_rbsp_hope_l3(utr0, probes = tprobe, type = 'mom')
    if size(hopemom,/type) ne 8 then continue
    
    uts = sfmepoch(hopemom.epoch_ele,'unix')
    store_data, pre0+'n', uts, hopemom.dens_e_200, $
        limits = {ytitle:'N!Ie!N!C(cm!E-3!N)', ylog:1, constant:1}
        
    store_data, pre0+'t', uts, hopemom.tperp_e_200, $
        limits = {ytitle:'T!Ie!I!N!C(eV)', ylog:1, $
        labels:'Tperp', constant:1000}

    
    ; reset time range.
    dt = 10*60  ; 10 min.
    if abs(utr[1]-utr[0]) le dt then utr = mean(utr)+[-1,1]*0.5*dt
    
    figfn = rootdir+'/psbl_de_detect_step_'+id+'.pdf'
    pos1 = [0.1,0.5,1,1]
    pos2 = [0.1,0,1,0.5]
    tutrs = []
    sgopen, figfn, xsize = 6, ysize = 8, /inch
    
    ; density.
    tvar = pre0+'n'
    psbl_de_detect_psbl, tvar, probe = tprobe, $
        trange = utr, id = id, filename = ofn, position = pos1, $
        steputr = tmp, infos = tinfo
    if n_elements(tmp) ne 0 then begin
        tutrs = [tutrs,tmp]
        infos = [infos,tinfo]
    endif
    
    ; temperature.
    tvar = pre0+'t'
    get_data, tvar, uts, dat
    idx = where(dat eq 1e20, cnt)
    if cnt ne 0 then dat[idx] = !values.d_nan
    store_data, tvar, uts, dat
    
    tvar = pre0+'t'
    psbl_de_detect_psbl, tvar, probe = tprobe, $
        trange = utr, id = id, filename = ofn, position = pos2, $
        steputr = tmp, infos = tinfo
    if n_elements(tmp) ne 0 then begin
        tutrs = [tutrs,tmp]
        infos = [infos,tinfo]
    endif

    sgclose
    if n_elements(tutrs) eq 0 then file_delete, figfn
endfor

rootdir = shomedir()+'/psbl_de_detect_step'
spc = '    '
save, infos, filename = rootdir+'/psbl_de_detect_psbl_info_'+tprobe+'.sav'
logfn = rootdir+'/psbl_de_detect_psbl_info_'+tprobe+'.log'
openw, lun, logfn, /get_lun
foreach tinfo, infos do $
    printf, lun, time_string(tinfo.utr[0])+spc+ $
        time_string(tinfo.utr[1])+spc+ $
        tinfo.var+spc+ $
        string(tinfo.jumpratio,format='(F5.1)')
free_lun, lun

end