

ilogfn = shomedir()+'/Google Drive/works/works/rbsp_de/'+$
    'dipolarization/list_large_de_round3.log'

infofn = shomedir()+'/Google Drive/works/data/rbsp_de/'+$
    'event_info.tplot'
pfmax0 = 0
densitymax0 = 3
betamax0 = 10

if file_test(infofn) ne 0 then tplot_restore, filename = infofn
get_data, 'psbl_de_info', tmp, infos

if size(infos,/type) ne 8 then begin
    psbl_de_update_info, ilogfn, infofn
    get_data, 'psbl_de_info', tmp, infos
endif

ninfo = n_elements(infos)

mlats = dblarr(ninfo)
mlts = dblarr(ninfo)
pfmaxs = dblarr(ninfo)
pfratios = dblarr(ninfo)
pfasyms = dblarr(ninfo)
densityjumps = dblarr(ninfo)
betajumps = dblarr(ninfo)

for i = 0, ninfo-1 do begin
    mlats[i] = infos[i].fptmlat
    mlts[i] = infos[i].mlt
    pfmaxs[i] = infos[i].maxpfb
    pfratios[i] = infos[i].pfbratio
    pfasyms[i] = infos[i].pfasym
    densityjumps[i] = infos[i].jumpindensity
    betajumps[i] = infos[i].jumpinbeta
endfor
mlts[where(mlts gt 12)]-= 24


; the conditions.
idx = where(abs(mlts) le 4 and pfasyms gt 0.5 and pfmaxs gt pfmax0, cnt)
idx = where(abs(mlts) le 4 and pfratios gt 0.8 and pfasyms gt 0.5 and pfmaxs gt 30, cnt)
idx = where(abs(mlts) le 4 and pfratios gt 0.8 and densityjumps ge densitymax0 and pfmaxs gt pfmax0, cnt)


if cnt eq 0 then message, 'no event found ...'
mlats = mlats[idx]
mlts = mlts[idx]
infos = infos[idx]
ninfo = n_elements(infos)



types = ['e','e_dot0']
maxv = 3d/10   ; deg/mV/m.
pos = [0.15,0.15,0.9,0.9]
xr = [-6,6]*15
yr = [-1,1]*90


foreach type, types do begin
    
    udcnt = intarr(2,2,2)   ; mlt, mlat, n/s.
    
    ofn = shomedir()+'/fig_'+type+'_vec_pflux_gt_'+sgnum2str(pfmax0)+'.pdf'
;    ofn = 0
    sgopen, ofn, xsize = 5, ysize = 5, /inch
    
    ; set the box and labels.
    plot, xr/15, yr, xrange = xr/15, yrange = yr, /nodata, position = pos, $
        noerase = noerase, $
        xstyle = 1, ystyle = 1, xtitle = 'MLT (hr)', ytitle = 'Mapped MLat (deg)', $
        title = 'RBSP dE at PSBL, '+strupcase(type)+', Max S!D||!N>'+sgnum2str(pfmax0)+' (mW/m!U2!N)'
        
    ; set the data coord.
    plot, xr, yr, xrange = xr, yrange = yr, /nodata, /noerase, position = pos, $
        xstyle = 5, ystyle = 5
    plots, mean(xr)+[0,0], yr, linestyle = 1
    plots, xr, 65+[0,0], linestyle = 1
    plots, xr,-65+[0,0], linestyle = 1
    
    for i = 0, ninfo-1 do begin
        v0 = [mlts[i]*15,mlats[i]]
        if type eq 'e' then begin
            v1 = infos[i].efac.data[[1,2]]
        endif else begin
            v1 = infos[i].edot0fac.data[[1,2]]
        endelse
        udcnt[mlts[i] gt 0,mlats[i] gt 0,v1[1] gt 0]++
        v1 = v0+v1*maxv
        arrow, v0[0],v0[1], v1[0],v1[1], /data, color = sgcolor('blue'), hsize = 50, /solid
    endfor
    
    plots, mean(xr)+maxv*0.5*[-1,1]*50, mean(yr)+[0,0], color = sgcolor('blue')
    xyouts, mean(xr), mean(yr)-10, /data, '50 mV/m', alignment = 0.5
    
    ; lower left.
    xyouts,-3*15,-40, /data, alignment = 0.5, $
        sgnum2str(udcnt[0,0,0])+'D:'+sgnum2str(udcnt[0,0,1])+'U'
    ; lower right.
    xyouts, 3*15,-40, /data, alignment = 0.5, $
        sgnum2str(udcnt[1,0,0])+'D:'+sgnum2str(udcnt[1,0,1])+'U'
    ; upper left.
    xyouts,-3*15, 40, /data, alignment = 0.5, $
        sgnum2str(udcnt[0,1,0])+'D:'+sgnum2str(udcnt[0,1,1])+'U'
    ; upper right.
    xyouts, 3*15, 40, /data, alignment = 0.5, $
        sgnum2str(udcnt[1,1,0])+'D:'+sgnum2str(udcnt[1,1,1])+'U'
        
    sgclose

endforeach

end
