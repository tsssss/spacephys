
; settings.
utr0 = time_double('2013-06-07/04:56'+[':02.00',':03.00'])  ; total data.
utr1 = time_double('2013-06-07/04:56'+[':02.20',':02.80'])  ; zoom in for wave packet.


tprobe='a'
pre0 = 'rbsp'+tprobe+'_'
dr0 = 1d/4096
tpad = 200*dr0 ; sec.
p0 = 0.02   ; sec.
utr = utr0+[-1,1]*tpad
timespan, utr[0], utr[1]-utr[0], /second


re = 6378d & re1 = 1d/re
rgb = sgcolor(['red','green','blue'])
uvw = ['U','V','W']


tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'num_lab_min', 5
tplot_options, 'xgridstyle', 1
tplot_options, 'xticklen', 1
tplot_options, 'version', 1
tplot_options, 'ygridstyle', 1
tplot_options, 'yticklen', 1



_2013_0607_0456_load_burst_data



; **** mat settings.
order = 3
sclinfo = [2d-3,40d-3,40]
ytickvs = 1000d/[30,60,120,240]
yrng = [max(ytickvs),min(ytickvs)]
yticks = n_elements(ytickvs)-1
ytickns = strarr(yticks+1) & for i = 0, yticks do ytickns[i] = sgnum2str(1000d/ytickvs[i])
tscales = smkgmtrc(sclinfo[0],sclinfo[1],sclinfo[2],'n')
rscales = suniq(floor(tscales/dr0))
tscales = rscales*dr0


vars = pre0+['eb1_uvw','mb1_uvw']
ztits = ['(mV/m)','10!U-3!N (nT)']
zrngs = [6,6]
nvar = n_elements(vars)
for j = 0, nvar-1 do begin
    tvar = vars[j]
    get_data, tvar, uts, dat
    for i = 0, 2 do begin
        f0 = dat[*,i]
        f0mat = swvmat(f0, order, scale = rscales)
        if tvar eq pre0+'mb1_uvw' then f0mat *= 1e3
        store_data, tvar+'_mat'+string(i+1,format='(I0)'), uts, f0mat, tscales*1e3, $
            limits = {ytitle:'(Hz)', spec:1, no_interp:1, $
            ylog:1, ystyle:1, yrange:yrng, $
            yticks:yticks, yminor:5, ytickv:ytickvs, ytickname:ytickns, $
            ztitle:ztits[j], zrange:[-1,1]*zrngs[j], zticks:4}
    endfor
endfor


; **** plot spectrograms.
    pos1 = [0,0.5,1,1]
    pos2 = [0,0,1,0.5]
    
    ofn = 0
    sgopen, ofn, xsize = 8.5, ysize = 11, /inch
    device, decomposed = 0
    loadct2, 66
    
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size


    vars = pre0+'eb1_uvw_mat'+['1','2','3']
    figlabs = ['a','b','c']+'. dE '+['U','V','W']
    
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, region=pos1, lmargin=20, rmargin=15)
    tplot, vars, trange = utr1, position = poss, /noerase
    for i = 0, nvar-1 do xyouts, $
        poss[0,i]-xchsz*12, poss[3,i]-ychsz*1, /normal, alignment = 0, figlabs[i]



    vars = pre0+'mb1_uvw_mat'+['1','2','3']
    figlabs = ['d','e','f']+'. dB '+['U','V','W']
    
    nvar = n_elements(vars)
    poss = sgcalcpos(nvar, region=pos2, lmargin=20, rmargin=15)
    tplot, vars, trange = utr1, position = poss, /noerase
;    for i = 0, nvar-1 do begin
;        plot, [0,1], yrng, /noerase, position = poss[*,i], /nodata, $
;             xstyle = 5, ystyle = 5, ylog = 1
;        for j = 0, yticks-1 do plots, [0,1], ytickvs[j]+[0,0], linestyle = 1
;    endfor
    for i = 0, nvar-1 do xyouts, $
        poss[0,i]-xchsz*12, poss[3,i]-ychsz*1, /normal, alignment = 0, figlabs[i]

    sgclose




; **** filter the bands at 50 and 70 Hz.
filts = [6.04,14.02,22.02]
nfilt = n_elements(filts)-1
vars = pre0+['eb1','mb1']
for i = 0, nfilt-1 do begin
    fid = string(i+1,format='(I0)')
    foreach tvar, vars do begin
        for j = 0, 2 do begin
            id = string(j+1,format='(I0)')
            ttvar = tvar+'_uvw_mat'+id
            get_data, ttvar, uts, dat, val
            idx = where(val ge filts[i] and val lt filts[i+1])
            f0 = total(dat[*,idx],2)
            ; if tvar eq pre0+'mb1' then f0 *= 1e-3
            ttvar = tvar+'_f'+fid+'_'+id
            store_data, ttvar, uts, f0
        endfor
        tvars = tnames(tvar+'_f'+fid+'_?')
        ttvar = tvar+'_f'+fid
        lab = (tvar eq pre0+'eb1')? 'dE ': 'dB '
        ytit = (tvar eq pre0+'eb1')? '(mV/m)': '10!U3!N (nT)'
        stplot_merge, tvars, newname = ttvar, ytitle = ytit, colors = rgb, $
            labels = lab+uvw, /delete
    endforeach
endfor



end
