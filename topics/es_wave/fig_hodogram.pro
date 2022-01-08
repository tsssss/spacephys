
; settings.
utr0 = time_double('2013-06-07/04:56'+[':02.00',':03.00'])  ; total data.
utr1 = time_double('2013-06-07/04:56'+[':02.20',':02.80'])  ; zoom in for wave packet.
utr2 = time_double('2013-06-07/04:56'+[':02.36',':02.41'])  ; the largest sub-wave packet.
utr2 = time_double('2013-06-07/04:56'+[':02.35',':02.40'])  ; the largest sub-wave packet.
utr2 = time_double('2013-06-07/04:56'+[':02.35',':02.40'])  ; the largest sub-wave packet.
utr3 = time_double('2013-06-07/04:56'+[':02.40',':02.47'])  ; the largest sub-wave packet.
;utr2 = time_double('2013-06-07/04:56'+[':02.37',':02.40'])  ; the largest sub-wave packet.


tprobe='a'
pre0 = 'rbsp'+tprobe+'_'
dr0 = 1d/4096
tpad = 200*dr0 ; sec.
p0 = 0.02   ; sec.
utr = utr0+[-1,1]*tpad
timespan, utr[0], utr[1]-utr[0], /second


re = 6378d & re1 = 1d/re
rgb = sgcolor(['red','green','blue'])
fac = ['para','west','north']
lmn = ['max','med','min']
wfs = [24d,48,74]   ; wave frequency.


tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'num_lab_min', 5
tplot_options, 'xgridstyle', 1
tplot_options, 'xticklen', 0
tplot_options, 'version', 1
tplot_options, 'ygridstyle', 0
tplot_options, 'yticklen', 0



_2013_0607_0456_load_burst_data


get_data, 'morlet_info', tmp, info
filters = info.filters
nband = n_elements(filters)-1
bandids = string(findgen(nband)+1,format='(I0)')
periods = info.periods


tinfo = {l:dblarr(3),m:dblarr(3),n:dblarr(3),eigenvals:dblarr(3)}
lmninfos = replicate(tinfo,nband)


ofn = shomedir()+'/es_wave/fig_hodogram_de.pdf'
;ofn = 0
sgopen, ofn, xsize = 3, ysize = 6, /inch

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

ys = [0.7,0.4,0.1]
xs = dblarr(3)+0.7

for i = 0, nband-1 do begin
    
    tutr = (i eq 0)? utr3: utr2
    utidx = where(uts ge tutr[0] and uts le tutr[1])
    tuts = uts[utidx]
    tnrec = n_elements(tuts)
       
    
    tvar = pre0+'eb1_fac_f'+bandids[i]
    get_data, tvar, uts, dat
    
    tuts = uts[utidx]
    defac = dat[utidx,*]
    
    eigvs = smva(defac, bvecs)
    delmn = defac # bvecs  ; [max,median,min].

    datlmn = delmn
    
    y0 = ys[i]
    y1 = ys[i]+0.22
    tx = xs[i]
    xpad = xchsz*0.5
    
    tmp = (y1-y0)
    tpos = [tx-xpad-tmp/0.5d,y0,tx-xpad,y1]
    
    
    case i of
        0: begin
            xr = [-1,1]*40
            yr = [-1,1]*40
            xtkv = [-1,1]*30
            ytkv = [-1,1]*30 & end
        1: begin
            xr = [-1,1]*40
            yr = [-1,1]*40
            xtkv = [-1,1]*30
            ytkv = [-1,1]*30 & end
        2: begin
            xr = [-1,1]*80
            yr = [-1,1]*80
            xtkv = [-1,1]*60
            ytkv = [-1,1]*60 & end
    endcase
    xr = [-1,1]*80
    yr = [-1,1]*80
    xtkv = [-1,1]*60
    ytkv = [-1,1]*60
    xtkfmt = (i eq 2)? '(I0)': '(A1)'
    xtitl = (i eq 2)? 'E!DM!N (nT)': ''
    ytitl = 'E!DL!N (nT)'
    
    plot, datlmn[*,1], datlmn[*,0], position = tpos, /noerase, $
        xstyle = 1, xrange = xr, xticks = 1, xminor = 6, xtickv = xtkv, $
        xtickformat = xtkfmt, xtitle = xtitl, $
        ystyle = 1, yrange = yr, yticks = 1, yminor = 6, ytitle = '', ytickv = ytkv
    plots, !x.crange, [0,0], linestyle = 1
    plots, [0,0], !y.crange, linestyle = 1
    plots, datlmn[0,1], datlmn[0,0], psym = 1, color = sgcolor('red')
    
    xyouts, tpos[0,0], y1+ychsz*0.5, /normal, alignment = 0, $
        'Wave '+bandids[i]+'  '+$
        time_string(tutr[0],tformat='hh:mm:ss.ff')+' to '+$
        time_string(tutr[1],tformat='ss.ff')+' UT'
    xyouts, tpos[0]-ychsz*5, (tpos[1]+tpos[3])*0.5, /normal, alignment = 0.5, $
        orientation = 90, ytitl
        
    xyouts, tpos[0]+xchsz, tpos[3]-ychsz*1.2, /normal, bandids[i]+'a.'
    
    
    
    tmp = (y1-y0)
    tpos = [tx+xpad,y0,tx+xpad+tmp/1d,y1]
    

    case i of
        0: begin
            xr = [-1,1]*15
            yr = [-1,1]*30
            xtkv = [-1,1]*15
            ytkv = [-1,1]*15 & end
        1: begin
            xr = [-1,1]*20
            yr = [-1,1]*40
            xtkv = [-1,1]*15
            ytkv = [-1,1]*30 & end
        2: begin
            xr = [-1,1]*40
            yr = [-1,1]*80
            xtkv = [-1,1]*30
            ytkv = [-1,1]*60 & end
    endcase
    xr = [-1,1]*40
    yr = [-1,1]*80
    xtkv = [-1,1]*30
    ytkv = [-1,1]*60
    xtitl = (i eq 2)? 'E!DN!N (nT)': ''
    
    plot, datlmn[*,2], datlmn[*,0], position = tpos, /noerase, $
        xstyle = 1, xrange = xr, xticks = 1, xminor = 3, xtickv = xtkv, $
        xtickformat = xtkfmt, xtitle = xtitl, $
        ystyle = 1, yrange = yr, yticks = 1, yminor = 6, ytickv = ytkv, $
        ytickformat = '(A1)'
    plots, !x.crange, [0,0], linestyle = 1
    plots, [0,0], !y.crange, linestyle = 1
    
    
    xyouts, tpos[0]+xchsz, tpos[3]-ychsz*1.2, /normal, bandids[i]+'b.'


endfor

sgclose

stop





ofn = shomedir()+'/es_wave/fig_hodogram_db.pdf'
; ofn = 0
sgopen, ofn, xsize = 3, ysize = 6, /inch

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

ys = [0.7,0.4,0.1]
xs = dblarr(3)+0.7

for i = 0, nband-1 do begin

    tutr = (i eq 0)? utr3: utr2
    utidx = where(uts ge tutr[0] and uts le tutr[1])
    tuts = uts[utidx]
    tnrec = n_elements(tuts)
    
    
    tvar = pre0+'mb1_fac_f'+bandids[i]
    get_data, tvar, uts, dat
    
    tuts = uts[utidx]
    dbfac = dat[utidx,*]
    
    eigvs = smva(dbfac, bvecs)
    dblmn = dbfac # bvecs  ; [max,median,min].
    lmninfos[i].l = bvecs[*,0]
    lmninfos[i].m = bvecs[*,1]
    lmninfos[i].n = bvecs[*,2]
    lmninfos[i].eigenvals = eigvs
    
    tvar = pre0+'mb1_lmn_f'+bandids[i]
    store_data, tvar, tuts, dblmn, limits = {ytitle:'(nT)',labels:lmn, colors:rgb}
    
    tvar = pre0+'eb1_fac_f'+bandids[i]
    get_data, tvar, uts, dat
    tuts = uts[utidx]
    defac = dat[utidx,*]
    delmn = defac # bvecs
    tvar = pre0+'eb1_lmn_f'+bandids[i]
    store_data, tvar, tuts, delmn, limits = {ytitle:'(mV/m)',labels:lmn, colors:rgb}
    
    
    datlmn = dblmn
    
    y0 = ys[i]
    y1 = ys[i]+0.22
    tx = xs[i]
    xpad = xchsz*0.5
    
    tmp = (y1-y0)
    tpos = [tx-xpad-tmp/0.5d,y0,tx-xpad,y1]
    
    xr = [-1,1]*0.08
    yr = [-1,1]*0.08
    xtkv = [-1,1]*0.06
    ytkv = [-1,1]*0.06
    xtkfmt = (i eq 2)? '(F5.2)': '(A1)'
    xtitl = (i eq 2)? 'B!DM!N (nT)': ''
    ytitl = 'B!DL!N (nT)'
    
    plot, datlmn[*,1], datlmn[*,0], position = tpos, /noerase, $
        xstyle = 1, xrange = xr, xticks = 1, xminor = 6, xtickv = xtkv, $
        xtickformat = xtkfmt, xtitle = xtitl, $
        ystyle = 1, yrange = yr, yticks = 1, yminor = 6, ytitle = '', ytickv = ytkv
    plots, !x.crange, [0,0], linestyle = 1
    plots, [0,0], !y.crange, linestyle = 1
    plots, datlmn[0,1], datlmn[0,0], psym = 1, color = sgcolor('red')
    
    xyouts, tpos[0,0], y1+ychsz*0.5, /normal, alignment = 0, $
        'Wave '+bandids[i]+'  '+$
        time_string(tutr[0],tformat='hh:mm:ss.ff')+' to '+$
        time_string(tutr[1],tformat='ss.ff')+' UT'
    xyouts, tpos[0]-ychsz*5, (tpos[1]+tpos[3])*0.5, /normal, alignment = 0.5, $
        orientation = 90, ytitl
        
    xyouts, tpos[0]+xchsz, tpos[3]-ychsz*1.2, /normal, bandids[i]+'a.'
    
    
    
    tmp = (y1-y0)
    tpos = [tx+xpad,y0,tx+xpad+tmp/1d,y1]
    
    xr = [-1,1]*0.04
    yr = [-1,1]*0.08
    xtkv = [-1,1]*0.03
    ytkv = [-1,1]*0.03
    xtitl = (i eq 2)? 'B!DN!N (nT)': ''
    
    plot, datlmn[*,2], datlmn[*,0], position = tpos, /noerase, $
        xstyle = 1, xrange = xr, xticks = 1, xminor = 3, xtickv = xtkv, $
        xtickformat = xtkfmt, xtitle = xtitl, $
        ystyle = 1, yrange = yr, yticks = 1, yminor = 6, ytickv = ytkv, $
        ytickformat = '(A1)'
    plots, !x.crange, [0,0], linestyle = 1
    plots, [0,0], !y.crange, linestyle = 1
    
    
    xyouts, tpos[0]+xchsz, tpos[3]-ychsz*1.2, /normal, bandids[i]+'b.'
    
    
endfor

sgclose




end
