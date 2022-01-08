
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
times = '!9'+string(180b)+'!X'
fac = ['para','west','north']
rgb = sgcolor(['red','green','blue'])



tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'num_lab_min', 5
tplot_options, 'xgridstyle', 1
tplot_options, 'xticklen', 1
tplot_options, 'version', 1
tplot_options, 'ygridstyle', 1
tplot_options, 'yticklen', 1



; rbspx_[vb1,eb1,mb1]. original E/B field in UVW.
; rbspx_[uvw]_gsm. UVW direction in GSM.
if tnames(pre0+'vb1') eq '' then _2013_0607_0456_load_burst_data



; load data.
tvar = pre0+'eb1'
get_data, tvar, uts, dat
tdat = dat[*,0]
nrec = n_elements(uts)


; **** morlet wavelet.
w0 = 6
s2p = 4*!dpi/(w0+sqrt(2+w0^2))

; coi.
tmid = 0.5*(uts[0]+uts[nrec-1])
coi = (tmid-uts[0]-abs(uts-tmid))/sqrt(2)

; prepare scales.
s0 = 16*dr0
dj = 0.125d
nj = 40


; spectrogram.
mor = wv_cwt(tdat, 'Morlet', w0, /pad, $
    start_scale = s0/dr0, dscale = dj, nscale = nj, scale = recscls)
timescls = recscls*dr0
periods = timescls*s2p


gen_plot = 0
if gen_plot then begin
    ofn = 0
    ;ofn = shomedir()+'/lh_test_morlet.pdf'
    sgopen, ofn, xsize = 6, ysize = 4, /inch
    
    poss = sgcalcpos(1)
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    
    zz = alog(abs(mor^2))
    xx = uts-uts[0]
    yy = periods
    tpos = poss[*,0]
    
    
    zrange = alog([1e-3,0.1])
    ztickv = [-1d31,smkarthm(zrange[0],zrange[1],10,'n')]
    zticks = n_elements(ztickv)
    colors = smkarthm(10,250,zticks,'n')
    tcolor = sgcolor(colors, ct = 43, file = 'ct2')
    
    xrange = minmax(xx)
    xtitle = 'Time (sec)'
    xticks = 6
    xminor = 4
    
    yrange = minmax(yy)
    ytitle = 'Frequency (Hz)'
    ytickn = ['12.5','25','50','100','200']
    ytickv = 1d/double(ytickn)
    yticks = n_elements(ytickv)-1
    ygrids = 1d/[20d,50,75]
    
    yfilts = 1d/[10,35,60,150]
    
    
    
    contour, zz, xx, yy, position = tpos, /noerase, $
        /fill, levels = ztickv, c_color = tcolor, $
        xstyle = 9, xrange = xrange, xlog = 0, xtitle = xtitle, $
        xticks = xticks, xminor = xminor, xticklen = -0.01, $
        ystyle = 9, yrange = yrange, ylog = 1, ytitle = ytitle, $
        ytickv = ytickv, yticks = yticks, ytickname = ytickn, yticklen = -0.01
        
    foreach ty, ygrids do begin
        plots, xrange, ty+[0,0], linestyle = 1, color = sgcolor('white')
        tmp = where(ytickv eq ty, cnt)
        if cnt ne 0 then continue
        tmp = convert_coord(xrange[0], ty, /data, /to_normal)
        xyouts, tmp[0]-xchsz*0.5, tmp[1]-ychsz*0.2, /normal, alignment = 1, $
            sgnum2str(1d/ty), charsize = 0.8
    endforeach
    foreach ty, yfilts do begin
        plots, xrange, ty+[0,0], linestyle = 0, color = sgcolor('red')
        tmp = where(ytickv eq ty, cnt)
        if cnt ne 0 then continue
        tmp = convert_coord(xrange[0], ty, /data, /to_normal)
        xyouts, tmp[0]-xchsz*0.5, tmp[1]-ychsz*0.2, /normal, alignment = 1, $
            sgnum2str(1d/ty), charsize = 0.8
    endforeach
    oplot, xx, coi, color = sgcolor('red')
    
    axis, yaxis = 1, ystyle = 1, yticklen = -0.01, $
        yrange = yrange*1e3, ylog = 1, ytitle = 'Period (msec)', $
        ytickv = ytickv*1e3, yticks = yticks
    axis, xaxis = 1, xstyle = 1, xticklen = 0, xtickformat = '(A1)', xticks = 1, xminor = 0
    
    tpos = [tpos[0],tpos[3]+ychsz*0.2,tpos[2],tpos[3]+ychsz*0.7]
    sgcolorbar, colors, zrange = zrange, position = tpos, /horizontal, $
        ztitle = 'Log!D10!N(mV/m)!U2!N', ct = 43, file = 'ct2'
        
    sgclose
endif




; reconstruction.
cdelta = 0.776d
gamma = 2.32
psi0 = !dpi^(-0.25)

yfilts = 1d/[10,35,60,150]      ; in period, determined by looking at spectrogram.
nfilt = n_elements(yfilts)-1

vars = pre0+['mb1','eb1']+'_fac'
foreach tvar, vars do begin     ; for each var.
    ytitle = (tvar eq pre0+'mb1_fac')? '(nT)': '(mV/m)'
    get_data, tvar, uts, dat
    
    fs = dblarr(nrec,3,nfilt)
    for i = 0, 2 do begin       ; for each component.
        tdat = dat[*,i]
        mor = wv_cwt(tdat, 'Morlet', w0, /pad, $
            start_scale = s0/dr0, dscale = dj, nscale = nj, scale = recscls)
        timescls = recscls*dr0
        periods = timescls*s2p
        
        for j = 0, nfilt-1 do begin ; for each band.
            filt = yfilts[j:j+1]
            filt = filt[sort(filt)]
            tys = dblarr(nrec)
            idx = where(periods ge filt[0] and periods le filt[1], cnt)
            for k = 0, cnt-1 do tys += real_part(mor[*,idx[k]]/sqrt(timescls[idx[k]]))
            tys *= (dj*sqrt(dr0))/cdelta/psi0
            fs[*,i,j] = tys
        endfor
    endfor
    
    for j = 0, nfilt-1 do begin
        id = string(j+1,format='(I0)')
        store_data, tvar+'_f'+id, uts, reform(fs[*,*,j]), limits = $
            {ytitle:ytitle, colors:rgb, labels:fac}
    endfor
endforeach





end
