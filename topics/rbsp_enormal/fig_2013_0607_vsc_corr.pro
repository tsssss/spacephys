
; plot vsc correlation.

_2013_0607_load_data


; settings.
autr = time_double(['2013-06-07/04:46','2013-06-07/05:05'])
butr = time_double(['2013-06-07/04:56','2013-06-07/05:00'])
;butr = time_double(['2013-06-07/04:54','2013-06-07/05:03'])
probes = ['a','b']

tvar = 'vsc'



device, decomposed = 0
loadct2, 43
tplot_options, 'version', 3
tplot_options, 'num_lab_min', 10
tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'yticklen', 1
tplot_options, 'ygridstyle', 1



; constants.
deg = 180d/!dpi
rad = !dpi/180
re = 6378d & re1 = 1d/re
r0 = 100d/re+1



get_data, 'rbspa_'+tvar, uts, avsc
get_data, 'rbspb_'+tvar, uts, bvsc
idx = where(uts ge butr[0] and uts le butr[1], tnrec)
tmp = bvsc[idx]
tmp = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
bvsc = tmp
buts = uts[idx]

;avsc = smooth(avsc,100,/edge_mirror)
;bvsc = smooth(bvsc,100,/edge_mirror)


; settings for time shift.
dr0 = sdatarate(uts)*10

dtrg = [-200,200]           ; the time range for time shift.
;dtrg = [-200,200]
;dtrg = [50,100]
dnrg = ceil(dtrg/dr0)       ; corresponding # of records.
dns = smkarthm(dnrg[0],dnrg[1],1,'dx')  ; each shift in # of record.
dts = dns*dr0
ndt = n_elements(dts)
corrs = fltarr(ndt)         ; the cross correlation at each shift.


;; this is another way of doing the cross correlation.
;dtrg = [-150,150]           ; the time range for time shift.
;get_data, 'rbspa_vsc', uts, avsc
;idx = where(uts ge butr[0]-150 and uts le butr[1]-150)
;avsc = avsc[idx]
;
;for i = 0, ndt-1 do begin
;    corrs[i] = c_correlate(avsc, bvsc, dns[i])
;endfor


for i = 0, ndt-1 do begin
    get_data, 'rbspa_'+tvar, uts, avsc
    tmp = avsc[where(uts ge butr[0]+dts[i] and uts le butr[1]+dts[i])]
    tmp = tmp-(tmp[0]+findgen(tnrec)/(tnrec-1)*(tmp[tnrec-1]-tmp[0]))
    avsc = tmp
    corrs[i] = c_correlate(avsc, bvsc, 0)
    corrs[i]*= stddev(avsc)
    
;    ofn = shomedir()+'/vsc_corr/dt_'+sgnum2str(dts[i])+'.pdf'
;    sgopen, ofn, xsize = 5, ysize = 3, /inch
;    device, decomposed = 0
;    loadct2, 43
;    plot, bvsc-bvsc[0], title = sgnum2str(corrs[i])
;    oplot, avsc-avsc[0], color = 6
;    sgclose
endfor

get_data, 'rbspa_'+tvar, uts, avsc
corrs *= 1d/stddev(bvsc)

plot, dts, corrs, xstyle = 1, ystyle = 1, $
    xrange = dtrg, yrange = [0,1], $
    xtitle = 'Time shift on RBSP-B (sec)', ytitle = 'Cross correlation'

maxcorr = max(corrs, idx)
maxcorrdt = dts[idx]
sigcorr = 2/sqrt(n_elements(bvsc))




ofn = shomedir()+'/fig_vsc_correlation.pdf'
;ofn = 0
sgopen, ofn, xsize = 6, ysize = 3, /inch

device, decomposed = 0
loadct2, 43

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

xr = [200,2000]
yr = [-50,-0]
tpos = [0.1,0.2, 0.5,0.9]
dvsc = -20  ; shift -B by -20 V.

tutr = butr+[-1.5,0.5]*300
get_data, 'rbspa_'+tvar, uts, avsc
idx = where(uts ge tutr[0] and uts le tutr[1])
txs = uts[idx]-uts[idx[0]]
tys = avsc[idx]

plot, txs, tys, xstyle = 1, ystyle = 9, yrange = yr, $
    ytitle = 'RBSP-A!CVsc (V)', position = tpos, $
    xticks = 3, xminor = 6, xticklen = 1, xgridstyle = 1, $
    xtitle = 'Second from '+time_string(uts[idx[0]])
axis, max(txs), /yaxis, yrange = yr-dvsc, ytitle = 'RBSP-B!CVsc (V)!Cshifted -20 V'

get_data, 'rbspb_'+tvar, uts, bvsc
idx = where(uts ge tutr[0]-maxcorrdt and uts le tutr[1]-maxcorrdt)
tys = bvsc[idx]+dvsc
oplot, txs, tys, color = 6

xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*0.5-ychsz*1.1*1, /normal, 'RBSP-A'
xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*0.5-ychsz*1.1*2, /normal, $
    'RBSP-B, dT = '+sgnum2str(maxcorrdt,nsgn=4)+' sec', color = 6


tpos = [0.72,0.2, 0.95,0.9]
tx = dts
ty = corrs
plot, tx, ty, position = tpos, /noerase, /ynozero, $
    ytitle = 'Cross Correlation', xtitle = 'Time shift on RBSP-B (sec)', yrange = [0,1]
    
xr = !y.crange
yr = !y.crange
plots, maxcorrdt+[0,0], yr, linestyle = 1
plots, !x.crange, sigcorr+[0,0], linestyle = 2

xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*0.5-ychsz*1.1*1, /normal, $
    'Max Correlation = '+sgnum2str(maxcorr,nsgn=2)
xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*0.5-ychsz*1.1*2, /normal, $
    'dT = '+sgnum2str(maxcorrdt,nsgn=4)+' sec'
xyouts, tpos[0]+xchsz*1, tpos[3]-ychsz*0.5-ychsz*1.1*3, /normal, $
    'Significant level = '+sgnum2str(sigcorr,nsgn=1)


sgclose




end
