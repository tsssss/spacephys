
; log file contains basic info for the events.
ilogfn = shomedir()+'/Google Drive/works/works/rbsp_de/'+$
    'dipolarization/list_large_de_round3.log'

; key info for all events stored as tplot var.
infofn = shomedir()+'/Google Drive/works/data/rbsp_de/'+$
    'event_info.tplot'

; load key info
if file_test(infofn) ne 0 then tplot_restore, filename = infofn
get_data, 'psbl_de_info', tmp, infos

if size(infos,/type) ne 8 then begin
    psbl_de_update_info, ilogfn, infofn
    get_data, 'psbl_de_info', tmp, infos
endif

rad = !dpi/180


; pick out events with unipolar E.
idx = where(infos.efac.bipolar eq 0, nevent)
infos = infos[idx]


; E data, [3,n], [parallel,east-west,normal].
edat = infos.efac.data

ofn = 0
ofn = shomedir()+'/fig_edir.pdf'
sgopen, ofn, xsize = 5, ysize = 5, /inch

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size

tpos = [.2,.2,.9,.9]
emax = 80
thick = 8
xr = emax*[-1,1]
yr = emax*[-1,1]
plot, xr, yr, position = tpos, xstyle = 1, ystyle = 1, /nodata, $
    xtitle = 'dEy east-west (mV/m)', ytitle = 'dEz normal (mV/m)'
plots, xr, [0,0], linestyle = 2
plots, [0,0], yr, linestyle = 2
for i = 0, nevent-1 do begin
    plots, [0,edat[1,i]], [0,edat[2,i]], psym = 1, symsize = 0.8
endfor

idx = where(edat[2,*] gt 0)
tx = mean(edat[1,idx])
ty = mean(edat[2,idx])
arrow, 0, 0, tx, ty, /data, /solid, thick = thick, hsize = 0
xyouts, tpos[0]+2*xchsz, 0.5*(tpos[1]+tpos[3])+ychsz, /normal, $
    '<'+string(tx,format='(F4.1)')+', '+string(ty,format='(I3)')+' (mV/m)>'
    

idx = where(edat[2,*] le 0)
tx = mean(edat[1,idx])
ty = mean(edat[2,idx])
arrow, 0, 0, tx, ty, /data, /solid, thick = thick, hsize = 0
xyouts, tpos[0]+2*xchsz, 0.5*(tpos[1]+tpos[3])-1.8*ychsz, /normal, $
    '<'+string(tx,format='(F4.1)')+', '+string(ty,format='(I3)')+' (mV/m)>'

sgclose

end
