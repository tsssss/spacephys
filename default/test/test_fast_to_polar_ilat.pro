
fn = sdiskdir('Works')+'/confs/seminar/seminar_2014_1021/data/seminar_dat.tplot'
tplot_restore, filename = fn
ofn = shomedir()+'/agu_2014_fig_eflux.eps'
eventid = '1998_1001_02'
logfile = sdiskdir('Works')+'/works/cusp/cusp_list_of_conjun.log'
info = cusp_read_conjun_list(logfile, event = eventid)
potr = info.polar.plot_time
fatr = info.fast.plot_time
potrcusp = info.polar.cusp_time
fatrcusp = info.fast.cusp_time

get_data, 'po_ilat', pot0, poilat
get_data, 'fa_ilat', fat0, failat
idx = where(fat0 ge fatr[0] and fat0 le fatr[1])
failat = failat[idx] & fat0 = fat0[idx]

podum = findgen(n_elements(poilat))
fadum = interpol(podum,poilat,failat)

nmajor = 6
tickv = findgen(nmajor)/(nmajor-1)*n_elements(poilat)
tickv = interpol(failat,fadum,tickv)
plot, fadum, fat0

; get faeflux.
get_data, 'fa_ion_keflux_map', tmp0, faikeflux0
tfadum = interpol(fadum,fat0,tmp0)
faikeflux1 = interpol(faikeflux0,tfadum,findgen(n_elements(tmp0)))
store_data, 'fa_ion_keflux_map_ilat', tmp0, faikeflux1
plot, tfadum, faikeflux0, xticks = nmajor-1, xminor = 10, xtickname = snum2str(tickv), xrange = minmax(podum), xstyle = 1
stop

vf = [1,1.2,3.5,2.3]  ; var x, Fast data.
latf = [59,64,70,76]    ; Fast lat.
vp = [6,5,7.5,6.3]  ; var x, Polar data.
latp = [60,62,70,78]    ; Polar lat.
yr = [0,9]
xr = [55,80]

window, 0
plot, xr, yr, /nodata, xstyle = 1, ystyle = 1
plot, latf, vf, xrange = xr, yrange = yr, /noerase, xstyle = 5, ystyle = 5
plot, latp, vp, xrange = xr, yrange = yr, /noerase, xstyle = 5, ystyle = 5

dump = [1d,2,3,4]          ; Polar virtual x-axis coord.
dumf = interpol(dump,latp,latf)
xr = [0,4.5]
tickv = [0d,1,2,3,4,5]
tickv = interpol(latf,dumf,tickv)

window, 1
plot, xr, yr, /nodata, xstyle = 1, ystyle = 1, xtickname = snum2str(tickv)
plot, dumf, vf, xrange = xr, yrange = yr, /noerase, xstyle = 5, ystyle = 5
plot, dump, vp, xrange = xr, yrange = yr, /noerase, xstyle = 5, ystyle = 5

dump = [50d,60,70,80]      ; dummy Polar.
dumf = interpol(dump,latp,latf)
xr = [40,90]
tickv = [40d,50,60,70,80,90]
tickv = interpol(latf,dumf,tickv)

window, 2
plot, xr, yr, /nodata, xstyle = 1, ystyle = 1, xtickname = snum2str(tickv)
plot, dumf, vf, xrange = xr, yrange = yr, /noerase, xstyle = 5, ystyle = 5
plot, dump, vp, xrange = xr, yrange = yr, /noerase, xstyle = 5, ystyle = 5

end