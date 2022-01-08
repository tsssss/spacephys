;+
; Calculate the phase diviation from linear in time.
;-

    ; 2012-10-06.
    fix_var = 'rbspa_e_fix'
    get_data, fix_var, times, efix
    phase = atan(efix[*,0], efix[*,1])
    
    spin_times = list()
    ntime = n_elements(times)
    for ii=0, ntime-2 do begin
        if phase[ii]-phase[ii+1] gt !dpi*1.8 then begin
            spin_times.add, times[ii]
        endif
    endfor
    spin_times = spin_times.toarray()

;---To get a circle.    
    ; 2012-10-06/00:02:00 2012-10-06/00:02:11.
    ; plot, shift(efix[index,0],12), efix[index,1]*0.8, psym=1, /iso
    ; 2012-10-06/00:02:11 2012-10-06/00:02:22
    ; plot, shift(efix[index,0], 0), efix[index,1]*0.8, psym=1, /iso
    ; 2012-10-06/00:02:22 2012-10-06/00:02:33
    ; plot, shift(efix[index,0],-5), efix[index,1]*0.86, psym=1, /iso
    
    
    nspin_time = n_elements(spin_times)
    for ii=0, nspin_time-2 do begin
        index = lazy_where(times, spin_times[ii:ii+1], count=count)
        xx = times[index]
        yy = phase[index]
        for jj=1, count-1 do if yy[jj-1] ge yy[jj] then yy[jj:*] += 2*!dpi
        xx = xx-xx[0]
        yy = yy-yy[0]
        fit = linfit(xx,yy, yfit=yfit)
        dy = yy-yfit
        stop
    endfor
    
    xx = times-times[0]
    yy = phase
    fit_res = linfit(xx, yy)
    stop

end