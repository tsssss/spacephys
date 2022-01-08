;+
; Load the flux of pixel 1 and 5 to see the net flux along the w-axis.
;-

pro test_hope_velocity_direction_l2, time, probe=probe

;---Settings.
    the_species = 'o'
    species_name = 'O!U+!N'
    zrange = [10,1e5]

;---Load level 2 data for ion.
    foreach id, 'l2%'+['ion','misc'] do begin
        rbsp_read_hope, time, id=id, probe=probe, errmsg=errmsg
        if errmsg ne '' then return
    endforeach

    flux_var = 'f'+the_species+'du'
    get_data, flux_var, times, fdat
    fdat = fdat*1e-3    ; from 1/s-cm^2-sr-keV to 1/s-cm^2-sr-eV.
    nrec = n_elements(times)
    dims = size(reform(fdat[0,*,*,*]),/dimensions)

;---Pick out pixel 1 and 5:
;   Pixel 1 points to the positive w-axis.
;   Pixel 5 points to the negative w-axis.
;   Notice that particles flow into the pixels, therefore
;   the net flux is positive along the w-axis,
;   if the flux of pixel 5 is larger than pxiel 1.
    w_flux = reform(fdat[*,*,*,4]-fdat[*,*,*,0])
    w_flux = total(w_flux, 3)   ; To be simple, sum over all sectors.

;---Save data.
    pre0 = 'rbsp'+probe+'_'
    pre1 = pre0+the_species+'_'

    the_type = (the_species eq 'e')? 'ele': 'ion'
    en_var = 'hope_energy_'+the_type
    energy_bins = get_var_data(en_var)     ; in eV.

    var = pre1+'w_flux'
    store_data, var, times, w_flux, energy_bins
    add_setting, var, /smart, {$
        display_type: 'spec', $
        unit: '1/cm!U2!N-s-sr-keV', $
        zrange: zrange, $
        species_name: species_name, $
        ytitle: species_name+' Energy (eV)', $
        short_name: ''}

;---Plot data.
    vars = pre1+'w_flux'
    nvar = n_elements(vars)

    file = join_path([shomedir(),'test_hope_velocity_direction_l2_rbsp'+probe+'.pdf'])
    if keyword_set(test) then file=test
    sgopen, file, xsize=8.5, ysize=3

    tpos = sgcalcpos(nvar, lmargin=16, tmargin=5, rmargin=12, xchsz=xchsz, ychsz=ychsz)

    device, decomposed=0
    loadct2, 43

    tplot, vars, position=tpos, /noerase, /novtitle, trange=time
    labels = ['a. L2 Px5 - Px1']
    tx = xchsz*2
    ty = tpos[3]-ychsz*0.8
    xyouts, tx,ty,/normal, labels

    tx = tpos[0]
    ty = tpos[3]+ychsz*0.5
    xyouts, tx,ty,/normal, 'Test HOPE L2 Px5 and Px1 for '+strupcase(the_species)+$
        ', RBSP-'+strupcase(probe)+', Date='+time_string(time[0],tformat='YYYY-MM-DD')

    if keyword_set(test) then stop
    sgclose

end

time = time_double(['2013-06-07/04:50','2013-06-07/05:02'])
probes = ['a','b']
foreach probe, probes do test_hope_velocity_direction_l2, time, probe=probe
end
