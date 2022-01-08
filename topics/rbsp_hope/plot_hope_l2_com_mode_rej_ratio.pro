;+
; Type: procedure.
; Purpose: Plot common mode rejection ratio for pixel 2-2, 1-3, 0-4, 
;   and the number flux for the 5 pixels at selected energy.
; Parameters:
;   ion, in, strng, opt. 'proton','oxygen','helium'. Default is 'proton'.
; Keywords:
;   probes, in, string, opt. 'a','b'. Default is 'a'.
;   tr, in, dblarr[2], req. Time range in unix time.
; Notes: none.
; Dependence: slib.
; History:
;   2016-02-11, Sheng Tian, create.
;-

pro plot_hope_l2_com_mode_rej_ratio, ion, probes = probe, tr = utr

    if n_elements(utr) ne 2 then message, 'no time range ...'
    if n_elements(ion) eq 0 then ion = 'proton'

    vars = ['FPDU','Epoch_Ion','HOPE_ENERGY_Ion', $
        'Epoch_Ion_DELTA', 'L_Ion','MLT_Ion', $
        'Epoch','Sector_Collapse_Cntr','Energy_Collapsed']

    case strlowcase(ion) of
        'proton': vars[0] = ['FPDU']
        'oxygen': vars[0] = ['FODU']
        'helium': vars[0] = ['FHEDU']
    endcase

    hope = sread_rbsp_hope(utr, probes = probe, vars = vars, type = 'level2')



    ; assume that ion time is the start of the ion mode.
    ; reconstruct the times for each sector.

    ; decompose the sectors.
    ; out: 'fxdu_px[0,1,2,3,4]
    npixel = 5
    nenergy = 72
    for i = 0, npixel-1 do rbsp_hope_l2_decompose_sector, hope, i, type = 'fpdu'



    ; energy based operation.
    enidx = indgen(nenergy)
    enidx = [5,55]
    
    for i = 0, n_elements(enidx)-1 do begin
        
        ; calc the flux on opposite pixels for all energy bins.
        rbsp_hope_l2_pixel_flux, hope, [2,2], type = 'fpdu', index = enidx[i]
        rbsp_hope_l2_pixel_flux, hope, [1,3], type = 'fpdu', index = enidx[i]
        rbsp_hope_l2_pixel_flux, hope, [0,4], type = 'fpdu', index = enidx[i]

        ; decompose the pixels.
        vars = 'fpdu'+'_en'+string(enidx[i],format='(I02)')+ $
            '_px'+string(indgen(npixel),format='(I0)')
        for j = 0, npixel-1 do begin
            tvar = 'fpdu'+'_px'+string(j,format='(I0)')
            get_data, tvar, t0, d0, e0      
            store_data, vars[j], t0, d0[*,enidx[i]], $
                limits = {ylog:1, yrange:[1e5,1e11], $
                labels:'Pixel '+string(j,format='(I0)')+'!C  '+$
                    sgnum2str(e0[0,enidx[i]],nsgn=3)+' eV', $
                ytitle:'Flux!C(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)'}
        endfor
        
    endfor

    sgindexcolor, 43
    tplot_options, 'constant', 0
    tplot_options, 'labflag', -1
    
    vars = 'fpdu_en??_df??'
    options, vars, 'ystyle', 1
    options, vars, 'yrange', [-1,1]

    vars = 'fpdu_en55_*'
    titl = 'RBSP-'+strupcase(probe)+' '+ion
    tplot, vars, title = titl

end

utr = time_double(['2013-04-14/00:00','2013-04-14/00:10'])
probe = 'a'

plot_hope_l2_com_mode_rej_ratio, tr = utr, probes = probe

end
