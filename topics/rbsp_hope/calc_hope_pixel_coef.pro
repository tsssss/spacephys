;+
; Type: function.
; Purpose: Caclulate the coefficients that balance opposite pixels, 
;   for ion species, energy bins, and pixels.
; Parameters:
;   ion, in, strng, opt. 'proton','oxygen','helium'. Default is 'proton'.
; Keywords:
;   probes, in, string, opt. 'a','b'. Default is 'a'.
;   tr, in, dblarr[2], req. Time range in unix time.
; Return: [nenerngy,npixel]. The coef for pixel 0,1,2 are 1.
; Notes: none.
; Dependence: slib.
; History:
;   2016-02-11, Sheng Tian, create.
;-


function calc_hope_pixel_coef, ion, probes = probe, tr = utr

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

    nenergy = 72
    npixel = 5


    ; assume that ion time is the start of the ion mode.
    ; reconstruct the times for each sector.

    ; decompose the sectors.
    ; out: 'fxdu_px[0,1,2,3,4]
    npixel = 5
    nenergy = 72
    for i = 0, npixel-1 do rbsp_hope_l2_decompose_sector, hope, i, type = 'fpdu'



    ; energy based operation.
    enidx = indgen(nenergy)
;    enidx = [55]
    nenergy = n_elements(enidx)

    coefs = dblarr(nenergy,npixel)+1

    
    for i = 0, n_elements(enidx)-1 do begin
        
        ; calc the flux on opposite pixels for all energy bins.
        rbsp_hope_l2_balance_pixel_flux, hope, [2,2], type = 'fpdu', index = enidx[i], coef & coefs[i,2] = coef[1]
        rbsp_hope_l2_balance_pixel_flux, hope, [1,3], type = 'fpdu', index = enidx[i], coef & coefs[i,3] = coef[1]
        rbsp_hope_l2_balance_pixel_flux, hope, [0,4], type = 'fpdu', index = enidx[i], coef & coefs[i,4] = coef[1]

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
    
    ; cannot be this large.
    idx = where(coefs ge 1.5, cnt)
    if cnt ne 0 then coefs[idx] = 1
    idx = where(coefs le 1d/1.5, cnt)
    if cnt ne 0 then coefs[idx] = 1
    
    return, coefs

end


date = '2013-05-14'
;date = '2014-04-01'

timespan, date, 1, /day
utr = time_double(date)+[0,86400d]

probe = 'a'

coefs = calc_hope_pixel_coef(tr = utr, probes = probe)


end
