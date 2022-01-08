;+
; Generate survey plot for FAST, include (1) survey of ele and ion,
; (2) ele spectrogram, density, pressure, etc. (3) ion the same quantities.
; Determine FAST cusp time, using ctime.
;-
pro cusp_fast_esa_survey_figure, fn, trange = tr, id = id, noplot = noplot

    device, decomposed = 0 & loadct2, 43
    !p.font = 1
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', 1
    tplot_options, 'zcharsize', 0.9
    time_stamp, /off
    red = 6

    tmp = file_basename(fn)
    strput, tmp, 'fld', strpos(tmp,'esa')
    tmp = file_dirname(fn)+'/'+tmp
    tplot_restore, filename = tmp
    tplot_restore, filename = fn
    labs = 'fa_'+['ilat','mlt','dis']
    get_data, 'fa_dis', t0
    options, ['ion_n','ion_p','ele_n','ele_p'], 'ylog', 1
    vars = ['ion_p','ele_p']
    for i = 0, n_elements(vars)-1 do begin
        get_data, vars[i], t0, tmp
        store_data, vars[i], t0, tmp*1e-8
        options, vars[i], 'ytitle', '(ba)!Ccgs'
    endfor
    
    coef = 1.6e-9   ; convert # cm^-2 s^-1 to uA m^-2.
    vars = 'ion_'+['para_spec','perp_spec','anti_spec',$
        'j','eflux','n','en_spec','pa_spec']
    vars = [vars,'ele_'+['j','eflux','n','en_spec','pa_spec']]
    vars = [vars,'E_ALONG_V','E_NEAR_B']
    titl = 'FAST ESA and E field, '+time_string(t0[0],tformat='YYYY-MM-DD')

    window, xsize = 700, ysize = 960
    tplot, vars, trange = tr, var_label = labs, title = titl
    ctime, cusptr
    print, time_string(cusptr)
    
    pre = 'ion_'
    vars = pre+['para_spec','perp_spec','anti_spec',$
        'j','eflux','n','p','en_spec','pa_spec']
    titl = 'FAST ESA ion, '+time_string(t0[0],tformat='YYYY-MM-DD')
    fn = shomedir()+'/'+id+'_esa_ion.eps'
    sgpsopen, fn, xsize = 6, ysize = 8, /inch
    tplot, vars, trange = tr, var_label = labs, title = titl
    timebar, cusptr, color = red, thick = 2
    if keyword_set(noplot) then begin
        sgpsclose
        file_delete, fn
    endif else sgpsclose, /pdf
    
    pre = 'ele_'
    vars = pre+['para_spec','perp_spec','anti_spec',$
        'j','eflux','n','p','en_spec','pa_spec']
    titl = 'FAST ESA ele, '+time_string(t0[0],tformat='YYYY-MM-DD')

    device, decomposed = 0 & loadct2, 43
    tplot, vars, trange = tr, var_label = labs, title = titl
    fn = shomedir()+'/'+id+'_esa_ele.eps'
    sgpsopen, fn, xsize = 6, ysize = 8, /inch
    tplot, vars, trange = tr, var_label = labs, title = titl
    timebar, cusptr, color = red, thick = 2
    if keyword_set(noplot) then begin
        sgpsclose
        file_delete, fn
    endif else sgpsclose, /pdf
        
    ; combine to one plot.
    mu = '!9'+string(109b)+'!X'
    gamma = '!9'+string(71b)+'!X'
    para = '||'
    perp = '!9'+string(94b)+'!X'
    
    titl = 'FAST ESA, ' + time_string(t0[0],tformat='YYYY-MM-DD')
    fn = shomedir()+'/'+id+'_esa.eps'
    get_data, 'ele_n', t0, elen
    get_data, 'ion_n', tmp, ionn
    ionn = interpol(ionn, tmp, t0)
    store_data, 'density', t0, [[ionn],[elen]], limits = $
        {labels:['Ni','Ne'], colors:[-1,6], ytitle:'(cm!U-3!N)',ylog:1}
    get_data, 'ele_j', t0, elej
    get_data, 'ion_j', tmp, ionj
    ionj = interpol(ionj, tmp, t0)
    store_data, 'nflux', t0, [[ionj],[elej]], limits = $
        {labels:['Ji','Je'], colors:[-1,6], ytitle:'(!9'+mu+'A/m!U2!N)'}
    vars = ['ion_en_spec','ion_pa_spec','nflux','ion_eflux','ele_eflux', $
        'density','ele_en_spec']
    sgpsopen, fn, xsize = 6, ysize = 8, /inch
    tplot, vars, trange = tr, var_label = labs, title = titl
    timebar, cusptr, color = red, thick = 2
    if keyword_set(noplot) then begin
        sgpsclose
        file_delete, fn
    endif else sgpsclose, /pdf
    
    print, time_string(cusptr)
end

rootdir = shomedir()+'/Google Drive/works'
eventid = '1998_0925_06'

fn = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
info = cusp_read_conjun_list(fn, event = eventid)
tr = info.fast.plot_time
orb = string(info.fast.orbit,format='(I05)')
fn = rootdir+'/data/cusp/fa_sdt_esa_'+eventid+'_'+orb+'.tplot'
cusp_fast_esa_survey_figure, fn, trange = tr, id = eventid, /noplot
end