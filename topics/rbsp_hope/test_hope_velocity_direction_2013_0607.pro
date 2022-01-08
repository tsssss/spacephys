;+
; Test the direction of HOPE velocity:
;   1. During O+ outflow events.
;-

pro test_hope_velocity_direction, time_range, probe=probe, $
    ion_energy_range=ion_energy_range, electron_energy_range=electron_energy_range


;---Check inputs.
    if n_elements(time_range) ne 2 then message, 'Invalid time range ...'
    if n_elements(probe) eq 0 then message, 'Invalid probe ...'

;---Constants.
    errmsg = ''
    test = 0

    time = time_range
    pre0 = 'rbsp'+probe+'_'
    species = 'o'
    mass = 1.67e-27*16
    q0 = 1.6e-19
    pre1 = pre0+species+'_'
    rbsp_read_hope_moments, time, probe=probe, species=species;, /renew_file
;    rbsp_read_hope_moments, time, probe=probe, species=species, energy_range=[1e3,1e7];, /renew_file


;---Load magnetic field data.
    rbsp_read_bfield, time_range, probe=probe, errmsg=errmsg
    if errmsg ne '' then message, errmsg
    bvar = pre0+'b_gsm'
    interp_time, bvar, to=pre1+'density'

    rbsp_read_orbit, time_range, probe=probe, errmsg=errmsg
    if errmsg ne '' then message, errmsg
    rvar = pre0+'r_gsm'
    interp_time, rvar, to=pre1+'density'

    define_fac, bvar, rvar, time_var=pre1+'vbulk'
    vars = pre1+['vbulk','number_flux','energy_flux','enthalpy']
    foreach var, vars do to_fac, var, to=var+'_fac'
    
;---Load V_ExB.
    if keyword_set(load_vexb) then begin
        id = '2013_0607_event_info'
        if tnames(id) eq '' then load = 1 else load = 0
        if load then _2013_0607_load_data3
        get_data, 'rb'+probe+'_vexb_fac', times, vexb, limits=lims
        store_data, pre0+'vexb_fac', times, -vexb
        interp_time, pre0+'vexb_fac', to=pre0+'o_vbulk'
        add_setting, pre0+'vexb_fac', /smart, {$
            coord: 'FAC', $
            coord_labels: ['||','west','outward'], $
            short_name: 'ExB V', $
            display_type: 'vector', $
            colors: sgcolor(['red','green','blue']), $
            unit: 'km/s'}
    endif

;---Load L3 pitch angle data.
    pitch_angles = [162.,144,36,18]
    npitch_angle = n_elements(pitch_angles)
    foreach pitch_angle, pitch_angles do rbsp_read_en_spec, time_range, probe=probe, pitch_angle=pitch_angle


;---Configuration.
    vars = pre0+[species+'_vbulk_fac','vexb_fac']
    foreach var, vars do add_setting, var, /smart, {coord_labels:['||','west','north']}
    vars = pre0+species+'_enspec_'+string(pitch_angles,format='(I0)')+'deg'
    energy_range =  [1,4e4]
    options, vars, 'yrange', energy_range

    vars = pre0+[species+'_'+['enspec_'+string(pitch_angles,format='(I0)')+'deg','vbulk'+['','_fac']]]
    if keyword_set(load_vexb) then vars = pre0+[species+'_'+['enspec_'+string(pitch_angles,format='(I0)')+'deg','vbulk'+['','_fac']],'vexb_fac']
    nvar = n_elements(vars)
    file = shomedir()+'/test_hope_v_parallel_in_o_outflow_'+time_string(time_range[0],tformat='YYYY_MMDD')+'_rbsp'+probe+'.pdf'
    if keyword_set(test) then file = 0
    sgopen, file, xsize=8.5, ysize=11
    
    poss = sgcalcpos(nvar, lmargin=16, tmargin=5, rmargin=12, xchsz=xchsz, ychsz=ychsz)
    
    device, decomposed=0
    loadct2, 43
    tpos = poss[*,0:npitch_angle-1]
    tplot, vars[0:npitch_angle-1], position=tpos, /noerase, /nouttick

    device, decomposed=1
    tpos = poss[*,npitch_angle:*]
    options, vars[npitch_angle:*], 'constant', 0
    tplot, vars[npitch_angle:*], position=tpos, /noerase, /novtitle
    
;---Add labels.
    labels = ['a','b','c','d','e','f','g']+'. '+$
        ['L3 PA='+string(pitch_angles,format='(I0)')+' deg','L2 V '+['GSM','FAC'],'-V ExB FAC']
    for ii=0, nvar-1 do begin
        tpos = poss[*,ii]
        tx = xchsz*2
        ty = tpos[3]-ychsz*0.8
        xyouts, tx,ty,/normal, labels[ii]
    endfor
    
    get_data, pre1+'vbulk_fac', times, vfac
    vfac_en = 0.5*mass*(vfac*1e3)^2/q0
    tpos = poss[*,0]
    plot, time_range, energy_range, xstyle=5, ystyle=5, ylog=1, xlog=1, $
        /nodata, /noerase, position=tpos
    oplot, times, vfac_en, color=sgcolor('red')
    txs = tpos[0]+[1,6]*xchsz
    tys = tpos[1]+[1,1]*ychsz
    plots, txs, tys, /normal, color=sgcolor('red')
    tx = txs[1]+xchsz
    ty = tys[0]-0.3*ychsz
    xyouts, tx,ty,/normal, 'Energy corresponds to m!D'+species+'!N|V!D||!N|!U2!N/2'
    
    tpos = poss[*,0]
    tx = tpos[0]
    ty = tpos[3]+ychsz*0.5
    xyouts, tx,ty,/normal, 'Test HOPE L2 V parallel for '+strupcase(species)+$
        ' against L3 data, RBSP-'+strupcase(probe)+', Date='+time_string(time_range[0],tformat='YYYY-MM-DD')
        
    tpos = poss[*,nvar-1]
    tx = tpos[0]+xchsz*1
    ty = tpos[3]-ychsz*1
    xyouts, tx,ty,/normal, 'Positive V!D||!N is earthward in the northerm hemisphere'
    if keyword_set(test) then stop
    sgclose

end

time_ranges = list()
time_ranges.add, time_double(['2013-06-07/04:50','2013-06-07/05:02'])
;time_ranges.add, time_double(['2013-06-07/04:50','2013-06-07/05:02'])
;time_ranges.add, time_double(['2013-05-01/07:30','2013-05-01/07:45'])
; Gkioulidou+2019
;time_range = time_double(['2015-06-23/03:00','2015-06-23/07:00'])
probes = ['a','b']
foreach time_range, time_ranges do foreach probe, probes do test_hope_velocity_direction, time_range, probe=probe

;time_range = time_double(['2013-05-01/07:30','2013-05-01/07:45'])
;probe = 'b'


;probe = 'b'


end
