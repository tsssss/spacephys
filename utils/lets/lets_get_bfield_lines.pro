;+
;-
function lets_get_bfield_lines, input_times, r_gsms, $
    external_model=external_model, internal_model=internal_model, t89_use_kp=t89_use_kp, h0=h0, refine=refine

    if n_elements(external_model) eq 0 then external_model = 't89'
    if n_elements(internal_model) eq 0 then internal_model = 'igrf'

    npoint = n_elements(r_gsms[*,0])
    if n_elements(input_times) eq npoint then times = input_times else times = dblarr(npoint)+input_times[0]
    time_range = minmax(times)

    if n_elements(h0) eq 0 then h0 = 100d
    r0 = h0/constant('re')+1
    if n_elements(refine) eq 0 then refine = 1

;---Prepare external model parameters.
    has_external_model = 0
    avail_external_models = ['t89','t96','t01','t04s','ts04']
    index = where(strlowcase(external_model) eq avail_external_models, count)
    if count ne 0 then has_external_model = 1
    if has_external_model then begin
        t89_par = keyword_set(t89_use_kp)? !null: 2d
        par_time_range = time_range+[-1,1]*300
        par_var = geopack_read_par(par_time_range, model=external_model, t89_par=t89_par)
        pars = get_var_data(par_var, at=times)
    endif

    tmp = geopack_resolve_model(external_model)
    t89 = tmp.t89
    t96 = tmp.t96
    t01 = tmp.t01
    ts04 = tmp.ts04
    storm = tmp.storm

    if external_model eq 'igrf' then igrf = 1
    if internal_model eq 'dip' or internal_model eq 'dipole' then igrf = 0
    if internal_model eq 'igrf' then igrf = 1
    if n_elements(igrf) eq 0 then igrf = 0


    flines = list()
    last_time = !null
    foreach time, times, time_id do begin
        if last_time ne time then begin
            ps = geopack_recalc(time)
        endif
        last_time = time
        rx = r_gsms[time_id,0]
        ry = r_gsms[time_id,1]
        rz = r_gsms[time_id,2]

        my_fline = []
        ; -1 for northern hemisphere, 1 for southern hemispher.
        foreach trace_dir, [-1,1] do begin
            geopack_trace, rx,ry,rz, trace_dir, reform(pars[time_id,*]), $
                xf,yf,zf, r0=r0, refine=refine, ionosphere=1, fline=the_fline, $
                t89=t89, t96=t96, t01=t01, ts04=ts04, storm=storm, igrf=igrf
            if n_elements(the_fline) le 3 then continue ; failed to trace.
            if trace_dir eq -1 then my_fline = [my_fline,the_fline] else my_fline = [reverse(the_fline,1),my_fline]
        endforeach
        flines.add, my_fline
    endforeach

    return, flines

end