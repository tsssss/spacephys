;+
; A wrapper of cotran.
;-

function stplot_cotran, in_var, msg, input=input, output=output, mission_probe=mission_probe, errmsg=errmsg, _extra=ex

    errmsg = ''
    retval = ''
    if n_elements(in_var) eq 0 then begin
        errmsg = 'No in_var ...'
        return, retval
    endif
    if tnames(in_var) eq '' then begin
        errmsg = 'in_var: '+in_var+' does not exist ...'
        return, retval
    endif

    if n_elements(msg) ne 0 then begin
        strs = strsplit(msg, '2', extract=1)
        in_coord = strs[0]
        out_coord = strs[1]
    endif else begin
        in_coord = strlowcase(get_setting(in_var, 'coord'))
        if n_elements(input) ne 0 then in_coord = input
        out_coord = output
    endelse
    in_coord = strlowcase(in_coord)
    out_coord = strlowcase(out_coord)
    the_msg = in_coord+'2'+out_coord

    index = strpos(in_var, in_coord)
    out_var = strmid(in_var,0,index)+out_coord+strmid(in_var,index+strlen(in_coord))
    stop

    get_data, in_var, times, vec_in, limit=lim
    vec_out = cotran(vec_in, times, the_msg, _extra=ex)
    store_data, out_var, times, vec_out, limit=lim
    add_setting, out_var, smart=1, dictionary($
        'display_type', 'vector', $
        'short_name', lim.short_name, $
        'unit', lim.unit, $
        'coord', strupcase(out_coord), $
        'coord_labels', lim.coord_labels )
    return, out_var

end