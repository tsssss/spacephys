;+
; Load RBSP injections for a given probe and time_range.
;-

function rbsp_themis_injection_read_rbsp_injection, input_time_range, probe=input_probe

    pinfo = rbsp_themis_injection_load_project()
    data_dir = pinfo['data_dir']
    injection_file = join_path([data_dir,'rbsp_injection_list_v01.txt'])
    if file_test(injection_file) eq 0 then begin
        injection_file = rbsp_themis_injection_gen_rbsp_injection_list(filename=injection_file)
    endif

    lines = read_all_lines(injection_file)
    nline = n_elements(lines)
    nheader = 2
    ninjection = nline-nheader
    the_probes = strarr(ninjection)
    time_ranges = dblarr(ninjection,2)
    r_gses = fltarr(ninjection,6)
    lshells = fltarr(ninjection,2)
    mlts = fltarr(ninjection,2)
    mlats = fltarr(ninjection,2)
    lines = lines[nheader:*]
    foreach line, lines, ii do begin
        infos = strsplit(line,' ,',extract=1)
        the_probes[ii] = infos[0]
        time_ranges[ii,*] = time_double(infos[1:2])
        r_gses[ii,*] = float(infos[3:8])
        lshells[ii,*] = float(infos[9:10])
        mlts[ii,*] = float(infos[11:12])
        mlats[ii,*] = float(infos[13:14])
    endforeach

    injection_info = dictionary()
    if n_elements(input_probe) eq 0 then begin
        probes = ['a','b']
    endif else begin
        probes = input_probe
    endelse
    foreach probe, probes do begin
        index = where(the_probes eq probe, nrec)
        the_info = dictionary($
            'time_ranges', time_ranges[index,*], $
            'r_gses', reform(r_gses[index,*],[nrec,3,2]), $
            'lshells', lshells[index,*], $
            'mlts', mlts[index,*], $
            'mlats', mlats[index,*] )
        
        ; Injections are identified based on natural days. This will cut one injection into two.
        ; Here we merge the cut ones into one.
        the_time_ranges = the_info['time_ranges']
        index = where(the_time_ranges[0:nrec-2,1] eq the_time_ranges[1:nrec-1,0], count)
        if count ne 0 then begin
            foreach key, the_info.keys() do begin
                val = the_info[key]
                if key eq 'r_gses' then begin
                    val[index,*,1] = val[index+1,*,1]
                    val[index+1,*,0] = val[index,*,0]
                endif else begin
                    val[index,1] = val[index+1,1]
                    val[index+1,0] = val[index,0]
                endelse
                the_info[key] = val
            endforeach
        endif
        

        ; Remove duplicated time ranges.
        the_time_ranges = the_info['time_ranges']
        index = uniq(the_time_ranges[*,0], sort(the_time_ranges[*,0]))
        foreach key, the_info.keys() do begin
            val = the_info[key]
            if key eq 'r_gses' then begin
                val = val[index,*,*]
            endif else begin
                val = val[index,*]
            endelse
            the_info[key] = val
        endforeach
        
        ; Apply time range.
        if n_elements(input_time_range) eq 2 then begin
            time_range = time_double(input_time_range)
            the_time_ranges = the_info['time_ranges']
            index = where(the_time_ranges[*,0] ge time_range[0] and the_time_ranges[*,1] le time_range[1], count)
            if count eq 0 then begin
                foreach key, the_info.keys() do the_info[key] = !null
            endif else begin
                foreach key, the_info.keys() do begin
                    val = the_info[key]
                    if key eq 'r_gses' then begin
                        val = val[index,*,*]
                    endif else begin
                        val = val[index,*]
                    endelse
                    the_info[key] = val
                endforeach
            endelse
        endif
        
        ; Save the results.
        injection_info[probe] = the_info
    endforeach

    return, injection_info

end


tmp = rbsp_themis_injection_read_rbsp_injection(['2014','2015'])
end