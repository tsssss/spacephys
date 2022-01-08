
pro rbsp_plot_aurora_efw, tr0, prob, type = type

    on_error, 2
    compile_opt idl2
    
    ; prepare.
    if n_elements(probe) eq 0 then prob = 'a'
    if n_elements(type) eq 0 then type = 'ast'
    re = 6374.4D    ; km.
    re1 = 1D/re
    tr = tr0 + [-1,1]*600   ; add 10 min pad.
    rbsp_efw_init
    timespan, tr0[0], tr0[1]-tr0[0], /second
    rbsp_load_spice_kernels
    no_spice_load = 1
    
    ; load efw.
    rbsp_load_efw_iono, tr, prob
    
    ; plot aurora and efw.
    dr = 11     ; sec.
    nrec = (tr0[1]-tr0[0])/dr
    for i = 0, nrec-1 do begin
        ut0 = tr0[0]+i*dr
        et0 = 1000D*ut0+62167219200000D
        wid = read_thm_asi(et0, type = type, /plot)
    endfor

end