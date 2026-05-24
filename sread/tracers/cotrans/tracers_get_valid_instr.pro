;+
; :Returns: string array, valid TRACERS instrument names.
;-

function tracers_get_valid_instr

    compile_opt idl2
    return, ['mag','efi','aci','ace','magic','msc']

end

