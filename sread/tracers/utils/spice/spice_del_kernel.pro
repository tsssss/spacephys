;+
; :Purpose: Delete specified kernels.
; :Parameters:
;   kernels: in, optional, string array of kernels to be deleted.
;-

pro spice_del_kernel, kernels
    compile_opt idl2

    if n_elements(kernels) eq 0 then begin
        cspice_kclear
        return
    endif

    loaded_kernels = spice_get_loaded_kernels()
    foreach kernel, kernels do begin
        index = where(loaded_kernels eq kernel, count)
        if count eq 0 then begin
            print, 'Kernel '+kernel+' is not loaded ...'
            continue
        end
        lprmsg, 'Unloading kernel '+kernel+' ...'
        cspice_unload, kernel
    endforeach
    

end