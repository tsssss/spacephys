;+
; :Returns: A string array of loaded SPICE kernels.
; :Purpose: Get loaded SPICE kernels.
;-

function spice_get_loaded_kernel
    compile_opt idl2

  
    ; Check currently loaded kernels.
    cspice_ktotal, 'ALL', nkernel
    kernels = list()
    for ii=0,nkernel-1 do begin
        cspice_kdata, ii, 'ALL', file, file_type, source, handle, found
        kernels.add, file
    end

    return, kernels.toarray()
end

compile_opt idl2
kernels = spice_get_loaded_kernel()
foreach kernel, kernels do print, kernel
end
