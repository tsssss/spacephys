function tracers_get_local_root
    compile_opt idl2
    return, join_path([default_local_root(),'tracers'])
end
