;+
; Test to use quaternion in rotating between FAC and other coords.
;-

time = time_double(['2013-06-07','2013-06-08'])
probe = 'a'

rbsp_read_orbit, time, probe=probe
rbsp_read_bfield, time, probe=probe

prefix = 'rbsp'+probe+'_'
bvar = prefix+'b_gsm'
rvar = prefix+'r_gsm'
components = ['b','w','o']

    ; The low-res version.
    define_fac, bvar, rvar, time_var=rvar
    r_fac_var = prefix+'r_fac'
    to_fac, rvar, to=r_fac_var
    vars = prefix+components+'hat_gsm'
    foreach var, vars, ii do begin
        newvar = var+'_lowres'
        rename_var, var, to=newvar
    endforeach

    ; The high-res version.
    define_fac, bvar, rvar, time_var=bvar
    vars = prefix+components+'hat_gsm'
    foreach var, vars, ii do begin
        newvar = var+'_highres'
        rename_var, var, to=newvar
    endforeach
    
    
    ; The quaternion.
    vars = prefix+components+'hat_gsm_lowres'
    get_data, vars[0], times
    ntime = n_elements(times)
    ndim = 3
    m_xxx2fac = fltarr(ntime,ndim,ndim)
    for ii=0,ndim-1 do m_xxx2fac[*,ii,*] = get_var_data(vars[ii])
    q_xxx2fac = mtoq(m_xxx2fac)
    xxx = get_setting(bvar, 'coord')
    xxx2fac = strlowcase(xxx)+'2fac'
    q_var = prefix+'q_'+xxx2fac
    store_data, q_var, times, q_xxx2fac
    add_setting, q_var, {$
        in_coord: xxx, $
        out_coord: 'FAC', $
        in_coord_labels: ['x','y','z'], $
        out_coord_labels: components}
    m_fac2xxx = fltarr(ntime,ndim,ndim)
    for ii=0,ntime-1 do m_fac2xxx[ii,*,*] = transpose(m_xxx2fac[ii,*,*])
    
    get_data, r_fac_var, times, r_fac
    r_gsm1 = rotate_vector(r_fac, m_fac2xxx)
    get_data, rvar, limits=lim
    store_data, prefix+'r_gsm1', times, r_gsm1, limits=lim
stop
end
