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


    define_fac_quaternion, b_var=bvar, r_var=rvar
    r_fac_var = prefix+'r_fac'
    to_fac, rvar, to=r_fac_var

    ; The quaternion.
    define_fac_quaternion, r_var=rvar, b_var=bvar
    qvar = prefix+'q_gsm2fac'
    newvar2 = prefix+'r_gsm1'
    from_fac, r_fac_var, to=newvar2, q_var=qvar

end