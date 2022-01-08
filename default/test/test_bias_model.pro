; bonnell's paper on themis e field.

n0 = 3e-1   ; density, in cc.
te = 600    ; electron temperature, in eV.
th = 7*te   ; ion temperature, in eV.
i0 = 240    ; photo current, in nA.

mh = 1.67e-27   ; kg.   Proton mass.
me = 0.91e-30   ; kg.   Electron mass.
kb = 1.38e-23   ; J/K.  Boltzmann constant.
qe = 1.6e-19    ; in C.
qh = 1.6e-19    ; in C.
n0 = n0*1e6     ; in m^-3.
te = te*abs(qe) ; in J.
th = th*abs(qe) ; in J.

i_e0 = n0*abs(qe)*sqrt(te/me/!dpi/2)    ; in A.
i_h0 = n0*abs(qh)*sqrt(th/mh/!dpi/2)    ; in A.
i_p0 = i0*1e-9
v_p0 = 1                                ; in V.

vfs = smkarthm(-10,150,1,'dx')

i_e = i_e0*exp(abs(qe)*vfs/te)
i_h = i_h0*exp(abs(qh)*vfs/ti)
i_p = i_p0*exp(-vfs/v_p0)

i_e = i_e*1e9
i_h = i_h*1e9
i_p = i_p*1e9

end
