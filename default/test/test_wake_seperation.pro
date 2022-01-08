fn = '/Users/Sheng/Downloads/po_sdt_wake_1998_0914.sdt'

sdt = ssdtread(fn)

tmp = sdt.var.polar_v12l_c.depend_0
dr = 0.05   ; 20 Hz.
t0 = smkarthm(min(tmp), max(tmp), dr, 'dx')
nrec = n_elements(t0)
v12 = interpol(sdt.var.polar_v12l_c.value, tmp, t0)
v34 = interpol(sdt.var.polar_v34l_c.value, tmp, t0)

order = 2

d0 = v12
d1 = smooth(d0,3/dr)        ; de-wake.
w1 = d0-d1                  ; wake.
d2 = d1+smooth(w1,3/dr)     ;
w2 = d0-d2

store_data, 'v12', t0, d0
store_data, 'v12_dewake', t0, d2
store_data, 'v12_wake', t0, w2


d0 = v34
d1 = smooth(d0,3/dr)        ; de-wake.
w1 = d0-d1                  ; wake.
d2 = d1+smooth(w1,3/dr)     ;
w2 = d0-d2

store_data, 'v34', t0, d0
store_data, 'v34_dewake', t0, d2
store_data, 'v34_wake', t0, w2

tplot_options, 'constant', 0

tplot, ['v12','v12_dewake','v12_wake', $
    'v34','v34_dewake','v34_wake']

end
