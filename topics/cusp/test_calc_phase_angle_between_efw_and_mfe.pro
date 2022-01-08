

fn = '/Users/Sheng/po_rawfld_1998_0914_03.sdt'
sdt = ssdtread(fn)

; uniform time.
dr = 0.125
tmp = sdt.var.magmfex.depend_0
uts = smkarthm(min(tmp),max(tmp),dr,'dx')

; B field. 0.12 sec, ~8 Hz.
tmp = sdt.var.magmfex.depend_0
bx = interpol(sdt.var.magmfex.value, tmp, uts)
by = interpol(sdt.var.magmfey.value, tmp, uts)

store_data, 'po_bx', uts, bx
store_data, 'po_by', uts, by

; E field. 0.025 sec, 40 Hz.
tmp = sdt.var.polar_v12l_c.depend_0
v12 = interpol(sdt.var.polar_v12l_c.value, tmp, uts)
v34 = interpol(sdt.var.polar_v34l_c.value, tmp, uts)

store_data, 'po_v12', uts, v12
store_data, 'po_v34', uts, v34


end
