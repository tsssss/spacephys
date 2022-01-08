
tutr = time_double(['2012-09-25','2012-09-26'])
tuts = smkarthm(tutr[0], tutr[1], 60, 'dx')
nrec = n_elements(tuts)
tets = stoepoch(tuts, 'unix')
rgsm = dblarr(nrec,3)
rgsm[*,0] = smkarthm(1,2, nrec, 'n')
rgsm[*,1] = smkarthm(-1,-2, nrec, 'n')
rgsm[*,2] = smkarthm(1,-1, nrec, 'n')

store_data, 'rgsm0', tuts, rgsm

cotrans, 'rgsm0', 'rgse0', /gsm2gse
cotrans, 'rgse0', 'rgei0', /gse2gei

rgse = sgsm2gse(rgsm, tets)
store_data, 'rgse1', tuts, rgse
rgei = sgse2gei(rgse, tets)
store_data, 'rgei1', tuts, rgei


device, decomposed=0
loadct2, 43
tplot_options, 'labflag', -1

vars = ['rgsm0','rgse0','rgei0','rgse1','rgei1']
options, vars, 'colors', [6,4,2]
options, vars, 'labels', ['x','y','z']
tplot, vars, trange=tutr

end