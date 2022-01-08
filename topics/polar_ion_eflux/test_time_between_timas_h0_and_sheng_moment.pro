
utr0 = time_double(['1996-03-17','1998-12-08'])
secofday = 86400d

nday = (utr0[1]-utr0[0])/secofday
ut1s = smkarthm(utr0[0], utr0[1], nday+1, 'n')
ut2s = ut1s+secofday

for i=0, nday-1 do begin
    utr = [ut1s[i],ut2s[i]]
    dat1 = sread_polar_timas(utr, type='h0', /local_only)
    dat2 = sread_polar_timas_sheng_moment(utr, /local_only)
    
    ut1s = sfmepoch(dat1.epoch_h,'unix')
    ut2s = dat2.h_ut_sec
    
    print, sdatarate(ut1s)
    print, sdatarate(ut2s)
    print, minmax(ut1s-ut2s)
endfor


end