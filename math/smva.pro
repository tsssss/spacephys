
function smva, bs, bvecs

    mb = dblarr(3,3)
    for i = 0, 2 do $
        for j = 0, i do $
            mb[i,j] = mean(bs[*,i]*bs[*,j])-mean(bs[*,i])*mean(bs[*,j])
    mb[0,1] = mb[1,0]
    mb[0,2] = mb[2,0]
    mb[1,2] = mb[2,1]

    vals = eigenql(mb, eigenvectors = bvecs, residual = bres)
    ; bvecs[0,*] maximam,
    ; bvecs[1,*] intermdeiate,
    ; bvecs[2,*] minimum variance.

    return, sunitvec(vals)

end

; load some data.
fn = shomedir()+'/Google Drive/works/data/rbsp_de/psbl_de_2013_0501_0735_b.tplot'
tplot_restore, filename = fn

stop
utr = time_double(['2013-05-01/07:33','2013-05-01/07:42'])
tvar = 'rbspb_b_gse'
get_data, tvar, uts, dat
idx = where(uts ge utr[0] and uts le utr[1])

tmp = smva(dat[idx,*])


end
