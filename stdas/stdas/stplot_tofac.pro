;+
; Convert a vector in xxx coord to fac, using a bvar and posvar in xxx coord.
; The vector, bvar, and posvar must be on the same times and in the same coord.
;-

pro stplot_tofac, tvar, posvar=rvar, bvar=bvar, newname=newname, labels=label, $
    addto=varlist

    get_data, tvar, uts, vec0, limits=lim
    get_data, bvar, uts, bvec
    get_data, rvar, uts, rvec

    ; b-para, p-perp in west, n-perp2 in north.
    ; in a tail-like configuration, [b,p,n] is close to gsm [x,y,z].
    bhat = sunitvec(bvec)
    rhat = sunitvec(rvec)
    phat = sunitvec(scross(rhat,bhat))
    nhat = scross(bhat,phat)
    vec1 = [[sdot(vec0,bhat)],[sdot(vec0,phat)],[sdot(vec0,nhat)]]

    if not keyword_set(newname) then newname = tvar+'_fac'
    store_data, newname, uts, vec1, limits=lim
    options, newname, 'labels', label+['b','w','n']
    
    if not keyword_set(varlist) then varlist=[]
    varlist = [varlist, newname]

end
