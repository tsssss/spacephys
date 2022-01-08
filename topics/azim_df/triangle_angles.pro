;+
; Calculate the angles of a triangle.
;
; vecs. An array in [ntime,3 (xyz),3 (npoint)].
;-

function triangle_angles, vecs

    ndim = 3.
    nrec = n_elements(vecs)
    if nrec eq 0 then return, !null

    angles = fltarr(nrec/ndim^2,ndim)
    index = indgen(ndim)
    for ii=0,ndim-1 do begin
        index = shift(index,1)
        angles[*,ii] = sang($
            vecs[*,*,index[0]]-vecs[*,*,index[1]], $
            vecs[*,*,index[0]]-vecs[*,*,index[2]], /deg)
    endfor

    return, angles

end