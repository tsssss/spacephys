
pro test_sg_routine_array

sgindexcolor, ct = 43

nrec = 100
xx = findgen(nrec)*2*!pi/(nrec-1)
yy = sin(xx)

kws = {xstyle:1,yrange:[-1.2,1.2],position:[0.1,0.1,0.45,0.85],xticks:3}
call_procedure, 'plot', xx, yy, _extra = kws

kws = {normal:1,charsize:2,color:6}
call_procedure, 'xyouts', 0.1, 0.9, 'Plot', _extra = kws

data = dist(nrec,nrec)
kws = {xstyle:1,ystyle:1,position:[0.55,0.1,0.9,0.85],noerase:1}
tcall = {routine:'contour', data:ptr_new(data), keywords:kws}
calls = ptr_new(tcall)
kws = {normal:1,charsize:2,color:4}
tcall = {routine:'xyouts', data:{x:0.55,y:0.9,v:'Contour'}, keywords:kws}
calls = [calls, ptr_new(tcall)]
sgcallexe, calls
sgcallexport, calls

end