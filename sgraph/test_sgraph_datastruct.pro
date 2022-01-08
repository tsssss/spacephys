pro test_sgraph_datastruct, dat
    compile_opt idl2
    
    x0 = findgen(101)/100*2*!const.pi
    dat = !null
    
    tdat = {name:'rbspa_vb1', $
        x:{data:ptr_new(x0),axis:ptr_new(!x)}, $
        y:{data:ptr_new(1.5*sin(x0)),axis:ptr_new(replicate(!y,2))}, $
        z:{data:ptr_new(),axis:ptr_new(!z)}}
    dat = [dat, tdat]
    
    tdat = {name:'rbspb_vb1', $
        x:{data:ptr_new(x0),axis:ptr_new(!x)}, $
        y:{data:ptr_new(0.7*cos(x0)),axis:ptr_new(!y)}, $
        z:{data:ptr_new(),axis:ptr_new(!z)}}
    dat = [dat, tdat]

end

test_sgraph_datastruct, dat
end