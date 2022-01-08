
_2014_0828_load_data
get_data, 'mos', uts, mos, pxidx
get_data, 'mos_info', 0, mosinfo
imgsz = mosinfo.imgsz


;---Settings.
    id = '2014_0828_10'
    rootdir = sparentdir(srootdir())
    datadir = rootdir+'/data' & if file_test(datadir,/directory) eq 0 then file_mkdir, datadir
    plotdir = rootdir+'/plot' & if file_test(plotdir,/directory) eq 0 then file_mkdir, plotdir
    
;---Gen mosaic movie.
    tfn = plotdir+'/asf_img.png'
    vfn = datadir+'/asf_img.mp4'
    fps = 10d
    ovid = idlffvideowrite(vfn)
    vidstream = ovid.AddVideoStream(imgsz[0],imgsz[1],fps)
    mlonrng = [-1,1]*90
    mlatrng = [minlat,90]
    mlonticks = [-1,0,1]*90
    mlatticks = [60,70,80]
    asirng = [50,200]
    tpos = [0,0,1,1]
    white = sgcolor('white')
    
    for i=0, nrec-1 do begin
        timg = fltarr(imgsz)
        timg[pxidx] = mos[i,*]
        tmp = bytscl(timg, min=asirng[0], max=asirng[1], top=254)
        
        sgopen, tfn, xsize=imgsz[0], ysize=imgsz[1]
        
        sgtv, tmp, position=tpos, ct=1
        sgset_map, position=tpos, color=white, $
            xrange=mlonrng, $
            yrange=mlatrng, ytickv=mlatticks, yminor=1, $
            minorgrid=1
        xyouts, 10,10, /device, time_string(uts[i],tformat='YYYY-MM-DD/hh:mm:ss'), color=white
        ; plot s/c footpoint on?
        
        sgclose
        
        read_png, tfn, timg
        time = ovid.put(vidstream, timg)
    endfor
    ovid.Cleanup
    file_delete, tfn

end
