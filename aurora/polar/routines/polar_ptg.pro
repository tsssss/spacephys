;to get help: idl> ptg,/help


; ancillary routines --------------------------------------------

function atan2d,x1,x2
   return,double(atan(x1,x2)/!dtor)
end

function dtand,x
    return,double(tan(x*!dtor))
end

function datand,x
    return,double(atan(x)/!dtor)
end

function fgeodep,a,b,v1x,v1y,v1z,v2x,v2y,v2z
    return,v1x*v2x + v1y*v2y + v1z*v2z * a*a/(b*b)
end

function dfmag,x1,x2,x3
    return,double(sqrt(x1*x1+x2*x2+x3*x3))
end

;---------------------------------------------------------------

pro vector_to_ra_decp,x,y,z,ra,dec


    fill_value = -1.d31
    ndx = where(z ne 0,count)
    if(count gt 0) then dec(ndx) = 90.*z(ndx)/abs(z(ndx))

    tmp = sqrt(x*x + y*y)
    ndx = where(tmp ne 0,count)
    if (count gt 0) then begin
      dec(ndx) = atan2d(z(ndx),tmp(ndx))
      ra(ndx)  = atan2d(y(ndx),x(ndx))
    endif

    ndx = where((ra lt 0) and (ra ne fill_value),count)
    if (count gt 0) then ra(ndx) = ra(ndx) + 360.

end

;---------------------------------------------------------------

pro drtollp,x,y,z,lat,lon,r

      lat = atan2d(z,sqrt(x*x + y*y))
      lon = atan2d(y,x)
      r   = sqrt(x*x + y*y + z*z)

      tmp = where(x eq y and x eq 0)
      if ((size(tmp))(0) ne 0) then begin
         lat(tmp)  = 90.d * z(tmp)/abs(z(tmp))
         lon(tmp) = 0.d
         r = 6371.d
      endif

      tmp2 = where(lon lt 0)
      if ((size(tmp2))(0) ne 0) then begin
         lon(tmp2) = lon(tmp2) + 360.d
      endif

end

;---------------------------------------------------------------

pro get_scalarp,ox,oy,oz,lx,ly,lz,emis_hgt,ncols,nrows,s,f, $
                off_axis,num_off_axis

;...  equatoral radius (km) and polar flattening of the earth
;     ref: table 15.4, 'explanatory supplement to the
;          astronomical almanac,' k. seidelmann, ed. (1992).
      re_eq = 6378.136d
      inv_f = 298.257d

;...  initialize output
      s =  dblarr(ncols,nrows)
      s1 = dblarr(ncols,nrows)
      s2 = dblarr(ncols,nrows)

;...  get polar radius
      re_po = re_eq*(1.d - 1.d / inv_f)

;...  get radii to assumed emission height
      ree = re_eq + emis_hgt
      rep = re_po + emis_hgt

;...  get flattening factor based on new radii
      f = (ree - rep)/ree

;...  get elements of quadratic formula
      a = fgeodep(ree,rep,lx,ly,lz,lx,ly,lz)
      b = fgeodep(ree,rep,lx,ly,lz,ox,oy,oz) * 2.d
      c = fgeodep(ree,rep,ox,oy,oz,ox,oy,oz) - ree*ree

;...  check solutions to quadratic formula
      determinant = b*b - 4.d * a*c
;...  remove points off the earth
      determinant = determinant > 0.
      off_axis = where(determinant eq 0.,num_off_axis)
      if(num_off_axis gt 0) then b(off_axis) = 0.d
;...  solve quadratic formula (choose smallest solution)
      s1 = ( -b + sqrt(determinant) ) / ( 2.d * a )
      s2 = ( -b - sqrt(determinant) ) / ( 2.d * a )

      s = s1<s2

end

;-------------------------------------------------------------------------
;  routine:	ptg
;-------------------------------------------------------------------------

pro polar_ptg, $
    et, emis_hgt, att, orb, system, l0, geodetic = geodetic, $  ; input vars.
    gclat, gclon    ; output vars.
    ; getra=getra,ra=ra,dec=dec,s=s,lpixx=lpixx,lpixy=lpixy,lpixz=lpixz, $
    ; posx=posx,posy=posy,posz=posz

     ncols = 200
     nrows = 228
     zrot  = dblarr(ncols,nrows)
     yrot  = dblarr(ncols,nrows)
     gclat = dblarr(ncols,nrows)
     gclon = dblarr(ncols,nrows)
     ra    = dblarr(ncols,nrows)
     dec   = dblarr(ncols,nrows)

     primary   = 1
     secondary = 2
     fill_value = -1.d31

;... define orthonormal coordinate axes
     xax = l0/dfmag(l0(0),l0(1),l0(2))
     yax = crossp(att,l0)
     yax = yax/dfmag(yax(0),yax(1),yax(2))
     zax = crossp(xax,yax)

;... single pixel angular resolution
; old code
;     pr = 0.03449d       ; 9-bin mean primary detector 9/26/97 (pyth)
;     pc = 0.03983d       ; same

; new code  7/28/00
case system of
  1: begin
       pc = 0.03987d  ; w. swift 10/97
       pr = 0.03380d  ; w. swift 10/97
     end
  2: begin
       pc = 0.04178d  ; w. swift 10/97
       pr = 0.03511d  ; w. swift 10/97
;       pc = 0.04159d ; w. peria 7/28/00
;       pr = 0.03384d ; w. peria 7/28/00
     end
endcase

; new approach per k. clark 8/4/00
pc = tan(pc*!dtor)  ; m. brittnacher 8/4/00
pr = tan(pr*!dtor)  ; m. brittnacher 8/4/00

;... initialize output arrays to default
     gclat(*,*) = fill_value
     gclon(*,*) = fill_value
        ra(*,*) = fill_value
       dec(*,*) = fill_value

;... find rotation angles for each pixel
     if (system eq secondary) then begin
       a = (findgen(200)-99.5)*pc
       b = replicate(1.,228)
       zrot = a#b
       c = (findgen(228)-113.5)*pr
       d = replicate(1.,200)
       yrot = d#c
     endif else begin
       if (system eq primary) then begin
         a = (findgen(200)-99.5)*pc
         b = replicate(1.,228)
         zrot = a#b
         c = -(findgen(228)-113.5)*pr
         d = replicate(1.,200)
         yrot = d#c
       endif else begin
         ;  error trap
         return
       endelse
     endelse

;... determine lpix
;     tanz = tan(zrot*!dtor)
;     tany = tan(yrot*!dtor)
; new approach per k. clark 8/4/00
     tanz = zrot
     tany = yrot

     lpx = 1.d /sqrt(1.d + tany*tany + tanz*tanz)
     lpy = lpx*tanz
     lpz = lpx*tany

     lpixx = lpx*xax(0) + lpy*yax(0) + lpz*zax(0)
     lpixy = lpx*xax(1) + lpy*yax(1) + lpz*zax(1)
     lpixz = lpx*xax(2) + lpy*yax(2) + lpz*zax(2)

;... convert lpix to a unit vector
     mag = dfmag(lpixx,lpixy,lpixz)

     lpixx = lpixx/mag
     lpixy = lpixy/mag
     lpixz = lpixz/mag

;    calculate right ascension and declination
     if(keyword_set(getra)) then $
        vector_to_ra_decp,lpixx,lpixy,lpixz,ra,dec

;... find scalar (s) such that s*l0 points to
;    the imaged emission source.  if the line of
;    sight does not intersect the earth s=0.0
     ox = orb(0)
     oy = orb(1)
     oz = orb(2)
     get_scalarp,ox,oy,oz,lpixx,lpixy,lpixz,emis_hgt,ncols,nrows,s,f, $
                 off_axis,num_off_axis

     posx = ox + s*lpixx
     posy = oy + s*lpixy
     posz = oz + s*lpixz

;... convert from gci to geo coordinates.  rotm is the
;    rotation matrix.
     cdf_epoch, et, yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
     time = long([yr*1000+stodoy(yr,mo,dy), 1000d*(sc+60d*(mi+60d*hr))+msc])
     ic_gci_to_geo,time,rotm

     p_geox = rotm(0,0)*posx + rotm(1,0)*posy + rotm(2,0)*posz
     p_geoy = rotm(0,1)*posx + rotm(1,1)*posy + rotm(2,1)*posz
     p_geoz = rotm(0,2)*posx + rotm(1,2)*posy + rotm(2,2)*posz

;... get geocentric lat/lon.  this converts from
;    a 3 element vector to two angles: lat & longitude
     drtollp,p_geox,p_geoy,p_geoz,gclat,gclon,r
     gclat = gclat < 90.

;... convert to geodetic lat/lon.  f is the flattening
;    factor of the earth.  see get_scalar for details.
;    ref: spacecraft attitude determination and control,
;    j.r. wertz, ed., 1991, p.821.
     if(keyword_set(geodetic)) then begin
        gdlat = 90.d + 0.d * gclat
        ndx = where(gclat lt 90.,count)
        if(count gt 0) then begin
           gdlat(ndx) = datand(dtand(gclat(ndx))/((1.d - f)*(1.d - f)))
        endif
        gclat = gdlat
     endif

     if (num_off_axis gt 0) then begin
       gclat(off_axis) = fill_value
       gclon(off_axis) = fill_value
     endif

end
