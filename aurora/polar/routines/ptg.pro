;to get help: IDL> ptg,/help


; ancillary routines --------------------------------------------

FUNCTION atan2d,x1,x2
   RETURN,DOUBLE(ATAN(x1,x2)/!DTOR)
END

FUNCTION dtand,x
    RETURN,DOUBLE(TAN(x*!DTOR))
END

FUNCTION datand,x
    RETURN,DOUBLE(ATAN(x)/!DTOR)
END

FUNCTION fgeodeP,a,b,v1x,v1y,v1z,v2x,v2y,v2z
    RETURN,v1x*v2x + v1y*v2y + v1z*v2z * a*a/(b*b)
END

FUNCTION dfmag,x1,x2,x3
    RETURN,DOUBLE(SQRT(x1*x1+x2*x2+x3*x3))
END

;---------------------------------------------------------------

PRO vector_to_ra_decP,x,y,z,ra,dec


    fill_value = -1.D31
    ndx = WHERE(z NE 0,count)
    IF(count GT 0) THEN dec(ndx) = 90.*z(ndx)/ABS(z(ndx))

    tmp = SQRT(x*x + y*y)
    ndx = WHERE(tmp NE 0,count)
    IF (count GT 0) THEN BEGIN
      dec(ndx) = atan2d(z(ndx),tmp(ndx))
      ra(ndx)  = atan2d(y(ndx),x(ndx))
    ENDIF

    ndx = WHERE((ra LT 0) AND (ra NE fill_value),count)
    IF (count GT 0) THEN ra(ndx) = ra(ndx) + 360.

END

;---------------------------------------------------------------

PRO drtollP,x,y,z,lat,lon,r

      lat = atan2d(z,SQRT(x*x + y*y))
      lon = atan2d(y,x)
      r   = SQRT(x*x + y*y + z*z)

      tmp = WHERE(x EQ Y AND x EQ 0)
      IF ((size(tmp))(0) NE 0) THEN BEGIN
         lat(tmp)  = 90.D * z(tmp)/ABS(z(tmp))
         lon(tmp) = 0.D
         r = 6371.D
      ENDIF

      tmp2 = WHERE(lon LT 0) 
      IF ((size(tmp2))(0) NE 0) THEN BEGIN
         lon(tmp2) = lon(tmp2) + 360.D
      ENDIF

END

;---------------------------------------------------------------

PRO get_scalarP,Ox,Oy,Oz,Lx,Ly,Lz,emis_hgt,ncols,nrows,s,f, $
                off_axis,num_off_axis

;...  Equatoral radius (km) and polar flattening of the earth
;     Ref: Table 15.4, 'Explanatory Supplement to the
;          Astronomical Almanac,' K. Seidelmann, ed. (1992).
      re_eq = 6378.136D
      inv_f = 298.257D

;...  initialize output
      s =  DBLARR(ncols,nrows)
      s1 = DBLARR(ncols,nrows)
      s2 = DBLARR(ncols,nrows)

;...  get polar radius
      re_po = re_eq*(1.D - 1.D / inv_f)

;...  get radii to assumed emission height
      ree = re_eq + emis_hgt
      rep = re_po + emis_hgt

;...  get flattening factor based on new radii
      f = (ree - rep)/ree

;...  get elements of quadratic formula
      a = fgeodeP(ree,rep,Lx,Ly,Lz,Lx,Ly,Lz)
      b = fgeodeP(ree,rep,Lx,Ly,Lz,Ox,Oy,Oz) * 2.D
      c = fgeodeP(ree,rep,Ox,Oy,Oz,Ox,Oy,Oz) - ree*ree

;...  check solutions to quadratic formula
      determinant = b*b - 4.D * a*c 
;...  remove points off the earth
      determinant = determinant > 0. 
      off_axis = WHERE(determinant EQ 0.,num_off_axis) 
      IF(num_off_axis GT 0) THEN b(off_axis) = 0.D
;...  solve quadratic formula (choose smallest solution) 
      s1 = ( -b + SQRT(determinant) ) / ( 2.D * a ) 
      s2 = ( -b - SQRT(determinant) ) / ( 2.D * a ) 

      s = s1<s2

END

;-------------------------------------------------------------------------
;  ROUTINE:	ptg
;-------------------------------------------------------------------------

PRO ptg,system,time,l0,att,orb,emis_hgt,gclat,gclon $
       ,geodetic=geodetic,getra=getra,ra=ra,dec=dec,s=s $
       ,LpixX=LpixX,LpixY=LpixY,LpixZ=LpixZ $
       ,posX=posX,posY=posY,posZ=posZ $
       ,versStr=versStr,help=help

    IF(KEYWORD_SET(help)) THEN BEGIN
       PRINT,''
       PRINT,' PRO ptg,system,time,l0,att,orb,emis_hgt,gclat,gclon
       PRINT,''
       PRINT,' Original base code:  UVIPTG'
       PRINT,' 7/31/95  Author:  G. Germany'
       PRINT,' Development into PTG: 01/15/98'
       PRINT,' Authors:  Mitch Brittnacher & John O''Meara'
       PRINT,''
       PRINT,' calculates geocentric lat,lon, for a complete image
       PRINT,' 
       PRINT,' input
       PRINT,'    system          =1 primary; =2 secondary
       PRINT,'    time            time(1)=yyyyddd, time(2)=msec of day 
       PRINT,'    L0              gci look direction (from uvilook)
       PRINT,'    att             gci attitude 
       PRINT,'    orb             gci position
       PRINT,'    emis_hgt        los altitude
       PRINT,'
       PRINT,' output
       PRINT,'    gclat           geocentric latitude
       PRINT,'    gclon           geocentric longitude
       PRINT,'
       PRINT,' keywords
       PRINT,'    geodetic        (set) returns geodetic values if set
       PRINT,'    getra           (set) calulates ra & dec if set
       PRINT,'       ra           (out) right ascension (deg)
       PRINT,'      dec           (out) declination (deg)
       PRINT,'        s           (out) scalar for lpix
       PRINT,'    lpixX           (out) x component of unit look direction
       PRINT,'    lpixY           (out) y component of unit look direction
       PRINT,'    lpixZ           (out) z component of unit look direction
       PRINT,'     posX           (out) x,y,z components of vector from
       PRINT,'     posY           (out)       earth center to emission
       PRINT,'     posZ           (out) 
       PRINT,'  versStr           (out) software version string
       PRINT,'
       PRINT,' external library routines required
       PRINT,'    ic_gci_to_geo
       PRINT,'
       PRINT,' NOTES:
       PRINT,'
       PRINT,' 1. Unlike UVIPTG, this routine returns latitude and longitude
       PRINT,'    for all pixels in an image.  It does the calculation in a
       PRINT,'    fraction of the time required by UVIPTG.
       PRINT,'
       PRINT,' 2. The default lat/lon values are in geocentric coordinates.
       PRINT,'    Geographic (geocentric) coordinates assume the earth is
       PRINT,'    a sphere and are defined from the center of the sphere.
       PRINT,'    For geodetic coordinates, the earth is assumed to be an 
       PRINT,'    ellipsoid of revolution.  See the routine fgeode for 
       PRINT,'    details.  
       PRINT,'    Geodetic coordinates are defined from the normal to the 
       PRINT,'    geode surface.  To enable geodetic calculations, set the 
       PRINT,'    keyword /GEODETIC.
       PRINT,' 
       PRINT,' 3. The look direction for a specified pixel (Lpix) is
       PRINT,'    calculated from the look direction of the center of the
       PRINT,'    UVI field of view (L0) by successive rotations in
       PRINT,'    row and column directions.  Each pixel is assumed to have
       PRINT,'    a fixed angular width.  The angular distance from the center
       PRINT,'    of the pixel to the center of the fov is calculated and then
       PRINT,'    L0 is rotated into Lpix.
       PRINT,'
       PRINT,'    Unlike UVIPTG, this routine explicitly calculates three
       PRINT,'    orthogonal axes whereas UVIPTG implicitly assumed the image
       PRINT,'    z-axis was given by the attitude vector.
       PRINT,'
       PRINT,' 4. The secondary and primary detectors have different 
       PRINT,'    orientations and require different rotations between L0 and 
       PRINT,'    Lpix.
       PRINT,'    
       PRINT,' 5. Geocentric lat/lon values are the intersection
       PRINT,'    of the look direction for the specified pixel (Lpix) and
       PRINT,'    the surface of the earth.  The geocentric values are then
       PRINT,'    transformed into geodetic values.  The vector from the
       PRINT,'    center of the earth to the intersection is pos so that
       PRINT,'    pos = orb + S*Lpix, where orb is the GCI orbit vector
       PRINT,'    and S is a scalar.
       PRINT,'
       PRINT,' 6. The intersection of Lpix and the earth is calculated first
       PRINT,'    in GCI coordinates and then converted to geographic 
       PRINT,'    coordinates.  The conversion is by means of ic_gci_to_geo.  
       PRINT,'    This routine and its supporting routines, was taken from 
       PRINT,'    the CDHF and is part of the ICSS_TRANSF_orb call.
       PRINT,'
       PRINT,' 7. The viewed emissions are assumed to originate emis_hgt km
       PRINT,'    above the surface of the earth.  See get_scalar for details.
       PRINT,'
       PRINT,'10. The keywords POS(xyz) are needed for LOS corrections.
       PRINT,'
       RETURN
     ENDIF   

     versStr = 'PTG v1.5  7/27/00'
     ncols = 200
     nrows = 228
     zrot  = DBLARR(ncols,nrows)
     yrot  = DBLARR(ncols,nrows)
     gclat = DBLARR(ncols,nrows)
     gclon = DBLARR(ncols,nrows)
     ra    = DBLARR(ncols,nrows)
     dec   = DBLARR(ncols,nrows)

     primary   = 1
     secondary = 2
     fill_value = -1.D31

;... Define orthonormal coordinate axes
     xax = l0/dfmag(l0(0),l0(1),l0(2))
     yax = CROSSP(att,l0)
     yax = yax/dfmag(yax(0),yax(1),yax(2))
     zax = CROSSP(xax,yax)

;... single pixel angular resolution
; old code
;     pr = 0.03449D       ; 9-bin mean primary detector 9/26/97 (Pyth)
;     pc = 0.03983D       ; same

; new code  7/28/00
CASE system OF
  1: BEGIN 
       pc = 0.03987D  ; W. Swift 10/97
       pr = 0.03380D  ; W. Swift 10/97
     END
  2: BEGIN
       pc = 0.04178D  ; W. Swift 10/97
       pr = 0.03511D  ; W. Swift 10/97
;       pc = 0.04159D ; W. Peria 7/28/00
;       pr = 0.03384D ; W. Peria 7/28/00
     END
ENDCASE

; new approach per K. Clark 8/4/00
pc = tan(pc*!DTOR)  ; M. Brittnacher 8/4/00
pr = tan(pr*!DTOR)  ; M. Brittnacher 8/4/00

;... initialize output arrays to default
     gclat(*,*) = fill_value
     gclon(*,*) = fill_value
        ra(*,*) = fill_value
       dec(*,*) = fill_value

;... find rotation angles for each pixel
     IF (system EQ secondary) THEN BEGIN
       a = (FINDGEN(200)-99.5)*pc
       b = REPLICATE(1.,228)
       zrot = a#b
       c = (FINDGEN(228)-113.5)*pr
       d = REPLICATE(1.,200)
       yrot = d#c 
     ENDIF ELSE BEGIN 
       IF (system EQ primary) THEN BEGIN
         a = (FINDGEN(200)-99.5)*pc
         b = REPLICATE(1.,228)
         zrot = a#b
         c = -(FINDGEN(228)-113.5)*pr
         d = REPLICATE(1.,200)
         yrot = d#c 
       ENDIF ELSE BEGIN 
         ;  error trap
         RETURN 
       ENDELSE 
     ENDELSE

;... Determine Lpix
;     tanz = tan(zrot*!DTOR)
;     tany = tan(yrot*!DTOR)
; new approach per K. Clark 8/4/00
     tanz = zrot
     tany = yrot

     lpx = 1.D /SQRT(1.D + tany*tany + tanz*tanz)
     lpy = lpx*tanz
     lpz = lpx*tany

     LpixX = lpx*xax(0) + lpy*yax(0) + lpz*zax(0)
     LpixY = lpx*xax(1) + lpy*yax(1) + lpz*zax(1)
     LpixZ = lpx*xax(2) + lpy*yax(2) + lpz*zax(2)

;... Convert Lpix to a Unit Vector
     mag = dfmag(LpixX,LpixY,LpixZ)

     LpixX = LpixX/mag
     LpixY = LpixY/mag
     LpixZ = LpixZ/mag

;    calculate right ascension and declination
     IF(KEYWORD_SET(getra)) THEN $
        vector_to_ra_decP,LpixX,LpixY,LpixZ,ra,dec

;... Find scalar (s) such that s*L0 points to
;    the imaged emission source.  If the line of
;    sight does not intersect the earth s=0.0
     Ox = orb(0)
     Oy = orb(1)
     Oz = orb(2)
     get_scalarP,Ox,Oy,Oz,LpixX,LpixY,LpixZ,emis_hgt,ncols,nrows,s,f, $
                 off_axis,num_off_axis

     posX = Ox + s*LpixX
     posY = Oy + s*LpixY
     posZ = Oz + s*LpixZ

;... Convert from GCI to GEO coordinates.  ROTM is the
;    rotation matrix.

     ic_gci_to_geo,time,rotm
     
     p_geoX = rotm(0,0)*posX + rotm(1,0)*posY + rotm(2,0)*posZ
     p_geoY = rotm(0,1)*posX + rotm(1,1)*posY + rotm(2,1)*posZ
     p_geoZ = rotm(0,2)*posX + rotm(1,2)*posY + rotm(2,2)*posZ

;... Get geocentric lat/lon.  this converts from
;    a 3 element vector to two angles: lat & longitude
     drtollP,p_geoX,p_geoY,p_geoZ,gclat,gclon,r
     gclat = gclat < 90.

;... Convert to geodetic lat/lon.  F is the flattening
;    factor of the Earth.  See get_scalar for details.
;    Ref: Spacecraft Attitude Determination and Control,
;    J.R. Wertz, ed., 1991, p.821.
     IF(KEYWORD_SET(geodetic)) THEN BEGIN
        gdlat = 90.D + 0.D * gclat
        ndx = WHERE(gclat LT 90.,count)
        IF(count GT 0) THEN BEGIN
           gdlat(ndx) = datand(dtand(gclat(ndx))/((1.D - f)*(1.D - f)))
        ENDIF
        gclat = gdlat
     ENDIF

     IF (num_off_axis GT 0) THEN BEGIN
       gclat(off_axis) = fill_value
       gclon(off_axis) = fill_value
     ENDIF

END
