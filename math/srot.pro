
; rotating the vector.

pro srot, vx, vy, vz, ang, axis, degree = degree

    t = keyword_set(degree)? !dtor*ang: ang
    cost = cos(t[0])
    sint = sin(t[0])
    
    if size(axis,/type) ne 7 then ax = axis[0] else $
    ax = where(strlowcase(axis[0]) eq ['x','y','z'])
    
    case ax of
        0: begin & vt = vy*cost + vz*sint
            vz = temporary(vz)*cost - temporary(vy)*sint
            vy = temporary(vt) & end
        1: begin & vt = vx*cost - vz*sint
            vz = temporary(vx)*sint + temporary(vz)*cost
            vx = temporary(vt) & end
        2: begin & vt = vx*cost + vy*sint
            vy = temporary(vy)*cost - temporary(vx)*sint
            vx = temporary(vt) & end
        else: message, 'invalid axis'
    endcase
end