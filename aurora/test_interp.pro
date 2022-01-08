
FUNCTION interp2p,x1,x2,x,f1,f2

    RETURN, f1 + (f2-f1)*(x-x1)/(x2-x1)
    
END

;--------------------------------------------------------------

FUNCTION favg3p,x1,x2,x3,f1,f2,f3

    RETURN, 0.5*( f2 + ( f1*(x2-x1) + f3*(x3-x2) )/(x3-x1) )
    
END

;--------------------------------------------------------------

FUNCTION favg4p,x1,x2,x3,x4,f1,f2,f3,f4

    RETURN, 0.5*( f1*(x2-x1) + f2*(x3-x1) + f3*(x4-x2) $
        + f4*(x4-x3) )/(x4-x1)
        
END

;--------------------------------------------------------------

FUNCTION favg5p,x1,x2,x3,x4,x5,f1,f2,f3,f4,f5

    RETURN, 0.5*( f1*(x2-x1) + f2*(x3-x1) + f3*(x4-x2) $
        + f4*(x5-x3) + f5*(x5-x4) )/(x5-x1)
        
END

;--------------------------------------------------------------

FUNCTION favg6p,x1,x2,x3,x4,x5,x6,f1,f2,f3,f4,f5,f6

    RETURN, 0.5*( f1*(x2-x1) + f2*(x3-x1) + f3*(x4-x2) $
        + f4*(x5-x3) + f5*(x6-x4) + f6*(x6-x5) )/(x6-x1)
        
END

;--------------------------------------------------------------

FUNCTION favg7p,x1,x2,x3,x4,x5,x6,x7,f1,f2,f3,f4,f5,f6,f7

    RETURN, 0.5*( f1*(x2-x1) + f2*(x3-x1) + f3*(x4-x2) $
        + f4*(x5-x3) + f5*(x6-x4) + f6*(x7-x5) $
        + f7*(x7-x6) )/(x7-x1)
        
END

x = findgen(7)
y = sin(x/100)
x0 = 3.4

print, y
print, spl_interp(x,y,spl_init(x,y),x0)
print, favg7p(x[0],x[1],x[2],x[3],x[4],x[5],x[6], $
    y[0],y[1],y[2],y[3],y[4],y[5],y[6])
end
