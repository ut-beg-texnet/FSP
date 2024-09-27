%interpolate pressure values from grid at desired locations


function pout = getP_interp(x,y,X,Y,Z)
pout = interp2(X,Y,Z,x,y) ;
end