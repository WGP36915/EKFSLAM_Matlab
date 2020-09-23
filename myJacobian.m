function Gt = myJacobian(angle, v, w, dt)
%     xUp = dt*v*cos(angle);
%     yUp = dt*v*sin(angle);
    
    xUp = -(v/w)*sin(angle) + (v/w)*sin(angle+w*dt);
    yUp = (v/w)*cos(angle) - (v/w)*cos(angle+w*dt);
    Gt = [1,0,-yUp; 0,1,xUp; 0,0,1];    
end