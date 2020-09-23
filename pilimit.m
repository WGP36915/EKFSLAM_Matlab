function angle=pilimit(angle)                    %½Ç¶ÈÇø·Ö
    %angle
    if angle > pi
        angle = angle - 2*pi;
    else if angle < -pi 
        angle = angle + 2*pi;
    end
    %angle
end
