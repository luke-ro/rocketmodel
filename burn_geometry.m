function [Ab,Vb,Vc] = burn_geometry(r,h,rb)

%  This burnback model only assumes burning with the cylindrical 
%  perforation of the grain.  Modify as needed for your rocket motor. 

    if rb >= r % motor is burnt out
        Ab = 0; % [m^2] 
    else % there is grain remaining
        %% BURN AREA
        Ab = 2*3.1416*(r + rb)*h; % [m^2] total burn area for cylindrically perforated grain

        %% VOLUME of PROPELLANT CONSUMED
        Vb = 3.1416*((r + rb)^2 - r^2)*h; % [m^3]
        
        %% CHAMBER VOLUME
        Vc = 3.1416*(r + rb)^2*h; % [m^3]  ignores the volume between the end of 
        % grain and nozzle entrance 
        
        end