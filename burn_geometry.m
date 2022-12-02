function [Ab,Vb,Vc] = burn_geometry(ri,ro,h,bd)

%  This burnback model only assumes burning with the cylindrical 
%  perforation of the grain.  Modify as needed for your rocket motor. 

    if bd >= ro-ri % motor is burnt out
        Ab = 0; % [m^2] 
        Vb = 3.1416*((ri + bd)^2 - ri^2)*h; % [m^3]
        Vc = 3.1416*(ri + bd)^2*h; % [m^3]  ignores the volume between the end of 
        
    else % there is grain remaining
        %% BURN AREA
        Ab = 2*3.1416*(ri + bd)*h; % [m^2] total burn area for cylindrically perforated grain

        %% VOLUME of PROPELLANT CONSUMED
        Vb = 3.1416*((ri + bd)^2 - ri^2)*h; % [m^3]
        
        %% CHAMBER VOLUME
        Vc = 3.1416*(ri + bd)^2*h; % [m^3]  ignores the volume between the end of 
        % grain and nozzle entrance        
    end