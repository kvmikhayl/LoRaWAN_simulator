%km_LoRaWAN_model_GW_distribution - generates the coordinates of the GWs
%for LoRaWAN simulation model
%
% Syntax:  [cood_cartesian,cood_polar] = km_LoRaWAN_model_GW_distribution(input_number_GWs,mode)
%
% Inputs:
%    input_number_GWs - number of GWs
%    mode - mode of operation, 0 - normal, other values - different test
%    cases
%
% Outputs:
%    cood_cartesian - cartesian 2D coordinates (distances are in meters)
%    cood_polar - polar 2D coordinates (distances are in meters)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author: Konstantin Mikhaylov, Dr. Tech, Wireless Communications
% University of Oulu, Oulu, Finland
% email address: konstantin.mikhaylov(at)oulu.fi
% Website: http://cc.oulu.fi/~kmikhayl/index.html
% March 2020; Last revision: 15-March-2020

function [cood_cartesian,cood_polar] = km_LoRaWAN_model_GW_distribution(input_number_GWs,mode)
    if(mode==0)%normal operation
        rmin=1*10^2;%minimum radius in meters
        rmax=3*10^3;%maximum radius in meters
        if input_number_GWs==1
            cood_cartesian(1,1)=0;%
            cood_cartesian(1,2)=0;% 
        elseif  input_number_GWs==2   
            cood_cartesian(1,1)=-rmax/2;%
            cood_cartesian(1,2)=0;%   
            cood_cartesian(2,1)=rmax/2;%
            cood_cartesian(2,2)=0;%
        elseif  input_number_GWs==4  
            cood_cartesian(1,1)=-rmax/2;%
            cood_cartesian(1,2)=-rmax/2;;%   
            cood_cartesian(2,1)=rmax/2;%
            cood_cartesian(2,2)=-rmax/2;;%    
            cood_cartesian(3,1)=-rmax/2;%
            cood_cartesian(3,2)=rmax/2;;%   
            cood_cartesian(4,1)=rmax/2;%
            cood_cartesian(4,2)=rmax/2;;%              
        else
            rmin=1*10^2;%minimum radius in meters
            rmax=4.79*10^3;%maximum radius in meters
            cood_cartesian=zeros(input_number_GWs,2);
            cood_polar=zeros(input_number_GWs,2);
            for i=1:input_number_GWs
                while 1
                   cood_cartesian(i,1)=(rand(1,1)-0.5)*2*rmax;%
                   cood_cartesian(i,2)=(rand(1,1)-0.5)*2*rmax;%
                   [cood_polar(i,1),cood_polar(i,2)]=cart2pol(cood_cartesian(i,1),cood_cartesian(i,2));
                   if (cood_polar(i,2)>=rmin) && (cood_polar(i,2)<=rmax)
                        break
                   end
                end
            end            
        end
        for i=1:input_number_GWs
            [cood_polar(i,1),cood_polar(i,2)]=cart2pol(cood_cartesian(i,1),cood_cartesian(i,2));
        end
    else%test case
        switch mode
            case num2cell(1:100)
                cood_cartesian(1,1)=-1000;%
                cood_cartesian(1,2)=0;%
                cood_cartesian(2,1)=1000;%
                cood_cartesian(2,2)=0;% 
            case num2cell(101:200)
                cood_cartesian(1,1)=-1000;%
                cood_cartesian(1,2)=0;%
                cood_cartesian(2,1)=1000;%
                cood_cartesian(2,2)=0;%    
            case num2cell(201:300)
                cood_cartesian(1,1)=-2000;%
                cood_cartesian(1,2)=0;%
                cood_cartesian(2,1)=2000;%
                cood_cartesian(2,2)=0;%    
            case num2cell(301:400)
                cood_cartesian(1,1)=0;%
                cood_cartesian(1,2)=0;%               
            case default
                
        end
        for i=1:input_number_GWs
            [cood_polar(i,1),cood_polar(i,2)]=cart2pol(cood_cartesian(i,1),cood_cartesian(i,2));
        end
    end
end

