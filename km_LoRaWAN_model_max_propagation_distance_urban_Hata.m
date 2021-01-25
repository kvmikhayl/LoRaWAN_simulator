%km_LoRaWAN_model_max_propagation_distance_urban_Hata - calculates maximum
%propagation distance in meters for specified maximum attenuation in dB
%using urban Hata model
%
% Syntax:  [output_distance_meters] = km_LoRaWAN_model_max_propagation_distance_urban_Hata(input_max_loss_dB)
%
% Inputs:
%    input_max_loss_dB - vector of attenuation values in dB
%
% Outputs:
%    output_distance_meters - vector of resulting maximum communication
%    distances in meters
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author: Konstantin Mikhaylov, Dr. Tech, Wireless Communications
% University of Oulu, Oulu, Finland
% email address: konstantin.mikhaylov(at)oulu.fi
% Website: http://cc.oulu.fi/~kmikhayl/index.html
% March 2020; Last revision: 15-March-2020
function [output_distance_meters] = km_LoRaWAN_model_max_propagation_distance_urban_Hata(input_max_loss_dB)
    %SRC: https://en.wikipedia.org/wiki/Hata_model
    CH=0.8+(1.1*log10(868)-0.7)*1-1.56*log10(868);
    output_distance_meters=(10.^((input_max_loss_dB-(69.55+26.16*log10(868)-13.82*log10(24)-CH))/(44.9-6.55*log10(24))))*1000;
end

