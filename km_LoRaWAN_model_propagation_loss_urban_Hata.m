%km_LoRaWAN_model_propagation_loss_urban_Hata - calculates the attenuation
%for the specified distance in meters using urban Hata model
%
% Syntax:  [output_propagation_loss_dB] = km_LoRaWAN_model_propagation_loss_urban_Hata(input_distance_meters)
%
% Inputs:
%    input_distance_meters - vector of distances in meters
%
% Outputs:
%    output_propagation_loss_dB - vector of resulting attenuation values in
%    dB
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author: Konstantin Mikhaylov, Dr. Tech, Wireless Communications
% University of Oulu, Oulu, Finland
% email address: konstantin.mikhaylov(at)oulu.fi
% Website: http://cc.oulu.fi/~kmikhayl/index.html
% March 2020; Last revision: 15-March-2020
function [output_propagation_loss_dB] = km_LoRaWAN_model_propagation_loss_urban_Hata(input_distance_meters)
    %SRC: https://en.wikipedia.org/wiki/Hata_model
    CH=0.8+(1.1*log10(868)-0.7)*1-1.56*log10(868);
    output_propagation_loss_dB=69.55+26.16*log10(868)-13.82*log10(24)-CH+(44.9-6.55*log10(24))*log10(input_distance_meters/1000);
end

