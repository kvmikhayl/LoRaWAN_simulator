%km_LoRaWAN_model_UL_transmissions - generates the matrix of UL
%transmission starts for LoRaWAN simulation model
%
% Syntax:  [UL_TX_start_milliseconds] = km_LoRaWAN_model_UL_transmissions(mode,num_EDs,timespan_seconds)
%
% Inputs:
%    mode - mode of operation, 0 - normal, other values - different test
%    cases
%    num_EDs - number of EDs
%    timespan_seconds - duration of the simulation time
%
% Outputs:
%    UL_TX_start_milliseconds - matrix of size (num_EDs,x) with the
%    timestamps of UL transmission start
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author: Konstantin Mikhaylov, Dr. Tech, Wireless Communications
% University of Oulu, Oulu, Finland
% email address: konstantin.mikhaylov(at)oulu.fi
% Website: http://cc.oulu.fi/~kmikhayl/index.html
% March 2020; Last revision: 15-March-2020

function [UL_TX_start_milliseconds] = km_LoRaWAN_model_UL_transmissions(mode,num_EDs,timespan_seconds, input_T)
    min_Period_s=5;
    rand_Period_s=input_T;
    %need to make sure that period between packets exceeds maximum time for
    %RX tranmission
    if(mode==0)%normal operation
        for i_ED=1:1:num_EDs
            i_packet=1;
            rand_start_ms=rand()*rand_Period_s*1000;
            timestamp_ms=rand_start_ms;
            while(timestamp_ms<timespan_seconds*1000)
                UL_TX_start_milliseconds(i_ED,i_packet)=timestamp_ms;
                timestamp_ms=timestamp_ms+min_Period_s*1000+rand()*rand_Period_s*1000;
                i_packet=i_packet+1;
            end
        end
    else%test case
        switch mode
            case 1 %1 packet per ED
                i_packet=1;
                for i_ED=1:1:num_EDs
                    rand_start_ms=rand()*rand_Period_s*1000;
                    timestamp_ms=rand_start_ms;
                    UL_TX_start_milliseconds(i_ED,i_packet)=timestamp_ms;
                end
    
                            
            case 2
                Period_s=10;
                
                rand_start_ms=rand()*Period_s*1000;  
                i_packet=1;
                timestamp_ms=rand_start_ms;
                while(timestamp_ms<timespan_seconds*1000)
                    for i_ED=1:1:num_EDs
                        UL_TX_start_milliseconds(i_ED,i_packet)=timestamp_ms;
                    end
                    timestamp_ms=timestamp_ms+Period_s*1000;
                    i_packet=i_packet+1;
                end  
          case 3
                Period_s=5;
                rand_start_ms=rand()*Period_s*1000;  
                i_packet=1;
                timestamp_ms=rand_start_ms;
                while(timestamp_ms<timespan_seconds*1000)
                    for i_ED=1:1:num_EDs
                        UL_TX_start_milliseconds(i_ED,i_packet)=timestamp_ms;
                    end
                    timestamp_ms=timestamp_ms+Period_s*1000;
                    i_packet=i_packet+1;
                end  
          case 101
                Period_s=10;
                %1s->ED1->10s->ED2->10s->ED1->...
                rand_start_ms=1*1000;  
                timestamp_ms=rand_start_ms;
                i_ED=1
                i_packet=1
                while(timestamp_ms<timespan_seconds*1000)
                    UL_TX_start_milliseconds(i_ED,i_packet)=timestamp_ms;
                    timestamp_ms=timestamp_ms+Period_s*1000;
                    if(i_ED<2)
                        i_ED=i_ED+1;
                    else
                        i_ED=1
                        i_packet=i_packet+1;
                    end
                end  
          case num2cell(102:104)
                Period_s=5;
                %1s->ED1->5s->ED2->5s->ED1->...
                rand_start_ms=1*1000;  
                timestamp_ms=rand_start_ms;
                i_ED=1
                i_packet=1
                while(timestamp_ms<timespan_seconds*1000)
                    UL_TX_start_milliseconds(i_ED,i_packet)=timestamp_ms;
                    timestamp_ms=timestamp_ms+Period_s*1000;
                    if(i_ED<2)
                        i_ED=i_ED+1;
                    else
                        i_ED=1
                        i_packet=i_packet+1;
                    end
                end       
          case num2cell(201:202)
                %ED1-ED2->10s->%ED1-ED2->10s->...
                Period_s=10;
                rand_start_ms=rand()*Period_s*1000;  
                i_packet=1;
                timestamp_ms=rand_start_ms;
                while(timestamp_ms<timespan_seconds*1000)
                    for i_ED=1:1:num_EDs
                        UL_TX_start_milliseconds(i_ED,i_packet)=timestamp_ms;
                    end
                    timestamp_ms=timestamp_ms+Period_s*1000;
                    i_packet=i_packet+1;
                end
            case num2cell(301:400)
                min_Period_s=5;
                rand_Period_s=10;
                for i_ED=1:1:num_EDs
                    i_packet=1;
                    rand_start_ms=rand()*rand_Period_s*1000;
                    timestamp_ms=rand_start_ms;
                    while(timestamp_ms<timespan_seconds*1000)
                        UL_TX_start_milliseconds(i_ED,i_packet)=timestamp_ms;
                        timestamp_ms=timestamp_ms+min_Period_s*1000+rand()*rand_Period_s*1000;
                        i_packet=i_packet+1;
                    end
                end             
        end
    end
    
end

