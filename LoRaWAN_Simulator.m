% Authors: Konstantin Mikhaylov, Dr. Tech, Wireless Communications, University of Oulu, Oulu, Finland
%          Riccardo Marini, PhD student, Radio Networks, DEI, University of Bologna, Bologna, Italy
% email address: konstantin.mikhaylov(at)oulu.fi
%                r.marini(at)unibo.it

clearvars;%clear all variables
clc;%clear command window

%SIMULATION PARAMETERS
input_simulation_duration_seconds=600; %duration of the simulated time window
input_draw_playground_map=false; %draw the "map" of the playground with GWs, EDs, and DR zones? true or false
input_all_nodes_use_same_DR=false; %all nodes use the same DR
input_energy_simulate=true; %simulate energy consumption
config_num_iterations=1; %number of simulation iterations, each iteration may have different distributions of ED,GWs and UL packets (depends on the respective files)

%CONFIGS FOR LOOPs
config_loop_GWs=[2]; %number of GWs
config_loop_EDs=[3]; %number of EDs
config_loop_packet_periodicity=[30]; %seconds
config_loop_R=[3]; %circular area radius
config_loop_DR_RX2=[0]; %data rate in RX2
config_loop_RX_priority=[1];

%INPUT PARAMETERS

%APPLICATION PARAMETERS
input_B_UL=20; % uplink packet dimension [bytes]
input_B_DL=20; % downlink payload dimension [bytes]

%TRANSMISSION PARAMETERS
input_sigma=3; %shadowing standard deviation, [dB]
input_UL_tx_power_dBm = [14 12 10 8 6 4 2]; %different power level, sorted in descending order, in dBm, https://lora-alliance.org/sites/default/files/2018-04/lorawantm_regional_parameters_v1.1rb_-_final.pdf  page 17
input_DL_tx_power_GW_dBm = 16; %[dBm]
input_G_GW_UL_dBI=5;%GW RX antenna gain in dBi
input_G_GW_DL_dBI=5;%GW TX antenna gain in dBi 
input_G_ED_UL_dBI=5;%ED TX antenna gain in dBi 
input_G_ED_DL_dBI=5;%ED RX antenna gain in dBi 
input_NF=6; %receiver noise figure, [dB]

%LoRa PARAMETERS
input_CR=1; %Coding rate = 4/CR+4;
input_DR=2; %if all nodes use the same SF
input_BW=125000; %Bandwidth, [Hz]
input_SF=[12 11 10 9 8 7]; %DR0...5 respectively
input_H=1; %0 or 1, presence of header (0 means enabled)
input_preamble_length=8; %symbols
input_SNR_min=[-20, -17.5, -15, -12.5, -10, -7.5]; %[dB] %https://www.rs-online.com/designspark/rel-assets/ds-assets/uploads/knowledge-items/application-notes-for-the-internet-of-things/LoRa%20Design%20Guide.pdf
%interference matrix for intra-/inter-SF interference
input_interference_matrix=[1 -23 -24 -25 -25 -25; -20 1 -20 -21 -22 -22; -18 -17 1 -17 -18 -19; -15 -14 -13 1 -13 -15; -13 -13 -12 -11 1 -11; -9 -9 -9 -9 -8 1];% SINR level to receive the packet, in dB, https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8430542                

%DOWNLINK CONFIGURATION
input_DL_scheduling_probability=1;%probability the NS will generate a downlink in response to an UL
input_RX1_delay_s=1;%delay between end of UL packet to start of RX1 in seconds
input_RX1_DR_shift=0;%difference between DRs used in UL and RX1
input_RX2_delay_s=1;%delay between start of RX1 and start of RX2 in seconds
input_RX2_DR_notless_UL_DR=true;%if set to true in case is RX2 DR is higher than UL, the DL packet is not schedulled [debugging:not checked]
input_DL_tranmission_in_one_RW_only=true;%if set to true in case if the packet is already schedulled for one RX, it will not be schedulled for another one

%DUTY CYCLE CONFIG: no packet will be sent within packet_duration_ms*(100/input_dc_limit_percent) from previous packet start
input_ED_dc_limit_percent=1;%duty cycle limit for ED UL
input_RX1_dc_limit_percent=1;%duty cycle limit for GW RX1 DL
input_RX2_dc_limit_percent=10;%duty cycle limit for GW RX2 DL
input_full_duplex_GW=true;%if set to true, a GW may receive and send packets in RX1 concurrently with no interferences, if false - GW can alternatively receive or send in RX1, and in parallel to this send in RX2

%ENERGY CONSUMPTION PARAMETERS
input_energy_voltage_V=3.3;
input_energy_current_receive_mA = 38; %https://www.readcube.com/articles/10.1007%2Fs10776-017-0341-8?author_access_token=bZ8Zu4O9p5KPvkkwnJr5Hve4RwlQNchNByi7wbcMAY7aktQGE1CpPht2FP-_od4xuWG9k1Yn7J9CKz0tVrPjSfrDAHTWdqoyHsuZreDzXSsYkVqyaoPPcT2zJnT7c9c1B-Yd0dKQdvr2_CG1bW9Ivg==
input_energy_current_transmit_per_DR_mA = [38 35.1 32.4 30 27.5 24.7 22.3]; %order of TX power levels matches that of input_UL_tx_power_dBm, http://ww1.microchip.com/downloads/en/DeviceDoc/RN2483-Low-Power-Long-Range-LoRa-Technology-Transceiver-Module-Data-Sheet-DS50002346D.pdf
input_energy_current_sleep_mA=0.0016; %http://ww1.microchip.com/downloads/en/DeviceDoc/RN2483-Low-Power-Long-Range-LoRa-Technology-Transceiver-Module-Data-Sheet-DS50002346D.pdf
input_energy_current_RXDelay_mA=27; %same for RX_DELAY 1 and 2 %https://www.mdpi.com/1424-8220/17/10/2364
input_energy_RX1_empty_duration_per_DR_ms = [262.14 131.07 98.30 49.15 24.58 12.29]; %order of DRs matches that of input_UL_sensitivity_dBm, %https://www.mdpi.com/1424-8220/17/10/2364
input_energy_RX2_empty_duration_per_DR_ms = [33.02 16.64 8.45 4.35 2.30 1.28]; %order of DRs matches that of input_UL_sensitivity_dBm, %https://www.mdpi.com/1424-8220/17/10/2364

%FUNCTIONS: %user-defined models
input_func_ED_distribution=@km_LoRaWAN_model_ED_distribution;%controls ED placement 
input_func_GW_distribution=@km_LoRaWAN_model_GW_distribution;%controls GW placement
input_func_UL_max_distance=@km_LoRaWAN_model_max_propagation_distance_urban_Hata;%determines max distance reachable in UL for each DR, in meters 
input_func_UL_propagation_model=@km_LoRaWAN_model_propagation_loss_urban_Hata;%propagation of UL signals [debugging:not checked]
input_func_DL_propagation_model=@km_LoRaWAN_model_propagation_loss_urban_Hata;%propagation of DL signals [debugging:not checked]
input_func_GW_to_GW_propagation_model=@km_LoRaWAN_model_propagation_loss_urban_Hata;%propagation of signals between GWs (interference)
input_func_ED_to_ED_propagation_model=@km_LoRaWAN_model_propagation_loss_urban_Hata;%propagation of signals between EDs (interference) 

config_testcase=0; %if non-zero defined - launch the respective test case
config_testcase_transmission=0; %case0 - multiple transmissions per ED, case1 - one transmission per ED

%Some computation fixed for all the simulation cycles
%SENSITIVITY
input_UL_sensitivity_dBm=-174+10*log10(input_BW)+input_NF+input_SNR_min; %Rx sensitivity, dBm
input_DL_sensitivity_dBm=input_UL_sensitivity_dBm;

%TIME ON AIR
Tsym=(2.^input_SF)/input_BW;
Tpreamble=(input_preamble_length+4.25)*Tsym;
%ToA_UL computation [s]
PayloadSymbNb_UL=8+ceil((8*input_B_UL-4*input_SF+28+16-20*input_H)/(4*input_SF))*(input_CR+4);
Tpayload_UL=PayloadSymbNb_UL*Tsym;
ToA_UL=Tpreamble+Tpayload_UL;
input_UL_packets_duration_ms=ToA_UL*1000;
%ToA_DL computation [s]
PayloadSymbNb_DL=8+ceil((8*input_B_DL-4*input_SF+28+16-20*input_H)/(4*input_SF))*(input_CR+4);
Tpayload_DL=PayloadSymbNb_DL*Tsym;
ToA_DL=Tpreamble+Tpayload_DL; %time on air of a packet received with SFi
input_DL_packets_duration_ms=ToA_DL*1000;

%interSF interference: row 446
for config_ctr_loop_RX_priority=1:1:size(config_loop_RX_priority,2)
    input_DL_RW_priority=config_loop_RX_priority(config_ctr_loop_RX_priority);
    for config_ctr_loop_DR_RX2=1:1:size(config_loop_DR_RX2,2)
        input_RX2_DR=config_loop_DR_RX2(config_ctr_loop_DR_RX2);
        for config_ctr_loop_EDs=1:1:size(config_loop_EDs,2)
            input_num_EDs=config_loop_EDs(1,config_ctr_loop_EDs);
            for config_ctr_loop_GWs=1:1:size(config_loop_GWs,2)
                input_num_GWs=config_loop_GWs(1,config_ctr_loop_GWs);
                for config_ctr_loop_R=1:1:size(config_loop_R,2)
                    input_R=config_loop_R(1,config_ctr_loop_R);
                    for config_ctr_loop_T=1:1:size(config_loop_packet_periodicity,2)
                        input_T=config_loop_packet_periodicity(1,config_ctr_loop_T);
                
                clearvars -except input_* config_* simulation
                config_seed=1234578;

                %service
                config_colors={[1 0 0] [1 0.65 0] [1 1 0] [0 1 0] [0 1 1] [0 0 1] [1 0 1]};% red orange yellow green cyan blue violet
                
                %Simulations start here
                config_timestamp_started=datetime('now');
                fprintf('Experiments started at %s\n',config_timestamp_started);
                
                input_class_c_interference_probability=0;
                input_ED_class=zeros(1,config_loop_EDs); %0 - class A /// 1 - class C

                result_iterations_log_DL_RX1_schedulled_by_ED=ones(input_num_EDs,config_num_iterations)*0;
                result_iterations_log_DL_RX1_delivered_by_ED=ones(input_num_EDs,config_num_iterations)*0;
                result_iterations_log_DL_RX2_schedulled_by_ED=ones(input_num_EDs,config_num_iterations)*0;
                result_iterations_log_DL_RX2_delivered_by_ED=ones(input_num_EDs,config_num_iterations)*0;
                
                %Energy consumption initialization
                if input_energy_simulate==true
                    result_iterations_log_consumption_ED_active_time_ms=zeros(input_num_EDs,config_num_iterations);
                    result_iterations_log_consumption_ED_sleep_time_ms=zeros(input_num_EDs,config_num_iterations);
                    result_iterations_log_consumption_ED_UL_consumption_mJ=zeros(input_num_EDs,config_num_iterations);
                    result_iterations_log_consumption_ED_DL_useful_RX_mJ=zeros(input_num_EDs,config_num_iterations);%i.e., leading to packet reception; for both RX1 and RX2
                    result_iterations_log_consumption_ED_DL_idle_RX_mJ=zeros(input_num_EDs,config_num_iterations);%i.e., not leading to a packet reception (for ANY reason: (i) UL not received, (ii) no DL traffic, (iii) DC restriction of GW(s), (iv) DL trasnmitted but not received due to interference, ; for both RX1 and RX2
                    result_iterations_log_consumption_ED_sleep_mJ=zeros(input_num_EDs,config_num_iterations);
                    result_iterations_log_consumption_ED_total_mJ=zeros(input_num_EDs,config_num_iterations);
                end
                
%%%%%%%%%%%%%%%%cycle over iterations
                for config_iteration=1:1:config_num_iterations
                    
                    input_s_ED_GW=randn(input_num_EDs,input_num_GWs)*input_sigma; %shadowing samples
                    input_s_GW_GW=randn(input_num_GWs,input_num_GWs)*input_sigma; %shadowing samples
                    input_s_ED_ED=randn(input_num_EDs,input_num_EDs)*input_sigma; %shadowing samples
                    
                    calc_simulation_end_ms=input_simulation_duration_seconds*1000;
                    if input_draw_playground_map==true
                        clf
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation

                    calc_num_DRs=size(input_UL_sensitivity_dBm,2);
                    calc_num_UL_TXPower_levels=size(input_UL_tx_power_dBm,2);
                    calc_max_UL_tx_power_dBm=max(input_UL_tx_power_dBm);
                    
                    %calculate the maximum distance per SF
                    calc_max_UL_budget_dB=-input_UL_sensitivity_dBm+input_G_GW_UL_dBI+input_G_ED_UL_dBI+calc_max_UL_tx_power_dBm;
                    calc_max_UL_range_meters=input_func_UL_max_distance(calc_max_UL_budget_dB);
                                    
%%%%%%%%%%%%%%%%%%%%step 1: distribute the EDs (based on the desired rule and parameters defined in the function to be called)
                    [calc_coord_GW_cartesian,calc_coord_GW_polar]=input_func_GW_distribution(input_num_GWs,config_testcase); 
                    [calc_coord_ED_cartesian,calc_coord_ED_polar]=input_func_ED_distribution(input_num_EDs,config_testcase,input_R); 
                    
%%%%%%%%%%%%%%%%%%%%step 2: calculate the distances 
                    %a: from each ED to each GW, as well as the minimum distance and the index of the GW nearest to an ED
                    calc_distance_ED_GW_m=zeros(input_num_EDs,input_num_GWs);
                    for temp_i_num_GW=1:input_num_GWs
                       for temp_i_num_ED=1:input_num_EDs 
                          calc_distance_ED_GW_m(temp_i_num_ED,temp_i_num_GW)=pdist2([calc_coord_GW_cartesian(temp_i_num_GW,1) calc_coord_GW_cartesian(temp_i_num_GW,2)],[calc_coord_ED_cartesian(temp_i_num_ED,1) calc_coord_ED_cartesian(temp_i_num_ED,2)],'euclidean');
                       end    
                    end
                    %specify the minimum distance between a GW and an ED
                    [calc_distance_min_ED_GW_m, calc_UL_nearest_GW_to_ED_index]=min(calc_distance_ED_GW_m,[],2);
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation
                    %b: calculate the inter-ED distances for downlink
                    clearvars calc_distance_interED_m 
                    calc_distance_interED_m = ones(input_num_EDs,input_num_EDs);
                    for i=1:input_num_EDs
                        for j=1:input_num_EDs
                            if j==i
                                calc_distance_interED_m(i,j)=-1;
                            else
                                %between EDs
                                clearvars temp_dist_vector temp_dist_d 
                                temp_dist_vector_m=[calc_coord_ED_cartesian(i,:);calc_coord_ED_cartesian(j,:)];
                                temp_pair_dist_m = pdist(temp_dist_vector_m,'euclidean');
                                calc_distance_interED_m(i,j)=temp_pair_dist_m;
                            end        
                        end
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation
                    %c: also calculate the distances between the GWs
                    calc_distance_GW_GW_m=zeros(input_num_GWs,input_num_GWs);
                    for temp_i_num_GW1=1:input_num_GWs
                       for temp_i_num_GW2=1:input_num_GWs
                          calc_distance_GW_GW_m(temp_i_num_GW1,temp_i_num_GW2)=pdist2([calc_coord_GW_cartesian(temp_i_num_GW1,1) calc_coord_GW_cartesian(temp_i_num_GW1,2)],[calc_coord_GW_cartesian(temp_i_num_GW2,1) calc_coord_GW_cartesian(temp_i_num_GW2,2)],'euclidean');
                       end    
                    end  
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation

%%%%%%%%%%%%%%%%%%%%%step 3: calculate the attenuation in UL between each ED and each
                    calc_attenuation_UL_ED_GW_dB=zeros(input_num_EDs,input_num_GWs);
                    for temp_i_num_GW=1:input_num_GWs
                       for temp_i_num_ED=1:input_num_EDs 
                           calc_attenuation_UL_ED_GW_dB(temp_i_num_ED,temp_i_num_GW)=input_func_UL_propagation_model(calc_distance_ED_GW_m(temp_i_num_ED,temp_i_num_GW))-input_G_GW_UL_dBI-input_G_ED_UL_dBI;
                       end    
                    end
                    %specify the minimum attenuation between a GW and an ED (indexes same as distance, but for future extension)
                    [calc_UL_attenuation_min_ED_GW_dB, calc_UL_attenuation_min_GW_to_ED_index]=min(calc_attenuation_UL_ED_GW_dB,[],2);    
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation

%%%%%%%%%%%%%%%%%%%%%step 4: standard ADR - optimize the DRs and TX power of EDs in UL
                    clearvars calc_DR_matrix calc_TX_Power_matrix calc_channel_UL 
                    calc_DRs_UL = zeros(input_num_EDs,1);
                    calc_UL_TX_power_matrix_dBm = ones(input_num_EDs,1)*-1;
                    
                    %compute maximum DR
                    for i=1:input_num_EDs
                        %maximum DR
                        for DRi=calc_num_DRs-1:-1:0
                           if abs(input_UL_sensitivity_dBm(DRi+1))>calc_UL_attenuation_min_ED_GW_dB(i)-calc_max_UL_tx_power_dBm
                                calc_DRs_UL_max(i)=DRi;
                                break;
                           end
                        end
                        calc_UL_TX_power_matrix_dBm(i,1)=max(input_UL_tx_power_dBm);
                    end
                    %ADR
                    for temp_i_num_ED=1:input_num_EDs
                        s_vector=randn(1,20)*input_sigma;
                        Pr_vector=calc_max_UL_tx_power_dBm-calc_UL_attenuation_min_ED_GW_dB(temp_i_num_ED)-s_vector;
                        SNR_vector=Pr_vector+174-10*log10(input_BW)-input_NF;
                        SNR_max=max(SNR_vector);
                        SNR_margin=SNR_max-input_SNR_min(1)-10; %margin=10dB
                        Nstep=round(SNR_margin/3);
                        
                        while (Nstep>0) && (calc_DRs_UL(temp_i_num_ED)<calc_DRs_UL_max(temp_i_num_ED))
                            calc_DRs_UL(temp_i_num_ED)=calc_DRs_UL(temp_i_num_ED)+1;
                            Nstep=Nstep-1;
                        end
                        while (Nstep>0) && (calc_UL_TX_power_matrix_dBm(temp_i_num_ED,1)>min(input_UL_tx_power_dBm))
                            calc_UL_TX_power_matrix_dBm(temp_i_num_ED,1)=calc_UL_TX_power_matrix_dBm(temp_i_num_ED,1)-2;
                            Nstep=Nstep-1;
                        end
                        while (Nstep<0) && (calc_UL_TX_power_matrix_dBm(temp_i_num_ED,1)<max(input_UL_tx_power_dBm))
                            calc_UL_TX_power_matrix_dBm(temp_i_num_ED,1)=calc_UL_TX_power_matrix_dBm(temp_i_num_ED,1)+2;
                            Nstep=Nstep+1;
                        end
                        
                        calc_DL_TX_power_matrix_dBm(temp_i_num_ED,1)=input_DL_tx_power_GW_dBm;
                        calc_DRs_RX2(temp_i_num_ED,1)=input_RX2_DR;
                    end
                    calc_DRs_RX1=calc_DRs_UL;

                    %match DR to SF
                    for i=1:input_num_EDs
                    %UL and RX1
                    if(calc_DRs_UL(i)==0)
                        calc_SF_UL(i)=12;
                        calc_SF_DL_RX1(i)=12;
                    elseif(calc_DRs_UL(i)==1)
                        calc_SF_UL(i)=11;
                        calc_SF_DL_RX1(i)=11;
                    elseif(calc_DRs_UL(i)==2)
                        calc_SF_UL(i)=10;
                        calc_SF_DL_RX1(i)=10;
                    elseif(calc_DRs_UL(i)==3)
                        calc_SF_UL(i)=9;
                        calc_SF_DL_RX1(i)=9;
                    elseif(calc_DRs_UL(i)==4)
                        calc_SF_UL(i)=8;
                        calc_SF_DL_RX1(i)=8;
                    elseif(calc_DRs_UL(i)==5)
                        calc_SF_UL(i)=7;
                        calc_SF_DL_RX1(i)=7;
                    end
                    
                    %RX2
                  	if(calc_DRs_RX2(i)==0)
                        calc_SF_DL_RX2(i)=12;
                    elseif(calc_DRs_RX2(i)==1)
                        calc_SF_DL_RX2(i)=11;
                    elseif(calc_DRs_RX2(i)==2)
                        calc_SF_DL_RX2(i)=10;
                    elseif(calc_DRs_RX2(i)==3)
                        calc_SF_DL_RX2(i)=9;
                    elseif(calc_DRs_UL(i)==4)
                        calc_SF_DL_RX2(i)=8;
                    elseif(calc_DRs_UL(i)==5)
                        calc_SF_DL_RX2(i)=7;
                    end
                    end

                    %everyone use the same SF - skip if false
                    if input_all_nodes_use_same_DR==true
                       for i=1:input_num_EDs
                           calc_DRs_UL(i)=input_DR;
                           calc_DRs_RX2(i)=input_RX2_DR;
                       end
                       calc_DRs_RX1=calc_DRs_UL;                
                    end

                    %calculate RSSI and SNR for UL signal 
                    for temp_i_num_GW=1:input_num_GWs
                       for temp_i_num_ED=1:input_num_EDs 
                           calc_RSSI_UL_ED_GW_dBm(temp_i_num_ED,temp_i_num_GW)=calc_UL_TX_power_matrix_dBm(temp_i_num_ED,1)-calc_attenuation_UL_ED_GW_dB(temp_i_num_ED,temp_i_num_GW)-input_s_ED_GW(temp_i_num_ED,temp_i_num_GW);
                           calc_SNR_UL(temp_i_num_ED,temp_i_num_GW)=calc_RSSI_UL_ED_GW_dBm(temp_i_num_ED,temp_i_num_GW)+174-10*log10(input_BW)-input_NF;
                       end    
                    end
                    %specify the maximum SNR between a GW and an ED (indexes same as distance, but for future extension)
                    [calc_UL_RSSI_max_ED_GW_dB, calc_UL_RSSI_min_GW_to_ED_index]=max(calc_RSSI_UL_ED_GW_dBm,[],2);  
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation    
                    %calculate the resulting SNR in dB for UL of each ED to each GW
                    for temp_i_num_GW=1:input_num_GWs
                       for temp_i_num_ED=1:input_num_EDs 
                           calc_sensitivity_gap_UL_ED_GW_dB(temp_i_num_ED,temp_i_num_GW)=calc_RSSI_UL_ED_GW_dBm(temp_i_num_ED,temp_i_num_GW)-input_UL_sensitivity_dBm(calc_DRs_UL(temp_i_num_ED,1)+1);
                       end    
                    end
                    %specify the minimum SNR between a GW and an ED (indexes same as distance, but for future extension)
                    [calc_UL_SNR_max_ED_GW_dB, calc_UL_SNR_min_GW_to_ED_index]=max(calc_sensitivity_gap_UL_ED_GW_dB,[],2);  
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation

%%%%%%%%%%%%%%%%%%%%%step 5: calculate the other attenuations
                    %a: downlink attenuation between each GW and each ED
                    for temp_i_num_GW=1:input_num_GWs
                       for temp_i_num_ED=1:input_num_EDs 
                           calc_attenuation_DL_ED_GW_dB(temp_i_num_ED,temp_i_num_GW)=input_func_DL_propagation_model(calc_distance_ED_GW_m(temp_i_num_ED,temp_i_num_GW))-input_G_GW_DL_dBI-input_G_ED_DL_dBI+input_s_ED_GW(temp_i_num_ED,temp_i_num_GW);
                       end    
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation
                    %b: downlink attenuation between the GWs
                    for temp_i_num_GW1=1:input_num_GWs
                       for temp_i_num_GW2=1:input_num_GWs 
                           if temp_i_num_GW1~=temp_i_num_GW2
                               calc_attenuation_GW_GW_dB(temp_i_num_GW1,temp_i_num_GW2)=input_func_GW_to_GW_propagation_model(calc_distance_GW_GW_m(temp_i_num_GW1,temp_i_num_GW2))-input_G_GW_DL_dBI-input_G_GW_UL_dBI+input_s_GW_GW(temp_i_num_GW1,temp_i_num_GW2);
                           else
                               calc_attenuation_GW_GW_dB(temp_i_num_GW1,temp_i_num_GW2)=0;
                           end
                       end    
                    end  
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation
                    %c: between the EDs
                    for temp_i_num_ED1=1:input_num_EDs 
                        for temp_i_num_ED2=1:input_num_EDs 
                            if temp_i_num_ED1~=temp_i_num_ED2
                               calc_attenuation_ED_ED_dB(temp_i_num_ED1,temp_i_num_ED2)=input_func_ED_to_ED_propagation_model(calc_distance_interED_m(temp_i_num_ED1,temp_i_num_ED2))-input_G_ED_UL_dBI-input_G_ED_DL_dBI+input_s_ED_ED(temp_i_num_ED1,temp_i_num_ED2);
                            else
                               calc_attenuation_ED_ED_dB(temp_i_num_ED1,temp_i_num_ED2)=0;
                            end
                        end
                    end 

%%%%%%%%%%%%%%%%%%%%%step 6: calculate RSSIs
                    %a: GW to ED
                    for temp_i_num_GW=1:input_num_GWs
                       for temp_i_num_ED=1:input_num_EDs 
                           calc_RSSI_DL_ED_GW_dBm(temp_i_num_ED,temp_i_num_GW)=-calc_attenuation_DL_ED_GW_dB(temp_i_num_ED,temp_i_num_GW)+calc_DL_TX_power_matrix_dBm(temp_i_num_ED,1);
                       end    
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_*  simulation    
                    %b: GW to GW
                    for temp_i_num_GW1_target=1:input_num_GWs
                       for temp_i_num_GW2_interferer=1:input_num_GWs 
                           if temp_i_num_GW1_target~=temp_i_num_GW2_interferer
                               calc_RSSI_GW_GW_dBm(temp_i_num_GW1_target,temp_i_num_GW2_interferer)=-calc_attenuation_GW_GW_dB(temp_i_num_GW1_target,temp_i_num_GW2_interferer)+input_DL_tx_power_GW_dBm;
                           else
                               calc_RSSI_GW_GW_dBm(temp_i_num_GW1_target,temp_i_num_GW2_interferer)=0;
                           end
                       end    
                    end  
                    clearvars -except calc_* input_* result_* DEF_* config_*  simulation  
                    %c: ED to ED
                    for temp_i_num_ED_target=1:input_num_EDs 
                        for temp_i_num_ED_interferer=1:input_num_EDs 
                            if temp_i_num_ED_target~=temp_i_num_ED_interferer
                               calc_RSSI_ED_ED_dB(temp_i_num_ED_target,temp_i_num_ED_interferer)=-calc_attenuation_ED_ED_dB(temp_i_num_ED_target,temp_i_num_ED_interferer)+calc_UL_TX_power_matrix_dBm(temp_i_num_ED_interferer);
                            else
                               calc_RSSI_ED_ED_dB(temp_i_num_ED_target,temp_i_num_ED_interferer)=0;
                            end
                        end
                    end  
                    
                    calc_sensitivity_gap_DL_ED_GW_RX1_dB=zeros(input_num_EDs,input_num_GWs);
                    calc_sensitivity_gap_DL_ED_GW_RX2_dB=zeros(input_num_EDs,input_num_GWs);
                    for temp_i_num_GW=1:input_num_GWs
                       for temp_i_num_ED=1:input_num_EDs 
                           calc_sensitivity_gap_DL_ED_GW_RX1_dB(temp_i_num_ED,temp_i_num_GW)=calc_RSSI_DL_ED_GW_dBm(temp_i_num_ED,temp_i_num_GW)-input_DL_sensitivity_dBm(calc_DRs_RX1(temp_i_num_ED,1)+1);
                           calc_sensitivity_gap_DL_ED_GW_RX2_dB(temp_i_num_ED,temp_i_num_GW)=calc_RSSI_DL_ED_GW_dBm(temp_i_num_ED,temp_i_num_GW)-input_DL_sensitivity_dBm(calc_DRs_RX2(temp_i_num_ED,1)+1);
                           calc_SNR_DL(temp_i_num_ED,temp_i_num_GW)=calc_RSSI_DL_ED_GW_dBm(temp_i_num_ED,temp_i_num_GW)+174-10*log10(input_BW)-input_NF;

                       end    
                    end  
                                                           
%%%%%%%%%%%%%%%%%%%%%step 7: make the interference matrixes
                    %a: UL vs UL
                    calc_UL_EDs_RX_levels=zeros(input_num_EDs,input_num_EDs,input_num_GWs);
                    calc_UL_EDs_blockers=zeros(input_num_EDs,input_num_EDs,input_num_GWs);
                    calc_UL_EDs_blockers_interSF=zeros(input_num_EDs,input_num_EDs,input_num_GWs);
                    calc_UL_EDs_blockers_intraSF=zeros(input_num_EDs,input_num_EDs,input_num_GWs);
                    for i_GW=1:input_num_GWs
                        for ED_target=1:input_num_EDs
                            for ED_interferer=1:input_num_EDs
                                if ED_interferer~=ED_target
                                    %check if this ED_interferer can prevent GW from receiving the
                                    %signal of ED_target
                                    ED_target_DR=calc_DRs_UL(ED_target,1);
                                    ED_interferer_DR=calc_DRs_UL(ED_interferer,1);
                                    GW_ED_target_RSSI_dBm=calc_RSSI_UL_ED_GW_dBm(ED_target,i_GW);
                                    GW_ED_interferer_RSSI_dBm=calc_RSSI_UL_ED_GW_dBm(ED_interferer,i_GW);
                                    calc_UL_EDs_RX_levels(ED_target,ED_interferer,i_GW)=GW_ED_target_RSSI_dBm-GW_ED_interferer_RSSI_dBm;
                                    if (GW_ED_target_RSSI_dBm-GW_ED_interferer_RSSI_dBm)<input_interference_matrix(ED_target_DR+1,ED_interferer_DR+1)
                                        %calc_UL_EDs_blockers(ED_target,ED_interferer,i_GW)=1;
                                        if(ED_target_DR==ED_interferer_DR)
                                            calc_UL_EDs_blockers_intraSF(ED_target,ED_interferer,i_GW)=1;
                                            calc_UL_EDs_blockers(ED_target,ED_interferer,i_GW)=1; %same SF
                                        else
                                            calc_UL_EDs_blockers_interSF(ED_target,ED_interferer,i_GW)=1;
                                            calc_UL_EDs_blockers(ED_target,ED_interferer,i_GW)=1; %different SF
                                        end
                                    end
                                end
                            end%
                            clearvars ED_interferer 
                        end%
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation
                    %b: showing which GW DL block GW from receiving UL
                    
                    calc_UL_GW_ED_GWint_SF_blockers=zeros(input_num_GWs,input_num_EDs,input_num_GWs,calc_num_DRs);%1:receiving GW,2:ED of whom UL trying to receive,3:interfereing GW,4:DR of interfering GW
                    for i_GW_target=1:input_num_GWs
                        for ED_target=1:input_num_EDs
                            for i_GW_interferer=1:input_num_GWs
                                if i_GW_target~=i_GW_interferer
                                    %define the attenuations
                                    ED_target_DR=calc_DRs_UL(ED_target,1);
                                    GW_ED_target_RSSI_dBm=calc_RSSI_UL_ED_GW_dBm(ED_target,i_GW_target);
                                    GW_GW_interferer_RSSI_dBm=calc_RSSI_GW_GW_dBm(i_GW_target,i_GW_interferer);
                                    temp_SNR_dB=GW_ED_target_RSSI_dBm-GW_GW_interferer_RSSI_dBm;
                                    for i_DR_interferer=0:calc_num_DRs-1
                                        if(temp_SNR_dB<input_interference_matrix(ED_target_DR+1,i_DR_interferer+1))
                                            calc_UL_GW_ED_GWint_SF_blockers(i_GW_target,ED_target,i_GW_interferer,i_DR_interferer+1)=1;
                                        end
                                    end
                                else
                                    for i_DR_interferer=0:calc_num_DRs-1
                                        calc_UL_GW_ED_GWint_SF_blockers(i_GW_target,ED_target,i_GW_interferer,i_DR_interferer+1)=1;
                                    end
                                end    
                            end
                        end
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation
                    %c: UL interferer vs DL target
                    calc_DL_GW_ED_EDint_SF_blockers=zeros(input_num_GWs,input_num_EDs,input_num_EDs,calc_num_DRs);%1:receiving GW,2:ED of whom UL trying to receive,3:interfereing GW,4:DR of interfering GW
                    calc_DL_GW_ED_GWint_SF_blockers=zeros(input_num_GWs,input_num_EDs,input_num_EDs,calc_num_DRs);
                    for i_GW=1:input_num_GWs
                        for ED_target=1:input_num_EDs
                            for ED_interferer=1:input_num_EDs
                                if ED_interferer~=ED_target
                                    %DRs
                                    ED_interferer_DR=calc_DRs_UL(ED_interferer,1);
                                    %RSSI
                                    ED_DL_RSSI_dBm=calc_RSSI_DL_ED_GW_dBm(ED_target,i_GW);
                                    ED_interference_RSSI_dBm=calc_RSSI_ED_ED_dB(ED_target,ED_interferer);
                                    temp_SNR_dB=ED_DL_RSSI_dBm-ED_interference_RSSI_dBm;
                                    for i_DR_target=0:calc_num_DRs-1
                                        if(temp_SNR_dB<input_interference_matrix(i_DR_target+1,ED_interferer_DR+1))
                                            calc_DL_GW_ED_GWint_SF_blockers(i_GW,ED_target,ED_interferer,i_DR_target+1)=1;
                                        end
                                    end                    
                                else
                                     for i_DR_target=0:calc_num_DRs-1
                                        calc_DL_GW_ED_GWint_SF_blockers(i_GW,ED_target,ED_interferer,i_DR_target+1)=1;
                                    end                   
                                end
                            end
                            clearvars ED_interferer 
                        end
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation
                    %d: DL interferer vs DL target in RX1
                    calc_DL_RX1_GWtarg_ED_GWint_SF_blockers=zeros(input_num_GWs,input_num_EDs,input_num_GWs,calc_num_DRs);%1:receiving GW,2:ED of whom UL trying to receive,3:interfereing GW,4:DR of interfering GW
                    for ED_target=1:input_num_EDs
                        ED_target_DR=calc_DRs_RX1(ED_target,1);
                        for i_GW_target=1:input_num_GWs
                            GW_DL_target_RSSI_dBm=calc_RSSI_DL_ED_GW_dBm(ED_target,i_GW_target);
                            for i_GW_interferer=1:input_num_GWs
                                if i_GW_interferer~=i_GW_target
                                    GW_DL_interferer_RSSI_dBm=calc_RSSI_DL_ED_GW_dBm(ED_target,i_GW_interferer);
                                    temp_SNR_dB=GW_DL_target_RSSI_dBm-GW_DL_interferer_RSSI_dBm;
                                    for GW_interferer_DR=0:calc_num_DRs-1
                                        if(temp_SNR_dB<input_interference_matrix(ED_target_DR+1,GW_interferer_DR+1))
                                            calc_DL_RX1_GWtarg_ED_GWint_SF_blockers(i_GW_target,ED_target,i_GW_interferer,GW_interferer_DR+1)=1;
                                        end                   
                                    end
                                else
                                     for GW_interferer_DR=0:calc_num_DRs-1
                                        calc_DL_RX1_GWtarg_ED_GWint_SF_blockers(i_GW_target,ED_target,i_GW_interferer,GW_interferer_DR+1)=1;                  
                                    end                   
                                end
                            end%
                        end
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_*  simulation   
                    %e: DL interferer vs DL target in RX2
                    calc_DL_RX2_GWtarg_ED_GWint_blockers=zeros(input_num_GWs,input_num_EDs,input_num_GWs);%1:receiving GW,2:ED of whom UL trying to receive,3:interfereing GW
                     for ED_target=1:input_num_EDs
                        ED_target_DR=calc_DRs_RX2(ED_target,1);
                        for i_GW_target=1:input_num_GWs
                            GW_DL_target_RSSI_dBm=calc_RSSI_DL_ED_GW_dBm(ED_target,i_GW_target);
                            for i_GW_interferer=1:input_num_GWs
                                if i_GW_interferer~=i_GW_target
                                    GW_DL_interferer_RSSI_dBm=calc_RSSI_DL_ED_GW_dBm(ED_target,i_GW_interferer);
                                    temp_SNR_dB=GW_DL_target_RSSI_dBm-GW_DL_interferer_RSSI_dBm;
                                    GW_interferer_DR=input_RX2_DR;
                                    if(temp_SNR_dB<input_interference_matrix(ED_target_DR+1,GW_interferer_DR+1))
                                        calc_DL_RX2_GWtarg_ED_GWint_blockers(i_GW_target,ED_target,i_GW_interferer)=1;
                                    end  
                                else
                                     for GW_interferer_DR=0:calc_num_DRs-1
                                        calc_DL_RX2_GWtarg_ED_GWint_blockers(i_GW_target,ED_target,i_GW_interferer)=1;                  
                                    end                   
                                end
                            end%
                        end
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation    

%%%%%%%%%%%%%%%%%%%%%step 8: fill the transmit matrix showing when each tranmission starts
                    %and ends, accounting for UL DC limitation
                    calc_UL_traffic_packets_start_ms=km_LoRaWAN_model_UL_transmissions(config_testcase_transmission,input_num_EDs,input_simulation_duration_seconds,input_T);
                    %now need to check duty cycle restrictions and specify the timestamp of
                    %packet end
                    for temp_ED_i=1:1:input_num_EDs
                        EDi_DR = calc_DRs_UL(temp_ED_i,1);
                        clearvars  temp_ED_UL_traffic_packets_start_ms_nonzeros simulation
                        temp_ED_UL_traffic_packets_start_ms_nonzeros=nonzeros(calc_UL_traffic_packets_start_ms(temp_ED_i,:)');
                        temp_release_DC_timestamp=-1;
                        temp_scheduled_packet=1;
                        temp_UL_packet_duration_ms=input_UL_packets_duration_ms(EDi_DR+1);
                        for temp_ED_packet_i=1:1:size(temp_ED_UL_traffic_packets_start_ms_nonzeros,1)
                            temp_UL_packets_schedulled_start_ms=temp_ED_UL_traffic_packets_start_ms_nonzeros(temp_ED_packet_i,1);
                            if((temp_UL_packets_schedulled_start_ms>=temp_release_DC_timestamp)&&(temp_UL_packets_schedulled_start_ms+temp_UL_packet_duration_ms<(input_simulation_duration_seconds*1000)))%DC check ok
                                calc_UL_packets_scheduled_start_ms(temp_ED_i,temp_scheduled_packet)=temp_UL_packets_schedulled_start_ms; %#ok<SAGROW>
                                calc_UL_packets_scheduled_end_ms(temp_ED_i,temp_scheduled_packet)=temp_UL_packets_schedulled_start_ms+temp_UL_packet_duration_ms; %#ok<SAGROW>
                                temp_release_DC_timestamp=temp_UL_packets_schedulled_start_ms+temp_UL_packet_duration_ms*(100/input_ED_dc_limit_percent);
                                temp_scheduled_packet=temp_scheduled_packet+1;
                            else
                                breakpoint=1;
                            end
                        end
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation

%%%%%%%%%%%%%%%%%%%%%step 9: crawl through schedulled UL tranmissions
                    %           make decision which UL are received
                    %           schedule downlinks whenever needed
                    %sort all scheduled UL based on TX start time
                    temp_UL_sorted_packet_starts=sort(nonzeros(calc_UL_packets_scheduled_start_ms(:)));
                    calc_UL_num_packets=size(temp_UL_sorted_packet_starts,1);
                    %UL reception result matrixes
                    %for each packet shows if it has been received by each GW
                    result_UL_packet_delivery_matrix=zeros(size(calc_UL_packets_scheduled_start_ms,1),size(calc_UL_packets_scheduled_start_ms,2),input_num_GWs);
                    %number of GWs receiving each packet
                    result_UL_packet_delivery_matrix_num_GWs=zeros(size(calc_UL_packets_scheduled_start_ms,1),size(calc_UL_packets_scheduled_start_ms,2));
                    %extra resuls: reasons for each packet not received by each GW
                    result_extra_UL_packet_undelivery_matrix_cause_SNR=zeros(size(calc_UL_packets_scheduled_start_ms,1),size(calc_UL_packets_scheduled_start_ms,2),input_num_GWs);
                    result_extra_UL_packet_undelivery_matrix_cause_UL_interference=zeros(size(calc_UL_packets_scheduled_start_ms,1),size(calc_UL_packets_scheduled_start_ms,2),input_num_GWs);
                    result_extra_UL_packet_undelivery_matrix_cause_DL_interference=zeros(size(calc_UL_packets_scheduled_start_ms,1),size(calc_UL_packets_scheduled_start_ms,2),input_num_GWs); 
                    %extra results: reasons for DL packet not been schedulled for each UL
                    result_extra_DL_packet_unscheduled_matrix_cause_no_traffic=zeros(size(calc_UL_packets_scheduled_start_ms,1),size(calc_UL_packets_scheduled_start_ms,2),input_num_GWs);
                    result_extra_DL_packet_unscheduled_matrix_cause_RX1_backoff=zeros(size(calc_UL_packets_scheduled_start_ms,1),size(calc_UL_packets_scheduled_start_ms,2),input_num_GWs);
                    result_extra_DL_packet_unscheduled_matrix_cause_RX2_backoff=zeros(size(calc_UL_packets_scheduled_start_ms,1),size(calc_UL_packets_scheduled_start_ms,2),input_num_GWs);
                    result_extra_DL_packet_unscheduled_matrix_cause_RX2_low_DR=zeros(size(calc_UL_packets_scheduled_start_ms,1),size(calc_UL_packets_scheduled_start_ms,2),input_num_GWs);
                    result_extra_DL_packet_unscheduled_matrix_cause_planned_already=zeros(size(calc_UL_packets_scheduled_start_ms,1),size(calc_UL_packets_scheduled_start_ms,2),input_num_GWs);
                    %crawl through the UL packets
                    temp_previous_packet_start_timestamp_ms=-1;
                    temp_duplicate=0;%handles the case when multiple packets start exactly at the same time
                    calc_DL_packets_RX1_scheduled_DR=[];%empty matrix
                    calc_DL_packets_RX1_scheduled_start_ms=[];%double.empty(1,1,0);%empty matrix
                    calc_DL_packets_RX1_scheduled_ED_target=[];%empty matrix
                    calc_DL_packets_RX1_scheduled_end_ms=[];%empty matrix
                    calc_DL_packets_RX2_scheduled_DR=[];%empty matrix
                    calc_DL_packets_RX2_scheduled_start_ms=[];%empty matrix
                    calc_DL_packets_RX2_scheduled_end_ms=[];%empty matrix
                    calc_DL_packets_RX2_scheduled_ED_target=[];%empty matrix
                    calc_DL_channel_release_timestamp_ms=ones(2,input_num_GWs)*(-1);%2 - RX1 and RX2
                    calc_DL_GW_packet_RX1_scheduled=zeros(input_num_GWs,1);
                    calc_DL_GW_packet_RX2_scheduled=zeros(input_num_GWs,1);
                    for temp_i_packet=1:1:calc_UL_num_packets%go through all UL packets
                        temp_packet_start_timestamp_ms=temp_UL_sorted_packet_starts(temp_i_packet,1);
                        if(temp_packet_start_timestamp_ms==temp_previous_packet_start_timestamp_ms)
                            temp_duplicate=temp_duplicate+1;
                        else
                            temp_duplicate=0;
                        end
                        clearvars temp_found_ED temp_found_ED_packet_index
                        clearvars temp_packets_UL_overlapping_partially1_indexes temp_packets_UL_overlapping_partially2_indexes temp_packets_UL_overlapping_fully_indexes 
                        clearvars temp_packets_UL_overlapping_all_indexes 
                        [temp_found_ED, temp_found_ED_packet_index] = find(calc_UL_packets_scheduled_start_ms==temp_packet_start_timestamp_ms);
                        temp_ED=temp_found_ED(temp_duplicate+1,1);
                        temp_ED_packet_index=temp_found_ED_packet_index(temp_duplicate+1,1);
                        temp_UL_packet_end_time_ms=calc_UL_packets_scheduled_end_ms(temp_ED,temp_ED_packet_index);

                        if 1%process UL interferences
                            %find all packets starting or ending between
                            %temp_packet_start_timestamp_ms and temp_UL_packet_end_time_ms  (NOTE: includes also self)
                            [temp_packets_UL_overlapping_partially1_indexes(:,1),temp_packets_UL_overlapping_partially1_indexes(:,2)] = find(calc_UL_packets_scheduled_start_ms >= temp_packet_start_timestamp_ms & calc_UL_packets_scheduled_start_ms <= temp_UL_packet_end_time_ms);
                            [temp_packets_UL_overlapping_partially2_indexes(:,1), temp_packets_UL_overlapping_partially2_indexes(:,2)] = find(calc_UL_packets_scheduled_end_ms >= temp_packet_start_timestamp_ms & calc_UL_packets_scheduled_end_ms <= temp_UL_packet_end_time_ms);
                            %find all packets starting before temp_packet_start_timestamp_ms 
                            %and ending after temp_UL_packet_end_time_ms (NOTE: includes also self)         
                            [temp_packets_UL_overlapping_fully_indexes(:,1), temp_packets_UL_overlapping_fully_indexes(:,2)] = find(calc_UL_packets_scheduled_start_ms <= temp_packet_start_timestamp_ms & calc_UL_packets_scheduled_end_ms >= temp_UL_packet_end_time_ms);
                            %combine all these and remove duplicates
                            temp_packets_UL_overlapping_all_indexes=setdiff(unique(vertcat(temp_packets_UL_overlapping_partially1_indexes,temp_packets_UL_overlapping_partially2_indexes,temp_packets_UL_overlapping_fully_indexes),'rows'),[temp_ED temp_ED_packet_index],'rows');
                            if(isempty(temp_packets_UL_overlapping_all_indexes)==true)
                                %no overlapping UL packets
                                for temp_i_num_GW=1:input_num_GWs%for each GW 
                                    %see if sensitivity level is exceeded
                                    if calc_sensitivity_gap_UL_ED_GW_dB(temp_ED,temp_i_num_GW)>0 %packet received
                                        result_UL_packet_delivery_matrix(temp_ED,temp_ED_packet_index,temp_i_num_GW)=true;
                                    else
                                        result_extra_UL_packet_undelivery_matrix_cause_SNR(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                    end
                                end
                            else
                                %doing otherwise around - first mark for all GWs as received,
                                %then unmark the ones which get blocked
                                for temp_i_num_GW=1:input_num_GWs%for each GW 
                                    %see if sensitivity level is exceeded                                   
                                    if calc_sensitivity_gap_UL_ED_GW_dB(temp_ED,temp_i_num_GW)>0 %packet received
                                        result_UL_packet_delivery_matrix(temp_ED,temp_ED_packet_index,temp_i_num_GW)=true;
                                    else
                                        result_extra_UL_packet_undelivery_matrix_cause_SNR(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                    end
                                end           
                                %check for the UL interferences and mark if any blocks the UL
                                for temp_i_interf_packet=1:1:size(temp_packets_UL_overlapping_all_indexes,1)
                                    temp_ED_interferer=temp_packets_UL_overlapping_all_indexes(temp_i_interf_packet,1);
                                    for temp_i_num_GW=1:input_num_GWs
                                        %next row should be without intraSF
                                        %to consider all the EDs
                                        if calc_UL_EDs_blockers_intraSF(temp_ED,temp_ED_interferer,temp_i_num_GW)==1
                                            result_UL_packet_delivery_matrix(temp_ED,temp_ED_packet_index,temp_i_num_GW)=false;
                                            result_extra_UL_packet_undelivery_matrix_cause_UL_interference(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                        end
                                    end
                                end
                            end
                        end
                        clearvars temp_packets_DL_overlapping_partially1_indexes temp_packets_DL_overlapping_partially2_indexes temp_packets_DL_overlapping_fully_indexes 
                        clearvars temp_packets_DL_overlapping_all_indexes 
                        if 1%process DL interferences
                             if(isempty(calc_DL_packets_RX1_scheduled_start_ms)==false)
                                %clearvars calc_DL_packets_RX1_scheduled_start_ms_2d
                                %calc_DL_packets_RX1_scheduled_start_ms_2d=reshape(calc_DL_packets_RX1_scheduled_start_ms,2,[])
                                %find all DL packets starting or ending between
                                %temp_packet_start_timestamp_ms and temp_UL_packet_end_time_ms  (NOTE: includes also self)
                                temp_packets_DL_overlapping_all_indexes=[];
                                try
                                    [temp_packets_DL_overlapping_partially1_indexes(:,1), temp_packets_DL_overlapping_partially1_indexes(:,2)] = find(calc_DL_packets_RX1_scheduled_start_ms >= temp_packet_start_timestamp_ms & calc_DL_packets_RX1_scheduled_start_ms <= temp_UL_packet_end_time_ms);
                                    if(isempty(temp_packets_DL_overlapping_partially1_indexes)==false)
                                       temp_packets_DL_overlapping_all_indexes=unique(vertcat(temp_packets_DL_overlapping_all_indexes,temp_packets_DL_overlapping_partially1_indexes),'rows'); 
                                    end
                                end
                                try
                                    [temp_packets_DL_overlapping_partially2_indexes(:,1), temp_packets_DL_overlapping_partially2_indexes(:,2)] = find(calc_DL_packets_RX1_scheduled_end_ms >= temp_packet_start_timestamp_ms & calc_DL_packets_RX1_scheduled_end_ms <= temp_UL_packet_end_time_ms);
                                    if(isempty(temp_packets_DL_overlapping_partially2_indexes)==false)
                                       temp_packets_DL_overlapping_all_indexes=unique(vertcat(temp_packets_DL_overlapping_all_indexes,temp_packets_DL_overlapping_partially2_indexes),'rows'); 
                                    end
                                %find all DL packets starting before temp_packet_start_timestamp_ms 
                                %and ending after temp_UL_packet_end_time_ms (NOTE: includes also self)     
                                end
                                try
                                    [temp_packets_DL_overlapping_fully_indexes(:,1), temp_packets_DL_overlapping_fully_indexes(:,2)] = find(calc_DL_packets_RX1_scheduled_start_ms <= temp_packet_start_timestamp_ms & calc_DL_packets_RX1_scheduled_end_ms >= temp_UL_packet_end_time_ms);
                                    if(isempty(temp_packets_DL_overlapping_fully_indexes)==false)
                                       temp_packets_DL_overlapping_all_indexes=unique(vertcat(temp_packets_DL_overlapping_all_indexes,temp_packets_DL_overlapping_fully_indexes),'rows'); 
                                    end
                                end 
                                if exist('temp_packets_DL_overlapping_all_indexes','var') == 1
                                    if(isempty(temp_packets_DL_overlapping_all_indexes)==false)
                                        for temp_i_interf_packet=1:1:size(temp_packets_DL_overlapping_all_indexes,1)
                                            temp_GW_interferer=temp_packets_DL_overlapping_all_indexes(temp_i_interf_packet,1);
                                            temp_GW_DL_packet=temp_packets_DL_overlapping_all_indexes(temp_i_interf_packet,2);
                                            temp_GW_DR=calc_DL_packets_RX1_scheduled_DR(temp_GW_interferer,temp_GW_DL_packet);
                                            %go through the DL packets and see if any blocks the UL
                                            for temp_i_num_GW=1:input_num_GWs%for each GW 
                                                %see if blocking
                                                if calc_UL_GW_ED_GWint_SF_blockers(temp_i_num_GW,temp_ED,temp_GW_interferer,temp_GW_DR+1)==1
                                                    if ((input_full_duplex_GW==true) && (temp_i_num_GW==temp_GW_interferer))
                                                       %ignore case of self-interference
                                                    else
                                                        result_UL_packet_delivery_matrix(temp_ED,temp_ED_packet_index,temp_i_num_GW)=false;
                                                        result_extra_UL_packet_undelivery_matrix_cause_DL_interference(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                                    end
                                                end
                                            end 
                                        end
                                    end
                                end
                             end
                             
                            %generate class c interference with probability
                            %P=input_class_c_interference_probability
                            for temp_i_num_GW=1:input_num_GWs%for each GW
                                for i_GW_interferer=1:input_num_GWs
                                    for i_DR_interferer=0:calc_num_DRs-1
                            %see if blocking
                            	if calc_UL_GW_ED_GWint_SF_blockers(temp_i_num_GW,temp_ED,i_GW_interferer,i_DR_interferer+1)==1
                                	if ((input_full_duplex_GW==true) && (temp_i_num_GW==i_GW_interferer))
                                    	%ignore case of self-interference
                                    else
                                    	if(input_class_c_interference_probability>=rand())
                                            result_UL_packet_delivery_matrix(temp_ED,temp_ED_packet_index,temp_i_num_GW)=false;
                                            result_extra_UL_packet_undelivery_matrix_cause_DL_interference(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                        end
                                    end
                                end
                                    end
                                end
                            end
                            
                        end
                        %some defines
                        temp_num_GWs_received=0;
                        temp_UL_RSSI_GWs_dBm=ones(input_num_GWs,1)*(-1000);     
                        %calculate how many GWs have received each packet and 
                        for temp_i_num_GW=1:input_num_GWs%for each GW 
                            if(result_UL_packet_delivery_matrix(temp_ED,temp_ED_packet_index,temp_i_num_GW)==1)
                                temp_num_GWs_received=temp_num_GWs_received+1;
                                temp_UL_RSSI_GWs_dBm(temp_i_num_GW,1)=calc_RSSI_UL_ED_GW_dBm(temp_ED,temp_i_num_GW);
                            end
                        end       
                        result_UL_packet_delivery_matrix_num_GWs(temp_ED,temp_ED_packet_index)=temp_num_GWs_received;
                        %schedule the DL for this UL (if packet has been received)
                        if temp_num_GWs_received>0
                            temp_rand=rand();
                            temp_flag_downlink_schedulled=false;
                            if(temp_rand<=input_DL_scheduling_probability)%try to schedule this packet

                                temp_timestamp_DL_RX_slot_start_ms=temp_UL_packet_end_time_ms;
                                temp_RX_packet_duration_ms=input_DL_packets_duration_ms(calc_DRs_RX1(temp_ED,1)+1);
                                temp_timestamp_DL_RX_end_ms=temp_timestamp_DL_RX_slot_start_ms+temp_RX_packet_duration_ms;
                                 
                                temp_timestamp_DL_RX1_slot_start_ms=temp_UL_packet_end_time_ms+input_RX1_delay_s*1000;
                                temp_RX1_packet_duration_ms=input_DL_packets_duration_ms(calc_DRs_RX1(temp_ED,1)+1);
                                temp_timestamp_DL_RX1_end_ms=temp_timestamp_DL_RX1_slot_start_ms+temp_RX1_packet_duration_ms;
                                temp_timestamp_DL_RX2_slot_start_ms=temp_timestamp_DL_RX1_slot_start_ms+input_RX2_delay_s*1000;
                                temp_RX2_packet_duration_ms=input_DL_packets_duration_ms(calc_DRs_RX2(temp_ED,1)+1);
                                temp_timestamp_DL_RX2_end_ms=temp_timestamp_DL_RX2_slot_start_ms+temp_RX2_packet_duration_ms;
                                %sort all GWs based on their UL RSSI level
                                [temp_GWs_received_RSSI_sorted_dBm,temp_GWs_received_index_sorted] = sort(temp_UL_RSSI_GWs_dBm,'descend');
                                %check which of the RWs has priorities
                                for temp_GW_receiver_ID=1:1:temp_num_GWs_received
                                    temp_i_num_GW=temp_GWs_received_index_sorted(temp_GW_receiver_ID,1);
                                    
%                                % if device works in class a - 0
                                 if(input_ED_class(temp_ED)==0)                                             
                                    %RX1 has the priority over RX2
                                    if(input_DL_RW_priority==1)
                                        %check if DC limitations allow schedulling to RX1
                                        if calc_simulation_end_ms>=temp_timestamp_DL_RX1_end_ms
                                            if(calc_DL_channel_release_timestamp_ms(1,temp_i_num_GW)<=temp_timestamp_DL_RX1_slot_start_ms)
                                               %allocate the packet to RX1 to this GW
                                               calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1)=calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1)+1;
                                               calc_DL_packets_RX1_scheduled_start_ms(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=temp_timestamp_DL_RX1_slot_start_ms;
                                               calc_DL_packets_RX1_scheduled_end_ms(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=temp_timestamp_DL_RX1_end_ms;
                                               calc_DL_packets_RX1_scheduled_ED_target(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=temp_ED;
                                               calc_DL_packets_RX1_scheduled_DR(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=calc_DRs_RX1(temp_ED,1);
                                               %update the duty cycle release timestamp
                                               calc_DL_channel_release_timestamp_ms(1,temp_i_num_GW)=temp_timestamp_DL_RX1_slot_start_ms+temp_RX1_packet_duration_ms*(100/input_RX1_dc_limit_percent);
                                               temp_flag_downlink_schedulled=true;
                                               if input_DL_tranmission_in_one_RW_only==true
                                                   break; %break the for cycle over GWs
                                               end
                                            else
                                                result_extra_DL_packet_unscheduled_matrix_cause_RX1_backoff(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                            end
                                        else

                                        end
                                        %check if RX2 config allows delivery of packet in
                                        %downlink
                                        if calc_simulation_end_ms>=temp_timestamp_DL_RX2_end_ms
                                            if((input_RX2_DR_notless_UL_DR==false)||(calc_DRs_RX2(temp_ED,1)<=calc_DRs_UL(temp_ED,1)))
                                                %check if DC limitations allow schedulling to RX2
                                                if(calc_DL_channel_release_timestamp_ms(2,temp_i_num_GW)<=temp_timestamp_DL_RX2_slot_start_ms)
                                                   %allocate the packet to RX1 to this GW
                                                   calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1)=calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1)+1;
                                                   calc_DL_packets_RX2_scheduled_start_ms(temp_i_num_GW,calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1))=temp_timestamp_DL_RX2_slot_start_ms;
                                                   calc_DL_packets_RX2_scheduled_end_ms(temp_i_num_GW,calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1))=temp_timestamp_DL_RX2_end_ms;
                                                   calc_DL_packets_RX2_scheduled_ED_target(temp_i_num_GW,calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1))=temp_ED;
                                                   calc_DL_packets_RX2_scheduled_DR(temp_i_num_GW,calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1))=calc_DRs_RX2(temp_ED,1);
                                                   %update the duty cycle release timestamp
                                                   calc_DL_channel_release_timestamp_ms(2,temp_i_num_GW)=temp_timestamp_DL_RX2_slot_start_ms+temp_RX2_packet_duration_ms*(100/input_RX2_dc_limit_percent);
                                                   temp_flag_downlink_schedulled=true;
                                                   if input_DL_tranmission_in_one_RW_only==true
                                                        break; %break the for cycle over GWs
                                                   end
                                                else
                                                    result_extra_DL_packet_unscheduled_matrix_cause_RX2_backoff(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                                end   
                                            else
                                                result_extra_DL_packet_unscheduled_matrix_cause_RX2_low_DR(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                            end
                                        else

                                        end
                                        if temp_flag_downlink_schedulled==true
                                            for temp_GW_number=temp_GW_receiver_ID:1:temp_num_GWs_received
                                                result_extra_DL_packet_unscheduled_matrix_cause_planned_already(temp_ED,temp_ED_packet_index,temp_GW_number)=1;
                                            end
                                            break; %break the for cycle over GWs
                                        end    
                                    %RX2 has the priority over RX1
                                    elseif(input_DL_RW_priority==2)
                                        %check if RX2 config allows delivery of packet in
                                        %downlink
                                        if calc_simulation_end_ms>=temp_timestamp_DL_RX2_end_ms
                                            if((input_RX2_DR_notless_UL_DR==false)||(calc_DRs_RX2(temp_ED,1)<=calc_DRs_UL(temp_ED,1)))
                                                %check if DC limitations allow schedulling to RX2
                                                if(calc_DL_channel_release_timestamp_ms(2,temp_i_num_GW)<=temp_timestamp_DL_RX2_slot_start_ms)
                                                   %allocate the packet to RX1 to this GW
                                                   calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1)=calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1)+1;
                                                   calc_DL_packets_RX2_scheduled_start_ms(temp_i_num_GW,calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1))=temp_timestamp_DL_RX2_slot_start_ms;
                                                   calc_DL_packets_RX2_scheduled_end_ms(temp_i_num_GW,calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1))=temp_timestamp_DL_RX2_end_ms;
                                                   calc_DL_packets_RX2_scheduled_ED_target(temp_i_num_GW,calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1))=temp_ED;
                                                   calc_DL_packets_RX2_scheduled_DR(temp_i_num_GW,calc_DL_GW_packet_RX2_scheduled(temp_i_num_GW,1))=calc_DRs_RX2(temp_ED,1);
                                                   %update the duty cycle release timestamp
                                                   calc_DL_channel_release_timestamp_ms(2,temp_i_num_GW)=temp_timestamp_DL_RX2_slot_start_ms+temp_RX2_packet_duration_ms*(100/input_RX2_dc_limit_percent);
                                                   temp_flag_downlink_schedulled=true;
                                                   if input_DL_tranmission_in_one_RW_only==true
                                                        break; %break the for cycle over GWs
                                                   end
                                                else
                                                    result_extra_DL_packet_unscheduled_matrix_cause_RX2_backoff(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                                end     
                                            else
                                                result_extra_DL_packet_unscheduled_matrix_cause_RX2_low_DR(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                            end
                                        else

                                        end
                                        %check if DC limitations allow schedulling to RX1
                                        if calc_simulation_end_ms>=temp_timestamp_DL_RX1_end_ms
                                            if(calc_DL_channel_release_timestamp_ms(1,temp_i_num_GW)<=temp_timestamp_DL_RX1_slot_start_ms)
                                               %allocate the packet to RX1 to this GW
                                               calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1)=calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1)+1;
                                               calc_DL_packets_RX1_scheduled_start_ms(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=temp_timestamp_DL_RX1_slot_start_ms;
                                               calc_DL_packets_RX1_scheduled_end_ms(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=temp_timestamp_DL_RX1_end_ms;
                                               calc_DL_packets_RX1_scheduled_ED_target(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=temp_ED;
                                               calc_DL_packets_RX1_scheduled_DR(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=calc_DRs_RX1(temp_ED,1);
                                               %update the duty cycle release timestamp
                                               calc_DL_channel_release_timestamp_ms(1,temp_i_num_GW)=temp_timestamp_DL_RX1_slot_start_ms+temp_RX1_packet_duration_ms*(100/input_RX1_dc_limit_percent);
                                               temp_flag_downlink_schedulled=true;
                                               if input_DL_tranmission_in_one_RW_only==true
                                                    break; %break the for cycle over GWs
                                               end
                                            else
                                                result_extra_DL_packet_unscheduled_matrix_cause_RX1_backoff(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                            end 
                                        end
                                        if temp_flag_downlink_schedulled==true
                                           for temp_GW_number=temp_GW_receiver_ID:1:temp_num_GWs_received
                                               result_extra_DL_packet_unscheduled_matrix_cause_planned_already(temp_ED,temp_ED_packet_index,temp_GW_number)=1;
                                           end
                                           break; %break the for cycle over GWs
                                        end                         
                                    else
                                    %ERROR 
                                    end
                                end
%                                
%                                     %if class c
                                    if (input_ED_class(temp_ED)==1)
                                    %if calc_simulation_end_ms>=temp_timestamp_DL_RX_end_ms
                                        if(calc_DL_channel_release_timestamp_ms(1,temp_i_num_GW)<=temp_timestamp_DL_RX_slot_start_ms)
                                             %allocate the packet to RX1 to this GW
                                               calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1)=calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1)+1;
                                               calc_DL_packets_RX1_scheduled_start_ms(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=temp_timestamp_DL_RX_slot_start_ms;
                                               calc_DL_packets_RX1_scheduled_end_ms(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=temp_timestamp_DL_RX_end_ms;
                                               calc_DL_packets_RX1_scheduled_ED_target(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=temp_ED;
                                               calc_DL_packets_RX1_scheduled_DR(temp_i_num_GW,calc_DL_GW_packet_RX1_scheduled(temp_i_num_GW,1))=calc_DRs_RX1(temp_ED,1);
                                               %update the duty cycle release timestamp
                                               calc_DL_channel_release_timestamp_ms(1,temp_i_num_GW)=temp_timestamp_DL_RX_slot_start_ms+temp_RX_packet_duration_ms*(100/input_RX1_dc_limit_percent);
                                               temp_flag_downlink_schedulled=true;
                                               if input_DL_tranmission_in_one_RW_only==true
                                                   break; %break the for cycle over GWs
                                               end
                                        end
                                    %end
                                    end
                                    
                                end 
                            else
                                for temp_i_num_GW=1:input_num_GWs%for each GW 
                                    result_extra_DL_packet_unscheduled_matrix_cause_no_traffic(temp_ED,temp_ED_packet_index,temp_i_num_GW)=1;
                                end 
                            end           
                        end
                        temp_previous_packet_start_timestamp_ms=temp_packet_start_timestamp_ms;
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation


                    %for temp_i_num_GW=1:input_num_GWs%for each GW 
                    %    result_UL_packet_delivery_matrix_num_GWs=result_UL_packet_delivery_matrix_num_GWs+result_UL_packet_delivery_matrix(:,:,temp_i_num_GW);
                    %end     
                    %clearvars -except calc_* input_* result_* DEF_* config_* 

                    %crawl through RX1 DL and see if packets have been received
                    result_DL_RX1_packet_delivery_matrix=zeros(size(calc_DL_packets_RX1_scheduled_start_ms,1),size(calc_DL_packets_RX1_scheduled_start_ms,2));
                    result_extra_DL_RX1_packet_undelivery_matrix_cause_SNR=zeros(size(calc_DL_packets_RX1_scheduled_start_ms,1),size(calc_DL_packets_RX1_scheduled_start_ms,2));
                    result_extra_DL_RX1_packet_undelivery_matrix_cause_UL_interference=zeros(size(calc_DL_packets_RX1_scheduled_start_ms,1),size(calc_DL_packets_RX1_scheduled_start_ms,2));
                    result_extra_DL_RX1_packet_undelivery_matrix_cause_DL_interference=zeros(size(calc_DL_packets_RX1_scheduled_start_ms,1),size(calc_DL_packets_RX1_scheduled_start_ms,2));
                    temp_DL_RX1_sorted_packet_starts=sort(nonzeros(calc_DL_packets_RX1_scheduled_start_ms(:)));
                    calc_DL_RX1_num_packets=size(temp_DL_RX1_sorted_packet_starts,1);
                    temp_duplicate=0;%handles the case when multiple packets start exactly at the same time
                    temp_DL_RX1_previous_packet_start_timestamp_ms=-1;
                    for temp_i_DL_RX1_packet=1:1:calc_DL_RX1_num_packets%go through all DL RX1 packets
                        temp_DL_RX1_packet_start_timestamp_ms=temp_DL_RX1_sorted_packet_starts(temp_i_DL_RX1_packet,1);
                        if(temp_DL_RX1_packet_start_timestamp_ms==temp_DL_RX1_previous_packet_start_timestamp_ms)
                            temp_duplicate=temp_duplicate+1;
                        else
                            temp_duplicate=0;
                        end    
                        clearvars temp_GW_originating temp_DL_RX1_GW_packet_index 
                        [temp_found_GW_originating, temp_found_DL_RX1_GW_packet_index] = find(calc_DL_packets_RX1_scheduled_start_ms==temp_DL_RX1_packet_start_timestamp_ms);
                        temp_DL_RX1_originating_GW=temp_found_GW_originating(temp_duplicate+1,1);
                        temp_DL_RX1_originating_GW_packet_index=temp_found_DL_RX1_GW_packet_index(temp_duplicate+1,1);
                        temp_DL_RX1_packet_DR=calc_DL_packets_RX1_scheduled_DR(temp_DL_RX1_originating_GW,temp_DL_RX1_originating_GW_packet_index);
                        temp_DL_RX1_packet_end_time_ms=calc_DL_packets_RX1_scheduled_end_ms(temp_DL_RX1_originating_GW,temp_DL_RX1_originating_GW_packet_index);
                        temp_DL_RX1_packet_ED_target=calc_DL_packets_RX1_scheduled_ED_target(temp_DL_RX1_originating_GW,temp_DL_RX1_originating_GW_packet_index);
                        %check if SNR is sufficient
                        if calc_sensitivity_gap_DL_ED_GW_RX1_dB(temp_DL_RX1_packet_ED_target,temp_DL_RX1_originating_GW)>0
                            temp_DL_RX1_packet_flag_received=true;
                        else
                            result_extra_DL_RX1_packet_undelivery_matrix_cause_SNR(temp_DL_RX1_originating_GW,temp_DL_RX1_originating_GW_packet_index)=1;
                            temp_DL_RX1_packet_flag_received=false;
%remove DL interference     temp_DL_RX1_packet_flag_received=true;
                        end
                        %check for infererences from UL
                        if 1
                            clearvars temp_DL_RX1_EDs_interfering;
                            clearvars temp_packets_DL_RX1_UL_overlapping_partially1_indexes temp_packets_DL_RX1_UL_overlapping_partially2_indexes;
                            clearvars temp_packets_DL_RX1_UL_overlapping_fully_indexes temp_packets_DL_RX1_UL_overlapping_all_indexes;
                            temp_DL_RX1_EDs_interfering=calc_DL_GW_ED_EDint_SF_blockers(temp_DL_RX1_originating_GW,temp_DL_RX1_packet_ED_target,:,temp_DL_RX1_packet_DR+1);
                            temp_DL_RX1_num_interfering_EDs=sum(temp_DL_RX1_EDs_interfering,3);
                            if temp_DL_RX1_num_interfering_EDs~=0 %there are potential interferers
                                %see if any UL packets overlap in time
                                temp_packets_DL_RX1_UL_overlapping_all_indexes=[];
                                try
                                    [temp_packets_DL_RX1_UL_overlapping_partially1_indexes(:,1), temp_packets_DL_RX1_UL_overlapping_partially1_indexes(:,2)] = find(calc_UL_packets_scheduled_start_ms >= temp_DL_RX1_packet_start_timestamp_ms & calc_UL_packets_scheduled_start_ms <= temp_DL_RX1_packet_end_time_ms);
                                    if(isempty(temp_packets_DL_RX1_UL_overlapping_partially1_indexes)==false)
                                       temp_packets_DL_RX1_UL_overlapping_all_indexes=unique(vertcat(temp_packets_DL_RX1_UL_overlapping_all_indexes,temp_packets_DL_RX1_UL_overlapping_partially1_indexes),'rows'); 
                                    end                  
                                end
                                try
                                    [temp_packets_DL_RX1_UL_overlapping_partially2_indexes(:,1), temp_packets_DL_RX1_UL_overlapping_partially2_indexes(:,2)] = find(calc_UL_packets_scheduled_end_ms >= temp_DL_RX1_packet_start_timestamp_ms & calc_UL_packets_scheduled_end_ms <= temp_DL_RX1_packet_end_time_ms);
                                    if(isempty(temp_packets_DL_RX1_UL_overlapping_partially2_indexes)==false)
                                       temp_packets_DL_RX1_UL_overlapping_all_indexes=unique(vertcat(temp_packets_DL_RX1_UL_overlapping_all_indexes,temp_packets_DL_RX1_UL_overlapping_partially2_indexes),'rows'); 
                                    end    
                                end            
                                try
                                    [temp_packets_DL_RX1_UL_overlapping_fully_indexes(:,1), temp_packets_DL_RX1_UL_overlapping_fully_indexes(:,2)] = find(calc_UL_packets_scheduled_start_ms <= temp_DL_RX1_packet_start_timestamp_ms & calc_UL_packets_scheduled_end_ms >= temp_DL_RX1_packet_end_time_ms);
                                    if(isempty(temp_packets_DL_RX1_UL_overlapping_fully_indexes)==false)
                                       temp_packets_DL_RX1_UL_overlapping_all_indexes=unique(vertcat(temp_packets_DL_RX1_UL_overlapping_all_indexes,temp_packets_DL_RX1_UL_overlapping_fully_indexes),'rows'); 
                                    end    
                                end   
                                %go through each overlapping in time packet and see if any of
                                %them blocks reception
                                if exist('temp_packets_DL_overlapping_all_indexes','var') == 1
                                    if(isempty(temp_packets_DL_RX1_UL_overlapping_all_indexes)==false)
                                        for temp_i_interf_packet=1:1:size(temp_packets_DL_RX1_UL_overlapping_all_indexes,1)
                                            temp_ED_interferer=temp_packets_DL_RX1_UL_overlapping_all_indexes(temp_i_interf_packet,1)
                                            if(calc_DL_GW_ED_EDint_SF_blockers(temp_DL_RX1_originating_GW,temp_DL_RX1_packet_ED_target,temp_ED_interferer,temp_DL_RX1_packet_DR+1)==1)
                                                  temp_DL_RX1_packet_flag_received=false;
                                                  result_extra_DL_RX1_packet_undelivery_matrix_cause_UL_interference(temp_DL_RX1_originating_GW,temp_DL_RX1_originating_GW_packet_index)=1;
                                                   break;%break the for cycle
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        %check for infererences from DL
                        if 1
                            clearvars temp_DL_RX1_GWs_interfering;
                            clearvars temp_packets_DL_RX1_DL_overlapping_partially1_indexes temp_packets_DL_RX1_DL_overlapping_partially2_indexes;
                            clearvars temp_packets_DL_RX1_DL_overlapping_fully_indexes temp_packets_DL_RX1_DL_overlapping_all_indexes;
                            %find all packets starting or ending between
                            %temp_DL_RX1_packet_start_timestamp_ms and temp_DL_RX1_packet_end_time_ms  (NOTE: includes also self)
                            [temp_packets_DL_RX1_DL_overlapping_partially1_indexes(:,1), temp_packets_DL_RX1_DL_overlapping_partially1_indexes(:,2)] = find(calc_DL_packets_RX1_scheduled_start_ms >= temp_DL_RX1_packet_start_timestamp_ms & calc_DL_packets_RX1_scheduled_start_ms <= temp_DL_RX1_packet_end_time_ms);
                            [temp_packets_DL_RX1_DL_overlapping_partially2_indexes(:,1), temp_packets_DL_RX1_DL_overlapping_partially2_indexes(:,2)] = find(calc_DL_packets_RX1_scheduled_end_ms >= temp_DL_RX1_packet_start_timestamp_ms & calc_DL_packets_RX1_scheduled_end_ms <= temp_DL_RX1_packet_end_time_ms);
                            %find all packets starting before temp_DL_RX1_packet_start_timestamp_ms 
                            %and ending after temp_DL_RX1_packet_end_time_ms (NOTE: includes also self)         
                            [temp_packets_DL_RX1_DL_overlapping_fully_indexes(:,1), temp_packets_DL_RX1_DL_overlapping_fully_indexes(:,2)] = find(calc_DL_packets_RX1_scheduled_start_ms <= temp_DL_RX1_packet_start_timestamp_ms & calc_DL_packets_RX1_scheduled_end_ms >= temp_DL_RX1_packet_end_time_ms);
                            %combine all these and remove duplicates
                            temp_packets_DL_RX1_DL_overlapping_all_indexes=setdiff(unique(vertcat(temp_packets_DL_RX1_DL_overlapping_partially1_indexes,temp_packets_DL_RX1_DL_overlapping_partially2_indexes,temp_packets_DL_RX1_DL_overlapping_fully_indexes),'rows'),[temp_DL_RX1_originating_GW temp_DL_RX1_originating_GW_packet_index],'rows');
                            if(isempty(temp_packets_DL_RX1_DL_overlapping_all_indexes)==false)%some overlapping transmissions exist
                                for temp_i_interf_packet=1:1:size(temp_packets_DL_RX1_DL_overlapping_all_indexes,1)
                                    temp_GW_interferer=temp_packets_DL_RX1_DL_overlapping_all_indexes(temp_i_interf_packet,1);
                                    temp_GW_interferer_packet_index=temp_packets_DL_RX1_DL_overlapping_all_indexes(temp_i_interf_packet,2);
                                    temp_GW_interferer_DR=calc_DL_packets_RX1_scheduled_DR(temp_GW_interferer,temp_GW_interferer_packet_index);
                                    if calc_DL_RX1_GWtarg_ED_GWint_SF_blockers(temp_DL_RX1_originating_GW,temp_DL_RX1_packet_ED_target,temp_GW_interferer,temp_GW_interferer_DR+1)==1
                                         temp_DL_RX1_packet_flag_received=false;
                                         result_extra_DL_RX1_packet_undelivery_matrix_cause_DL_interference(temp_DL_RX1_originating_GW,temp_DL_RX1_originating_GW_packet_index)=1;
                                         break;%break the for cycle
                                    end
                                end
                            end
                        end
                        %make the final decision is the packet has been received
                        if 1
                            if temp_DL_RX1_packet_flag_received==true
                                result_DL_RX1_packet_delivery_matrix(temp_DL_RX1_originating_GW,temp_DL_RX1_originating_GW_packet_index)=1;
                            end
                        end
                        temp_DL_RX1_previous_packet_start_timestamp_ms=temp_DL_RX1_packet_start_timestamp_ms;
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_*  simulation

                    %crawl through RX2 DL and see if packets have been received
                    result_DL_RX2_packet_delivery_matrix=zeros(size(calc_DL_packets_RX2_scheduled_start_ms,1),size(calc_DL_packets_RX2_scheduled_start_ms,2));
                    result_extra_DL_RX2_packet_undelivery_matrix_cause_SNR=zeros(size(calc_DL_packets_RX2_scheduled_start_ms,1),size(calc_DL_packets_RX2_scheduled_start_ms,2));
                    result_extra_DL_RX2_packet_undelivery_matrix_cause_DL_interference=zeros(size(calc_DL_packets_RX2_scheduled_start_ms,1),size(calc_DL_packets_RX2_scheduled_start_ms,2));
                    temp_DL_RX2_sorted_packet_starts=sort(nonzeros(calc_DL_packets_RX2_scheduled_start_ms(:)));
                    calc_DL_RX2_num_packets=size(temp_DL_RX2_sorted_packet_starts,1);
                    temp_duplicate=0;%handles the case when multiple packets start exactly at the same time
                    temp_DL_RX2_previous_packet_start_timestamp_ms=-1;
                    for temp_i_DL_RX2_packet=1:1:calc_DL_RX2_num_packets%go through all DL RX2 packets    
                        temp_DL_RX2_packet_start_timestamp_ms=temp_DL_RX2_sorted_packet_starts(temp_i_DL_RX2_packet,1);
                        if(temp_DL_RX2_packet_start_timestamp_ms==temp_DL_RX2_previous_packet_start_timestamp_ms)
                            temp_duplicate=temp_duplicate+1;
                        else
                            temp_duplicate=0;
                        end    
                        clearvars temp_GW_originating temp_DL_RX2_GW_packet_index 
                        [temp_found_GW_originating, temp_found_DL_RX2_GW_packet_index] = find(calc_DL_packets_RX2_scheduled_start_ms==temp_DL_RX2_packet_start_timestamp_ms);
                        temp_DL_RX2_originating_GW=temp_found_GW_originating(temp_duplicate+1,1);
                        temp_DL_RX2_originating_GW_packet_index=temp_found_DL_RX2_GW_packet_index(temp_duplicate+1,1);
                        temp_DL_RX2_packet_DR=calc_DL_packets_RX2_scheduled_DR(temp_DL_RX2_originating_GW,temp_DL_RX2_originating_GW_packet_index);
                        temp_DL_RX2_packet_end_time_ms=calc_DL_packets_RX2_scheduled_end_ms(temp_DL_RX2_originating_GW,temp_DL_RX2_originating_GW_packet_index);
                        temp_DL_RX2_packet_ED_target=calc_DL_packets_RX2_scheduled_ED_target(temp_DL_RX2_originating_GW,temp_DL_RX2_originating_GW_packet_index);
                        %check if SNR is sufficient
                        if calc_sensitivity_gap_DL_ED_GW_RX2_dB(temp_DL_RX2_packet_ED_target,temp_DL_RX2_originating_GW)>0
                           temp_DL_RX2_packet_flag_received=true;
                        else
                            result_extra_DL_RX2_packet_undelivery_matrix_cause_SNR(temp_DL_RX2_originating_GW,temp_DL_RX2_originating_GW_packet_index)=1;
                            temp_DL_RX2_packet_flag_received=false;
                        end

                        %check for infererences from DL
                        if 1
                            clearvars temp_DL_RX2_GWs_interfering;
                            clearvars temp_packets_DL_RX2_DL_overlapping_partially1_indexes temp_packets_DL_RX2_DL_overlapping_partially2_indexes;
                            clearvars temp_packets_DL_RX2_DL_overlapping_fully_indexes temp_packets_DL_RX2_DL_overlapping_all_indexes;
                            %find all packets starting or ending between
                            %temp_DL_RX2_packet_start_timestamp_ms and temp_DL_RX2_packet_end_time_ms  (NOTE: includes also self)
                            [temp_packets_DL_RX2_DL_overlapping_partially1_indexes(:,1) temp_packets_DL_RX2_DL_overlapping_partially1_indexes(:,2)] = find(calc_DL_packets_RX2_scheduled_start_ms >= temp_DL_RX2_packet_start_timestamp_ms & calc_DL_packets_RX2_scheduled_start_ms <= temp_DL_RX2_packet_end_time_ms);
                            [temp_packets_DL_RX2_DL_overlapping_partially2_indexes(:,1) temp_packets_DL_RX2_DL_overlapping_partially2_indexes(:,2)] = find(calc_DL_packets_RX2_scheduled_end_ms >= temp_DL_RX2_packet_start_timestamp_ms & calc_DL_packets_RX2_scheduled_end_ms <= temp_DL_RX2_packet_end_time_ms);
                            %find all packets starting before temp_DL_RX2_packet_start_timestamp_ms 
                            %and ending after temp_DL_RX2_packet_end_time_ms (NOTE: includes also self)         
                            [temp_packets_DL_RX2_DL_overlapping_fully_indexes(:,1) temp_packets_DL_RX2_DL_overlapping_fully_indexes(:,2)] = find(calc_DL_packets_RX2_scheduled_start_ms <= temp_DL_RX2_packet_start_timestamp_ms & calc_DL_packets_RX2_scheduled_end_ms >= temp_DL_RX2_packet_end_time_ms);
                            %combine all these and remove duplicates
                            temp_packets_DL_RX2_DL_overlapping_all_indexes=setdiff(unique(vertcat(temp_packets_DL_RX2_DL_overlapping_partially1_indexes,temp_packets_DL_RX2_DL_overlapping_partially2_indexes,temp_packets_DL_RX2_DL_overlapping_fully_indexes),'rows'),[temp_DL_RX2_originating_GW temp_DL_RX2_originating_GW_packet_index],'rows');
                            if(isempty(temp_packets_DL_RX2_DL_overlapping_all_indexes)==false)%some overlapping transmissions exist
                                for temp_i_interf_packet=1:1:size(temp_packets_DL_RX2_DL_overlapping_all_indexes,1)
                                    temp_GW_interferer=temp_packets_DL_RX2_DL_overlapping_all_indexes(temp_i_interf_packet,1);
                                    if calc_DL_RX2_GWtarg_ED_GWint_blockers(temp_DL_RX2_originating_GW,temp_DL_RX2_packet_ED_target,temp_GW_interferer)==1
                                         result_extra_DL_RX2_packet_undelivery_matrix_cause_DL_interference(temp_DL_RX2_originating_GW,temp_DL_RX2_originating_GW_packet_index)=1;
                                         temp_DL_RX2_packet_flag_received=false;
                                         break;%break the for cycle
                                    end
                                end
                            end
                        end       
                        %make the final decision if the packet has been received
                        if 1
                            if temp_DL_RX2_packet_flag_received==true
                                result_DL_RX2_packet_delivery_matrix(temp_DL_RX2_originating_GW,temp_DL_RX2_originating_GW_packet_index)=1;
                            end
                        end       
                        temp_DL_RX2_previous_packet_start_timestamp_ms=temp_DL_RX2_packet_start_timestamp_ms;
                    end
                    clearvars -except calc_* input_* result_* DEF_* config_* simulation
                    
                    %update the loggs of results
                    %"primary" results
                    %UL packets sent by each ED
                    result_iterations_log_UL_tranmissions_by_ED(:,config_iteration)=sum(calc_UL_packets_scheduled_start_ms~=0,2);
                    %UL packets received by GWs and NS
                    %temp_UL_packet_delivery_matrix_cumulative=zeros(size(calc_UL_packets_scheduled_start_ms,1),size(calc_UL_packets_scheduled_start_ms,2));
                    %for temp_GW_number=1:1:input_num_GWs
                    %    result_iterations_log_UL_deliveries_by_GW_and_ED(:,temp_GW_number,config_iteration)=sum(result_UL_packet_delivery_matrix_num_GWs~=0,2);
                    %    temp_UL_packet_delivery_matrix_cumulative=temp_UL_packet_delivery_matrix_cumulative+result_iterations_log_UL_deliveries_by_GW_and_ED(:,temp_GW_number,config_iteration);
                    %end 
                    %for temp_GW_number=1:1:input_num_GWs
                    %    result_iterations_log_UL_deliveries_by_GW_and_ED(:,temp_GW_number,config_iteration)=sum(result_UL_packet_delivery_matrix_num_GWs~=0,2);
                    %    temp_UL_packet_delivery_matrix_cumulative=temp_UL_packet_delivery_matrix_cumulative+result_iterations_log_UL_deliveries_by_GW_and_ED(:,temp_GW_number,config_iteration);
                    %end                     
                    
                    result_iterations_log_UL_deliveries_by_ED(:,config_iteration)=sum(result_UL_packet_delivery_matrix_num_GWs~=0,2);

                    %mean and standard deviation number of GWs receiving each packets of each ED
                    result_iterations_log_UL_mean_GW_EDs(:,config_iteration)=mean(result_UL_packet_delivery_matrix_num_GWs,2);

                    %RX1: convert to 1D array and find non-zeros
                    clearvars temp_nonzero_idx_DL_scheduled_index 
                    if(isempty(calc_DL_packets_RX1_scheduled_ED_target)==false)
                        [temp_nonzero_idx_DL_scheduled_index]=find(calc_DL_packets_RX1_scheduled_ED_target(:)~=0);

                        %DL packets scheduled in RX1 by DR
                        if(isempty(calc_DL_packets_RX1_scheduled_DR)==true)
                            for temp_DR_number=0:1:calc_num_DRs-1
                                result_iterations_log_DL_RX1_schedulled_by_DR(temp_DR_number+1,config_iteration)=0;
                            end
                        else
                            clearvars temp_vector_calc_DL_packets_RX1_scheduled_DR 
                            temp_vector_calc_DL_packets_RX1_scheduled_DR=calc_DL_packets_RX1_scheduled_DR(:);
                            for temp_DR_number=0:1:calc_num_DRs-1
                                result_iterations_log_DL_RX1_schedulled_by_DR(temp_DR_number+1,config_iteration)=sum(temp_vector_calc_DL_packets_RX1_scheduled_DR(temp_nonzero_idx_DL_scheduled_index)==temp_DR_number,'all');
                            end          
                        end

                        %DL packets scheduled in RX1 by ED
                        if(isempty(calc_DL_packets_RX1_scheduled_ED_target)==true)
                            for temp_ED_number=1:1:input_num_EDs
                                result_iterations_log_DL_RX1_schedulled_by_ED(temp_ED_number,config_iteration)=0;
                            end
                        else
                            clearvars temp_vector_calc_DL_packets_RX1_scheduled_ED 
                            temp_vector_calc_DL_packets_RX1_scheduled_ED=calc_DL_packets_RX1_scheduled_ED_target(:);
                            for temp_ED_number=1:1:input_num_EDs
                                result_iterations_log_DL_RX1_schedulled_by_ED(temp_ED_number,config_iteration)=sum(temp_vector_calc_DL_packets_RX1_scheduled_ED(temp_nonzero_idx_DL_scheduled_index)==temp_ED_number,'all');
                            end          
                        end

                        %DL packets scheduled in RX1 by GW
                        if(isempty(calc_DL_packets_RX1_scheduled_ED_target)==true)
                            for temp_GW_number=1:1:input_num_GWs
                                result_iterations_log_DL_RX1_schedulled_by_GW(temp_GW_number,config_iteration)=0;
                            end 
                        else
                            result_iterations_log_DL_RX1_schedulled_by_GW(:,config_iteration)=sum(calc_DL_packets_RX1_scheduled_ED_target~=0,2);
                        end

                        if(isempty(result_DL_RX1_packet_delivery_matrix)==false)
                            clearvars temp_nonzero_idx_DL_delivery_index temp_vector_result_DL_RX1_packet_delivery_matrix 
                            [temp_nonzero_idx_DL_RX1_delivery_index]=find(result_DL_RX1_packet_delivery_matrix(:)~=0);
                            temp_vector_result_DL_RX1_packet_delivery_matrix=result_DL_RX1_packet_delivery_matrix(:);

                            %DL packet delivery by DR
                            if(isempty(calc_DL_packets_RX1_scheduled_DR)==true)
                                for temp_DR_number=0:1:calc_num_DRs-1
                                    result_iterations_log_DL_RX1_delivered_by_DR(temp_DR_number+1,config_iteration)=0;
                                end
                            else
                                clearvars temp_vector_calc_DL_packets_RX1_scheduled_DR 
                                temp_vector_calc_DL_packets_RX1_scheduled_DR=calc_DL_packets_RX1_scheduled_DR(:);
                                for temp_DR_number=0:1:calc_num_DRs-1
                                    result_iterations_log_DL_RX1_delivered_by_DR(temp_DR_number+1,config_iteration)=sum(temp_vector_calc_DL_packets_RX1_scheduled_DR(temp_nonzero_idx_DL_RX1_delivery_index)==temp_DR_number,'all');
                                end          
                            end           

                            %DL packets delivery in RX1 by ED
                            if(isempty(calc_DL_packets_RX1_scheduled_ED_target)==true)
                                for temp_ED_number=1:1:input_num_EDs
                                    result_iterations_log_DL_RX1_delivered_by_ED(temp_ED_number,config_iteration)=0;
                                end
                            else
                                clearvars temp_vector_calc_DL_packets_RX1_scheduled_ED 
                                temp_vector_calc_DL_packets_RX1_scheduled_ED=calc_DL_packets_RX1_scheduled_ED_target(:);
                                for temp_ED_number=1:1:input_num_EDs
                                    result_iterations_log_DL_RX1_delivered_by_ED(temp_ED_number,config_iteration)=sum(temp_vector_calc_DL_packets_RX1_scheduled_ED(temp_nonzero_idx_DL_RX1_delivery_index)==temp_ED_number,'all');
                                end          
                            end
                            %DL packet delivery by GW
                            if(isempty(result_DL_RX1_packet_delivery_matrix)==true)
                                for temp_GW_number=1:1:input_num_GWs
                                    result_iterations_log_DL_RX1_delivered_by_GW(temp_GW_number,config_iteration)=0;
                                end 
                            else
                                result_iterations_log_DL_RX1_delivered_by_GW(:,config_iteration)=sum(result_DL_RX1_packet_delivery_matrix~=0,2);
                            end         
                        else
                            for temp_DR_number=0:1:calc_num_DRs-1
                                result_iterations_log_DL_RX1_delivered_by_DR(temp_DR_number+1,config_iteration)=0;
                            end             
                            for temp_ED_number=1:1:input_num_EDs
                                result_iterations_log_DL_RX1_by_ED(temp_ED_number,config_iteration)=0;
                            end          
                        end
                    else
                        for temp_GW_number=1:1:input_num_GWs
                            result_iterations_log_DL_RX1_schedulled_by_GW(temp_GW_number,config_iteration)=0;
                        end 
                        for temp_ED_number=1:1:input_num_EDs
                            result_iterations_log_DL_RX1_schedulled_by_ED(temp_ED_number,config_iteration)=0;
                        end    
                        for temp_DR_number=0:1:calc_num_DRs-1
                            result_iterations_log_DL_RX1_schedulled_by_DR(temp_DR_number+1,config_iteration)=0;
                        end            
                    end

                    %RX2
                    if(isempty(calc_DL_packets_RX2_scheduled_ED_target)==false)
                        [temp_nonzero_idx_DL_RX2_scheduled_index]=find(calc_DL_packets_RX2_scheduled_ED_target(:)~=0);

                        %DL packets scheduled in RX2 by DR
                        if(isempty(calc_DL_packets_RX2_scheduled_DR)==true)
                            for temp_DR_number=0:1:calc_num_DRs-1
                                result_iterations_log_DL_RX2_schedulled_by_DR(temp_DR_number+1,config_iteration)=0;
                            end
                        else
                            clearvars temp_vector_calc_DL_packets_RX2_scheduled_DR 
                            temp_vector_calc_DL_packets_RX2_scheduled_DR=calc_DL_packets_RX2_scheduled_DR(:);
                            for temp_DR_number=0:1:calc_num_DRs-1
                                result_iterations_log_DL_RX2_schedulled_by_DR(temp_DR_number+1,config_iteration)=sum(temp_vector_calc_DL_packets_RX2_scheduled_DR(temp_nonzero_idx_DL_RX2_scheduled_index)==temp_DR_number,'all');
                            end          
                        end

                        %DL packets scheduled in RX2 by ED
                        if(isempty(calc_DL_packets_RX2_scheduled_ED_target)==true)
                            for temp_ED_number=1:1:input_num_EDs
                                result_iterations_log_DL_RX2_schedulled_by_ED(temp_ED_number,config_iteration)=0;
                            end
                        else
                            clearvars temp_vector_calc_DL_packets_RX2_scheduled_ED 
                            temp_vector_calc_DL_packets_RX2_scheduled_ED=calc_DL_packets_RX2_scheduled_ED_target(:);
                            for temp_ED_number=1:1:input_num_EDs
                                result_iterations_log_DL_RX2_schedulled_by_ED(temp_ED_number,config_iteration)=sum(temp_vector_calc_DL_packets_RX2_scheduled_ED(temp_nonzero_idx_DL_RX2_scheduled_index)==temp_ED_number,'all');
                            end          
                        end

                        %DL packets scheduled in RX2 by GW
                        if(isempty(calc_DL_packets_RX2_scheduled_ED_target)==true)
                            for temp_GW_number=1:1:input_num_GWs
                                result_iterations_log_DL_RX2_schedulled_by_GW(temp_GW_number,config_iteration)=0;
                            end 
                        else
                            result_iterations_log_DL_RX2_schedulled_by_GW(:,config_iteration)=sum(calc_DL_packets_RX2_scheduled_ED_target~=0,2);
                        end

                        if(isempty(result_DL_RX2_packet_delivery_matrix)==false)
                            clearvars temp_nonzero_idx_DL_RX2_delivery_index temp_vector_result_DL_RX2_packet_delivery_matrix 
                            [temp_nonzero_idx_DL_RX2_delivery_index]=find(result_DL_RX2_packet_delivery_matrix(:)~=0);
                            temp_vector_result_DL_RX2_packet_delivery_matrix=result_DL_RX2_packet_delivery_matrix(:);

                            %DL packet delivery by DR
                            if(isempty(calc_DL_packets_RX2_scheduled_DR)==true)
                                for temp_DR_number=0:1:calc_num_DRs-1
                                    result_iterations_log_DL_RX2_delivered_by_DR(temp_DR_number+1,config_iteration)=0;
                                end
                            else
                                clearvars temp_vector_calc_DL_packets_RX2_scheduled_DR 
                                temp_vector_calc_DL_packets_RX2_scheduled_DR=calc_DL_packets_RX2_scheduled_DR(:);
                                for temp_DR_number=0:1:calc_num_DRs-1
                                    result_iterations_log_DL_RX2_delivered_by_DR(temp_DR_number+1,config_iteration)=sum(temp_vector_calc_DL_packets_RX2_scheduled_DR(temp_nonzero_idx_DL_RX2_delivery_index)==temp_DR_number,'all');
                                end          
                            end           

                            %DL packets delivery in RX2 by ED
                            if(isempty(calc_DL_packets_RX2_scheduled_ED_target)==true)
                                for temp_ED_number=1:1:input_num_EDs
                                    result_iterations_log_DL_RX2_delivered_by_ED(temp_ED_number,config_iteration)=0;
                                end
                            else
                                clearvars temp_vector_calc_DL_packets_RX2_scheduled_ED 
                                temp_vector_calc_DL_packets_RX2_scheduled_ED=calc_DL_packets_RX2_scheduled_ED_target(:);
                                for temp_ED_number=1:1:input_num_EDs
                                    result_iterations_log_DL_RX2_delivered_by_ED(temp_ED_number,config_iteration)=sum(temp_vector_calc_DL_packets_RX2_scheduled_ED(temp_nonzero_idx_DL_RX2_delivery_index)==temp_ED_number,'all');
                                end          
                            end
                            %DL packet delivery by GW
                            if(isempty(result_DL_RX2_packet_delivery_matrix)==true)
                                for temp_GW_number=1:1:input_num_GWs
                                    result_iterations_log_DL_RX2_delivered_by_GW(temp_GW_number,config_iteration)=0;
                                end 
                            else
                                result_iterations_log_DL_RX2_delivered_by_GW(:,config_iteration)=sum(result_DL_RX2_packet_delivery_matrix~=0,2);
                            end         
                        else
                            for temp_DR_number=0:1:calc_num_DRs-1
                                result_iterations_log_DL_RX2_delivered_by_DR(temp_DR_number+1,config_iteration)=0;
                            end             
                            for temp_ED_number=1:1:input_num_EDs
                                result_iterations_log_DL_RX2_by_ED(temp_ED_number,config_iteration)=0;
                            end          
                        end
                    else
                        for temp_GW_number=1:1:input_num_GWs
                            result_iterations_log_DL_RX2_schedulled_by_GW(temp_GW_number,config_iteration)=0;
                        end 
                        for temp_ED_number=1:1:input_num_EDs
                            result_iterations_log_DL_RX2_schedulled_by_ED(temp_ED_number,config_iteration)=0;
                        end    
                        for temp_DR_number=0:1:calc_num_DRs-1
                            result_iterations_log_DL_RX2_schedulled_by_DR(temp_DR_number+1,config_iteration)=0;
                        end  
                        for temp_GW_number=1:1:input_num_GWs
                            result_iterations_log_DL_RX2_delivered_by_GW(temp_GW_number,config_iteration)=0;
                        end 
                        for temp_ED_number=1:1:input_num_EDs
                            result_iterations_log_DL_RX2_delivered_by_ED(temp_ED_number,config_iteration)=0;
                        end    
                        for temp_DR_number=0:1:calc_num_DRs-1
                            result_iterations_log_DL_RX2_delivered_by_DR(temp_DR_number+1,config_iteration)=0;
                        end          
                    end

                    %"secondary" results
                    %coordinates
                    result_iterations_support_log_ED_coordinate_x(:,config_iteration)=calc_coord_ED_cartesian(:,1);
                    result_iterations_support_log_ED_coordinate_y(:,config_iteration)=calc_coord_ED_cartesian(:,2);
                    result_iterations_support_log_GW_coordinate_x(:,config_iteration)=calc_coord_GW_cartesian(:,1);
                    result_iterations_support_log_GW_coordinate_y(:,config_iteration)=calc_coord_GW_cartesian(:,2);
                    %DRs
                    result_iterations_support_log_ED_UL_DR(:,config_iteration)=calc_DRs_UL(:,1);
                    result_iterations_support_log_ED_DL_RX1_DR(:,config_iteration)=calc_DRs_RX1(:,1);
                    result_iterations_support_log_ED_DL_RX2_DR(:,config_iteration)=calc_DRs_RX2(:,1);

                    
                    %RMadd -> Comparison with analytical Ps [UL]
                    ED_DR_counter=zeros(1,6);
                    for i=1:1:input_num_EDs
                        tempDR=calc_DRs_UL(i)+1;
                        ED_DR_counter(tempDR)=ED_DR_counter(tempDR)+1;
                    end
                    ToA=input_UL_packets_duration_ms/1000;
                    for i=1:1:6
                        if(ED_DR_counter(i)==0)
                            result_analytical_Ps1_temp(i)=0;
                            result_analytical_Ps2_temp(i)=0;
                        else
                            result_analytical_Ps1_temp(i)=(1-(ToA(i)/input_T)).^(2*(ED_DR_counter(i)-1));
                            result_analytical_Ps2_temp(i)=exp(-2*ED_DR_counter(i)*ToA(i)/input_T);
                        end
                    end
                    
                    result_iterations_log_analytical_Ps1(:,config_iteration)=dot(ED_DR_counter,result_analytical_Ps1_temp)/sum(ED_DR_counter);
                    result_iterations_log_analytical_Ps2(:,config_iteration)=dot(ED_DR_counter,result_analytical_Ps2_temp)/sum(ED_DR_counter);
                               
                    %ENERGY MODEL
                    if input_energy_simulate==true
                        result_energy_consumption_EDs_total_mJ=zeros(input_num_EDs,1);
                        result_energy_consumption_EDs_UL_mJ=zeros(input_num_EDs,1);
                        result_energy_consumption_EDs_idleRWs_mJ=zeros(input_num_EDs,1);
                        result_energy_consumption_EDs_DL_mJ=zeros(input_num_EDs,1);
                        result_energy_consumption_EDs_sleep_mJ=zeros(input_num_EDs,1);
                        result_energy_consumption_EDs_active_time_ms=zeros(input_num_EDs,1);%all activities, except sleep
                        result_energy_consumption_EDs_sleep_time_ms=zeros(input_num_EDs,1); 
                        result_energy_consumption_EDs_idle_mJ=zeros(input_num_EDs,1);

                        for temp_ED_number=1:1:input_num_EDs
                            %basic
                            temp_num_UL_transmissions_delivered=result_iterations_log_UL_tranmissions_by_ED(temp_ED_number,config_iteration);
                            temp_UL_DR=calc_DRs_UL(temp_ED_number)+1;
                            temp_UL_TX_power_dBm=calc_UL_TX_power_matrix_dBm(temp_ED_number);
                            temp_UL_TX_power_index=find(input_UL_tx_power_dBm==temp_UL_TX_power_dBm);
                            temp_num_DL_received_in_RX1=result_iterations_log_DL_RX1_delivered_by_ED(temp_ED_number,config_iteration);
                            temp_DL_RX1_DR=calc_DRs_RX1(temp_ED_number)+1;
                            temp_num_DL_received_in_RX2=result_iterations_log_DL_RX2_delivered_by_ED(temp_ED_number,config_iteration);
                            temp_DL_RX2_DR=calc_DRs_RX2(temp_ED_number)+1;
                            
                            %cumulative durations of the different phases
                            temp_ED_UL_duration_ms=temp_num_UL_transmissions_delivered*input_UL_packets_duration_ms(temp_UL_DR);
                            
                            temp_ED_DL_useful_receiving_RX1_duration_ms=temp_num_DL_received_in_RX1*input_DL_packets_duration_ms(temp_DL_RX1_DR);
                            temp_ED_DL_idle_receiving_RX1_duration_ms=(temp_num_UL_transmissions_delivered-temp_num_DL_received_in_RX1)*input_energy_RX1_empty_duration_per_DR_ms(temp_DL_RX1_DR);          
                            
                            temp_ED_DL_useful_receiving_RX2_duration_ms=temp_num_DL_received_in_RX2*input_DL_packets_duration_ms(temp_DL_RX2_DR);
                            if input_DL_tranmission_in_one_RW_only==true %i.e., the ED skips RX2 if received in RX1
                                temp_ED_DL_idle_receiving_RX2_duration_ms=(temp_num_UL_transmissions_delivered-temp_num_DL_received_in_RX1-temp_num_DL_received_in_RX2)*input_energy_RX1_empty_duration_per_DR_ms(temp_DL_RX2_DR);
                            else%otherwise
                                temp_ED_DL_idle_receiving_RX2_duration_ms=(temp_num_UL_transmissions_delivered-temp_num_DL_received_in_RX2)*input_energy_RX1_empty_duration_per_DR_ms(temp_DL_RX2_DR);                   
                            end
                            temp_ED_sleep_duration_ms=input_simulation_duration_seconds*1000-temp_ED_UL_duration_ms-temp_ED_DL_useful_receiving_RX1_duration_ms-temp_ED_DL_idle_receiving_RX1_duration_ms-temp_ED_DL_useful_receiving_RX2_duration_ms-temp_ED_DL_idle_receiving_RX2_duration_ms;
                            %RM add -< energy consumption in RX_DELAY
                            temp_ED_idle_duration_ms=temp_num_DL_received_in_RX1*input_RX1_delay_s*1000+temp_num_DL_received_in_RX2*(input_RX1_delay_s+input_RX2_delay_s)*1000;
                            
                            %cumulative energy consumption of the different phases
                            temp_ED_UL_consumption_mJ=temp_ED_UL_duration_ms*input_energy_current_transmit_per_DR_mA(temp_UL_TX_power_index)*input_energy_voltage_V/1000;% WRONG ms*mA*V*1000=uJ*1000=mJ
                            temp_ED_DL_consumption_RX1_useful_mJ=temp_ED_DL_useful_receiving_RX1_duration_ms*input_energy_current_receive_mA*input_energy_voltage_V/1000;% CORRECT ms*mA*V/1000=mJ
                            temp_ED_DL_consumption_RX1_idle_mJ=temp_ED_DL_idle_receiving_RX1_duration_ms*input_energy_current_receive_mA*input_energy_voltage_V/1000;
                            temp_ED_DL_consumption_RX2_useful_mJ=temp_ED_DL_useful_receiving_RX2_duration_ms*input_energy_current_receive_mA*input_energy_voltage_V/1000;
                            temp_ED_DL_consumption_RX2_idle_mJ=temp_ED_DL_idle_receiving_RX2_duration_ms*input_energy_current_receive_mA*input_energy_voltage_V/1000;
                            temp_ED_sleep_consumption_mJ=temp_ED_sleep_duration_ms*input_energy_current_sleep_mA*input_energy_voltage_V/1000;
                            temp_ED_idle_consumption_mJ=temp_ED_idle_duration_ms*input_energy_current_RXDelay_mA*input_energy_voltage_V/1000;
                            
                            %store to results
                            result_energy_consumption_EDs_active_time_ms(temp_ED_number,1)=input_simulation_duration_seconds*1000-temp_ED_sleep_duration_ms;
                            result_energy_consumption_EDs_sleep_time_ms(temp_ED_number,1)=temp_ED_sleep_duration_ms;
                            result_energy_consumption_EDs_UL_mJ(temp_ED_number,1)=temp_ED_UL_consumption_mJ;
                            result_energy_consumption_EDs_idleRWs_mJ(temp_ED_number,1)=temp_ED_DL_consumption_RX1_idle_mJ+temp_ED_DL_consumption_RX2_idle_mJ;
                            result_energy_consumption_EDs_DL_mJ(temp_ED_number,1)=temp_ED_DL_consumption_RX1_useful_mJ+temp_ED_DL_consumption_RX2_useful_mJ;
                            result_energy_consumption_EDs_sleep_mJ(temp_ED_number,1)=temp_ED_sleep_consumption_mJ;
       
                            tempEC_RX_Delay1=input_RX1_delay_s*1000*input_energy_current_RXDelay_mA*input_energy_voltage_V/1000;

                            result_energy_consumption_EDs_total_mJ(temp_ED_number,1)=tempEC_RX_Delay1+temp_ED_UL_consumption_mJ+temp_ED_DL_consumption_RX1_idle_mJ+temp_ED_DL_consumption_RX2_idle_mJ+temp_ED_DL_consumption_RX1_useful_mJ+temp_ED_DL_consumption_RX2_useful_mJ+temp_ED_sleep_consumption_mJ;
                            result_energy_consumption_EDs_total_mJ_nosleep(temp_ED_number,1)=tempEC_RX_Delay1+temp_ED_UL_consumption_mJ+temp_ED_DL_consumption_RX1_idle_mJ+temp_ED_DL_consumption_RX1_useful_mJ;                                                       
                            result_EC_EDs_single_mJ(temp_ED_number,1)=result_energy_consumption_EDs_total_mJ_nosleep(temp_ED_number,1)/temp_num_UL_transmissions_delivered;                            
                            tempEC_UL=input_UL_packets_duration_ms(temp_UL_DR)*input_energy_current_transmit_per_DR_mA(temp_UL_TX_power_index)*input_energy_voltage_V/1000;
                            tempEC_DL_full=input_DL_packets_duration_ms(temp_DL_RX1_DR)*input_energy_current_receive_mA*input_energy_voltage_V/1000;
                            tempEC_DL_empty=input_energy_RX1_empty_duration_per_DR_ms(temp_DL_RX1_DR)*input_energy_current_receive_mA*input_energy_voltage_V/1000; 
                            temp_delivery_rate_UL=result_iterations_log_UL_deliveries_by_ED(temp_ED_number, config_iteration)/result_iterations_log_UL_tranmissions_by_ED(temp_ED_number, config_iteration);
                            temp_delivery_rate_DL_RX1=result_iterations_log_DL_RX1_delivered_by_ED(temp_ED_number,config_iteration)/result_iterations_log_DL_RX1_schedulled_by_ED(temp_ED_number,config_iteration);
                            
                            %Energy model
                            result_analytical_EC_Ps1(temp_ED_number,1)=tempEC_UL+tempEC_RX_Delay1+result_analytical_Ps1_temp(temp_UL_DR)*tempEC_DL_full+(1-result_analytical_Ps1_temp(temp_UL_DR))*tempEC_DL_empty;
                            result_analytical_EC_Ps2(temp_ED_number,1)=tempEC_UL+tempEC_RX_Delay1+result_analytical_Ps2_temp(temp_UL_DR)*tempEC_DL_full+(1-result_analytical_Ps2_temp(temp_UL_DR))*tempEC_DL_empty;
                        end     
                    
                        result_iterations_log_consumption_ED_active_time_ms(:,config_iteration)=result_energy_consumption_EDs_active_time_ms(:,1);
                        result_iterations_log_consumption_ED_sleep_time_ms(:,config_iteration)=result_energy_consumption_EDs_sleep_time_ms(:,1);
                        result_iterations_log_consumption_ED_UL_consumption_mJ(:,config_iteration)=result_energy_consumption_EDs_UL_mJ(:,1);
                        result_iterations_log_consumption_ED_DL_useful_RX_mJ(:,config_iteration)=result_energy_consumption_EDs_DL_mJ(:,1);
                        result_iterations_log_consumption_ED_DL_idle_RX_mJ(:,config_iteration)=result_energy_consumption_EDs_idleRWs_mJ(:,1);
                        result_iterations_log_consumption_ED_sleep_mJ(:,config_iteration)=result_energy_consumption_EDs_sleep_mJ(:,1);
                        
                        result_iterations_log_consumption_ED_idle_mJ(:, config_iteration)=result_energy_consumption_EDs_idle_mJ(:,1);
                        result_iterations_log_consumption_ED_total_mJ(:,config_iteration)=result_energy_consumption_EDs_total_mJ(:,1);                      
                        result_iterations_log_consumption_ED_single_mJ(:,config_iteration)=result_EC_EDs_single_mJ(:,1);
                        result_iterations_log_consumption_ED_analytical_Ps1_mJ(:,config_iteration)=result_analytical_EC_Ps1(:,1);
                        result_iterations_log_consumption_ED_analytical_Ps2_mJ(:,config_iteration)=result_analytical_EC_Ps2(:,1);
                    end
                                      
                    %notification
                    calc_timestamp_interation_ended=datetime('now');
                    fprintf('Iteration %i of %i finished at %s\n', config_iteration,config_num_iterations,calc_timestamp_interation_ended); 
                    calc_time_since_start_of_experiment=calc_timestamp_interation_ended-config_timestamp_started;
                    calc_time_average_per_iteration=calc_time_since_start_of_experiment/config_iteration;
                    calc_time_left=calc_time_average_per_iteration*(config_num_iterations-config_iteration);
                    calc_timestamp_experiment_end_estimation=datetime('now')+calc_time_left;
                    fprintf('Experiments estimated to end at %s\n', calc_timestamp_experiment_end_estimation);     
                end%loop over iterations
                
                %ALL ITERATIONS FINISHED - STATISTICAL PROCESSING OF RESULTS
                
                if(isempty(result_iterations_log_UL_deliveries_by_ED)==false && isempty(result_iterations_log_UL_tranmissions_by_ED)==false)
                    results_calc_all_iterations_UL_delivery_rate=sum(result_iterations_log_UL_deliveries_by_ED,'all')/sum(result_iterations_log_UL_tranmissions_by_ED,'all');
                else 
                end
                if(isempty(result_iterations_log_UL_tranmissions_by_ED)==false)
                    results_calc_all_iterations_UL_mean_transmissions=mean(result_iterations_log_UL_tranmissions_by_ED,'all');
                else
                    results_calc_all_iterations_UL_mean_transmissions=0;
                end                
                if(isempty(result_iterations_log_DL_RX1_delivered_by_ED)==false && isempty(result_iterations_log_DL_RX1_schedulled_by_ED)==false)
                    results_calc_all_iterations_DL_RX1_delivery_rate=sum(result_iterations_log_DL_RX1_delivered_by_ED,'all')/sum(result_iterations_log_DL_RX1_schedulled_by_ED,'all');
                else
                    results_calc_all_iterations_DL_RX1_delivery_rate=-1;
                end                
                if(isempty(result_iterations_log_DL_RX1_schedulled_by_ED)==false)
                    results_calc_all_iterations_DL_RX1_mean_transmissions=mean(result_iterations_log_DL_RX1_schedulled_by_ED,'all');
                else
                    results_calc_all_iterations_DL_RX1_mean_transmissions=0;
                end                
                if(isempty(result_iterations_log_DL_RX2_delivered_by_ED)==false && isempty(result_iterations_log_DL_RX2_schedulled_by_ED)==false)
                    results_calc_all_iterations_DL_RX2_delivery_rate=sum(result_iterations_log_DL_RX2_delivered_by_ED,'all')/sum(result_iterations_log_DL_RX2_schedulled_by_ED,'all');
                else
                    results_calc_all_iterations_DL_RX2_delivery_rate=-1;
                end                
                if(isempty(result_iterations_log_DL_RX2_schedulled_by_ED)==false)
                    results_calc_all_iterations_DL_RX2_mean_transmissions=mean(result_iterations_log_DL_RX2_schedulled_by_ED,'all');
                else
                    results_calc_all_iterations_DL_RX2_mean_transmissions=0;
                end                
                if(isempty(result_iterations_log_DL_RX1_delivered_by_ED)==false)
                    temp_RX1_deliveries=sum(result_iterations_log_DL_RX1_delivered_by_ED,'all');
                else
                    temp_RX1_deliveries=0;
                end
                if(isempty(result_iterations_log_DL_RX2_delivered_by_ED)==false)
                    temp_RX2_deliveries=sum(result_iterations_log_DL_RX2_delivered_by_ED,'all');
                else
                    temp_RX2_deliveries=0;
                end                
                if(isempty(result_iterations_log_UL_tranmissions_by_ED)==false)
                    results_calc_all_iterations_DL_response_rate=(temp_RX1_deliveries+temp_RX2_deliveries)/sum(result_iterations_log_UL_tranmissions_by_ED,'all');
                else
                    results_calc_all_iterations_DL_response_rate=-1;
                end
                results_calc_all_iterations_analytical_Ps1=mean(result_iterations_log_analytical_Ps1, 'all');
                results_calc_all_iterations_analytical_Ps2=mean(result_iterations_log_analytical_Ps2, 'all');
                %ENERGY
                if input_energy_simulate==true
                    results_calc_all_iterations_energy_mean_active_time_ms=mean(result_iterations_log_consumption_ED_active_time_ms,'all');
                    results_calc_all_iterations_energy_mean_sleep_time_ms=mean(result_iterations_log_consumption_ED_sleep_time_ms,'all');
                    results_calc_all_iterations_energy_mean_UL_consumption_mJ=mean(result_iterations_log_consumption_ED_UL_consumption_mJ,'all');
                    results_calc_all_iterations_energy_mean_DL_full_consumption_mJ=mean(result_iterations_log_consumption_ED_DL_useful_RX_mJ,'all');
                    results_calc_all_iterations_energy_mean_DL_idle_consumption_mJ=mean(result_iterations_log_consumption_ED_DL_idle_RX_mJ,'all');
                    results_calc_all_iterations_energy_mean_sleep_consumption_mJ=mean(result_iterations_log_consumption_ED_sleep_mJ,'all');
                    results_calc_all_iterations_energy_mean_total_consumption_mJ=mean(result_iterations_log_consumption_ED_total_mJ,'all');                 
                    results_calc_all_iterations_energy_mean_single_mJ=mean(result_iterations_log_consumption_ED_single_mJ, 'all');
                    results_calc_all_iterations_energy_mean_analytical_Ps1_mJ=mean(result_iterations_log_consumption_ED_analytical_Ps1_mJ, 'all');
                    results_calc_all_iterations_energy_mean_analytical_Ps2_mJ=mean(result_iterations_log_consumption_ED_analytical_Ps2_mJ, 'all');
                end

                    %save workspace
                    temp_filename = sprintf('LoRaWANSim_GW_%i_ED_%i_R_%i_T_%i_CR_%i_iter_%i.mat',input_num_GWs,input_num_EDs,input_R,input_T,input_CR,config_num_iterations);
                    save(temp_filename);    
                    end
                end
            end
        end
    end
end

%optionally: draw the map
if input_draw_playground_map==1
   set(gcf, 'Color', 'w');
   screen=get(0,'ScreenSize');
   min_screen_dimension=min(screen(3:4));
   set(gcf, 'Units', 'Pixels', 'OuterPosition', [50 50 min_screen_dimension-50 min_screen_dimension-50]);
   ylabel('Coordinate Y, m', 'FontSize', 28)
   xlabel('Coordinate X, m', 'FontSize', 28)
   %GWs
   fig_scatter=scatter(calc_coord_GW_cartesian(:,1),calc_coord_GW_cartesian(:,2));
   fig_scatter.Marker='d';       
   for temp_i_num_GW=1:input_num_GWs 
       txt = ['GW' num2str(temp_i_num_GW)];
       text(calc_coord_GW_cartesian(temp_i_num_GW,1),calc_coord_GW_cartesian(temp_i_num_GW,2),txt)
   end
   hold on
   %EDs
   fig_scatter=scatter(calc_coord_ED_cartesian(:,1),calc_coord_ED_cartesian(:,2));
   fig_scatter.Marker='+';    
   for temp_i_num_ED=1:input_num_EDs 
       txt = ['ED' num2str(temp_i_num_ED)];
       text(calc_coord_ED_cartesian(temp_i_num_ED,1),calc_coord_ED_cartesian(temp_i_num_ED,2),txt)
   end
   hold on       
end

%print results
fprintf('All experiments finished!\n');
fprintf('Uplink Delivery Rate: %i\n', results_calc_all_iterations_UL_delivery_rate);
fprintf('Energy Consumption: %i [mJ]\n', results_calc_all_iterations_energy_mean_total_consumption_mJ);