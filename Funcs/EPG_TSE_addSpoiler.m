function [signals_steady,signals_TE0] = EPG_TSE_addSpoiler(T1,R2,R2_prime,param,show)
if nargin<5
    show=0;
end
%% simulation parameters
Nstates = 64;	% Number of states to simulate (easier to keep constant if showing state evolution.
T2 = 1/R2;	    % unit: s

P = zeros(3,Nstates);	% State matrix
P(3,1)=1;		% Equilibrium magnetization. 

np_SE_beforeTE = sum(param.TEs_SE <= param.Time_TR_TSE/2);
%%
signals_TE0 = zeros(param.N_sampling,1);
signals_SE = zeros(param.nechoSE,param.N_sampling);
signals_all=[];
num_sample = 1;
for k = 1:param.N_events     
        if param.mask_sampling(k)==0
            P(1:2,:) = 0;
            P = epg_rf(P, param.FA_array(k), param.RF_phase_array(k));
        else
            P = epg_grelax(P, T1, T2, param.Time_TR_TSE/2);
            P = epg_rf(P, param.FA_array(k), param.RF_phase_array(k));
            P = epg_grelax(P, T1, T2, param.Time_TR_TSE/2);
            
            signals_TE0(num_sample) = P(1,1)*exp(-1i*(param.RF_phase_array(k)));  	% Signal is F0 state.
            signals_SE(1:np_SE_beforeTE,num_sample) =signals_TE0(num_sample).*exp(-(param.TEs_SE(1:np_SE_beforeTE)-param.Time_TR_TSE/2).*(R2-R2_prime));	% Phase-Demodulated signal.
            signals_SE(np_SE_beforeTE+1:end,num_sample) = signals_TE0(num_sample).*exp(-(param.TEs_SE(np_SE_beforeTE+1:end)-param.Time_TR_TSE/2).*(R2_prime+R2));	% Phase-Demodulated signal.
            signals_all=cat(1,signals_all,signals_SE(:,num_sample));
            num_sample=num_sample+1;
            
            if k < param.N_events
                if param.mask_sampling(k+1)==0
                    P = epg_grelax(P,T1,T2,param.Time_array(k+1)-param.Time_array(k)-param.Time_TR_TSE/2,1,0,0,1);   % spoiler gradient, relaxation.
                    P(1:2,:) = 0;
                end 
            end
        end
end
%% Plot of all signals
% signals_all= real(signals_all);
% signals_all = abs(signals_all);  % let image recon to figure out phase 20191104
if show==1
    figure;
    plot(param.Time_sampling,signals_all(:),'b-');
%     scatter(param.Time_sampling,signals_all(:));

    hold on;
    plot_index_IR=param.Time_array(1:param.N_TSE+1:end);
    plot_index_IR=[plot_index_IR';plot_index_IR'];
    min_show = min(signals_all(:));
    max_show = max(signals_all(:));
    plot(plot_index_IR,[min_show, max_show],'r--','linewidth',2);
    xlabel('time (s)');
    ylabel('Signal at TE = 0');
    legend('EPG Signal','Pulses');
    title('TSE Evolution');
    grid on; 
end
signals_all=signals_all(:);
signals_steady=signals_all(end-param.N_signal+1:end);
end

