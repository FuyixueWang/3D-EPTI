function [signals_steady,signals_TE0] = EPG_IRSPGR(T1,R2,R2_prime,param,show)
if nargin<5
    show=0;
end
%% simulation parameters
Nstates = param.N_GRE;	% Number of states to simulate (easier to keep constant if showing state evolution.
T2 = 1/R2;	    % unit: s

P = zeros(3,Nstates);	% State matrix
P(3,1)=1;		% Equilibrium magnetization. 

%%
signals_TE0 = zeros(param.N_sampling,1);
signals_GE = zeros(param.nechoGE,param.N_sampling);
signals_all=[];
num_sample = 1;
for k = 1:param.N_events
      P = epg_rf(P,param.FA_array(k),param.RF_phase_array(k));		% RF pulse
      if param.mask_sampling(k) == 0
         P(1:2,:) = 0;
      end
      if param.mask_sampling(k)==1
          signals_TE0(num_sample) = P(1,1)*exp(-1i*(param.RF_phase_array(k)+pi/2));      			% Signal is F0 state.
          signals_GE(:,num_sample) = signals_TE0(num_sample).*exp(-param.TEs_GE.*(R2_prime+R2));	% Phase-Demodulated signal.
          signals_all=cat(1,signals_all,signals_GE(:,num_sample));
          num_sample=num_sample+1;
      end
      % -- Simulate relaxation and spoiler gradient
      if k<param.N_events
          P = epg_grelax(P,T1,T2,param.Time_array(k+1)-param.Time_array(k),1,0,1,1);   % spoiler gradient, relaxation.
      end
      P(1:2,:) = 0;
end
%% Plot of all signals
% signals_all=-real(signals_all);
% signals_TE0(:)=-signals_TE0(:);
if show==1
    figure;
    plot(param.Time_sampling,signals_all(:),'b-');
%     scatter(param.Time_sampling,signals_all(:));

    hold on;
    plot_index_IR=param.Time_array(1:param.N_GRE+1:end);
    plot_index_IR=[plot_index_IR';plot_index_IR'];
    min_show = min(signals_all(:));
    max_show = max(signals_all(:));
    plot(plot_index_IR,[min_show, max_show],'r--','linewidth',2);
    xlabel('time (s)');
    ylabel('Signal at TE = 0');
    legend('EPG Signal','IR pulses');
    title('GE Evolution');
    grid on; 
end
signals_all=signals_all(:);
signals_steady=signals_all(end-param.N_signal+1:end);

end

