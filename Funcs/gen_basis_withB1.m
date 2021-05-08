function [U, X0] = gen_basis_withB1(N1, N2, T2svals, T2vals, T1vals,B1_factor, param)
% Generate basis
% Inputs:
%  N1 -- maximum number of T2* signals to simulate
%  N2 -- maximum number of T2 signals to simulate
%  T2svals, T2vals, T1vals (s) -- array of T2* T2 T1 values to simulate
%  B1_factor -- array of B1 values to simulate
%  param -- acquisition parameters of the sequence
%
% Outputs:
%  U -- temporal basis based on PCA
%  X0 -- [T, L] matrix of simulated signals
% randomly choose T2 values if more than N are given

if length(T2svals) > N1
    idx = randperm(length(T2svals));
    T2svals = T2svals(idx(1:N1));
end
if length(T2vals) > N2
    idx = randperm(length(T2vals));
    T2vals = T2vals(idx(1:N2));
end

LT1 = length(T2svals);
LT2 = length(T2vals);
LT3 = length(T1vals);
LT4 = length(B1_factor);

X0 = zeros(param.N_signal, LT1*LT2*LT3*LT4);
num=1;
for nb1 = 1:LT4
    disp(['----- Processing B1: ',num2str(nb1),'/',num2str(LT4),' ------']);
    param_tmp = param;
    mask_nonIR = (param_tmp.param1.FA_array~=deg2rad(180));
    param_tmp.param1.FA_array(mask_nonIR) = param_tmp.param1.FA_array(mask_nonIR)*B1_factor(nb1);
    param_tmp.param2.FA_array = param_tmp.param2.FA_array*B1_factor(nb1);
    
    for ii=1:LT3
        disp(['Processing: ',num2str(ii),'/',num2str(LT3)]);
        T1=T1vals(ii);
        for jj=1:LT2
            R2 = 1/T2vals(jj);
            for kk=1:LT1
                R2s = 1/T2svals(kk);
                if (R2s>=R2) && (R2s/R2<6)
                    R2_p=R2s-R2;
                        [signal_a,~] = EPG_IRSPGR(T1,R2,R2_p,param_tmp.param1,0);
                        [signal_b,~] = EPG_TSE_addSpoiler(T1,R2,R2_p,param_tmp.param2,0);
                        signal1 = cat(1,-real(signal_a),real(signal_b));
                    X0(:,num) = signal1;
                    num=num+1;
                end
            end
        end
    end
end
X0=X0(:,1:num-1);

[U, ~, ~] = svd(X0, 'econ');

end