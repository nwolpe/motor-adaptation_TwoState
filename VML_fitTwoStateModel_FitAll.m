% a script to fit Smith et al (2006) model of two learning processes - fast
% and slow - to CamCAN VML data
clear

% cd ~/Documents/VML_TwoStateFit/
load /imaging/nw03/CamCAN/VML_TwoState/SubjectTrialTrajectoryError.mat

paramList_twoState = [];
GoF_twoState = [];
R2_twoState = [];
AIC_twoState = [];

n_trials = [24, 120, 48]; % [preexposure, exposure, postexposure]
ifcycle = 1;

% set up options for bads
nonbcon = @(x) (x(:,2) - x(:,1) <= 0) | (x(:,3) - x(:,4) <= 0);
lb = [0.001, 0.001, 0.001, 0.001];
ub = [0.999, 0.999, 0.999, 0.999];

options_bads = bads('defaults'); %options_bads.Display = 'final';

% run simulation on McDougle et al. data to generate plausible bounds
% rng default
% m=0.05; CI95=0.02;
% s=CI95*sqrt(9)/1.960;
% truncate = @(m,s,u,l) trandn((l-m)/s, (u-m)/s)*s + m; % m=mean, s=std
% numIter = 10000;
% for iteration=1:numIter
%     simX(iteration) = truncate(m,s,1,0);
% end
% plb_sim = min(simX); pub_sim = max(simX);

Af_plb = 0.3918; Af_pub = 0.999;
As_plb = 0.9327; As_pub = 0.999;
Bf_plb = 0.0963; Bf_pub = 0.7686;
Bs_plb = 0.001; Bs_pub = 0.1752;
plb = [Af_plb, As_plb, Bf_plb, Bs_plb];
pub = [Af_pub, As_pub, Bf_pub, Bs_pub];
    % parameter starting values from McDougle et al. JoN:
    Af = 0.85; As = 0.99; Bf = 0.44; Bs = 0.05; 

    % tunning parameters for fminsearch
%     options=optimset('fmincon');
%     options.Display='final';
%     options.MaxFunEvals=10000;
%     options.MaxIter=20000;
    
%     b = [0; 0];
%     A = [1, -1, 0, 0; 0, 0, -1, 1];
%     

rng default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start subject loop

for iSub = 1:length(SubjectTrialTrajectoryError)
    
    clearvars data2fit x1 x2 e subjectData dataAll error_test R2_test
    subjectData = SubjectTrialTrajectoryError(iSub,:);
    
    if ifcycle
        for iCycle = 1:length(subjectData)/4
            dataAll(iCycle) = nanmedian(subjectData((iCycle-1)*4+1:iCycle*4));
        end
        n = n_trials./4; % if we want it cycle-wise - divide by 4
        dataStruc = 'cycles';
    else
        dataAll = subjectData;
        dataStruc = 'trials';
        n = n_trials;
    end
    
    ntot=sum(n);
    data_train = dataAll(1:sum(n(1:2)));
    data_test = dataAll(sum(n(1:2))+1:ntot);
    
    p_train=[zeros(1,n(1)) 30.*ones(1,n(2))];
    p_test = zeros(1,n(3)); % perturbation=30 degrees in exposure, and otherwise zero
    
    x1(1)=0; x2(1)=0; e(1)=0;
%     [params, SSE] = fmincon('sumEvidenceTwoState', [Af, As, Bf, Bs],...
%         A, b, [], [], lb, ub, [], options, {p_train, data_train, 'train', dataStruc});
%     
        
        forDisplay = sprintf('now running subject %d / %d ...', iSub, length(SubjectTrialTrajectoryError));
        display(forDisplay)
        
        [params_bads, RMSE_bads] = bads('sumEvidenceTwoState', [Af, As, Bf, Bs],...
            lb, ub, plb, pub, nonbcon, options_bads, {[p_train, p_test], [data_train, data_test],...
            'all', dataStruc});

        Af_all = params_bads(1); As_all = params_bads(2);
        Bf_all = params_bads(3); Bs_all = params_bads(4);
        p_all = [p_train, p_test];
        for i=0:sum(n(1:3))-1
            if i==0
                x1(i+1)=Af_all*0 + Bf_all*0;
                x2(i+1)=As_all*0 + Bs_all*0;
                x(i+1)=x1(i+1)+x2(i+1);
                e(i+1)=p_all(i+1)-x(i+1);
            else
                x1(i+1)=Af_all*x1(i)+Bf_all*e(i);
                x2(i+1)=As_all*x2(i)+Bs_all*e(i);
                x(i+1)=x1(i+1)+x2(i+1);
                e(i+1)=p_all(i+1)-x(i+1);
            end
        end
        
        
%         AIC = n*log(RSS/n) + 2*k
        [R2_test, SSE] = feval(@Rsquared_TwoState, params_bads, {[p_train, p_test], [data_train, data_test],...
            'all', dataStruc});
        AIC_test = length(data_test)*log(SSE/length(data_test)) + 2*4;
        
        
        AIC_twoState = [AIC_twoState; AIC_test];
        GoF_twoState = [GoF_twoState, RMSE_bads];
        R2_twoState = [R2_twoState; R2_test];
        paramList_twoState = [paramList_twoState; params_bads];
        
end



%% plot parameters

scatter(paramList(:,1), paramList(:,2), '.', 'b')
% scatter(paramList(:,3), paramList(:,4), '.', 'r')
% plot(GoF)

%%     % plot subject data against model


for iSubject = 1:length(paramList)

    clf; subjectNum = iSubject;
    Af = paramList(subjectNum,1); As = paramList(subjectNum,2); 
    Bf = paramList(subjectNum,3); Bs = paramList(subjectNum,4);
    for i=1:(ntot-1)
        x1(i+1)=Af(1)*x1(i)+Bf(1)*e(i);
        x2(i+1)=As(1)*x2(i)+Bs(1)*e(i);
        x(i+1)=x1(i+1)+x2(i+1);
        e(i+1)=p(i+1)-x(i+1);
    end
    
    subjectData = SubjectTrialTrajectoryError(subjectNum,:);
    
    for iCycle = 1:length(subjectData)/4
        data_cycle(iCycle) = nanmean(subjectData((iCycle-1)*4+1:iCycle*4));
    end
    
    plot(e, 'b'); hold on; plot(data_cycle, 'r'); h=legend('Model', 'Data'); set(h,'Location','SouthWest');
    xlabel('Cycle'); ylabel('Trajectory Error (degrees)');
    
    pause
end
    
    %% Compare goodness of fit with one state model
    
%     A = Af;
%     B = Bf;
%     [params, SSE] = fminsearch('sumEvidenceOneState', [A, B],...
%         options, {p, x_data})


% save training and test set

