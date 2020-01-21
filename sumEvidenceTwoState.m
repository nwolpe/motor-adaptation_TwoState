function RMSE = sumEvidenceTwoState(paramsToFit, data)

Af = paramsToFit(1);
As = paramsToFit(2);
Bf = paramsToFit(3);
Bs = paramsToFit(4);

p = data{1};
subjectData = data{2};
fit_mode = data{3};

if strcmp(fit_mode, 'train') % trials for training set
    n_trials = [24 120];
    x1_start=0; x2_start=0; e_start=0;
elseif strcmp(fit_mode, 'test') % trials for test set - deadapt
    n_trials = 48;
    x1_start=data{5}(1); x2_start=data{5}(2); e_start=data{5}(3);
elseif strcmp(fit_mode, 'all') 
    n_trials = [24 120 48];
    x1_start=0; x2_start=0; e_start=0;
end

if strcmp(data{4}, 'cycles')
    n = n_trials./4;
elseif strcmp(data{4}, 'trials')
    n = n_trials;
end

ntot=sum(n);

% for i=1:(ntot-1)
%     x1(i+1)=Af(1)*x1(i)+Bf(1)*e(i);
%     x2(i+1)=As(1)*x2(i)+Bs(1)*e(i);
%     x(i+1)=x1(i+1)+x2(i+1);
%     e(i+1)=p(i+1)-x(i+1);
% end


for i=0:(ntot-1)
    if i==0
        x1(i+1)=Af*x1_start + Bf*e_start;
        x2(i+1)=As*x2_start + Bs*e_start;
        x(i+1)=x1(i+1)+x2(i+1);
        e(i+1)=p(i+1)-x(i+1);
    else
        x1(i+1)=Af(1)*x1(i)+Bf(1)*e(i);
        x2(i+1)=As(1)*x2(i)+Bs(1)*e(i);
        x(i+1)=x1(i+1)+x2(i+1);
        e(i+1)=p(i+1)-x(i+1);
    end
end


RMSE = sqrt(nanmean((subjectData - e).^2));

% SSE = nansum((subjectData - e).^2); %SSE

% evidence = -corr(e', subjectData', 'rows','complete');
