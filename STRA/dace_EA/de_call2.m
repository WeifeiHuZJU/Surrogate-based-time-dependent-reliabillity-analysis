function [bestvalue,bestind]=de_call2(lb,ub,par)

lu=[lb;ub];
n=size(lb,2);

% Main body
popsize = 30;
MaxIt=12;
% Initialize the main population
p = repmat(lu(1, :), popsize, 1) + rand(popsize, n) .* (repmat(lu(2, :) - lu(1, :), popsize, 1));

fit=zeros(popsize,1);
for i=1:popsize
%     [fit(i,:),~,~]=Kriging_likelihood(p(i,:),X_sample_likelihood,Y_sample_likelihood,X_norm_flag,Y_norm_flag);
    [fit(i,:),~]=objfunc(p(i,:), par);
end
save de_call2
% Record the number of function evaluations (FES)
FES = popsize;
gen=1;
[valBest, indBest] = sort(fit, 'ascend');
bestind=p(indBest(1),:);
bestvalue=valBest(1,:);
% Show Iteration Information
disp(['Iteration ' num2str(gen) ': Best Cost = ' num2str(bestvalue)]);
while gen < MaxIt
    
    pTemp = p;
    fitTemp = fit;
    
    % uSet: the set of trial vectors
    uSet = zeros(3 * popsize, n);
    
    for i = 1 : popsize
        
        % the three control parameter settings
        F    = [1.0 1.0 0.8];
        CR = [0.1 0.9 0.2];
        
        % Uniformly and randomly select one of the control
        % parameter settings for each trial vector generation strategy
        paraIndex = floor(rand(1, 3) * length(F)) + 1;
        
        % Generate the trail vectors
        u = generator(p, lu, i, F, CR, popsize, n, paraIndex,bestind);
        
        uSet(i * 3 - 2 : 3 * i, :) = u;
        
        FES = FES + 3;
        
    end
    
    % Evaluate the trial vectors
    fitSet = zeros(size(uSet,1),1);
    for i=1:size(uSet,1)
%         [fitSet(i,:),~,~]=Kriging_likelihood(uSet(i,:),X_sample_likelihood,Y_sample_likelihood,X_norm_flag,Y_norm_flag);
        [fitSet(i,:),~]=objfunc(uSet(i,:), par);
    end
    
    for i = 1 : popsize
        
        % Choose the best trial vector from the three trial vectors
        [minVal, minID] = min(fitSet(3 * i - 2 : 3 * i, :));
        bestInd = uSet(3 * (i - 1) + minID, :);
        bestIndFit = fitSet(3 * (i - 1) + minID, :);
        
        % Choose the better one between the trial vector and the
        % target vector
        if fit(i) >= bestIndFit
            pTemp(i, :) = bestInd;
            fitTemp(i, :) = bestIndFit; 
        end
        
    end
    
    
    p = pTemp;
    fit = fitTemp;
    [valBest, indBest] = sort(fit, 'ascend');
    bestind=p(indBest(1),:);
    bestvalue=valBest(1,:);
    gen=gen+1;
    % Show Iteration Information
    disp(['Iteration ' num2str(gen) ': Best Cost = ' num2str(bestvalue)]);
    
    
end

end
