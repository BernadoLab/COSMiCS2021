function [chi, fit] = compare2curves(curvaExp, curvaSint);

    % [chi, fit] = compare2curves(curvaExp, curvaSint)
    % The curves must have the same number of points
    % curvaExp(:,3) - S | I| E
    % curvaSint(:,1) or curvaSint(:,2). The last column is the I value, S
    % must have the same scale than the experimental curve

    EPS = 1e-50;


    % Trick to speed up calculation xi2
    weif(:,1) = 1./max(curvaExp(:,3).^2,EPS);

    % Scaling coefficient
    sum1 = 0;
    sum2=0;
    for ind=1:length(curvaExp)
        sum1 = sum1 + (curvaExp(ind,2) * curvaSint(ind,end) * weif(ind,1));
        sum1Parcial(ind,1) = sum1;
        sum2 = sum2 + (curvaSint(ind,end) *curvaSint(ind,end) * weif(ind,1));
        %sum2 = sum(sum2P);
    end
    ratio = sum1/(max(sum2,EPS));

    % Evaluate Chi-square
    chi_o = 0;
    for ind = 1:length(curvaExp)
        fit(ind,1) = curvaSint(ind,end) .* ratio;
        chi_o = chi_o + ( (curvaExp(ind,2)-fit(ind,1)).^2 .* weif(ind,1) );
    end
    chi = chi_o./(length(curvaExp)-1);
end
