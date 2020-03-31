function plot_score(lambdas, scores)
% plot the score function from cross validation
%
% INPUT
% lambdas - values of lambda tested
% scores - a matrix of score values for each trial and each value of lambda

% calculte min lambda and lambda with the +1SE rule
[~, minlambdaind] = min(mean(scores, 1));
minlambda = lambdas(minlambdaind);

minlambdaindse = find( mean(scores, 1) < (mean(scores(:,minlambdaind)) + std(scores(:,minlambdaind))/sqrt(size(scores, 1))), 1 );
minlambdase = lambdas(minlambdaindse);

s = mean(scores, 1);
min_s = min(s);
max_s = max(s);

figure
loglog(lambdas, s, 'k-', 'linewidth', 2);
hold on
loglog(lambdas, s, 'k.', 'markersize', 20);
plot([minlambda minlambda], [min_s max_s], 'b-');
plot([minlambdase minlambdase], [min_s max_s], 'r-');
xlabel('\lambda');
ylabel('Score');

end