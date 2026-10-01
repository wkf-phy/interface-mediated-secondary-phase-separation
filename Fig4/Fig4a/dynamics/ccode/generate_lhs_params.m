function [c1_rand, c2_rand] = generate_lhs_params(c1_range, c2_range, N)
%GENERATE_LHS_PARAMS Sample paired concentrations by Latin hypercube sampling.
%   c1_range, c2_range: [minimum, maximum] for each concentration.
%   N: number of samples.
%   c1_rand, c2_rand: paired N-by-1 sample vectors.

    c1_min = c1_range(1); c1_max = c1_range(2);
    c2_min = c2_range(1); c2_max = c2_range(2);

    % Draw one point from each interval in both dimensions.
    x1_base = ((0:N-1)' + rand(N, 1)) / N; 
    x2_base = ((0:N-1)' + rand(N, 1)) / N;

    % Randomly pair the samples.
    x1_shuffled = x1_base(randperm(N));
    x2_shuffled = x2_base(randperm(N));

     % Map samples to the requested concentration ranges.
    c1_rand = c1_min + (c1_max - c1_min) * x1_shuffled;
    c2_rand = c2_min + (c2_max - c2_min) * x2_shuffled;
end