classdef mapsi_options
%MAPSI_OPTIONS creates a data object used to set all of the options for
% the optimization.
%
% Usage: opts = mapsi_options([key1,value1,...])
%
% Valid options include (default options in parentheses):
% 'N' (42)
%       The number of vertices for the 2D mesh.
% 'folds' (20)
%       Number of folds to use in cross validation. Each fold consists of a
%       different random sample of the data. In general, more folds leads
%       to a more representative solution but takes more computational
%       effort.
% 'regularizer' ('ridge')
%       Type of regularization to use during optimization. Valid options
%       include penalizing the first derivative ('d1') or entropy
%       regularization ('entropy').
% 'lambdas' (logspace(-3, 4, 30))
%       A list of the regularization hyperparameters to test during cross
%       validation. The range logspace(-3, 4, 30) generally works well for 
%       regularizers. Note that this range is automatically divided by numbasis!
% 'cvsplit' (0.5)
%       The size of the training data for each sample of data drawn during
%       the Monte Carlo cross validation procedure as a fraction of the
%       total size.
% 'mctrials' (1e5)
%       The number of Monte Carlo trials used to estimate error from the
%       resulting posterior distribution.
% 'optimality' (1e-8)
%       The optimality tolerance used during optimization.
% 'steptol' (1e-12)
%       The step tolerance used during optimization.
% 'n_sub' (4)
%       The number of subdivisions of triangles for integration

properties
    N
    folds
    regularizer
    lambdas
    cvsplit
    mctrials
    optimality
    steptol
    n_sub
end

methods
    
    function opts = mapsi_options(varargin)
        % ensure arguments come in pairs
        if mod(length(varargin), 2) ~= 0
            error('Arguments must be supplied in pairs\n%s','Usage: opts = mapsi_options([key1,value1,...])')
        end
        
        % default values
        opts.N = 42;
        opts.folds = 20;
        opts.regularizer = 'ridge';
        opts.lambdas = logspace(-3, 4, 30);
        opts.cvsplit = 0.5;
        opts.mctrials = 1e5;
        opts.optimality = 1e-8;
        opts.steptol = 1e-12;
        opts.n_sub = 4;
        
        for i = 1:2:length(varargin)
            switch varargin{i}
                case 'N'
                    opts.N = varargin{i+1};
                case 'folds'
                    opts.folds = varargin{i+1};
                case 'regularizer'
                    opts.regularizer = varargin{i+1};
                case 'lambdas'
                    opts.lambdas = varargin{i+1};
                case 'cvsplit'
                    opts.cvsplit = varargin{i+1};
                case 'mctrials'
                    opts.mctrials = varargin{i+1};
                case 'optimality'
                    opts.optimality = varargin{i+1};
                case 'steptol'
                    opts.steptol = varargin{i+1};
                case 'n_sub'
                    opts.n_sub = varargin{i+1};
                otherwise
                    error('The key ''%s'' is not a valid option', string(varargin{i}))
            end
        end
        
    end
    
    % ------------------------- properties ------------------------------
    
    function opts = set.N(opts, N)
        if numel(N) ~= 1 || ( N < 3 )
            error('N must consist of 1 integers greater than 3')
        end
        opts.N = floor(N);
    end
    
    function opts = set.folds(opts, folds)
        if folds <= 0
            error('''folds'' must be a positive integer')
        end
        opts.folds = floor(folds);
    end
    
    function opts = set.regularizer(opts, reg)
        regoptions = {'ridge'};
        
        if ~ismember(reg, regoptions)
            error('''%s'' is a not a valid regularization option',string(reg))
        end
        
        opts.regularizer = reg;
    end
    
    function opts = set.lambdas(opts, ls)
        if isempty(ls)
            error('''lambdas'' vector must be nonempty')
        end
        % ensure that lambdas are sorted from largest to smallest
        opts.lambdas = sort(ls, 'descend');
    end
    
    function opts = set.cvsplit(opts, split)
        if split < 0 || split > 1
            error('Invalid ''cvsplit'' value')
        end
        
        opts.cvsplit = split;
    end
    
    function opts = set.mctrials(opts, mctrials)
        if mctrials <= 0
            error('''mctrials'' must be a positive integer')
        end
        
        opts.mctrials = floor(mctrials);
    end
    
    function opts = set.optimality(opts, optim)
        opts.optimality = optim;
    end
    
    function opts = set.steptol(opts, steptol)
        opts.steptol = steptol;
    end
    
end

end