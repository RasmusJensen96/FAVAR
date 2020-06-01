classdef swirlTraining < hgsetget
	properties
		params;
		G ;
		tp;
		tp_initial;
		ants_movements;
		best_net;
		HFC;
	end
	methods
		%% Set and get methods
		function obj = set.G(obj, G)
			obj.G = G;
		end
		function obj = set.tp(obj, tp)
			obj.tp = tp;
		end
		function obj = set.tp_initial(obj, tp_initial)
			obj.tp_initial = tp_initial;
		end
		function obj = set.ants_movements(obj, ants_movements)
			obj.ants_movements = ants_movements;
		end
		function obj =set.best_net(obj,best_net)
			obj.best_net = best_net;
		end
		function [best_net] = get.best_net(obj)
			[~, idx] = max(obj.G);
			obj.best_net = obj.tp{idx}.nets{obj.tp{idx}.bestIdx};
			best_net = obj.best_net;
		end
		%% Core methods
		function obj =  swirlTraining()
			addpath('swirl/')
			obj.params.max_neurons  = 10;
			obj.params.max_iter     = 20;
			obj.params.ants_factor  = 1;
			% Set parameters for the ants
			obj.params.antParams = [];
			%Set parameters nnet
			obj.params.netParams.trainClass = @gdNet;
			obj.params.netParams.number_of_nets = 20;
			obj.HFC = 1;
			obj.params.antParams.beta = 0.5;
			obj.params.antParams.alpha = 2.;
			obj.params.antParams.rho = .5;
			obj.params.factor_train= 0.7;
		end

		function train(obj, input, target)
			% normalize
			data = input;
			data = data - repmat(mean(data,2), 1, size(data, 2));
			data = data ./ repmat(std(data, 0, 2), 1, size(data,2));
			trainInd   =  1:2:length(data);
			valInd     =  2:2:length(data);
			obj.params.netParams.data       =  data;
			obj.params.netParams.target     =  target;
			indexes = randperm(length(data));
			num_train = round(length(indexes)*obj.params.factor_train);
			obj.params.netParams.trainInd  =  indexes(1:num_train);
			obj.params.netParams.valInd  =  indexes((num_train+1):end);
			addpath('swirl/')
			[ants_movements, tp, tp_initial, G] = ...
				SWIRL(obj.HFC, obj.params);
			obj.ants_movements = ants_movements;
			obj.tp = tp;
			obj.tp_initial = tp_initial;
			obj.G = G;
		end

	end
end%classdef
