classdef TimeSeriesDataQuality < hgsetget
	properties
		estimated_serie
		serie
		% Quality indicators
		mape
		nmse1
		nmse2
		r
		lag0
		lag1
		validation_corridor
	end

	methods
		%% Initializer
		function obj= TimeSeriesDataQuality(estimated_serie, serie)
			obj.estimated_serie = estimated_serie;
			obj.serie = serie;
		end
		%% Set and Get methods

		%% Core Methods
		function mape = calculateMAPE(obj)
			mape =1;
			obj.mape =mape;
		end

		function [nmse1, nmse2] = calculateNMSE(obj)
			nmse1 =1;
			nmse2 =1;
			obj.nmse1 = nmse1;
			obj.nmse2 = nmse2;
		end
		function calculateValidationCorridor(obj)
			[nSeries, nEvents] = size(obj.estimated_serie);
			vCorridor = struct();
			vCorridor.upper = cell(nSeries, 1);
			vCorridor.lower = cell(nSeries, 1);
			for s=1:nSeries
				upper = zeros(3,nEvents);
				lower = zeros(3,nEvents);
				delta =  - length(obj.estimated_serie(s,:)) + length(obj.serie(s,:))+1;
				err = obj.estimated_serie(s,:) - obj.serie(s,delta:end);

				std_dev_err = nanstd(err);
				for k=1:3
					upper(k,:) = obj.estimated_serie(s,:) + k * std_dev_err;
					lower(k) = obj.estimated_serie(s,:) - k * std_dev_err;
				end
				vCorridor.upper{s} = upper;
				vCorridor.lower{s} = lower;
			end
			obj.validation_corridor = vCorridor;
		end

		function plotSerie(obj, k)
			serie = obj.serie(k,:)';
			estimated_serie = obj.estimated_serie(k,:)';
			vCorridor_lower = obj.validation_corridor.lower{k};
			vCorridor_upper = obj.validation_corridor.upper{k};
			figure;
			hold on;
			jbfill(1:length(vCorridor_upper(3,:)), vCorridor_upper(3,:), vCorridor_lower(3,:), [0.3 0.3 0.3]);
			hold on;
			jbfill(1:length(vCorridor_upper(2,:)), vCorridor_upper(2,:), vCorridor_lower(2,:), [0.3 0.3 0.3]);
			hold on;
			jbfill(1:length(vCorridor_upper(1,:)), vCorridor_upper(1,:), vCorridor_lower(1,:), [0.1 0.1 0.1]);
			hold on;
			
			hold on;
			plot(serie, 'y-');
            plot(estimated_serie, 'w-');
		end

		function r = calculateR(obj)
			nmse1 =1;
			nmse2 =1;
			obj.nmse1 = nmse1;
			obj.nmse2 = nmse2;
		end
		function r = calculateLag0(obj)
			lag0 =1;
			obj.lag0 = lag0;
		end
		function lag1 = calculateLag1(obj)
			lag1 =1;
			obj.lag1 = lag1;
		end
	end
end
