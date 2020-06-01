% Author: Fernando Ferreira

% email: fferreira@lps.ufrj.br
% Oct 2011
% Matlab 2011a required

classdef TimeSeries < hgsetget
	properties
		ts_original              % Original set of series to analyze
		ts                       % Timeseries
		model = cell({})         % cell containing features from the series
		nnParams = cell({})      % Parameters for NN
		test_serie
		test_serie_residue
		test_serie_original
		input_test
		estimated_output
		output
		output_gaps
	end
	methods


		%% Core Methods
		function obj = TimeSeries(serie)
			obj.ts_original = serie;
			obj.ts = serie;
			nSeries =  size(serie,1);
			%% Temporary
			obj.nnParams.corr_lag  = 10;
			obj.nnParams.corr_nstd =  3; % 99%
		end
		%% Set and get methods
		function obj = set.test_serie(obj, serie)
			obj.test_serie=serie;
			obj.test_serie_original= serie;
		end
		function obj = set.test_serie_original(obj, serie)
			obj.test_serie_original= serie;
		end

		%% Pre-processing block
		function removeHeteroscedastic(obj,opt)
			if nargin == 1
				[nSeries,nEvents] =size(obj.ts);
				obj.model.has_heteroscedastic = false(1,nSeries);
				for k= 1:nSeries
					serie = obj.ts(k,:);
					plotSerie(serie);
					fprintf('\nRemove heteroscedastic using logarithmic function?');
					if yesno
						obj.model.has_heteroscedastic(k) = true;
					end
				end
				fprintf('\n');
				H = repmat(...
					(obj.model.has_heteroscedastic)'...
					, 1, nEvents) ./ ...
					repmat(exp((1:nEvents)/nEvents), nSeries, 1);
				H(H==0)=1;
				obj.ts = obj.ts .* H;
			elseif nargin == 2
				if ~strcmpi(opt, 'test')
					err = MException('TimeSeries:removeHeteroscedastic', ...
					'Usage: removeHeteroscedastic() or removeHeteroscedastic(''test'')');
					throw(err);
				end
				if size(obj.test_serie,2) == 0
					err = MException('TimeSeries:removeHeteroscedastic', ...
					'No serie for testing. Use: TimeSerie::test_serie(serie) ');
					throw(err);
				end
				serie = obj.test_serie;
				if isfield(obj.model, 'has_heteroscedastic')
					[nSeries,nEvents] =size(serie);
					if length(obj.model.has_heteroscedastic) ~= size(serie, 1)
						err = MException('TimeSeries:removeHeteroscedastic', ...
							'Model and given serie dimensions are not compatible.');
						throw(err);
					end
					H = repmat(...
						(obj.model.has_heteroscedastic)'...
						, 1, nEvents) ./ ...
						repmat(exp((1:nEvents)/nEvents), nSeries, 1);
					H(H==0)=1;
					obj.test_serie = serie .* H;
				else
					err = MException('TimeSeries:removeHeteroscedastic', ...
					'No model was created. You should probably run TimeSeries::preprocess().');
					throw(err);
				end
			end
			close all
		end

		function removeStochasticTrend(obj, opt)
			if nargin == 1
				[fs, n] = removeNaN(obj.ts);
				ndiff = ordint(fs);
				fprintf('\n Suggested number of root unit: ');
				disp(ndiff);
				for i=1:length(ndiff)
					ndiff(i) = input('Number of unitary root should be used? ');
				end
			elseif nargin == 2
				if ~strcmpi(opt, 'test')
					err = MException('TimeSeries:removeStochasticTrend', ...
					'Usage: removeStochasticTrend() or removeStochasticTrend(''test'')');
					throw(err);
				end
				if size(obj.test_serie,2) == 0
					err = MException('TimeSeries:removeStochasticTrend', ...
					'No serie for testing. Use: TimeSerie::test_serie(serie) ');
					throw(err);
				end
				if ~isfield(obj.model, 'n_unit_root')
					err = MException('TimeSeries:removeStochasticTrend', ...
					'No model was created. You should probably run TimeSeries::preprocess().');
					throw(err);
				end
				if length(obj.model.n_unit_root) ~= size(obj.test_serie, 1)
					err = MException('TimeSeries:removeStochasticTrend', ...
						'Model and given serie dimensions are not compatible.');
					throw(err);
				end
				[fs, n] = removeNaN(obj.test_serie);
				ndiff = obj.model.n_unit_root;
			end
			obj.model.s_trend.initial_condition= cell(length(ndiff),1);
			for k=1:length(ndiff)
				[fs_diff, X0] = diff2(fs(k,:), ndiff(k));
				nan_array = NaN(1, ndiff(k));
				fs(k,:) = [nan_array fs_diff];
				obj.model.s_trend.initial_condition{k} = X0;
			end
			fs = addNaN(fs, n);
			if nargin == 1
				obj.ts = fs;
				obj.model.n_unit_root = ndiff;
			elseif nargin == 2
				obj.test_serie = fs;
			end
		end

		function removeSeasonality(obj, opt)
			nSTD = 3;
			if nargin == 1
				[fs, n] = removeNaN(obj.ts);
			%   Visual Test
				[nSeries, nEvents] = size(fs);
				residue = NaN(nSeries, nEvents);
				obj.model.seasonality = cell(nSeries,1);
				for k=1:nSeries
					[acf, ~, bounds] = crosscorr(fs(k,:), fs(k,:), nEvents-2, nSTD);
					acf = acf(floor(size(acf,2)/2)+2:end);
					p = figure(); hold on; grid on;
					stem(1:size(acf,2), acf, 'b.');
					plot([1 nEvents-2] , [bounds(1) bounds(1)], 'k-');
					plot([1 nEvents-2] , [bounds(2) bounds(2)], 'k-');
					title(sprintf('Visual test for serie %d', k));
					hold off;

					period = input(sprintf('Seasonality period for serie %d: ', k));
					close(p);
					if period ~=0
						y = zeros(nEvents-period,1);
						for index=1:(nEvents - period)
							y(index) = fs(k,index+period) - fs(k,index);
						end
						obj.model.seasonality{k}.x0 = fs(k,1:period);
						residue(k,period+1:end) = y';
					else
						residue(k,1:end)=fs(k,:);
						obj.model.seasonality{k}.x0 = [];
					end
					obj.model.seasonality{k}.period = period;
				end
				obj.ts = addNaN(residue,n);
			elseif nargin == 2
				if ~strcmpi(opt, 'test')
					err = MException('TimeSeries:removeSeasonality', ...
					'Usage: removeSeasonality() or removeSeasonality(''test'')');
					throw(err);
				end
				if size(obj.test_serie,2) == 0
					err = MException('TimeSeries:removeSeasonality', ...
					'No serie for testing. Use: TimeSerie::test_serie(serie) ');
					throw(err);
				end
				if ~isfield(obj.model, 'seasonality')
					err = MException('TimeSeries:removeSeasonality', ...
					'No model was created. You should probably run TimeSeries::preprocess().');
					throw(err);
				end
				[fs, n] = removeNaN(obj.test_serie);
				[nSeries, nEvents] = size(fs);
				residue = NaN(nSeries, nEvents);
				for k=1:nSeries
					period = obj.model.seasonality{k}.period;
					y = zeros(1,nEvents-period);
					for index=1:(nEvents - period)
						y(index) = fs(k,index+period) - fs(k,index);
					end
					residue(k,period+1:end) = y;
				end
				obj.test_serie = addNaN(residue,n);
			end
		end

		function removeCyclesAndTrend(obj,opt)
			if nargin == 1
				[fs, n] = removeNaN(obj.ts);
				[nSeries, nEvents] = size(fs);
				max_degree = 3;
				threshold = 0.9;

				obj.model.trend.degree = cell(nSeries,1);
				obj.model.trend.coeffs = cell(nSeries,1);
				obj.model.trend.yp     = cell(nSeries,1);

				obj.model.cycles.components = cell(nSeries,1);
				obj.model.cycles.ccos = cell(nSeries,1);
				obj.model.cycles.csin = cell(nSeries,1);
				obj.model.cycles.norm = cell(nSeries,1);
				obj.model.cycles.w    = cell(nSeries,1);

				for k=1:nSeries
					% Remove deterministic trend
					yp = cell(max_degree,1);
					coeffs = cell(max_degree,1);
					r2 = zeros(max_degree,1);
					[serie, nNaN] = removeNaN(fs(k,:));
					for d=1:length(r2)
						coeffs{d}= polyfit(1:length(serie), serie, d);
						yp{d} = polyval(coeffs{d}, 1:length(serie));
						r2(d) = 1 - norm(serie - yp{d},2)^2/...
							((length(serie)-1)*std(serie)^2);
					end
					degree = find(r2 > threshold ,1 );
					p =  figure();
					subplot(2,1,1);
					plot(serie, 'b-');
					if ~isempty(degree)
						hold on;
						plot(1:nEvents, yp{degree}, 'r--');
					else
						degree = 0;
					end
					grid on;
					subplot(2,1,2);
					bar(1:max_degree, r2);
					grid on;
					hold off;
					fprintf('\n The time serie can be described by a  degree %d polynomial?'...
						, degree);
					if ~yesno
						degree = inf;
						while (degree >=  max_degree)
							degree = input(sprintf('\nUse (degree < %d): ', max_degree));
						end
					end
					close(p)
					obj.model.trend.degree{k} = degree;
					if degree ~= 0
						obj.model.trend.coeffs{k} = coeffs{degree};
						obj.model.trend.yp{k}     = yp{degree};
						fs(k,:) = fs(k,:) -  yp{degree};
					else
						obj.model.trend.coeffs{k} = [];
						obj.model.trend.yp{k}     = [];
					end
					% Remove Cycles
					Y = serie';
					L = length(Y);
					fft_Y = fft(Y, length(Y));
					ccos = real(fft_Y'/(floor(L/2)));
					csin = imag(fft_Y'/(floor(L/2)));

					% Plot signal spectrum
					figure;
					stem(1:floor(L/2)+1, 2*abs(fft_Y(1:floor(L/2)+1)));
					xlabel('component');
					ylabel('|Y|')
					grid on;

					nComponents = input('\n Number of meaningful components: ');
					%status_fig = close(p);
					components = zeros(nComponents,1);
					%Extract Components
					for nc=1:nComponents
						[~ , index_max] = max(abs(fft_Y((1:floor(L/2)+1))));
						components(nc) = index_max;
						fft_Y(index_max) = 0;
						if index_max ~= 1
							fft_Y(L - index_max +2) = 0;
						else
							fprintf(...
							'\nWarning: The first component is not a cycle, is DC.');
							fprintf('\nAre you sure the pre-processing is correct?');
						end
					end

					% Reconstruct signal size
					serie = ifft(fft_Y)';
					% Store parameters
					obj.model.cycles.ccos{k} = ccos;
					obj.model.cycles.csin{k} = csin;
					obj.model.cycles.components{k} = components;
					obj.model.cycles.w{k} = 2*pi*(0:L-1)/L;
					fs(k,:) = addNaN(serie,nNaN);
				end
				obj.ts = addNaN(fs,n);
			elseif nargin == 2
				if ~strcmpi(opt, 'test')
					err = MException('TimeSeries:removeCycles', ...
					'Usage: removeCycles() or removeCycles(''test'')');
					throw(err);
				end
				if size(obj.test_serie,2) == 0
					err = MException('TimeSeries:removeCycles', ...
					'No serie for testing. Use: TimeSerie::test_serie(serie) ');
					throw(err);
				end
				if ~isfield(obj.model, 'cycles')
					err = MException('TimeSeries:removeCycles', ...
					'No model was created. You should probably run TimeSeries::preprocess().');
					throw(err);
				end
				[fs, n] = removeNaN(obj.test_serie);
				for k=1:size(fs, 1)
					if obj.model.trend.degree{k}
						% Remove Trend
						yp = polyval(obj.model.coeffs{d}, 1:size(fs));
						fs(k,:) = fs(k,:) - yp;
					end
					% Remove Cycles
					component =  obj.model.cycles.components{k};
					ccos = obj.model.cycles.ccos{k};
					csin = obj.model.cycles.csin{k};
					w    = obj.model.cycles.w{k};
					t    = 0:(length(fs(k,:)) - 1);
					for nc=1:length(component)
						c = component(nc);
						fs(k,:) = fs(k,:) - (...
						ccos(c)*cos(w(c)*t) +...
						csin(c)*sin(w(c)*t));
					end
				end
				obj.test_serie = addNaN(fs,n);
			end
			close all;
		end

		function preprocess(obj)
			%Gathers all functions concerned to preprocessing stage
			%  Heteroscedastic
			obj.removeHeteroscedastic()
			% Remove Trends
			obj.removeStochasticTrend()
			% Remove Seasons
			obj.removeSeasonality()
			% Remove cycles
			obj.removeCyclesAndTrend()
		end
		%%
		%% Pos-processing block
		function addCyclesAndTrend(obj)
			series = obj.output;
			[nSeries, nEvents] = size(series);
			for k=1:nSeries
				fs = series(k,:);
				% Add Cycles
				component =  obj.model.cycles.components{k};
				ccos = obj.model.cycles.ccos{k};
				csin = obj.model.cycles.csin{k};
				w    = obj.model.cycles.w{k};
				t    = 0:(length(fs) - 1);
				residue = zeros(1, nEvents);
				for nc=1:length(component)
					c = component(nc);
					residue =  residue + (...
						ccos(c)*cos(w(c)*t) +...
						csin(c)*sin(w(c)*t));
				end
				fs = fs +residue;
				% Add Trend
				if obj.model.trend.degree{k}
					yp = polyval(obj.model.trend.coeffs{k}, 1:size(fs));
					fs = fs + yp;
					figure;
					plot(yp);
				end
				obj.output(k,:) = fs;
			end
		end

		function addSeasonality(obj)
			series = obj.output;
			[nSeries, nEvents] = size(series);
			periods = zeros(nSeries,1);
			for k=1:nSeries
				periods(k) =  obj.model.seasonality{k}.period;
			end
			new_output = zeros(nSeries, nEvents + max(periods));
			for k=1:nSeries
				period = periods(k);
				if period ~= 0
					[serie, n] = removeNaN(series(k,:));
					% Ensure we have a whole number of periods in the current serie
					nPeriods = ceil(nEvents/period);
					x = [ serie NaN(1, nPeriods*period - nEvents) ];
					% Chop the serie into periods and add x0 to the beginning
					y = [ obj.model.seasonality{k}.x0; reshape(x, [period nPeriods])'];
					% Integrate on periods
					z = cumsum(y);
					% Get back to a single row
					fs = reshape(z', [1 numel(z)]);
					new_output(k,:)=   removeNaN(fs);
				else
					new_output(k,:) = fs;
				end
			end
			obj.output = new_output;
		end

		function addStochasticTrend(obj)
			series = obj.output;
			[nSeries, nEvents] = size(series);
			nRoot = obj.model.n_unit_root;
			series = [series NaN(nSeries, max(nRoot))];
			for k=1:nSeries
				series(k,:) = integrate(...
					series(k,:), ...
					obj.model.s_trend.initial_condition(k));
			end
			obj.output = series;
		end

		function addHeteroscedastic(obj)
			series = obj.output;
			[nSeries, nEvents] = size(series);
			gaps = obj.output_gaps;
			H = zeros(nSeries, nEvents);
			for index=1:nSeries
				H(index,:) = obj.model.has_heteroscedastic(index)...
					./ exp((gaps(index)+1:nEvents+gaps(index))...
					/(nEvents+gaps(index)));
			end

			H(H==0)=1;
			fs = series ./ H;
			obj.output = fs;
		end

		%% Estimator block
		function assembleData(obj)
			close all;
			%% Study xcorrelations
			fs = removeNaN(obj.ts);
			nSerie = size(fs,1);
			obj.model.estimator.used_lags = cell(nSerie);
			obj.model.estimator.input = cell(nSerie,1);
			obj.model.estimator.target = cell(nSerie,1);
			for index=1:nSerie % create one estimator for each serie
				for k=1:nSerie
					% Find correlation in the time window
					[xcf, ~, bounds] = crosscorr(fs(k,:)', fs(index,:)',...
						obj.nnParams.corr_lag, obj.nnParams.corr_nstd);
					xcf = xcf(floor(size(xcf,1)/2)+2:end);
					used_lags = find(abs(xcf) > bounds(1));
					obj.model.estimator.used_lags{index,k} = used_lags;
				end %end for k
			end %end for i
			no_corr = cellfun(@isempty, obj.model.estimator.used_lags);
			obj.model.estimator.use_random_walk = ...
				logical(sum(no_corr,2) == nSerie)';
			%% Assemble the estimator input and target dataset
			nNodes = sum(cellfun(@length, obj.model.estimator.used_lags),2);
			for index=1:nSerie
				if (~no_corr(index))
					cumindex = 0;
					nInputNodes = nNodes(index);
					used_lags = obj.model.estimator.used_lags(index, :);
					used_lags(cellfun(@isempty, used_lags))= {0};
					events_to_ignore  = 1 + max(cellfun(@(x) max(x),used_lags));
					input = zeros(nInputNodes, size(fs,2) - events_to_ignore);
					for k=1:nSerie
						n = length(used_lags{k});
						for event=events_to_ignore:size(fs,2)
							input(cumindex+1:cumindex+n, event+1-...
								events_to_ignore) = fs(k, event - used_lags{k});
						end
						cumindex = cumindex + n;
					end
				else
					input = fs(index, events_to_ignore-1 :end-1);
				end
				obj.model.estimator.input{index}  = input;
				obj.model.estimator.target{index} = fs(index, events_to_ignore:end);
			end
		end %assembleData

		function assembleDataForTest(obj)
			fs = removeNaN(obj.test_serie);
			nSerie = size(fs,1);
			obj.input_test = cell(nSerie,1);
			nNodes = sum(cellfun(@length, obj.model.estimator.used_lags),2);
			for index = 1:nSerie
				if ~obj.model.estimator.use_random_walk
					cumindex = 0;
					nInputNodes = nNodes(index);
					used_lags = obj.model.estimator.used_lags(index, :);
					used_lags(cellfun(@isempty, used_lags))= {0};
					events_to_ignore  = 1 + max(cellfun(@(x) max(x),used_lags));
					input = zeros(nInputNodes, size(fs,2) - events_to_ignore);
					for k=1:nSerie
						n = length(used_lags{k});
						for event=events_to_ignore:size(fs,2)
							input(cumindex+1:cumindex+n, event+1-...
								events_to_ignore) = fs(k, event - used_lags{k});
						end
						cumindex = cumindex + n;
					end
				else
					input = fs(index, events_to_ignore-1 :end-1);
				end
				obj.input_test{index} = input;
			end
		end

		function createEstimator(obj, nnobj)
			obj.model.estimator.net = cell(size(obj.ts,1),1);
			for index=1:size(obj.ts,1)
				if ~obj.model.estimator.use_random_walk
					input  = obj.model.estimator.input{index};
					target = obj.model.estimator.target{index};
					nnobj.train(input, target);
					obj.model.estimator.net{index} = nnobj.best_net;
				end
			end
		end%createEstimator

		%% Applying model

		function applyModel(obj, test_serie)
			obj.test_serie = test_serie;
			% Preprocess test dataset
			obj.removeHeteroscedastic('test');
			obj.removeStochasticTrend('test');
			obj.removeSeasonality('test');
			obj.removeCyclesAndTrend('test');
			obj.test_serie_residue = obj.test_serie;
			close all;
			% Format estimator input
			obj.assembleDataForTest();
			obj.estimated_output = cell(size(test_serie,1),1);
			for k=1:length(obj.model.estimator.use_random_walk)
				if obj.model.estimator.use_random_walk
					% The best estimator is the current event
					obj.estimated_output{k} = obj.input_test{k};
				else
					y = sim(obj.model.estimator.net{k}.net, obj.input_test{k});
					obj.estimated_output{k} = y;
				end
			end
			obj.model.estimator.output = obj.estimated_output;
			[obj.output, obj.output_gaps] = cropSeries(obj.estimated_output);
			obj.estimated_output = obj.output;
			obj.addCyclesAndTrend();
			obj.addSeasonality();
			obj.addStochasticTrend();
			obj.addHeteroscedastic();
		end %applyModel
	end%methods
end

function [Y, X0] = diff2(X, n)
	X0 = zeros(n,1);
	Y = X;
	for k=1:n
		X0(k) = Y(1);
		Y=diff(Y,1);
	end
end

function ndiff = ordint(serie)
	nSerie = size(serie,1);
	ndiff = zeros(1, nSerie);
	for i=1:nSerie
		s = serie(1,:);
		has_root = true;
		DWbound = 5e-2;
		while has_root,
			[ADF, ~, ~, ~] = unitroot (s);
			if(ADF(3,4) <= 1e-1 && (abs(ADF(1,2)-2) < DWbound))
				has_root = false;
			elseif ADF(3,1) <= 0.1
				has_root = false;
			else
				ndiff(i) = ndiff(i) + 1;
				s = diff(s);
			end
		end
	end
end

function [Y] = integrate(X, X0)
	if iscell(X)
		X = cell2mat(X);
	end
	if iscell(X0)
		X0 = cell2mat(X0);
	end
	[Y, n] = removeNaN(X);
	for k=length(X0):-1:1
		Y = [X0(k) cumsum(Y)+X0(k)];
	end
	Y = [Y NaN(1, n-length(X0))];
end

function [A, gaps] = cropSeries(series)
	%nSeries = length(series);
	smaller_serie_length = min(cellfun(@(x) size(x,2), series));
	A = cell2mat(cellfun(@(x) x(end-smaller_serie_length+1:end), series,...
		'UniformOutput', false));
	gaps = cell2mat(cellfun(@(x) length(x)-smaller_serie_length, ...
		series, 'UniformOutput', false));
end

function plotSerie(serie)
	figure;
	hold on;
	plot(serie);
	grid on;
end

function fs = addNaN(serie, n)
	if n ~= 0
		n_serie = size(serie,1);
		nan_matrix = NaN(n_serie, n);
		fs = [nan_matrix serie];
	else
		fs = serie;
	end
end

function [fs,n] = removeNaN(serie)
	fs = serie;
	indexes = find(sum(isnan(fs),1)==size(fs,1));
	fs(:,indexes) = [];
	n = length(indexes);
end

