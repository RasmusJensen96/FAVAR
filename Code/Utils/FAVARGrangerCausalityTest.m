function [TablePval]  = FAVARGrangerCausalityTest(FY, numLags, vnames)
%{ 
 Performs Granger Causality test on the system:
--------------------------------Inputs-------------------------------------
        FY        = Endogenous observations; [F, Y]
        numLags   = Number of lags to consider
        vnames    = Names of endogenous [F1, F2,...] variables
--------------------------------Outputs------------------------------------
        TablePval = table containing p-values, Predictor: rows, predicted:
                    cols.  diagonal defined as 0.

Rasmus M. Jensen 2020
%}


[T, N] = size(FY);

idxpre = 1:numLags;
idxpost = (numLags+1):T;

FY0 = FY(idxpre,:);
FY  = FY(idxpost,:);


for sr = 1:N
    Cause = [FY0((end-numLags+1):end,sr); FY(:,sr)];
    nameE   = vnames{sr};
    for i = 1:N
        if i == sr
            continue
        end
        Effect = [FY0((end-numLags+1):end,i); FY(:,i)];
        nameC   = vnames{i};
        [~,pvalue] = gctest(Cause,Effect,"NumLags",numLags);
        pval(i,sr) = pvalue;
    end
end

TablePval = array2table(pval);
[TablePval.Properties.VariableNames, TablePval.Properties.RowNames]= deal(string(vnames));







end