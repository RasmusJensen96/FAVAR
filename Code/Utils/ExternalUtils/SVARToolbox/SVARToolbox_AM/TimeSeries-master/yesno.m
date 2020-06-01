function isyes = yesno(m) 
% YESNO - returns true for y-key, false for n-key
%   ISYES = YESNO waits for a keypress of either the Y-key or N-key and
%   returns a logical one (true) if the Y-key was pressed, and a logical
%   zero (false) if the N-key was pressed. 
%
%   Example:
%   fprintf('\nHere is a number for you: %5.0f', 1000*rand)
%   while(1)
%       fprintf('\nDo you want another number (y/n)?  ') ;
%       if yesno, fprintf(' %5.0f', 1000*rand) ;
%       else      fprintf('\nBye!\n') ;
%                 break
%       end
%   end
%
%  See also INPUT, GETKEY (on the FEX)

% for Matlab 6.5 and upwards
% version 1.0 (apr 2008)
% author : Jos van der Geest
% email  : jos@jasen.nl
%
% History
% 1.0 (apr 2008) creation, simplified version of GETKEY

% Set up the figure. May be the position property  should be individually
% tweaked to avoid visibility
fh = figure('keypressfcn','uiresume', ...
    'windowstyle','modal',...    
    'position',[-10 -10 1 1],... 
    'Name','YESNO' ) ; 
try
    while(1)        
        uiwait(fh) ; % Wait for a keypress
        ch = lower(get(fh,'Currentcharacter')) ;
        if isequal(ch,'y') || isequal(ch,'n'),
            break
        end
    end
    isyes = ch == 'y' ;
catch
    % Something went wrong, return false and a warning.
    warning('Call to YESNO failed.') ;
    isyes = false ;
end

delete(fh) ;

