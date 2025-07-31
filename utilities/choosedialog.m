%Author: LWeissinger, 18.07.2025
function choice = choosedialog
answer = questdlg('We advise to pick atleast three points!', ...
	'Warning', ...
	'OK, pick more','Ignore');
% Handle response
switch answer
    case 'OK, pick more'
        choice = 1;
    case 'Ignore'
        choice = 0;
end