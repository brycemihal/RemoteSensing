function[] = generateReport(rptName)
%% Function to write figures to word document
% Bryce Mihalevich
% 8/21/18
% code addapted from save2word.m
%% cd to report directory
% rptPath = fullfile(pwd,'Results');
% if ~exist(rptPath,'dir')
%     mkdir(rptPath)
% end
% cd(rptPath);

%% open word and creat/open file
rptName = fullfile(pwd,rptName);
wordApp = actxserver('Word.Application');

if ~exist(rptName,'file')
  op = invoke(wordApp.Documents,'Add'); % Create new doc
else
  op = invoke(wordApp.Documents,'Open',rptName); % Open existing doc
end

%% Get list of open figures
figs = findobj('Type','Figure');

%% Loop through each open figure and add it to the word doc
for i = 1:length(figs)
    figure(i);
    % Find end of document and make it the insertion point
    end_of_doc = get(wordApp.ActiveDocument.Content,'end');
    set(wordApp.Application.Selection,'Start',end_of_doc);
    set(wordApp.Application.Selection,'End',end_of_doc);

    print('-clipboard','-dmeta') % copy current figure
    invoke(wordApp.Selection,'Paste'); %paste current figure
end

%% Save file
if ~exist(rptName,'file')
  % Save file as new:
  invoke(op,'SaveAs',rptName,1);
else
  invoke(op,'Save'); % Save existing file
end

invoke(op,'Close');
invoke(wordApp,'Quit');
delete(wordApp)