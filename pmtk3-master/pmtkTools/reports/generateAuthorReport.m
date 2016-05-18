function generateAuthorReport()
%% Generate the html contributing author report
% PMTKneedsMatlab 

% This file is from pmtk3.googlecode.com

dest            = fullfile(pmtk3Root(), 'docs', 'authors');
outputFile      = fullfile(dest, 'fileAuthors.html');
excludedAuthors = tokenize(getConfigValue('PMTKauthors'), ',');
colNames        = {'FILE NAME', 'AUTHOR(S)', 'SOURCE URL', 'DATE'};
pmtkRed         = getConfigValue('PMTKred');
R               = pmtkTagReport(); 
authored        = cellfun(@(a)~isempty(setdiff(a, excludedAuthors)),  R.authors);
fnames          = fnameOnly(R.files(authored)); 
authors         = cellfuncell(@(a)catString(a, ', '),   R.authors(authored)); 
links           = cellfuncell(@googleCodeLink, fnames); 
dates           = gatherTagText(R, 'PMTKdate', find(authored)); 
sourceUrls      = gatherTagText(R, 'PMTKurl', find(authored)); 
sourceUrls      = cellfunNonEmpty(@(u)sprintf('<a href="%s"> website </a>', u), sourceUrls); 
data            = [links(:), authors(:), sourceUrls(:), dates(:)]; 
perm            = sortidx(lower(fnames)); 
data            = data(perm, :); 
header = formatHtmlText({
    '<font align="left" style="color:%s"><h2>Contributed Files</h2></font>'
    ''
    'Revision Date: %s'
    ''
    'Auto-generated by %s.m'
    ''
    ''}, pmtkRed, date(), mfilename()); 
htmlTable('data'           , data                                  , ...
    'colNames'             , colNames                              , ...
    'doSave'               , true                                  , ...
    'filename'             , outputFile                            , ...
    'colNameColors'        , {pmtkRed, pmtkRed, pmtkRed, pmtkRed}  , ...
    'header'               , header                                , ...
    'dataAlign'            , 'left'                                , ...
    'caption'              , '<br> <br>'                           , ...
    'captionLoc'           , 'bottom'                              , ...
    'doshow'               , false);
end
