function [X,W,squic_info]=SQUIC(Y,lambda,X_pattern,max_iter_in,drop_tol_in,term_tol_in,verbose_in)

% This is the ignore text for optional inputs
IGNORE='NULL';

% Set Working directly as SQUIC_TEMP_{UNIX_TIMESTAMP}
CWD=strcat(pwd,'/SQUIC_TEMP',num2str(posixtime(datetime('now')) * 1e6 +randi(100)));
mkdir(CWD);

% SQUIC_CMD : p n Y_loc lambda LambdaMatrix_loc max_iter drop_tol term_tol X0_loc W0_loc index_offset verbose X_loc W_loc " << std::endl;
% [integer:p] Number of random variables
% [integer:n] Number of samples.
[n,p]=size(Y); % data must be saved as nxp matrix

save(strcat(CWD,'/Y.dat'),'Y','-ascii');

% [string:Y_loc] Input loc/of/file/Y.dat
Y_loc=strcat(CWD,'/Y.dat');

% [double:lambda] Scalar lambda value.
% [string:LambdaMatrix_loc] loc/of/file/LambdaMatrix.dat
lambda=lambda;

if(numel(X_pattern)==0)
    LambdaMatrix_loc=IGNORE;
else
    LambdaMatrix_loc=strcat(CWD,'/lambda_matrix.dat');
    [row,col,v] = find(X_pattern);
    dlmwrite(LambdaMatrix_loc,[col,row,v], 'delimiter', '\t');
end

% [integer:max_iter] Max iteration.
max_iter=max_iter_in;

% [double:drop_tol] Drop out tolderance.
drop_tol=drop_tol_in;

% [double:term_tol] Terminal Tolerence
term_tol=term_tol_in;

% [string:X0_loc] Input loc/of/file/X0.dat
% [string:W0_loc] Input loc/of/file/W0.dat
X0_loc=IGNORE;
W0_loc=IGNORE;

% [integer:index_offset] Offset of indexing
index_offset=1;

% [integer:verbose] Verbosity;
verbose=verbose_in;

% [string:X_loc] Output loc/of/file/X.dat
% [string:W_loc] Output loc/of/file/W.dat
X_loc=strcat(CWD,'/X.dat');
W_loc=strcat(CWD,'/W.dat');

% [string:log_loc] Output loc/of/file/info.log
log_loc=strcat(CWD,'/info.log');

% Set call string
cmd=[
'./SQUIC_CMD ', ...
num2str(p),' ', ...
num2str(n),' ', ...
Y_loc,' ', ...
num2str(lambda),' ', ...
LambdaMatrix_loc,' ', ...
num2str(max_iter),' ',...
num2str(drop_tol),' ',...
num2str(term_tol),' ',...
X0_loc,' ',...
W0_loc,' ', ...
num2str(index_offset),' ', ...
num2str(verbose),' ',...
X_loc,' ', ...
W_loc,' ',...
log_loc, ...
' >', CWD,'/ouput.txt' % here i just dump all the output out...
];

%Execute call
system(cmd);

% Load save data & Convert to sparse matrix
load(strcat(CWD,'/X.dat'),'Y','-ascii');
load(strcat(CWD,'/W.dat'),'W','-ascii');
X=spconvert(X);
W=spconvert(W);

squic_info=jsondecode(fileread(log_loc));

%Remove the temp folder
rmdir(CWD, 's'); 

end