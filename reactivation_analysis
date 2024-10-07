function reactivation_analysis
    
% Main script

ISC=1; % 1 for between-participant analysis, 0 for within-participant analysis

% remove diagonals by analysis type
if ISC==1
    numDiag=1;
    suffix='ISC';
else
    numDiag=10;
    suffix='intra';
end

% Define the base directory
base_dir = '/path/to/your/base/directory/';

% Define the subdirectory and filename
sub_dir = 'PUNC';
filename = '';

% Construct the full path to the file
full_file_path = fullfile(base_dir, sub_dir, filename);

% Display the full file path for verification
disp(['full file path: ', full_file_path]);

% Check if the file exists
if exist(full_file_path, 'file')
    disp('File exists.');
    load(full_file_path);  % Load the file
else
    error(['File does not exist: ', full_file_path]);
end

% List the contents of the directory for verification
disp('Listing contents of the directory:');
dir(fullfile(base_dir, sub_dir));

for iSub=1:size(EB,3)
    if ISC~=1 % within-participant analysis
        boundSceneE(:,:,iSub)=corr(squeeze(EB(:,:,iSub)),squeeze(scene(:,:,iSub)),'rows','pairwise');
        boundSceneP(:,:,iSub)=corr(squeeze(CTRp(:,:,iSub)),squeeze(scene(:,:,iSub)),'rows','pairwise');
        boundSceneF(:,:,iSub)=corr(squeeze(CTRf(:,:,iSub)),squeeze(scene(:,:,iSub)),'rows','pairwise');
    else % between-participant analysis
        ss=squeeze(scene(:,:,iSub)); % single-subject scene representations
        group=EB; % group boundary representations
        group(:,:,iSub)=[]; % exclude the subject whose scene representations will be used
        boundSceneE(:,:,iSub)=corr(squeeze(nanmean(group,3)),ss,'rows','pairwise');
        
        group=CTRp;
        group(:,:,iSub)=[];
        boundSceneP(:,:,iSub)=corr(squeeze(nanmean(group,3)),ss,'rows','pairwise');
        
        group=CTRf;
        group(:,:,iSub)=[];
        boundSceneF(:,:,iSub)=corr(squeeze(nanmean(group,3)),ss,'rows','pairwise');
    end
    
    % To make the matrix symmetrical, we must remove the rows (and not only the columns) of invalid representations
    idx=find(isnan(nanmean(boundSceneE(:,:,iSub))));
    boundSceneE(idx,:,iSub)=nan;
    boundSceneP(idx,:,iSub)=nan;
    boundSceneF(idx,:,iSub)=nan;
    
    numVals(iSub)=length(find(~isnan(boundSceneE(:,:,iSub))));
end

% remove diagonals
m=size(boundSceneE,1);
if numDiag>1
    X = full(spdiags(bsxfun(@times,ones(m,1),nan(1,numDiag)),[(numDiag/2-1)*-1:numDiag/2],m,m)); %create matrix of zeros with nans at the diagonals
else
    X = eye(m);
    X(X==1)=nan;
end

boundSceneE=boundSceneE+X;
boundSceneP=boundSceneP+X;
boundSceneF=boundSceneF+X;

% Fisher correction
boundSceneE=0.5*(log(1+boundSceneE)-log(1-boundSceneE));
boundSceneP=0.5*(log(1+boundSceneP)-log(1-boundSceneP));
boundSceneF=0.5*(log(1+boundSceneF)-log(1-boundSceneF));

% calculate the reactivation index
for iSub=1:size(EB,3)
    indexE(iSub)=nanmean(getTriangular(boundSceneE(:,:,iSub),1))-nanmean(getTriangular(boundSceneE(:,:,iSub),0));
    indexP(iSub)=nanmean(getTriangular(boundSceneP(:,:,iSub),1))-nanmean(getTriangular(boundSceneP(:,:,iSub),0));
    indexF(iSub)=nanmean(getTriangular(boundSceneF(:,:,iSub),1))-nanmean(getTriangular(boundSceneF(:,:,iSub),0));
end

% permutation tests

% uncontrolled effect
Total_vals=sum(numVals);
Stat=mean(indexE);

for itr=1:10000
    flip_idx=randi([0, 1], 1,size(indexE,2)); %sign-flip index
    flip_idx(flip_idx==0)=-1;
    rand_stat(itr)=mean(indexE.*flip_idx);
end

pval=1-sum(Stat>=rand_stat)/length(rand_stat);

% controlled effect
indexE_ctr=(indexE-(0.5*(indexP+indexF)));
Stat_ctr=mean(indexE_ctr);

for itr=1:10000
    flip_idx=randi([0, 1], 1,size(indexE,2)); %sign-flip index
    flip_idx(flip_idx==0)=-1;
    rand_stat_ctr(itr)=mean(indexE_ctr.*flip_idx);
end

pval_ctr=1-sum(Stat_ctr>=rand_stat_ctr)/length(rand_stat_ctr);

% Save the results
save_filename = ['results_' filename];
save(fullfile(base_dir, save_filename),'Stat','rand_stat','Stat_ctr','rand_stat_ctr','boundSceneE', 'boundSceneP', 'boundSceneF', 'indexE', 'indexP', 'indexF', 'pval', 'pval_ctr','indexE_ctr');
end
