function [ DataCell ] = readin_data(nodeFileName, linkFileName)
%Read in the rate file from nodeAdd_rate and linkAdd_rate generated
%     by two generators.
%     Input Format: 
% 	         nodeAdd_rate: t1 n(t1) t2 n(t2) ...
% 	         linkAdd_rate: t1 e(t1) t2 e(t2) ...
%     Output Format:
%            DataCell  #instances * 4
%            DataCell{i,1} time of node of i^th instance
%            DataCell{i,2} rate of node of i^th instance
%            DataCell{i,3} time of link of i^th instance
%            DataCell{i,4} rate of link of i^th instance
         
% fn = fopen('nodeAdd_rate', 'r');
% fe = fopen('linkAdd_rate', 'r');
fn = fopen(nodeFileName, 'r');
fe = fopen(linkFileName, 'r');


if fn < 0
    fprintf('Open %s error!\n', nodeFileName);
    return 
end

if fe < 0
    fprintf('Open %s error!\n', linkFileName);
    return
end


DataCell = [];
while ~feof(fn)
    ln = fgetl(fn);
    le = fgetl(fe);

    ntokes = regexp(ln, '\s+', 'split');
    etokes = regexp(le, '\s+', 'split');
    
    ns = sprintf('%s ', ntokes{:});
    ndata = sscanf(ns, '%f');
    
    nT = ndata(1:2:end);
    nRate = ndata(2:2:end);
        
    es = sprintf('%s ', etokes{:});
    edata = sscanf(es, '%f');
    eT = edata(1:2:end);
    eRate = edata(2:2:end);  
    DataCell = [DataCell; {nT,nRate,eT,eRate}];    
end

end

