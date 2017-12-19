
function contours =read_contours(path2cont)
    contours_list = dir([path2cont,'/*.txt']); 
    numFiles = size(contours_list,1);    

    %contours=struct; 
    for k=1:1:numFiles
        cname_i=contours_list(k).name;
        
        % path to text file 
        path2txt=[path2cont,'/',cname_i]
        
        % read endo contours files
        contours{k}=load(path2txt);
    end

