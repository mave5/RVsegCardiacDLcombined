%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RV Segmentation Challenge in Cardiac MRI
% Performance assessment of automatic algorithms vs manual ground truth
% main_perf_assessment.m: main program
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program assumes that:
% - the TrainingSetPrefix directory contains all training data as unzipped.
% - the TeamPrefix directory contains automatic contours, ie it shoud contain
%   NB_PATIENTS (16) individual directories (P%%contours-auto),
%   each one containing text files, eg P13-0040-icontour-auto.txt
% - both of these directories are in the current directory (where main_perf_assessment is).
%
% Warning! This program does not check for missing files.
% It lists files which are present in the automatic contour directories (P%%contours-auto).
% 
% INPUTS:
% Please see the Configuration section (l.39)
%
% OUTPUTS:
% - One txt file named TeamPrefix_ImageResults.txt, containing technical parameters for all images of all patients;
% - One txt file named TeamPrefix_PatientResults.txt, containing averaged technical parameters and clinical parameters for all patients;
% - One txt file named TeamPrefix_StatsResults.txt, containing technical parameters averaged over
%   patients and clinical parameters compared to ground truth and averaged over patients.
% - Optional plotting of manual vs. automatic data (area, ED and ES volumes, EF, mass if relevant) and linear regression fitting.
%   This option may be turned on/off in the corrcoef_regression.m function.
% For more details about what is in each txt file, please refer to the "Evaluation code" document.
%
% This program calls the following functions:
% - compare_rv_contours.m: returns Dice metric (DM) and Hausdorff distance (HD), and contour areas from a pair of contours.
% - corrcoef_regression.m: returns correlation coefficient and linear regression fit parameters from two measurements.

%   Copyright: LITIS EA 4108, Université de Rouen, France
%   Author: Caroline Petitjean (caroline.petitjean@univ-rouen.fr)
%   Revision: 2.0 - Date: June 25th, 2012
%  (bugs l.109/110 corrected from Rev.1.0 thanks to WB/SW) 

close all
clear all

%%%%%%%%%%%%%%%% Configuration section %%%%%%%%%%%%%%%%%%%%%%%
%This should be modified to match your configuration
EPI = 0; %EPI = 1 if your algorithm segments epicardium contours, EPI = 0 otherwise
%TeamPrefix = 'FakeAutomaticContours';
TeamPrefix = 'MyAutomaticContours';
%TeamPrefix = 'LITISAutomaticContours';
TrainingSetPrefix = 'TrainingSet';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abspathdir = cd;
manualdir = [abspathdir '\' TrainingSetPrefix '\'];
autodir = [abspathdir '\' TeamPrefix '\'];

%Number of patients in the training set
NB_PATIENTS = 5;

image_results_file = fopen([TeamPrefix '_ImageResults.txt'],'w');
patient_results_file = fopen([TeamPrefix '_PatientResults.txt'],'w');

%Array initialization
all_endoED=[]; all_endoES=[]; all_epiED=[]; all_epiES=[];
all_volmlED=[]; all_volmlES=[]; all_ef=[]; all_vm=[];
all_endo_area=[]; all_epi_area=[];
mean_epiED=[-1; -1; -1]; mean_epiES=[-1; -1 ;-1]; vm=[-1;-1];

for j=1:NB_PATIENTS
    
    fprintf('Processing patient #%d...\n', j)
    if j < 10    PatientPrefix = ['0' int2str(j)];
    else  PatientPrefix = int2str(j); end
    
    %BEGIN - get patient dicom info
    %NP = 20 and SpaceBetweenSlices = 8.4
    hdr1 = dicominfo([manualdir 'patient' PatientPrefix '\P' PatientPrefix 'dicom\P' PatientPrefix '-0000.dcm']);
    hdr2 = dicominfo([manualdir 'patient' PatientPrefix '\P' PatientPrefix 'dicom\P' PatientPrefix '-0020.dcm']);
    SpaceBetweenSlices = abs(hdr1.SliceLocation-hdr2.SliceLocation);
    ps = hdr1.PixelSpacing;
    NP = hdr1.CardiacNumberOfImages;
    M = double(hdr1.Rows);    N = double(hdr1.Columns);
    %END - get patient dicom info
    
    manualcontdir =  ['patient' PatientPrefix '\P' PatientPrefix 'contours-manual\']
    autocontdir = ['P' PatientPrefix 'contours-auto\']
    manualcontdir_list = dir([manualdir manualcontdir]);
    autocontdir_list = dir([autodir autocontdir]);
    
    volED=[0;0]; volES=[0;0];  epivolED=[0;0];
    endoED=[]; endoES=[]; epiED=[]; epiES=[]; 
    
    %Get patient images (do not process the first two entries of autocontdir_list (. and ..))
    for i=3:length(autocontdir_list)
        
        %Get image name
        manualcontname = manualcontdir_list(2*i-3).name;
        autocontname = autocontdir_list(i).name;
        manualcont = textread([manualdir manualcontdir manualcontname]);
        autocont = textread([autodir autocontdir autocontname]);   
        
        %Get phase number to know if the image is ED or ES
        %Filename example: P01-0080-icontour-auto.txt
        phasenb = str2num(autocontname(6:8));
        io = autocontname(10);
        
        
        %Compare manual vs automatic contour
        %area(1): auto, area(2): manual,
        %TechParam(1):DM, (2):HD
        
        [area(1), area(2), TechParam(1), TechParam(2)] = compare_rv_contours ( autocont, manualcont, N, M );
         TechParam(2) = TechParam(2)*ps(1);
        fprintf(image_results_file,'P%s-0%s-%s %.4f %.2f\n', PatientPrefix, autocontname(6:8), io, TechParam(1), TechParam(2));
        
        %Concatenate data for later volume, EF and mass computation
        if(io == 'i')
            all_endo_area = [all_endo_area area'*ps(1)*ps(2)*.01];
            if mod(phasenb,NP) == 0 %it's ED 
                endoED = [endoED TechParam'];
                volED  = volED + area';
            else %it's ES
                endoES = [endoES TechParam'];
                volES  = volES + area';
            end
        elseif(io == 'o')
            all_epi_area = [all_epi_area area'*ps(1)*ps(2)*.01];
            if mod(phasenb,NP) == 0 %it's ED 
                epiED = [epiED TechParam'];
                epivolED = epivolED + area';
            else
                epiES = [epiES TechParam'];
            end
        end
    end
    
    %Compute ED and ES volumes, EF (for both automatic and manual contours)
    volmlED = SpaceBetweenSlices*volED*ps(1)*ps(2)/1000;
    volmlES = SpaceBetweenSlices*volES*ps(1)*ps(2)/1000;
    ef = (volmlED-volmlES)./volmlED;
    
    if EPI
        mean_epiED = mean(epiED,2);
        mean_epiES = mean(epiES,2);      
        epivolmlED = SpaceBetweenSlices*epivolED*ps(1)*ps(2)/1000;
        
        % ventricular mass (g), with density : 1.05g/cm^3 = 0.001050g/mm^3
        vm = (epivolmlED-volmlED)*1.050;
    end
    
    fprintf (patient_results_file,'P%s-DM %.4f %.4f %.4f %.4f\n', PatientPrefix, mean(endoED(1,:)), mean(endoES(1,:)), mean_epiED(1), mean_epiES(1));
    fprintf (patient_results_file,'P%s-HD %.2f %.2f %.2f %.2f\n', PatientPrefix, mean(endoED(2,:)), mean(endoES(2,:)), mean_epiED(2), mean_epiES(2));
    fprintf (patient_results_file,'P%s %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n', PatientPrefix, volmlED(1), volmlED(2), volmlES(1), volmlES(2), ef(1), ef(2), vm(1), vm(2));
    
    %Concatenate data for later statistics computation
    all_endoED = [all_endoED endoED];     all_endoES = [all_endoES endoES];
    all_volmlED = [all_volmlED volmlED];  all_volmlES = [all_volmlES volmlES];
    all_ef = [all_ef ef];
    if EPI
        all_epiED = [all_epiED epiED];    all_epiES = [all_epiES epiES];
        all_vm = [all_vm vm];
    end
end

fclose(image_results_file);
fclose(patient_results_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Statistical results: averaging over patients, linear regression studies for clinical parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Computing statistical results...\n');
stats_results_file = fopen([TeamPrefix '_StatsResults.txt'],'w');

%Mean and std DM, PC, HD for endocardium (separately for ED and ES phases)
fprintf (stats_results_file,'Mean (std) endo DM: ED: %.4f (%.2f), ES: %.4f (%.2f)\n', mean(all_endoED(1,:)), std(all_endoED(1,:)), mean(all_endoES(1,:)), std(all_endoES(1,:)));
fprintf (stats_results_file,'Mean (std) endo HD: ED: %.2f (%.2f), ES: %.2f (%.2f)\n', mean(all_endoED(2,:)), std(all_endoED(2,:)), mean(all_endoES(2,:)), std(all_endoES(2,:)));

%Mean and std DM, PC, HD for epicardium (separately for ED and ES phases)
if EPI
    fprintf (stats_results_file,'Mean (std) epi DM: ED: %.4f (%.2f), ES: %.4f (%.2f)\n', mean(all_epiED(1,:)), std(all_epiED(1,:)), mean(all_epiES(1,:)), std(all_epiES(1,:)));
    fprintf (stats_results_file,'Mean (std) epi HD: ED: %.2f (%.2f), ES: %.2f (%.2f)\n', mean(all_epiED(2,:)), std(all_epiED(2,:)), mean(all_epiES(2,:)), std(all_epiES(2,:)));
else
    fprintf (stats_results_file,'Mean (std) epi DM: ED: -1.0000 (-1.00), ES: -1.0000 (-1.00)\n');
    fprintf (stats_results_file,'Mean (std) epi HD: ED: -1.00 (-1.00), ES: -1.00 (-1.00)\n');
end

%Mean and std DM, PC, HD for endocardium
fprintf (stats_results_file,'\n');
fprintf (stats_results_file,'Total mean (std) endo DM: %.4f (%.2f)\n', mean([all_endoED(1,:) all_endoES(1,:)]), std([all_endoED(1,:) all_endoES(1,:)]));
fprintf (stats_results_file,'Total mean (std) endo HD: %.2f (%.2f)\n', mean([all_endoED(2,:) all_endoES(2,:)]), std([all_endoED(2,:) all_endoES(2,:)]));

%Mean and std DM, PC, HD for epicardium
if EPI
    fprintf (stats_results_file,'Total mean (std) epi DM: %.4f (%.2f)\n', mean([all_epiED(1,:) all_epiES(1,:)]), std([all_epiED(1,:) all_epiES(1,:)]));
    fprintf (stats_results_file,'Total mean (std) epi HD: %.2f (%.2f)\n', mean([all_epiED(2,:) all_epiES(2,:)]), std([all_epiED(2,:) all_epiES(2,:)]));
else
    fprintf (stats_results_file,'Total mean (std) epi DM: -1.0000 (-1.00)\n');
    fprintf (stats_results_file,'Total mean (std) epi HD: -1.00 (-1.00)\n');
end

%Comparison between manual vs automatic areas, volumes (ED ans ES), EF, mass
endo_area_stats = corrcoef_regression(all_endo_area(2,:),all_endo_area(1,:),'Endocardium area');
volED_stats     = corrcoef_regression(all_volmlED(2,:), all_volmlED(1,:),'ED volume');
volES_stats     = corrcoef_regression(all_volmlES(2,:), all_volmlES(1,:),'ES volume');
ef_stats        = corrcoef_regression(all_ef(2,:), all_ef(1,:), 'Ejection fraction');
ef_error        = abs(all_ef(2,:)-all_ef(1,:))./all_ef(1,:);

%Variables set to -1 in case there is no epicardium segmentation
vm_error = -1;
vm_stats = struct('R',-1,'a',-1,'b',-1);
epi_area_stats = struct('R',-1,'a',-1,'b',-1);
mean_vm_error = -1;
std_vm_error = -1;
if EPI
    epi_area_stats = corrcoef_regression(all_epi_area(2,:),all_epi_area(1,:),'Epicardium area');
    vm_stats = corrcoef_regression(all_vm(2,:), all_vm(1,:), 'Ventricular mass');
    vm_error = abs(all_vm(2,:)-all_vm(1,:))./all_vm(1,:);
    mean_vm_error = mean(vm_error); std_vm_error = std(vm_error);
end

fprintf (stats_results_file,'\n');
fprintf (stats_results_file,'Endo area: R = %.4f, a = %.4f, b = %.4f\n', endo_area_stats.R, endo_area_stats.a, endo_area_stats.b);
fprintf (stats_results_file,'Epi area: R = %.4f, a = %.4f, b = %.4f\n', epi_area_stats.R, epi_area_stats.a, epi_area_stats.b);
fprintf (stats_results_file,'ED vol: R = %.4f, a = %.4f, b = %.4f\n', volED_stats.R, volED_stats.a, volED_stats.b);
fprintf (stats_results_file,'ES vol: R = %.4f, a = %.4f, b = %.4f\n', volES_stats.R, volES_stats.a, volES_stats.b);
fprintf (stats_results_file,'EF: R = %.4f, a = %.4f, b = %.4f, mean (std) error = %.4f (%.2f)\n', ef_stats.R, ef_stats.a, ef_stats.b, mean(ef_error), std(ef_error));
fprintf (stats_results_file,'vm: R = %.4f, a = %.4f, b = %.4f, mean (std) error = %.4f (%.2f)\n', vm_stats.R, vm_stats.a, vm_stats.b, mean_vm_error, std_vm_error);

fclose(stats_results_file);
fprintf('Done.\n');
