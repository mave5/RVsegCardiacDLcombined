% open figures and put them into subplots
clc
close all
clear all
addpath functions
%%
patient=29;

imset='Test1Set';

% get current folder
curr_dir=pwd;
figfolder=['Results/',imset,'/patient',num2str(patient),'/figs/'];
cd (figfolder);
figs_list= dir('*.fig');


% patient 29
ES=[2,4,6,8,10,12,14];

% patient 18
%ED=[1,2,3,5,7,9,11,13,15];
%ES=[4,6,8,10,12,14,16];


% patient 33
%ED=[1,2,3,5,7,9,11,13,15,17];
%ES=[4,6,8,10,12,14,16];

% patient 42
%ED=[];
%ES=[4,6,8,10,12,14,16,18,20,22];

Contours=ES;

k1=0;
for k=1:length(Contours)
k1=Contours(k);
fignk=figs_list(k1).name;
h1 = openfig(fignk,'reuse'); % open figure
ax(k) = gca; % get handle to axes of figure
end

% return to current folder
cd(curr_dir);

h3 = figure; %create new figure
set(gca,'YDir','reverse');

nf_row=2;
nf_col=4;
subplot1(nf_row, nf_col, 'Gap', [-.02 -.02]);
for k=1:length(Contours)
%s1 = subplottight(2,5,k); %create and get handle to the subplot axes
s1 = subplot1(k); %create and get handle to the subplot axes
fig1 = get(ax(k),'children'); %get handle to all the children in the figure
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
set(gca,'YDir','reverse');
%set(gca,'XDir','reverse');
axis off
%colormap(gray)
colormap(h3,gray)
brighten(.3)
end



