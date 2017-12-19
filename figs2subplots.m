% open figures and put them into subplots

clc
close all
clear all
addpath functions
%%

patient=2;
num_of_slices=15;
slice_num=1:num_of_slices;


for k1=1:length(slice_num)
k=slice_num(k1);
if k<10
snk=[num2str(floor(patient/10)),num2str(rem(k,10))];
else
    snk=k;
end

h1 = openfig(['figs/patient02/p02_slice',num2str(snk),'.fig'],'reuse'); % open figure
ax(k) = gca; % get handle to axes of figure
end


h3 = figure; %create new figure
set(gca,'YDir','reverse');


nf_row=3;
nf_col=5;
subplot1(nf_row, nf_col, 'Gap', [-.02 -.02]);
for k=1:num_of_slices
%s1 = subplottight(2,5,k); %create and get handle to the subplot axes
s1 = subplot1(k); %create and get handle to the subplot axes
fig1 = get(ax(k),'children'); %get handle to all the children in the figure
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
set(gca,'YDir','reverse');
%set(gca,'XDir','reverse');
axis off
colormap(gray)
end



