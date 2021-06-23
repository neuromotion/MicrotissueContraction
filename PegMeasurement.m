addpath("/Volumes/Emily_Project/Data/Project Data")
%this code is to measure the size of the pegs in your image set. You can
%then use that peg radius in the PegAlignment code
clear all
close all

PhaseSample1 = 'Phase_6_2_1_12'; %change this to file name

FileType = '.jpg';
PhaseFile = strcat(PhaseSample1,FileType);

pegs = 1; % set the number of pegs you want to measure based on how many you can see in your images. You want to be able to see enough of the peg to tell if the radius esimation is close.


I = imread(PhaseFile); %open image file
R1= []; % this is the matrix of averaged radius values measured from each peg
xcyc1=[]; %this is the calculated center point of each measured peg
for i = 1:pegs
figure(1)
imshow(I);
    circlePoints =[];
    triplePoints = [];
    title(strcat('Measurement 1 - PEG', num2str(i), ' - Select 3 points'))
    a1 = drawpoint('Color', 'w');
    triplePoints = cat(1, triplePoints,a1.Position);
    a2 = drawpoint('Color', 'w');
    triplePoints = cat(1, triplePoints,a2.Position);
    a3 = drawpoint('Color', 'w');
    triplePoints = cat(1, triplePoints,a3.Position);
    circlePoints = cat(2, circlePoints, triplePoints)
    [R,xcyc] = fit_circle_through_3_points(circlePoints);
    r = drawcircle('Center', xcyc(:,1)', 'Radius',R(:,1),'Color', 'm');
RA= R;
XA = xcyc;
figure(2)
imshow(I);

circlePoints =[];
triplePoints = [];
 title(strcat('Measurement 2 - PEG', num2str(i), ' - Select 3 points'))
a1 = drawpoint('Color', 'w');
    triplePoints = cat(1, triplePoints,a1.Position);
    a2 = drawpoint('Color', 'w');
    triplePoints = cat(1, triplePoints,a2.Position);
    a3 = drawpoint('Color', 'w');
    triplePoints = cat(1, triplePoints,a3.Position);
    circlePoints = cat(2, circlePoints, triplePoints)
    [R,xcyc] = fit_circle_through_3_points(circlePoints);
    r = drawcircle('Center', xcyc(:,1)', 'Radius',R(:,1),'Color', 'm');
RB = R
XB = xcyc;
figure(3)
imshow(I);
circlePoints =[];
triplePoints = [];
 title(strcat('Measurement 3 - PEG', num2str(i), ' - Select 3 points'))
a1 = drawpoint('Color', 'w');
    triplePoints = cat(1, triplePoints,a1.Position);
    a2 = drawpoint('Color', 'w');
    triplePoints = cat(1, triplePoints,a2.Position);
    a3 = drawpoint('Color', 'w');
    triplePoints = cat(1, triplePoints,a3.Position);
    circlePoints = cat(2, circlePoints, triplePoints)
    [R,xcyc] = fit_circle_through_3_points(circlePoints);
    r = drawcircle('Center', xcyc(:,1)', 'Radius',R(:,1),'Color', 'm');
RC = R
XC = xcyc;

R = (RA + RB +RC)/3
XX = (XA(1,1) + XB(1,1) +XC(1,1))/3;
XY = (XA(2,1) + XB(2,1) +XC(2,1))/3;
xcyc = cat(1,XX,XY)
R1 = cat(2,R1,R);
xcyc1 = cat(2,xcyc1,xcyc)
end 
'Peg radii'
R1
