clear all;
close all;

Path = '/Volumes/Emily_Project/Phase Images/Exp 14/';
FileName='14_2_1_34';
FileType = '.png';


File = strcat(Path,FileName,FileType);
rgb = double(imread(File))/255;
figure
imshow(rgb)
gray=rgb2gray(rgb);
BW=imbinarize(gray);
figure
imshow(BW)
    
img=rgb2gray(rgb); 
N_CIRCLES=4; %interested in finding 4 pegs

if ~exist(strcat('CircleData/',FileName,'.mat'),'file') %if location of pegs already exists, you don't need to refind them
bad = true;
figure
     imshow(img)
Rectangles=zeros(N_CIRCLES,4);
%draw 4 rectangles around pegs of interest
    for i=1:N_CIRCLES
        r=drawrectangle;
        Rectangles(i,3)=r.Position(1);
        Rectangles(i,1)=r.Position(2);
        Rectangles(i,4)=r.Position(3)+r.Position(1);
        Rectangles(i,2)=r.Position(4)+r.Position(2);
    end
    Rectangles=round(Rectangles);

[BW,thresh]=edge(img,'Canny',[0.11,0.12]);%find circles w Canny edge detection
  
    N_ITERATIONS=10000;
    Rs=zeros(1,N_CIRCLES);
    xs=zeros(1,N_CIRCLES);
    ys=zeros(1,N_CIRCLES);
    figure
    for i=1:N_CIRCLES %N_CIRCLES is goal
        bad = true;
        while bad
            hold off
        BW1=BW(Rectangles(i,1):Rectangles(i,2),Rectangles(i,3):Rectangles(i,4));
        white_pixels=find(BW1);
        best=0;
        for j=1:N_ITERATIONS
            while 1
                y=datasample(white_pixels,3,'Replace',false);
                [r1,c1]=ind2sub(size(BW1),y(1));
                [r2,c2]=ind2sub(size(BW1),y(2));
                [r3,c3]=ind2sub(size(BW1),y(3));

                [R,xcyc] = fit_circle_through_3_points([r1,c1;r2,c2;r3,c3]);
                if all(~isnan(xcyc))
                    break
                end
            end
            xcyc=round(xcyc);
            R=round(R);
            [xc, yc] = getMidpointCircle(xcyc(1), xcyc(2), R);
            xc2=xc(xc>0&xc<size(BW1,1)&yc>0&yc<size(BW1,2));
            yc2=yc(xc>0&xc<size(BW1,1)&yc>0&yc<size(BW1,2));
            goodness=sum(BW1(sub2ind(size(BW1),xc2,yc2)))/length(xc2);
            if goodness > best && R > 220 && R < 250 %input peg radius based on PegMeasurement.m code 
                best = goodness;
                xs(i)=xcyc(1);
                ys(i)=xcyc(2);
                Rs(i)=R;
            end
        end
        [v,goodness]=fminsearch(@(x) fitCircle(x(1),x(2),x(3),BW1),[xs(i),ys(i),Rs(i)]);
        xs(i)=v(1);
        ys(i)=v(2);
        Rs(i)=v(3);
        [xc, yc] = getMidpointCircle(xs(i), ys(i), Rs(i));
        imshow(img)
        hold on
        plot(yc+Rectangles(i,3),xc+Rectangles(i,1),'o')
        rectangle('Position',[Rectangles(i,3), Rectangles(i,1), Rectangles(i,4)-Rectangles(i,3), Rectangles(i,2)-Rectangles(i,1)])
        bad=~input("Circle good (1=good, 0=bad)?");%if circle properly finds peg type 1, if not correct type 0
        end
    end
    figure
    imshow(img)
    hold on
    for i=1:N_CIRCLES
        [xc, yc] = getMidpointCircle(xs(i), ys(i), Rs(i));
        plot(yc+Rectangles(i,3),xc+Rectangles(i,1),'o')
        rectangle('Position',[Rectangles(i,3), Rectangles(i,1), Rectangles(i,4)-Rectangles(i,3), Rectangles(i,2)-Rectangles(i,1)])
        text((Rectangles(i,3)+Rectangles(i,4))/2,(Rectangles(i,1)+Rectangles(i,2))/2, num2str(i), 'color', 'w', 'fontsize', 20);
    end
    save(strcat('/Users/emily/Desktop/Thesis/Thesis Images/CircleData/',FileName),'Rs','xs','ys','Rectangles');
else
    load(strcat('/Users/emily/Desktop/Thesis/Thesis Images/CircleData/',FileName))
end

BW2=imbinarize(img); %binarize image
figure
imshow(BW2)

%need to translate the coordinates to be for the whole image
for i=1:N_CIRCLES
xt(i)=xs(i)+Rectangles(i,1);
yt(i)=ys(i)+Rectangles(i,3);
end

p1 = drawline('Position', [yt(1) xt(1);yt(2) xt(2)]); %draws line between pegs 1 and 2
p2 = drawline('Position', [yt(3) xt(3);yt(4) xt(4)]); %draws line between pegs 3 and 4

[xi, yi]= polyxpoly([xt(1),xt(2)], [yt(1),yt(2)], [xt(3),xt(4)], [yt(3),yt(4)]);

centerROI = [xi, yi];
rROI = 450 %ROI radius based on PegMeasurement.m code

figure
imshow(img)
f = drawcircle('Center', flip(centerROI), 'Radius', rROI, 'Color', 'c');
ROI = createMask(f);
fpeg = drawcircle('Center', [ys(1)+Rectangles(1,3),xs(1)+Rectangles(1,1)], 'Radius', Rs(1), 'Color', 'c');
f2 = drawcircle('Center', flip(centerROI), 'Radius', rROI/2, 'Color', 'c');
TMask=createMask(f2);
tissueVals=img(TMask);
PegMask=createMask(fpeg);
pegVals=img(PegMask);
figure
histogram(tissueVals)
hold on
histogram(pegVals)

thresh=input("Input a threshold between the two distributions:"); % clean extra white pixels w threshold
tempImg=bwareaopen(bwmorph(imbinarize(img.*ROI,thresh),"open",10),1000); % clean extra black pixels w open function
invTemp=bwareaopen((~tempImg).*ROI,100);
masked=ROI-invTemp; % combine
figure
imshow(masked)
invert=input("Invert binary?");%want white to represent tissue and black to represent background
if invert
    masked=~masked.*ROI;
end
imshow(masked)

result=imfuse(img,masked); %overlay binary on top of phase

cleaned=masked; %so you don't override masked
figure
imshow(result);
title('Please Draw ROI');
hold on;
%determine how many areas need hand curation to remove extra pixels
N_Pen=input('How Many Areas Need to be Hand Curated?');
for i=1:N_Pen
    h = drawfreehand('Closed')
    binaryMask = createMask(h);
    cleaned=~binaryMask.*cleaned; %update each time
%define imshow as variable then multiply by mask for each loop 
end 
imshow(cleaned)

contraction = (nnz(cleaned)/nnz(ROI))*100

%Centroid Code for Asymmetry
s = regionprops(cleaned,'centroid');
centroids = cat(1,s.Centroid);
X=[centroids(:,1),centroids(:,2);centerROI(:,2),centerROI(:,1)];
d = pdist(X,'euclidean')%this is the euclidean distance between center ROI and center tissue
imshow(cleaned)
hold on
plot(centroids(:,1),centroids(:,2),'b*')
hold on
plot(centerROI(:,2),centerROI(:,1), 'r*')