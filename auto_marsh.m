%--------------------------------------------------------------------------
%   Author: Nicoletta Leonardi email: nicleona@bu.edu 
%   Copyrigth 2015 (C) Nicoletta Leonardi
%--------------------------------------------------------------------------
%   START of copying permission statement
%--------------------------------------------------------------------------
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   contact: Nicoletta Leonardi 
%   nicleona@bu.edu 
%   nicoletta_leonardi@hotmail.it
%--------------------------------------------------------------------------
%   END of copying permission statement
%--------------------------------------------------------------------------
%   These manuscripts are for an overview, information and/or background 
%   for the program, and citations to these would be appreciated:
%   Leonardi, N., and S. Fagherazzi (2014), How waves shape salt marshes, Geology , doi:10.1130/G35751.1.
%   http://geology.gsapubs.org/content/early/2014/08/28/G35751.1.abstract
%   Leonardi, N., and S. Fagherazzi (2015), Local variability in erosional resistance affects large scale 
%   morphodynamic response of salt marshes to wind waves, 
%   Geophysical Research Letters, 2015GL064730, doi:10.1002/2015GL064730.
%   http://onlinelibrary.wiley.com/doi/10.1002/2015GL064730/abstract
%--------------------------------------------------------------------------
function [eroded_x, eroded_y, countstored,randustar,t,time,N_bin,interval1,interval2,count_if_power_data]=...
auto_marsh(H,Power,bins,shottoplot,variability,ext,value_data,Power_data,H_data)
h=figure;
%   INPUT
%   H = number of grid cells (HxH=Domain size)
%   POWER=mean wave power value
%   bins = number of intervals you want to devide your total computational
%   time (dt=run_time/bins)
%   shottoplot = updata the figure every "shottoplot" time steps.
%	variability= is the standard deviaion of soil resistance. note that
%	resistance value cannot go below zero. which means
%   ( e.g. of the critical height (Terzaghi) of the scarp)
%   value data=is the average resistance (e.g. average critical height of
%   the marsh scarp)
%   Power_data= This is if you have real data of wave power to use for the
%   test;
%   H_data=This is if you have real data of wave height to test;

%   OUTPUT
%   eroded_x and eroded_y = store whata you have eroded. to plot again the
%   profile do plot(eroded_x,eroded_y,'bo')
%   r1 = erosion rate distribution
%   t = the time vector. at t(i) at least one cell is eroded (consider
%   t(1:length(countstored to do the plots)
%   time = when I stop to run the simulation. the time required to
%   erode half of the domain cells
%   countstored = How many cell I have eroded
%   N_bin = number of cells thet are eroded every dt with(dt=run_time/bins)
%   intervallo = time steps interval where I take the N_bin
%
%
%   %RUN EXAMPLES (uncomment lines below to run)
% 
%	%run HIGH wave power
%	H=100;
%	Power=100;
%	value_data=1;                                                          % mean critical height values
%	variability=value_data./3;                                             % standard deviation of critical height values; resistance values follow normal distribution; resistance values should not go below zero                                                                         % note value_data-3.*variability>0 otherwise resistance go below zero. 
%	ext=0;                                                                 % extreme events probability
%	Power_n=Power*ones(100000,1);
%	shottoplot=50;                                                             
%	bins=200;
%	cg=0.5;
%	E=Power_n./cg;
%	H_data=sqrt(8.*E./(1000.*9.81));
% 
%	[eroded_x, eroded_y, countstored,randustar,t,time,N_bin,interval1,...
%     interval2,count_if_power_data]=...
%     auto_marsh(H,Power,...
%     bins,shottoplot,variability,ext,value_data,Power,H_data);
% 
%	figure
%	plot(N_bin)                                                            % erosion events in time
%	ylabel('erosion events')
%	figure
%	hist(N_bin,bins);                                                      % see frequency magnitude distribution erosion events
%	ylabel('istogram erosion events');
%    
% 	%run LOW wave power
%	H=100;
%	Power=0.0025;
%	value_data=1;                                                          % mean critical height values
%	variability=value_data./3;                                             % standard deviation of critical height values; resistance values follow normal distribution; resistance values should not go below zero                                                                         % note value_data-3.*variability>0 otherwise resistance go below zero. 
%	ext=0;                                                                 % extreme events probability
%	Power_n=Power*ones(100000,1);
%	shottoplot=50;                                                             
%	bins=200;
%	cg=0.5;
%	E=Power_n./cg;
%	H_data=sqrt(8.*E./(1000.*9.81));
% 
%	[eroded_x, eroded_y, countstored,randustar,t,time,N_bin,interval1,...
%     interval2,count_if_power_data]=...
%     auto_marsh(H,Power,...
%     bins,shottoplot,variability,ext,value_data,Power,H_data);
% 
%	figure
%	plot(N_bin)                                                            % erosion events in time
%	ylabel('erosion events')
%	figure
%	hist(N_bin,bins);                                                      % see frequency magnitude distribution erosion events
%   ylabel('istogram erosion events');
%    

rr=length(Power_data);
if rr==1;
    P=Power;
    if value_data./H_data(1)>=500
        error('Please keep ratio 10^-40< resistance/wave heigh < 500');
%         it doesn't work well with values which are toot small, too many
%         0/Inf numbers');
    end
    
else
    tt=find(value_data./H_data>=500);
    H_data(tt)=value_data./500;
end


%	Initialize the system
% Rstar=ones(H,H+1);
% %    Initialize the matrix of resistance values cellulino
% randustar=value_data'.*((1-2.*variability).*rand(H,H-1)+variability);                         %    create a matrix of random values
% Rstar(:,2:(end-1))=0.35*P.^1.1*exp(-randustar/P);                                  %	give random resistance values to each grid cell
% r1star=Rstar./(P.^1.1);                                                            %    calculate the erosion rate dependent on cell resistance and wave power
% Rstar(:,1)=zeros(H,1);                                                         %	Set to zero startin
% Rstar(:,end)=zeros(H,1);

front_x=2*ones(H,1);                                                       %	The x coordinate of the eroded front 
front_y=(1:H)';                                                            %    The y_coordinate of the eroded front                                                                      %    Store eroded cells
eroded_x=ones(H,1);                                                        %    eroded_x is the x coordinate of eroded material
eroded_y=[(1:H)';(1:H)'];                                                  %    eroded_y is the y coordinates of eroded material
eroded_x=[eroded_x;(H+1)*ones(H,1)];

%    Now decide how much the simualtion should run and set a vecor size for
%    the dead particles equal to one half of the particle of the domain
%    (i.e. run the simulation until half fo the domain cells are eroded)

stopsimulation=4;
nn=floor(H*(H-1)/stopsimulation)+1;  
eroded_x=[eroded_x;ones(nn,1)];
eroded_y=[eroded_y;ones(nn,1)];

%   start the computation
ind=0;
t=zeros(nn,1);
k=1;
countremoved=0;
interval1=[];
interval2=[];

count_if_power_data=0;

R=ones(H,H+1);                                                             %    Initialize the matrix of resistance values
% randu=value_data'.*((1-2.*variability).*rand(H,H-1)+variability);
randu=value_data'+variability.*randn(H);
randu(randu<0)=value_data;
randu(:,end)=[];                                                           %    note: set a variability value such that randu doesn't go below zero.
randustar=randu;

% nl=10^50;

while countremoved<nn;
    count_if_power_data=count_if_power_data+1;
    extreme=(100-1).*rand(1,1) + 1;                                        %   random number between 1 and 100. ext give the probability of having an extreme event.
    
    if extreme<=ext ;
        P=100*Power;
        R(:,2:(end-1))=(0.35*P.^1.1*1);                                    %	give random resistance values to each grid cell
        r1=R./(P.^1.1);
        
    elseif rr>1
        P=Power_data(count_if_power_data);
        HSTAR=H_data(count_if_power_data);
        R(:,2:(end-1))=(0.35*P.^1.1*exp(-randu./HSTAR));                   %	give random resistance values to each grid cell
        r1=R./(P.^1.1);
    else
        HSTAR=nanmean(H_data);
        R(:,2:(end-1))=(0.35*P.^1.1*exp(-randu./HSTAR));                   %	give random resistance values to each grid cell
        r1=R./(P.^1.1);                                                    %    calculate the erosion rate dependent on cell resistance and wave power    
    end;
    
    
    k=k+1;
    eroded_k=(k-1)+2*H;                                                    %	There is the 2*H because I removed the first column at the beginning
    %======================================================================
    %---Find the amount of time required to erode the next marsh portion
    %   This time depend on the resistance of the front
    I=sub2ind([H,H],front_y,front_x);
    r=R(I);                                                                %    This is the resistance of the front
    %   This is the amount of time necessary because one cell is eroded
    %   with probability 1. The sites are eroded based on fron resistance.
    delta_t=sum(r)^-1;
    t(k)=t(k-1)+delta_t;
    %======================================================================
    %---Now decide which cell will be eroded
    %	Pick a random cell. For each cell the probability of being selected
    %   is random such that the probability of selecting it is r1/sum(r)
    %   where sum(r) is the resistance of the front
    p=rand(1)/delta_t;                                                     %	Choose a number at random between 0 and sum(r) which is the inverse of delta_t
    c=cumsum(r);                                                           %    cumulative probability
    pt=c-p;
    index=sum(pt<0);                                                       %	this is the index of the element that I am going to remove at random with probability p from the front
    index=index+1;
    [ROW,COL]=ind2sub([H,H+1],I(index));                                   %    ind2sub calculate the equivalent matrix indiecs for a certain element number
    %======================================================================
    %---find the neighbours of the previously eroded cell
    neigh_row=[ROW+1;ROW-1;ROW;ROW];
    neigh_col=[COL;COL;COL+1;COL-1];
    %   Once you have found the neighbours remove those neighbours which
    %   are external to the domain
    delete_this=find(neigh_row<1 | neigh_row>H);                           %    Find neighbours which are external to domain
    neigh_col(delete_this)=[];                                             % 	Delete the columns of these neighbours
    neigh_row(delete_this)=[];                                             %    Delete the row of these neighbours
    pos_i=sub2ind([H,H+1],neigh_row,neigh_col);                            %    Index of the neighbours of the eroded cells                                                                           %	remove the point wth probability p from the front
    %======================================================================
    front_x(index)=[];                                                     %    Here is where I remove a new cell
    front_y(index)=[];                                                     %    Here is where I remove a new cell
    eroded_x(eroded_k)=COL;                                                %	The column of the eroded cell (the cell is the "eroded_k" eroded cell
    eroded_y(eroded_k)=ROW;                                                %	The row of the eroded cell (the cell is the "eroded_k" eroded cell
    
    %======================================================================
    %	Erode when you are isolated
    eroded_test_x=eroded_x;
    eroded_test_y=eroded_y;
    eroded_test=[eroded_test_x eroded_test_y];
    eroded_test_toadd=[];
    %   find the maximum x coordinate of the dead particles for each row
    dd=[eroded_y, eroded_x];
    dd=sortrows(dd,1);
    dd1=dd;
    dd(:,1)=dd1(:,2);
    dd(:,2)=dd1(:,1);
    
    for kk=1:H;
        coor=find(dd(:,2)==kk);                                            %    choose  acertain column number
        doveprende=dd(coor,1);                                             %    doveprende="where to look" for the maximum take a certain row of eroded vectors
        doveprende(find(doveprende==H+1))=[];                              %    doveprende="where to look" for the maximum
        max_x=max(doveprende);                                             %    find the maximum
        eroded_x_toadd=[1:max_x]';                                         %    things that you need to erode because are isolated
        eroded_y_toadd=kk*ones(max_x,1);                                   %    things that you need to erode because are isolated
        eroded_test_toadd=[eroded_test_toadd; eroded_x_toadd eroded_y_toadd];
        
    end;
    %======================================================================
    %	Clean the situation. Remove from the vectors things repeated twice
    eroded_x_toadd=eroded_test_toadd(:,1);
    eroded_y_toadd=eroded_test_toadd(:,2);
    if isempty(eroded_x_toadd)==0;
        for jj=1:length(eroded_x_toadd);
            ff=find(eroded_x==eroded_x_toadd(jj));
            for gg=1:length(ff);
                if eroded_y(ff(gg))==eroded_y_toadd(jj);
                    eroded_x_toadd(jj)=0;
                    eroded_y_toadd(jj)=0;
                end
            end
        end
    end
    %======================================================================
    eroded_x_toadd(find(eroded_x_toadd==0))=[];
    eroded_y_toadd(find(eroded_y_toadd==0))=[];
    eroded_test_toadd=[eroded_x_toadd eroded_y_toadd];
    eroded_test=[eroded_test; eroded_test_toadd];
    
    eroded_x=eroded_test(:,1);
    eroded_y=eroded_test(:,2);
    eroded_i=sub2ind([H,H+1],eroded_y(1:eroded_k),eroded_x(1:eroded_k));
    removed_i=[eroded_i;I];                                                %	these are eroded plus front
    
    
    for j=1:length(pos_i);
        %	Do the loop for the 4 neighbours, at this point if they are not
        %	already eroded or part of the front they are added to the front.
        %	Here I update the front
        
        if sum(pos_i(j)==removed_i)==0;
            front_y=[front_y;neigh_row(j)];
            front_x=[front_x;neigh_col(j)];
            
            %==============================================================
            %	Clean from the front particles which remained behind
            front=[front_x front_y];
            front_delete=[];
            for kk=1:H
                coor=find(front(:,2)==kk);
                doveprende=front(coor,1);
                doveprende(find(doveprende==H+1))=[];
                
                %   Where to take the right neighbour
                coor_r=find(front(:,2)==kk-1);
                doveprende_r=front(coor_r,1);
                doveprende_r(find(doveprende_r==H+1))=[];
                
                %   where to take the left neighbour
                coor_l=find(front(:,2)==kk+1);
                doveprende_l=front(coor_l,1);
                doveprende_l(find(doveprende_l==H+1))=[];
                
                max_x=min(max(doveprende_r),min(max(doveprende),max(doveprende_l)));
                front_x_delete=[1:max_x-1]';
                front_y_delete=kk*ones(max_x-1,1);
                front_delete=[front_delete; front_x_delete front_y_delete];
            end;
            
            %   find rows present in both front_delete and in the front=[front_x front_y];
            %   I delete them because it means that they are behind the
            %   maximum fron point
            
            [~,indx]=ismember(front_delete,front,'rows');
            indx(find(indx==0))=[];
            front(indx,:)=[];
            
            %             plot(eroded_x,eroded_y,'go');
            %             hold on;
            %             plot(front(:,1),front(:,2),'ro');
            %==============================================================
            
        end
    end
    
    if floor(ind/shottoplot)==ind/shottoplot;
        plot(eroded_x,eroded_y,'wd','MarkerSize',5,'MarkerFaceColor','w')
        set(gca,'Color','k');
        hold on;
        axis([0,H,0,H]);
        drawnow;
    end
    ind=ind+1;
    countremoved=length(union([1 1],[eroded_x eroded_y],'rows'))-1;
    countstored(k)=countremoved;
end

% start some statistics
% bins=200;
t(t==0)=[];
t=[0 t'];
t_bin=linspace(0,max(t),bins+1);
time=max(t);
N_bin=ones(1,bins);
intervallo=[];
for i=2:bins-1;
    
    interval1_i=find(t<t_bin(i));
    interval1_i=interval1_i(end);
    interval1=[interval1 interval1_i];
    interval2_i=find(t>t_bin(i+1));
    if isempty(interval2_i)==1
        interval2_i=nan;
        interval2=[interval2 interval2_i];
        a=[interval1_i interval2_i];
        intervallo=[intervallo a'];
        N_bin(i)=nan;
    else
        interval2_i=interval2_i(1);
        interval2=[interval2 interval2_i];
        
        a=[interval1_i interval2_i];
        intervallo=[intervallo a'];
        N_bin(i)=countstored(interval2_i)-countstored(interval1_i);
    end
end
% plot(N_bin);
% end

%to save
% name=num2str(P);
% saveas(h,name,'fig');
end






