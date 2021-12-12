clear
clc

%swarm
Nswarm = 50;
sightRange = 5;
smellRange = 7;
prefDist = [ 3 5 ];
maxVel = 6;
minVel = 3;

%predator
Vp = 8; %note: this value should be greater then the max velocity of the swarm
tEnter = 6;
pInitPos = [ -1 -1 -1 ];
pSmellRange = 10;
first = true;
sizeRad = 0.5;

%random swarm
initPosR = 8;
initVelR = 3;
sPos = (rand([Nswarm,3])*2-1)*initPosR;
sVel = (rand([Nswarm,3])*2-1)*initVelR;
rad = 4*prefDist(1)*ones(Nswarm+1,1);
rad(length(rad)) = 2*rad(length(rad));

%examination area
AR = 10;
area = [ -AR AR -AR AR -AR AR ];
pPos = [ 2*AR 2*AR 2*AR ];
pVel = zeros(1,3);

%time
Tf = 25;
Ts = 0.05;
t = 0:Ts:Tf;


TimeC = 0;
while (TimeC < Tf)
    %find nearby
    nearBy = zeros(Nswarm);
    %creature examined
    for i = 1:Nswarm
        count = 1;
        %cerature compared
        for k = 1:Nswarm
            if (norm(sPos(i,:) - sPos(k,:)) < sightRange && k ~= i )
                nearBy(i,count) = k;
                count = count + 1;
            end
        end
    end
    
    maxNear = zeros(Nswarm,1);
    for i = 1:Nswarm
        [~,I] = min(nearBy(i,:));
        maxNear(i) = I - 1;
    end
    
    %adjust Velocity
    %for each creature
    for i = 1:Nswarm
        vadj = [ 0 0 0 ];
        near = nearBy(i,1:maxNear(i));
        if (~isempty(near))
            for j = 1:length(near)
                displ = sPos(near(j),:) - sPos(i,:);
                if( norm(displ) < prefDist(1) || norm(displ) > prefDist(2) )
                    adj = displ/norm(displ).*(norm(displ) - prefDist(1));
                else
                    adj = [ 0 0 0 ];
                end
                vadj = vadj + adj;
            end
            
            %tendancy to center
            tendancy = 0.8;
            center = -tendancy * sPos(i,:)/norm(sPos(i,:));

            vAvg = mean(sVel(near,:),1);
            
%             avgDir = ( vAvg/norm(vAvg) + sVel(i,:)/norm(sVel(i,:)) )/2;
%             scale = norm(sVel(i,:)) + ( norm(vAvg) - norm(sVel(i,:)) )/2;
%             sVel(i,:) = avgDir*scale + vadj;
            
            sVel(i,:) = mean([sVel(i,:);vAvg]) + vadj + center;
        end
        
        distP = pPos-sPos(i,:);
        
        if( norm(distP) < smellRange )
            sVel(i,:) = -distP/norm(distP)*norm(sVel(i,:))*1.5;
        end
        
        if ( norm(sVel(i,:)) > maxVel )
            sVel(i,:) = sVel(i,:)/norm(sVel(i,:))*maxVel;
        end
        if ( norm(sVel(i,:)) < minVel )
            sVel(i,:) = sVel(i,:)/norm(sVel(i,:))*minVel;
        end
    end
   
    
    %bounce if on edge
    for j = 1:Nswarm
        if sPos(j,1) < area(1)
            sVel(j,1) = -1*sVel(j,1);
            sPos(j,1) = area(1)*0.98;
            
        elseif sPos(j,1) > area(2)
            sVel(j,1) = -1*sVel(j,1);
            sPos(j,1) = area(2)*0.98;
            
        elseif sPos(j,2) < area(3)
            sVel(j,2) = -1*sVel(j,2);
            sPos(j,2) = area(3)*0.98;
            
        elseif sPos(j,2) > area(4)
            sVel(j,2) = -1*sVel(j,2);
            sPos(j,2) = area(4)*0.98;
            
       elseif sPos(j,3) < area(5)
            sVel(j,3) = -1*sVel(j,3);
            sPos(j,3) = area(5)*0.98;
            
        elseif sPos(j,3) > area(6)
            sVel(j,3) = -1*sVel(j,3);
            sPos(j,3) = area(6)*0.98;
        end
    end
    
    
    %determine preadator actions
    if( TimeC > tEnter )
        nearP = zeros(Nswarm,1);
        count = 1;
        
%         for i = 1:Nswarm
%             if( abs(norm( sPos(i,:) - pPos)) < pSmellRange )
%                 nearP(count) = i;
%                 count = count + 1;
%             end
%         end
%         
%         [~,index] = min(nearP);
%         nearP = nearP(1:index-1);
%         
%         pVelAdj = [ 0 0 0 ];
%         
%         for i = 1:length(nearP)
%             pVelAdj = pVelAdj + (sPos(i,:) - pPos)*norm(sPos(i,:) - pPos)/10;
%         end
%         
%         pVel = pVel + pVelAdj;
%         pVel = pVel/norm(pVel)*Vp;
        
        if(first)
            pPos = pInitPos/norm(pInitPos)*1.5*AR;
            pVel = -pInitPos/norm(pInitPos)*Vp;
            first = false;
        else
            pPos = pPos + Ts*pVel;
        end
    end
    
    colors = ones(Nswarm+1,3);
    
    for i = 1:Nswarm
        colors(i,:) = [ 0 0 256 ];
    end
    colors(length(colors),:) = [ 256 0 0 ];
    
    
    view(3) 
    hold on
    [x,y,z] = sphere(6);
    cla
    for i = 1:Nswarm
        surf( x/2-sPos(i,1), y/2-sPos(i,2), z/2-sPos(i,3) )
    end

    surf( x-pPos(1), y-pPos(2), z-pPos(3) );

    %plot swarm
    sPos = sPos + sVel*Ts;
    
%     scatter3( [ sPos(:,1)' pPos(1) ], [ sPos(:,2)' pPos(2) ], [ sPos(:,3)' pPos(3) ], rad , colors )
    axis(area)
    grid on
    pause(Ts)
    TimeC = TimeC + Ts;
end

% [x,y,z]=sphere;
% surf(x,y,z)
% surf(x/10,y/10,z/10)
% surf(x/10-5,y/10,z/10)
%     