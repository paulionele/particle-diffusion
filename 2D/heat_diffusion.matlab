%Temperature distribution in a 2D rectangular plate.

%Initialization of initial temperature distribution
clear;clf
a = zeros(21,21);
a(1,1:10)=100;
ITD = a; %ITD: Initial Temperature Distribution

time = input('Enter the time (as an integer) for desired temperature distribution: ');
assert(time>0,'The time you have entered is not valid.')

%Calculating the temperature distribution as time advances.

for q = 0:time
    
    for x = 1:20
        for y = 1:20
            
            if (x<=10) && (y==1)%Dirichlet
                ITD(y,x) = 100;
            
            elseif (x<=20) && (y==20)%Dirichlet
                ITD(y,x) = 0;
            
            elseif (x==1) && (y>=2) && (y<=19)%Insulted boundary left
                ITD(y,x) = 0.25*( ITD(y-1,x)+ITD(y+1,x)+2*ITD(y,x+1) );
            
            elseif (x==20) && (y>=2) && (y<=19)%Insulated boundary right
                ITD(y,x) = 0.25*( ITD(y-1,x)+ITD(y+1,x)+2*ITD(y,x-1) );
            
            elseif (y==1) && (x>=11) && (x<=19)%Insulted boundary top
                ITD(y,x) = 0.25*( ITD(y,x-1)+ITD(y,x+1)+2*ITD(y+1,x) );
            
            elseif (x>=2) && (x<=19) && (y>=2) && (y<=19)%For all else
                ITD(y,x) = 0.25*( ITD(y,x-1)+ITD(y,x+1)+ITD(y-1,x)+ITD(y+1,x) );
            
            else
                
            end
        end
    end
    

surf(ITD)
    title(sprintf('$u(x,y,t=%5.0f)$', q),'Interpreter', 'latex', 'FontSize', 14);
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 12)
    zlabel('$temperature$', 'Interpreter', 'latex', 'FontSize', 12)
    
    view(140,20);
    axis([1,20,1,20])
    shg
    pause(0.1)

end