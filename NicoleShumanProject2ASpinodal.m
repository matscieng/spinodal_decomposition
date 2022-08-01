%Nicole Shuman                        %
%EMS 164                              %     
%Project 2A: Spinodal Decomposition   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: MATLAB help and MathWorks was used to ask questions and conceptually understand
% how to graph in 2-D. 
% I looked at this page (below) for help in understanding 2D graphing and the Cahn-Hilliard method; however, no code was utilized
% https://www.math.utah.edu/~eyre/computing/matlab-intro/ch.txt,

%Defining Variables
W=1; %Constant (Given)
epsilon=0.1; %Constant (Given)
M=1; %Constant (Given)
Cavg=0.5; %Initial Concentration (Given)
deltaT=0.0002;%Time Step (Given)
delX=0.1; %Step across simulation domain (Given)
time=(0):deltaT:(10); %Time
Nx=0:0.1:10; %Simulation Region with Step of 0.1 (100x100) - 10/.1=100
Ny=0:0.1:10; %Simulation Region with step of 0.1 (100x100) 
C=zeros(length(Nx),length(Ny),length(time)); %Concentration array placeholder(Zeros)
C(:,:,1)=Cavg+(rand(length(Nx),length(Ny))-0.5).*0.2; %Concentration indexing with respect to the axis arrays (Given)

%For Loop calculates for concentration for(x-dir, y-dir, time) for given equations 
for t_i=2:length(time) %Gives a specific point in time
    t=time(t_i);
    for N_i=1:length(Nx) %Gives a specific position along x-direction
        I=Nx(N_i);
        for N_j=1:length(Ny) %Gives a specific position along y-direction
            J=Ny(N_j);
            i=N_i;
            j=N_j;
            highi=N_i+1; %Prevents values from exceeding boundary conditions
            lowi=N_i-1; %Prevents values from exceeding boundary conditions
            highj=N_j+1; %Prevents values from exceeding boundary conditions
            lowj=N_j-1; %Prevents values from exceeding boundary conditions
            dfdc=(W./2).*C(i,j,t_i-1).*((1-C(i,j,t_i-1)).*(1-2.*(C(i,j,t_i-1)))); %df/dc equation used for laplacian function (input value)
            if i>=length(Nx) %Primary boundary conditions in x-direction (greater than)
               highi=1;
            end
            if i<=1 %Primary boundary conditions in x-direction (less than)
               lowi=length(Nx);
            end               
            if j>=length(Ny) %Primary boundary conditions in y-direction (greater than)
               highj=1;
            end
            if j<=1 %Primary boundary conditions in y-direction (less than)
               lowj=length(Ny);
            end
            laplace2=(C(highi,j,t_i-1)+(C(lowi,j,t_i-1))+(C(i,highj,t_i-1))+(C(i,lowj,t_i-1))-4.*(C(i,j,t_i-1)))./(delX.^2); %Curvature Equation
            Q(N_i,N_j)=(dfdc-epsilon.^2.*laplace2); %Defining input equation w/ variables for first laplacian function
       
        end
    end    
    
    for N_i=1:length(Nx) %Iterates over the x-array for each position of x
        N_x=Nx(N_i);
       
        for N_j=1:length(Ny) %Iterates over the y-array for each position of y
            N_y=Ny(N_j);
            i=N_i; %Ensures values do not exceed boundary conditions (x-dir)
            j=N_j; %Ensures values do not exceed boundary conditions (y-dir)
            highi=N_i+1; %Prevents values from exceeding boundary conditions
            lowi=N_i-1; %Prevents values from exceeding boundary conditions
            highj=N_j+1; %Prevents values from exceeding boundary conditions
            lowj=N_j-1; %Prevents values from exceeding boundary conditions
            
            if i>=length(Nx) %Primary boundary conditions in x-direction
               highi=1;
            end
            if i<=1 %Primary boundary conditions in x-direction
               lowi=length(Nx);
            end               
            if j>=length(Ny) %Primary boundary conditions in y-direction
               highj=1;
            end
            if j<=1 %Primary boundary conditions in y-direction
               lowj=length(Ny);
            end
            
            laplace=(Q(highi,j)+(Q(lowi,j))+(Q(i,highj))+(Q(i,lowj))-4.*(Q(i,j)))./((delX).^2); %Equation to define second laplacian function inside dc/dt equation
            dcdt=M.*laplace; %Concentration of position over time (Given)
            val=C(N_i,N_j,t_i-1)+dcdt.*deltaT; %Euler equation (predicts concentration in x,y over time)
            if val>1
               val=1;
            end
            if val<0
               val=0;
            end   
            C(N_i,N_j,t_i)=val;
            
        end    
    end
    %Plot (Code provided by Dr. Gentry)
    if  mod(t_i,250)==0
        pcolor(C(:,:,t_i))
        caxis([0 1])
        colorbar
        axis equal
        shading flat
        pause (0.1)
    end
end

