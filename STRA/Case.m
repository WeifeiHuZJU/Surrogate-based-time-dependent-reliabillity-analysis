function [X_mean,X_std,x_mcs,z_mcs,Stochasticpro,X_initial]=Case(problem,total_Nmcs,N_initial)
%没有等价转换
switch problem % two mathematical examples       
    case 1
        X_mean=[3.5 3.5];X_std=[0.25 0.25];
        Nnode=31;
        F_Nz=3;%number of terms in the expansion
        F_mu=@(t) 0;%mean value of the random field
        F_std=@(t) 1;%std of the random field
        corr_func=@(x1,x2) exp(-(x1-x2).^2);
        T=[0 1];
        tnode=linspace(T(1),T(2),Nnode);
        Stochasticpro.T=T;
        Stochasticpro.Y_Nz=F_Nz;
        Stochasticpro.Nnode=Nnode;
        Stochasticpro.Y_mu=F_mu;
        Stochasticpro.Y_std=F_std;
        Stochasticpro.corr_func=corr_func;
        Stochasticpro.tnode=tnode;
               
        sample_initial=lhsdesign(N_initial,(length(X_mean)+F_Nz+1),'criterion','maximin','iterations',100);        
        x(:,1:length(X_mean))=sample_initial(:,1:length(X_mean))*8.*repmat(X_std,N_initial,1)+X_mean-4*X_std;
        z(:,1:F_Nz)=sample_initial(:,length(X_mean)+1:length(X_mean)+F_Nz)*8-4;
        time(:,1)=sample_initial(:,length(X_mean)+F_Nz+1)*(T(2)-T(1))+T(1);        
        for i=1:N_initial
            Ft(i,1)=EOLE_out(time(i),z(i,:),Stochasticpro);
        end
        X_initial=[x,Ft,time];         
        for h=1:length(X_mean)
            x_mcs(:,h)=normrnd(X_mean(h),X_std(h),total_Nmcs,1);
        end
        for m=1:F_Nz
            z_mcs(:,m)=normrnd(0,1,total_Nmcs,1);
        end
    

    case 2
        X_mean=[0.22 0.025 0.019 0.025];
        X_std=[2.2e-3 2.5e-4 1.9e-4 2.5e-4];
        am=[3.815 2.528 1.176 -0.07856];
        bm=[0.2895 0.5887 0.7619 2.183];
        cm=[-0.2668 0.9651 3.116 -3.161];
        as=[0.7382 1.013 1.875 1.283];
        bs=[6.456 4.075 9.913 1.035];
        cs=[0.9193 1.561 6.959 2.237];
        Nnode=200;
        F_Nz=2;%number of terms in the expansion
        F_mu=@(t) sum(am.*sin(bm*t+cm));%mean value of the random field
        F_std=@(t) sum(as.*exp(-((t-bs)./cs).^2));%std of the random field
        corr_func=@(x1,x2) cos(2*pi*(x2-x1));
        T=[0 12];
        tnode=linspace(T(1),T(2),Nnode);        
        Stochasticpro.T=T;
        Stochasticpro.Y_Nz=F_Nz;
        Stochasticpro.Nnode=Nnode;
        Stochasticpro.Y_mu=F_mu;
        Stochasticpro.Y_std=F_std;
        Stochasticpro.corr_func=corr_func;
        Stochasticpro.tnode=tnode;
        sample_initial=lhsdesign(N_initial,(length(X_mean)+F_Nz+1),'criterion','maximin','iterations',100);
        x(:,1:length(X_mean))=sample_initial(:,1:length(X_mean))*10.*repmat(X_std,N_initial,1)+X_mean-5*X_std;
        z(:,1:F_Nz)=sample_initial(:,length(X_mean)+1:length(X_mean)+F_Nz)*10-5;
        time(:,1)=sample_initial(:,length(X_mean)+F_Nz+1)*(T(2)-T(1))+T(1);           
        for i=1:N_initial
            Ft(i,1)=EOLE_out(time(i),z(i,:),Stochasticpro);
        end
        X_initial=[x,Ft];
        for h=1:length(X_mean)
            x_mcs(:,h)=normrnd(X_mean(h),X_std(h),total_Nmcs,1);
        end        
        for m=1:F_Nz
            z_mcs(:,m)=normrnd(0,1,total_Nmcs,1);
        end
       
end

end