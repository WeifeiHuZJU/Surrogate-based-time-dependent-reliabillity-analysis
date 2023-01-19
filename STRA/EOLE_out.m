function [Y_t,eigenvalues_truncate,eigenvectors_truncate,CovVector]=EOLE_out(t_i,Z_Y,Stochasticpro)
% function Y_t1=EOLE_out(t_i,Z_Y,Stochasticpro)
% Expansion optimal linear estimation method used for the discretization of stochastic process

    T=Stochasticpro.T;
    Nnode=Stochasticpro.Nnode;
    Y_Nz=Stochasticpro.Y_Nz;
    Y_mu=Stochasticpro.Y_mu;
    Y_std=Stochasticpro.Y_std;
    corr_func=Stochasticpro.corr_func;

    tnode=linspace(T(1),T(2),Nnode);
    T1=repmat(tnode',1,Nnode);
    T2=repmat(tnode,Nnode,1);
    
    for i=1:Nnode
        for j=1:Nnode
            CovMatrix(j,i)=Y_std(T1(j,i)).*Y_std(T2(j,i)).*corr_func(T1(j,i),T2(j,i));
%             CorrMatrix(j,i)=corr_func([T1(j,i),T2(j,i)]);
        end
    end
    
    [eigenvectors eigenvalues]=eig(CovMatrix);
    eigenvalues = diag(eigenvalues);
    [eigenvalues, I] = sort(eigenvalues, 'descend');
    eigenvalues_truncate=eigenvalues(1 : Y_Nz);
    eigenvectors_truncate=eigenvectors(:, I(1 : Y_Nz));
    
%     [eigenvectors1 eigenvalues1]=eig(CorrMatrix);
%     eigenvalues1 = diag(eigenvalues1);
%     [eigenvalues1, I1] = sort(eigenvalues1, 'descend');
%     eigenvalues_truncate1=eigenvalues1(1 : Y_Nz);
%     eigenvectors_truncate1=eigenvectors1(:, I1(1 : Y_Nz));

    t_Vector=[repmat(t_i,Nnode,1),tnode'];
    for i=1:Nnode
        CovVector(i,1)=Y_std(t_Vector(i,1))*Y_std(t_Vector(i,2))*corr_func(t_Vector(i,1),t_Vector(i,2));
    end

%     CorrVector=corr_func([tnode',repmat(t_i,Nnode,1)]);
%     CorrVector=corr_func([repmat(t_i,Nnode,1),tnode']);

    Yrep = zeros(size(Z_Y,1),1);
%     Yrep1 = zeros(size(Z_Y,1),1);

    for j=1:Y_Nz%支配特征根数量
        Yrep=Yrep+(1/sqrt(eigenvalues_truncate(j)).*eigenvectors_truncate(:,j))'*CovVector*Z_Y(:,j);
%         Yrep1=Yrep1+Y_std(t_i)*(1/sqrt(eigenvalues_truncate1(j)).*eigenvectors_truncate1(:,j))'*CorrVector*Z_Y(:,j);
    end
    Y_t=Y_mu(t_i)+Yrep;
%     Y_t1=Y_mu(t_i)+Yrep1;
end