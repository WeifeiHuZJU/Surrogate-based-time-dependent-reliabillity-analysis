function y=LSF(problem,x)
% g<0±íÊ¾Ê§Ð§
switch problem % two mathematical examples

    case 1 
        Yt=x(:,3);t=x(:,4);
        y=x(:,1).^2.*x(:,2)-5*x(:,1).*(1+Yt).*t+(x(:,2)+1).*t.^2-20;
       
    case 2
        l1=x(:,1);t1=x(:,2);t2=x(:,3);s_allow=x(:,4);v_t=x(:,5);
        I=l1.*(t1.^3-t2.^3)*2/3;
        M=0.5*1e3*v_t.^2.*0.3422;
        y=s_allow-M.*t1./(14e9*I);
       
end

end