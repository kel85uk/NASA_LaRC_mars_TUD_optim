function [X,fval,exitflag,optim_output] = optim_ez(X0,optionsga,options,Xmin,Xmax,mu,mu_1,mu_2,stay_return,dt_stay,Isp)

Nx = length(X0);
A = [eye(Nx,Nx);-eye(Nx,Nx)];
b = [Xmax;-Xmin];

func = @(X) transfer_ellipse_optim_obj(X,mu,mu_1,mu_2,stay_return,dt_stay,Isp);

[X0,fval,exitflag,ga_output,population] = ga(func,length(X0),A,b,[],[],Xmin,Xmax,[],optionsga);

[X,fval,exitflag,optim_output] = fmincon(func,X0,A,b,[],[],[],[],[],options);

end
