function [u_nmpc, uopt, xopt, sopt, lopt, vopt, ocp_res] = ...
    switching_nmpc(xhat, xtrg, ukm1, xkm1, skm1, lkm1, vkm1, args)
    %Nonlinear MPC with similar implementation to nmpc.m.
    %This scheme switches Q and Qf matrices depending on proximity to 
    %target bearing
    %Obtaining initial values
    x0 = xhat; %Initial state vector
    %Initial guesses for z, inequality and equality duals;
    u = ukm1; x = xkm1; s = skm1; l = lkm1; v = vkm1;
    nx = 4; %State dimension
    nu = 1; %Input dimension
    N = args.N; %Finite horizon
    n = nx*N; %Primal state dimension
    m = nu*N; %Primar input dimension
    nlp = args.nlp; % Functions required for OCP implementation
    z = [u;x;s]; %Getting initial z matrix (not used for program resolution)
    
    for i = 1:args.max_sqp_iters
       
       %Computing h (inequality function), g(equality function) and
       %dJ (Cost function gradient)
       if cos(abs(x0(3)-xtrg(3)))<0.90
           dJ = full(nlp.laplace_grad(x0,x,u,s,l,v,args.xub,args.xlb,args.uub,args.ulb,...
               args.Q,args.R,args.Qf,xtrg,args.gamma));
       else
           dJ = full(nlp.laplace_grad(x0,x,u,s,l,v,args.xub,args.xlb,args.uub,args.ulb,...
               args.Q2,args.R,args.Qf2,xtrg,args.gamma));
       end
       h = full(nlp.in_reg(x,u,s,args.xub,args.xlb,args.uub,args.ulb));
       g = full(nlp.eq_reg(x0,x,u));
       
       %Evaluating FNR to determine whether it is below the given tolerance
       FNR = [dJ;g;min(-h,v)];
       ocp_res = norm (FNR);
       if ocp_res - args.tol <= 0
           break;
       end
        
       %Computing G (equality gradient), C (inequality gradient)and
       %H (Cost function Hessian)
       if cos(abs(x0(3)-xtrg(3)))<0.90
           H = full(nlp.cost_hess(args.Q,args.R,args.Qf,x,xtrg));
           d = full(nlp.cost_grad(x,u,args.Q,args.R,args.Qf,xtrg,args.gamma));
       else
           H = full(nlp.cost_hess(args.Q2,args.R,args.Qf2,x,xtrg));
           d = full(nlp.cost_grad(x,u,args.Q2,args.R,args.Qf2,xtrg,args.gamma));
       end
       C = nlp.in_grad();
       C = full(C.hz);
       G = full(nlp.eq_grad(x,u,x0));
       
       %Quadrtic programming implementation
       qpopts = optimoptions('quadprog','display','off');
       [dz,~,~,~,duals] = quadprog(H,d,C,-h,G,-g, [],[],[],qpopts);
       %If quadratic program does not converge, have the program keep going
       if isempty(dz)
           dz = zeros(size(z));
       else
           %If program converges to solution, get duals
           l = duals.eqlin;
           v = duals.ineqlin;
       end
       %Increase u, x and s
       u = u + dz(1:m);
       x = x + dz(m+1:m+n);
       s = s + dz(m+n+1:end);
       
    end
    
    %Obtain the final (optimal) values calculated
    uopt = u;
    xopt = x;
    sopt = s;
    lopt = l;
    vopt = v;
    
    u_nmpc = u(1);

end