function dx = zika_model_nocontrol(t,x,b,betaH,muH,rho,r,omega,delta,n,K,nu,betaM,muM,gamma,e,gamma_D,gamma_B,lamdaP,lag)

dx = zeros(11,1);

M = n*(x(8,1) + x(9,1));
f = M*(1-M/K);

star = gamma_D*(x(3,1)+x(4,1)) + gamma_B*(x(8,1)+x(9,1));

function new = check_neg(old_value, change)
% check if a variable is about to be lower than zero
    if  old_value + change < 0
        new = -old_value;
    else
        new = change;
    end
end 

% Humans (Fully Susceptible, Protected Susceptible, Infectious, Protected Infectious, Recovered, Severe Outcomes (D))
dx(1,1) = b*(x(1,1)+x(2,1)+x(3,1)+x(4,1)+x(5,1)) - betaH*x(9,1)*x(1,1) - star*x(1,1) - muH*x(1,1) + lamdaP*x(2,1);
dx(2,1) = star*x(1,1) - rho*betaH*x(9,1)*x(2,1) - muH*x(2,1) - lamdaP*x(2,1);
dx(3,1) = betaH*x(9,1)*x(1,1) - r*x(3,1) - muH*x(3,1) - omega*x(3,1) - star*x(3,1);
dx(4,1) = star*x(3,1) - lamdaP*x(4,1) - muH*x(4,1) + rho*betaH*x(9,1)*x(2,1);
dx(5,1) = r*x(3,1) - muH*x(5,1);
dx(6,1) = delta*x(3,1);

% Mosquitos (Larvae, Susceptible adults, Infectious adults)
%dx(7,1) = f - nu*x(7,1);
tmp = f - nu*x(7,1);
dx(7,1) = check_neg(x(7,1),tmp);

dx(8,1) = nu*x(7,1) - betaM*(x(3,1)+rho*x(4,1))*x(8,1) - muM*x(8,1);
dx(9,1) = betaM*(x(3,1)+rho*x(4,1))*x(8,1) - muM*x(9,1);

dx(10,1) = f;

end
