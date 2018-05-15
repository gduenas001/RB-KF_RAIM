
function [f_u_max]= worst_case_fault(n_F,m_F,m,n_faults,Hp_1, DeltaI_DS,S,alpha)

% total number of sensor measurements
n= n_F * m_F;

% Computes all possible combinations
combinations= nchoosek(1:n_F,n_faults);
n_comb= length(combinations);


% Initialize the extraction matrix depending on previous faults
if Hp_1
    E_0= zeros(m_F*n_faults + m,n+m);
    E_0(end-(m-1):end,end-(m-1):end)= eye(m);
else
    E_0= zeros(m_F*n_faults,n+m);
end
    
% Find the fault with worst slope
g2_F_max= 0;
for i= 1:n_comb
    E= E_0;
    i_F= combinations(i,:);
    for j= 1:n_faults
        indcol= (m_F*i_F(j) - (m_F-1)):(m_F*i_F(j));
        indrow= (j*m_F - (m_F-1)):j*m_F;
        E(indrow,indcol)= eye(m_F);
    end
    f_u= E' * ( (E * (DeltaI_DS) * E' ) \ E * S' * alpha);
    g2_F= alpha'* S * f_u;
    if g2_F > g2_F_max
        g2_F_max= g2_F;
        f_u_max= f_u;
        E_max= E;
    end
end



