% CRM.m
% 
% Consumer resource model, adapted from Marsland et al., 2020.
%
% Alan R. Pacheco, 04/23/2020

function DX=CRM(t,X)

global S M C D g w l m k d kR monodTerm sequential initConcs

%% Define system of equations



N = X(1:S);
R = X(S+1:end);

if sequential
    A = zeros(length(N),length(R));
    
    for s = 1:S
        
        [~,preferenceVec] = sort(C(s,:),'descend');
        
        A(s,find(preferenceVec == 1)) = 1; % Allow consumption of first preferred nutrient always
        for i = 1:length(R) % Now go through the rest of the nutrients in order of preference
            if R(find(preferenceVec == i)) < initConcs(find(preferenceVec == i))*0.99 % If a preferred nutrient is exhausted,
                A(s,find(preferenceVec == i+1)) = 1; % Turn on the next preferred one
            end
        end
    end
else
    A = ones(length(N),length(R));
end

for i = 1:S
    
    summation = 0;
    for a = 1:length(R)
        if monodTerm
            summation = summation + A(i,a)*w(a)*(1-l(a))*C(i,a)*(R(a)/(R(a)+kR(a)));
        else
            summation = summation + A(i,a)*w(a)*(1-l(a))*C(i,a)*R(a);
        end
    end
    
    DN(i) = g(i)*N(i)*(summation-m(i))- (1/d)*N(i);
end

for a = 1:length(R)
    
    [summation1,summation2] = deal(0);
    for i = 1:S
        
        if monodTerm
            summation1 = summation1 + A(i,a)*C(i,a)*N(i)*(R(a)/(R(a)+kR(a)));
            for b = 1:length(R)
                summation2 = summation2 + A(i,b)*D(a,b,i)*(w(b)/w(a))*l(b)*C(i,b)*N(i)*(R(b)/(R(b)+kR(b)));
            end
        else
            summation1 = summation1 + A(i,a)*C(i,a)*N(i)*R(a);
            for b = 1:length(R)
                summation2 = summation2 + A(i,b)*D(a,b,i)*(w(b)/w(a))*l(b)*C(i,b)*N(i)*R(b);
            end
        end
    end
    DR(a) = k(a) - (1/d)*R(a) - summation1 + summation2;
end

%% Prepare for solver
DX = [];
for i = 1:S
    DX = vertcat(DX,DN(i));
end
for a = 1:length(R)
    DX = vertcat(DX,DR(a));
end