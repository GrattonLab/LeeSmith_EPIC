function L=listgen1(n, r)
% makes a list of elements so that all first and second order associations are matched.
%
% listgen1(n, r)
% n = number of different element types
% r = number of repetitions of each element



State = zeros(n);
Goal = ones(n);
p = 1;

while ~min(min(State == Goal))
	State = zeros(n);
	L(1) = p;
	Goal = State + r;
   for i=1:n*n*r
        o = State(:,p);
        m = min(o);
        o = find(o==m);
        j = floor(rand*length(o))+1;
        v = o(j);
        L(i+1) = v;
        State(v,p) = State(v,p)+1;
        p = v;
   end
end


% State
% Goal
% L=Goal == State
