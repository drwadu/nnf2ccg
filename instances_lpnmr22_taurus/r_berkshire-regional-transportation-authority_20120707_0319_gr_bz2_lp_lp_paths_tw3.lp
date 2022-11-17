%start is reachable
reach(X) :- X=#min{ Y: edge(Y,_); Y: edge(_,Y); Y: link(Y,_); Y: link(_,Y) }.

%goal
goal(X) :- X=#max{ Y: edge(Y,_); Y: edge(_,Y); Y: link(Y,_); Y: link(_,Y) }.

%guess select
select(X,Y) :- reach(X), edge(X,Y), not nots(X,Y).
select(X,Y) :- reach(X), edge(Y,X), not nots(X,Y).
nots(X,Y) :- reach(X), edge(Y,X), not select(X,Y).
nots(X,Y) :- reach(X), edge(X,Y), not select(X,Y).

reach(Y) :- select(X,Y).

:- reach(X), select(X,Y), select(X,Y2), Y < Y2. 
:- goal(X), select(X,_). 
:- not reach(X), goal(X).

% Graphs generated from publicly available GTFS transit feeds by Johannes K. Fichte
%
% References
% [1] https://en.wikipedia.org/wiki/General_Transit_Feed_Specification or
% [2] https://developers.google.com/transit/gtfs/
% [3] https://github.com/daajoe/transit_graphs/blob/master/transitfeeds-tw.pdf
% [4] https://github.com/daajoe/gtfs2graphs
% [5] https://github.com/daajoe/transit_graphs/converters/pace2lp.py
%
% GTFS feeds extracted using gtfs2graphs [4]
% PACE format to lparse converter [5]
% Author: hecher (2021-03-09 19:18:48)
%
edge(1,16).
edge(1,24).
edge(1,18).
edge(2,18).
edge(2,12).
edge(3,25).
edge(3,46).
edge(3,62).
edge(4,58).
edge(4,35).
edge(4,54).
edge(5,34).
edge(6,24).
edge(6,26).
edge(6,61).
edge(7,8).
edge(7,19).
edge(8,24).
edge(8,10).
edge(8,44).
edge(9,40).
edge(9,41).
edge(9,54).
edge(10,44).
edge(11,24).
edge(11,23).
edge(12,37).
edge(12,42).
edge(12,49).
edge(12,21).
edge(12,55).
edge(13,47).
edge(14,32).
edge(14,24).
edge(15,53).
edge(15,21).
edge(16,24).
edge(16,28).
edge(16,49).
edge(17,36).
edge(19,33).
edge(20,24).
edge(20,44).
edge(21,53).
edge(21,30).
edge(22,52).
edge(22,62).
edge(23,33).
edge(24,33).
edge(24,38).
edge(24,44).
edge(24,45).
edge(24,53).
edge(24,27).
edge(25,57).
edge(25,39).
edge(26,61).
edge(27,37).
edge(28,43).
edge(29,48).
edge(29,43).
edge(30,53).
edge(31,57).
edge(31,59).
edge(34,38).
edge(35,44).
edge(36,42).
edge(37,42).
edge(39,57).
edge(40,54).
edge(40,41).
edge(41,54).
edge(42,56).
edge(44,54).
edge(44,58).
edge(45,47).
edge(46,57).
edge(48,51).
edge(50,57).
edge(51,56).
edge(52,55).
edge(54,58).
edge(59,60).
