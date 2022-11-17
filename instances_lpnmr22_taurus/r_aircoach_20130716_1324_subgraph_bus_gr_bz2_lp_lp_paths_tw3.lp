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
% Author: hecher (2021-03-09 19:19:08)
%
edge(1,45).
edge(1,46).
edge(2,58).
edge(2,26).
edge(3,33).
edge(3,4).
edge(4,11).
edge(4,16).
edge(4,63).
edge(5,56).
edge(5,36).
edge(5,43).
edge(5,20).
edge(6,32).
edge(6,14).
edge(7,35).
edge(7,11).
edge(8,65).
edge(8,41).
edge(9,57).
edge(9,45).
edge(10,65).
edge(10,17).
edge(12,56).
edge(12,36).
edge(13,28).
edge(13,31).
edge(14,57).
edge(15,36).
edge(15,39).
edge(16,30).
edge(17,42).
edge(17,50).
edge(17,18).
edge(18,47).
edge(19,25).
edge(19,58).
edge(20,23).
edge(21,40).
edge(21,47).
edge(22,29).
edge(22,53).
edge(23,58).
edge(24,28).
edge(24,39).
edge(25,56).
edge(25,34).
edge(26,42).
edge(26,58).
edge(27,33).
edge(27,50).
edge(29,49).
edge(30,35).
edge(31,52).
edge(31,66).
edge(31,60).
edge(32,64).
edge(34,41).
edge(36,36).
edge(36,38).
edge(36,48).
edge(37,54).
edge(38,56).
edge(40,63).
edge(43,48).
edge(44,62).
edge(44,55).
edge(46,51).
edge(49,66).
edge(51,59).
edge(52,61).
edge(53,54).
edge(55,57).
edge(60,61).
edge(61,64).
edge(62,64).
