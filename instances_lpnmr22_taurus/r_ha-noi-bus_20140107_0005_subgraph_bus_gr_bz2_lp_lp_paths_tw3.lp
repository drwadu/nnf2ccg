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
% Author: hecher (2021-03-09 19:19:37)
%
edge(1,22).
edge(1,39).
edge(2,16).
edge(2,27).
edge(3,9).
edge(3,75).
edge(3,37).
edge(4,28).
edge(4,76).
edge(5,36).
edge(5,39).
edge(6,35).
edge(6,47).
edge(7,44).
edge(7,22).
edge(8,20).
edge(8,68).
edge(9,10).
edge(10,28).
edge(11,62).
edge(11,23).
edge(12,65).
edge(12,18).
edge(13,26).
edge(13,30).
edge(14,27).
edge(14,21).
edge(15,74).
edge(15,34).
edge(16,45).
edge(17,19).
edge(17,79).
edge(18,67).
edge(19,62).
edge(20,78).
edge(21,44).
edge(23,33).
edge(24,73).
edge(24,55).
edge(25,45).
edge(25,71).
edge(26,32).
edge(27,70).
edge(29,66).
edge(29,36).
edge(29,38).
edge(30,61).
edge(31,35).
edge(31,78).
edge(32,48).
edge(33,57).
edge(34,50).
edge(38,42).
edge(40,55).
edge(41,65).
edge(41,53).
edge(42,77).
edge(43,56).
edge(43,64).
edge(46,52).
edge(46,61).
edge(47,67).
edge(48,58).
edge(49,60).
edge(49,52).
edge(50,59).
edge(51,75).
edge(51,68).
edge(53,69).
edge(54,77).
edge(54,79).
edge(56,57).
edge(58,73).
edge(59,71).
edge(60,67).
edge(63,72).
edge(63,76).
edge(66,69).
edge(72,74).
