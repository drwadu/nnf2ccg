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
% Author: hecher (2021-03-09 19:19:11)
%
edge(1,10).
edge(1,3).
edge(1,14).
edge(2,35).
edge(2,30).
edge(3,21).
edge(4,34).
edge(4,19).
edge(5,35).
edge(5,36).
edge(6,26).
edge(7,34).
edge(7,38).
edge(8,11).
edge(8,13).
edge(9,27).
edge(9,14).
edge(10,31).
edge(11,41).
edge(11,18).
edge(12,24).
edge(12,13).
edge(13,29).
edge(15,37).
edge(16,18).
edge(16,44).
edge(17,41).
edge(17,21).
edge(18,41).
edge(19,33).
edge(20,26).
edge(20,43).
edge(22,32).
edge(22,39).
edge(23,33).
edge(23,42).
edge(24,43).
edge(25,44).
edge(25,38).
edge(28,40).
edge(28,30).
edge(29,40).
edge(31,37).
edge(32,42).
edge(32,39).
