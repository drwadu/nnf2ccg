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
% Author: hecher (2021-03-09 19:19:13)
%
edge(1,3).
edge(1,5).
edge(2,20).
edge(2,15).
edge(3,10).
edge(3,5).
edge(4,24).
edge(4,10).
edge(4,11).
edge(6,8).
edge(6,20).
edge(6,14).
edge(7,25).
edge(8,9).
edge(9,12).
edge(10,24).
edge(10,17).
edge(11,25).
edge(11,21).
edge(12,16).
edge(12,20).
edge(13,18).
edge(13,19).
edge(14,18).
edge(16,22).
edge(17,24).
edge(17,23).
edge(19,22).
edge(21,25).
edge(21,23).
