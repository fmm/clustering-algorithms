#ifndef MINIMUM_COST_MAXIMUM_FLOW_H_
#define MINIMUM_COST_MAXIMUM_FLOW_H_

#include "lib.h"

struct MinimumCostMaximumFlow {
	vector<int> graph;
	vector<int> from, to, weight, cap, prev;
	void add(unsigned int u, unsigned int v, int w, int c, bool reverse = true) {
		while(u >= graph.size()) graph.push_back(-1);
		from.push_back(u);
		to.push_back(v);
		weight.push_back(w);
		cap.push_back(c);
		prev.push_back(graph[u]);
		graph[u] = from.size() - 1;
		if(reverse) add(v, u, -w, 0, false);
	}
	pair<int,int> process(unsigned int source, unsigned int sink) {
		pair<int,int> p(0,0);
		for(unsigned int n = graph.size();;) {
			vector<int> dist(n, inf), path(n, -1);
			dist[source] = 0;
			for(unsigned int repeat = 0; repeat < n; ++repeat) {
				for(unsigned int node = 0; node < n; ++node) for(int edge = graph[node]; edge >= 0; edge = prev[edge]) {
					if(dist[node] < inf and cap[edge] and dist[node] + weight[edge] < dist[to[edge]]) {
						dist[to[edge]] = dist[node] + weight[edge];
						path[to[edge]] = edge;
					}
				}
			}
			if(dist[sink] >= inf) break;
			int cost = 0, flow = inf;
			for(unsigned int node = sink; node != source; node = from[path[node]]) if(cap[path[node]] < flow) {
				flow = cap[path[node]];
			}
			for(unsigned int node = sink; node != source; node = from[path[node]]) {
				cap[path[node]] -= flow, cap[path[node]^1] += flow;
				cost += weight[path[node]] * flow;
			}
			p.first += cost, p.second += flow;
		}
		return p;
	}
};

#endif
