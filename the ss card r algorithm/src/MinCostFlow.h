/*
 * MinCostFlow.h
 *
 *  Created on: May 18, 2010
 *      Author: fmm
 */

#ifndef MINCOSTFLOW_H_
#define MINCOSTFLOW_H_

#include "Includes.h"

#define INF 0x3f3f3f3f
#define MAXV 110
#define MAXE 1010
#define LIM (1<<16)

struct MinCostFlow {

  typedef int T;

  // Heap
  struct Heap {
    T no, dist;
    bool operator<(const Heap& h) const {
      return dist > h.dist;
    }
  } heap[LIM], atual;
  int Q;

  // Edge
  struct Edge {
    T u, v, cap, cost, ant;
  } edge[MAXE];
  int E;

  // Graph
  T adj[MAXV], pai[MAXV], dist[MAXV], pot[MAXV];
  int V;

  void init(int n) {
    memset(adj,-1,sizeof(adj));
    V = n, E = 0;
  }

  void add_edge(int u, int v, T cap, T cost, bool rev = false) {
    edge[E] = (Edge){ u, v, cap, cost, adj[u]};
    adj[u] = E++;

    edge[E] = (Edge){ v, u, 0, -cost, adj[v]};
    adj[v] = E++;

    if(rev) {
      add_edge(v,u,cap,cost);
    }
  }

  Heap& top() {
    pop_heap(heap,heap+Q);
    return heap[--Q];
  }

  void update(Heap& h, int e = -1) {
    if(h.dist < dist[h.no]) {
      dist[h.no] = h.dist, pai[h.no] = e;
      heap[Q++] = h;
      push_heap(heap,heap+Q);
    }
  }

  bool bellman_ford(int s) {
    memset(dist,INF,sizeof(dist));
    dist[s] = 0;
    for(int k = 0; k < V; k++) {
      for(int i = 0; i < V; i++) {
        for(int j = adj[i]; j != -1; j = edge[j].ant) {
          if(edge[j].cap <= 0) continue;
          T u = edge[j].u, v = edge[j].v, cost = edge[j].cost;
          if(dist[u] + cost < dist[v]) {
            dist[v] = dist[u] + cost;
          }
        }
      }
    }
    for(int i = 0; i < V; i++) {
      for(int j = adj[i]; j != -1; j = edge[j].ant) {
        if(edge[j].cap <= 0) continue;
        T u = edge[j].u, v = edge[j].v, cost = edge[j].cost;
        if(dist[u] + cost < dist[v]) { // negative cycle
          return true;
        }
      }
    }
    return false;
  }

  bool dijkstra(int s, int t) {
    memset(dist,INF,sizeof(dist)), Q = 0;
    memset(pai,-1,sizeof(pai));
    Heap h = {s,0};
    update(h);

    while(Q) {
      atual = top();
      if(atual.dist > dist[atual.no]) continue;
      for(int i = adj[atual.no]; i != -1; i = edge[i].ant) {
        if(edge[i].cap <= 0) continue;
        T u = edge[i].u, v = edge[i].v, cost = edge[i].cost;
        Heap h = { v, atual.dist + cost + pot[u] - pot[v] };
        update(h, i);
      }
    }

    return dist[t] < INF;
  }

  pair<T,T> solve(int s, int t) {

    memset(pot,0,sizeof(pot));

    /* dealing with negative costs */
    if(bellman_ford(s)) {
      assert(0); // negative cycle
    }

    for(int i = 0; i < V; i++) {
      pot[i] += dist[i];
    }
    /* dealing with negative costs */

    T flow = 0, cost = 0;
    while( dijkstra(s,t) ) {
      T vmin = INF;
      for(int i = pai[t]; i != -1; i = pai[edge[i].u]) {
        vmin = min(vmin, edge[i].cap);
      }
      flow += vmin;
      for(int i = pai[t]; i != -1; i = pai[edge[i].u]) {
        cost += (edge[i].cost * vmin);
        edge[i].cap -= vmin;
        edge[i^1].cap += vmin;
      }
      for(int i = 0; i < V; i++) {
        pot[i] += dist[i];
      }
    }

    return make_pair(flow,cost);
  }
};

#undef INF
#undef MAXV
#undef MAXE
#undef LIM

#endif /* MINCOSTFLOW_H_ */
