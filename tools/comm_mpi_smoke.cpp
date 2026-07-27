// comm_mpi_smoke -- exercise the real-MPI commt backend: a ring of edge-strip exchanges plus a
// gather to rank 0, mirroring the harness's patterns. Run with `mpirun -np <N>`; each rank sends
// its rank-id vector to its right neighbour and to rank 0, and rank 0 checks the gather.
#include "tools/comm_mpi.hpp"

#include <cstdio>
#include <vector>

namespace c = commt;

int main(){
  c::CommInitMPI();
  const int r = c::CommRank(), n = c::CommSize();
  int failures = 0;

  // Neighbour exchange (send-then-recv; must not deadlock with non-blocking sends).
  if(n>1){
    const int right = (r+1)%n, left = (r-1+n)%n;
    c::CommSend(std::vector<int>{r, r*10}, right, 1);
    std::vector<int> got;
    c::CommRecv(got, left, 1);
    if(got != std::vector<int>{left, left*10}) failures++;
  }

  // Gather to rank 0.
  if(r==0){
    std::vector<int> seen{0};
    for(int t=1;t<n;t++){ std::vector<int> v; c::CommRecv(v, t, 2); seen.push_back(v.empty()?-1:v[0]); }
    for(int t=0;t<n;t++) if(seen[t]!=t) failures++;
    printf(failures==0 ? "COMM-MPI-SMOKE-OK (n=%d)\n" : "COMM-MPI-SMOKE-FAIL\n", n);
  } else {
    c::CommSend(std::vector<int>{r}, 0, 2);
  }

  c::CommBarrier();
  c::CommFinalizeMPI();
  return failures==0 ? 0 : 1;
}
