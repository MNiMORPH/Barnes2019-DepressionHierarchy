#include "tools/comm_thread.hpp"
#include <cstdio>
#include <vector>
#include <atomic>
namespace c=commt;
static std::atomic<int> failures{0};
static void rank_main(){
  const int r=c::CommRank();
  if(r==0){
    std::vector<int> strip{10,20,30};           // stand-in for a (elev/label) edge strip
    c::CommSend(strip, 1, 42);
  } else if(r==1){
    std::vector<int> got;
    c::CommRecv(got, /*from*/0, /*tag*/42);
    if(got!=std::vector<int>{10,20,30}) failures++;
    printf("rank 1 got %zu ints: %d %d %d\n", got.size(), got[0],got[1],got[2]);
  }
  c::CommBarrier();
  printf("rank %d/%d past barrier\n", r, c::CommSize());
}
int main(){
  c::CommInit(2, rank_main);
  printf(failures==0 ? "COMM-SMOKE-OK\n" : "COMM-SMOKE-FAIL\n");
  return failures==0 ? 0 : 1;
}
