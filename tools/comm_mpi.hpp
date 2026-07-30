// comm_mpi -- a real-MPI implementation of the same commt API as comm_thread.hpp, so the
// distributed DepressionHierarchy harness compiles to a cluster build with the algorithm code
// unchanged. The thread shim (comm_thread.hpp) validates in-process; this runs one rank per MPI
// process (launch with `mpirun -np <ntiles>`).
//
// API parity with comm_thread.hpp: CommRank/CommSize/CommSend<T>/CommRecv<T>/CommBarrier, plus
// CommStartup/CommShutdown (MPI_Init/Finalize here; no-ops in the shim, whose lifecycle is the
// thread spawn/join inside CommInit -- so the harness driver brackets the run with them
// unconditionally, no #ifdef). Messages are cereal-serialized, matching (from,tag) as the shim does.
//
// Sends are non-blocking (MPI_Isend), matching the shim's fire-and-forget queue semantics so the
// harness's send-then-recv exchange patterns cannot deadlock. Each send's serialized buffer is
// kept alive in a pending list and completed at CommShutdown -- validation-grade (buffers live
// to the end); a production build would reclaim them per phase.
#pragma once

#include <mpi.h>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>

#include <list>
#include <sstream>
#include <string>
#include <utility>

namespace commt {

inline std::list<std::pair<MPI_Request,std::string>> g_pending;   // in-flight Isend buffers

// Process lifecycle. Matched with comm_thread's no-op CommStartup/CommShutdown so the driver can
// bracket the run identically under either backend (no #ifdef around the lifecycle calls).
inline void CommStartup(){ MPI_Init(nullptr, nullptr); }
inline void CommShutdown(){
  for(auto &p : g_pending) MPI_Wait(&p.first, MPI_STATUS_IGNORE);
  g_pending.clear();
  MPI_Finalize();
}

inline int CommRank(){ int r; MPI_Comm_rank(MPI_COMM_WORLD, &r); return r; }
inline int CommSize(){ int s; MPI_Comm_size(MPI_COMM_WORLD, &s); return s; }

template<class T>
void CommSend(const T &obj, int dest, int tag){
  std::ostringstream ss(std::ios::binary);
  { cereal::BinaryOutputArchive ar(ss); ar(obj); }
  g_pending.emplace_back(MPI_Request(), ss.str());               // keep the buffer alive for Isend
  auto &slot = g_pending.back();
  MPI_Isend(slot.second.data(), (int)slot.second.size(), MPI_BYTE, dest, tag, MPI_COMM_WORLD, &slot.first);
}

template<class T>
void CommRecv(T &obj, int from, int tag){                        // blocking; matches (from,tag)
  MPI_Status st;
  MPI_Probe(from, tag, MPI_COMM_WORLD, &st);
  int cnt; MPI_Get_count(&st, MPI_BYTE, &cnt);
  std::string buf; buf.resize(cnt);
  MPI_Recv(&buf[0], cnt, MPI_BYTE, from, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::istringstream ss(buf, std::ios::binary);
  cereal::BinaryInputArchive ar(ss); ar(obj);
}

inline void CommBarrier(){ MPI_Barrier(MPI_COMM_WORLD); }

// Source-compatibility shim for the thread backend's CommInit(n, fn): under MPI the ranks are
// processes that already exist, so just run the per-rank body once (this process is one rank).
template<class Fn>
void CommInit(int n, Fn fn){ (void)n; fn(); }

} // namespace commt
