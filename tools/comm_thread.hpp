// comm_thread -- a minimal thread-backed message-passing shim for validating the
// distributed DepressionHierarchy in-process (no cluster), matching the shape of
// richdem's Comm* API so the eventual code ports to real MPI unchanged.
//
// (richdem ships a communication-threads.hpp, but it is an unfinished WIP guarded by
// `#error This test module is not yet ready for compilation.` -- correctly kept from
// use. This is our own small, tested equivalent; the real MPI backend comes later.)
//
// Ranks are threads; each runs the same fn(). Messages are cereal-serialized so any
// cereal-serializable payload (incl. richdem Array2D, vectors, maps) crosses a rank.

#pragma once

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>

#include <atomic>
#include <list>
#include <mutex>
#include <sstream>
#include <thread>
#include <vector>

namespace commt {

struct Msg { int tag; int from; std::vector<char> bytes; };

inline int                              g_size = 0;
inline std::mutex                       g_mtx;
inline std::vector<std::list<Msg>>      g_queues;   // one inbox per rank
inline std::atomic<int>                 g_barrier{0};
inline std::atomic<int>                 g_epoch{0};
inline thread_local int                 t_rank = 0;

inline int CommRank(){ return t_rank; }
inline int CommSize(){ return g_size; }

// Process lifecycle. No-ops for the thread backend (the threads' whole life is CommInit, which spawns
// and joins them) -- matched with comm_mpi's MPI_Init/Finalize so the driver brackets the run the same
// way under either backend.
inline void CommStartup(){}
inline void CommShutdown(){}

// Spawn `n` thread-ranks, each running fn(). Joins them (CommShutdown is a no-op; join is in CommInit).
template<class Fn>
void CommInit(int n, Fn fn){
  g_size = n;
  g_queues.assign(n, {});
  g_barrier = 0; g_epoch = 0;
  std::vector<std::thread> pool;
  for(int i=0;i<n;i++) pool.emplace_back([i,fn]{ t_rank=i; fn(); });
  for(auto &t : pool) t.join();
}

template<class T>
void CommSend(const T &obj, int dest, int tag){
  std::ostringstream ss(std::ios::binary);
  { cereal::BinaryOutputArchive ar(ss); ar(obj); }
  const std::string s = ss.str();
  std::lock_guard<std::mutex> lock(g_mtx);
  g_queues[dest].push_back({tag, t_rank, std::vector<char>(s.begin(), s.end())});
}

template<class T>
void CommRecv(T &obj, int from, int tag){   // blocking; matches (from,tag)
  while(true){
    {
      std::lock_guard<std::mutex> lock(g_mtx);
      auto &q = g_queues[t_rank];
      for(auto it=q.begin(); it!=q.end(); ++it){
        if(it->from!=from || it->tag!=tag) continue;
        const std::string s(it->bytes.begin(), it->bytes.end());
        q.erase(it);
        std::istringstream ss(s, std::ios::binary);
        cereal::BinaryInputArchive ar(ss); ar(obj);
        return;
      }
    }
    std::this_thread::yield();
  }
}

inline void CommBarrier(){
  const int e = g_epoch.load();
  if(g_barrier.fetch_add(1)+1 == g_size){    // last arrival releases the rest
    g_barrier.store(0);
    g_epoch.fetch_add(1);
  } else {
    while(g_epoch.load()==e) std::this_thread::yield();
  }
}

} // namespace commt
