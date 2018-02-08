
#ifndef ADUPROP_PERF_HPP_
#define ADUPROP_PERF_HPP_

#include <map>
#include <iostream>

class cycle_counter
{
public:
  uint64_t StartCycleCount = 0;
  uint64_t CycleCount = 0;
  uint32_t HitCount = 0;
};

class perf {
  std::map<std::string, cycle_counter*> blocks;
  std::map<std::string, cycle_counter*>::iterator it;
  std::string name;
  
public:
  perf(std::string name_) : name(name_)  {}
  ~perf() {
    for(it = blocks.begin() ; it != blocks.end() ; ++it) {
      delete it->second;  
    }
  }
  
  void activate(std::string name) {
    if(blocks.find(name) != blocks.end()) return;
    cycle_counter *tmp = new cycle_counter;
    std::cout << "Activating " << name << " timer." << std::endl;
    blocks[name] = tmp;
  }
  
  void deactivate(std::string name) {
    it=blocks.find(name);
    if(it == blocks.end()) return;
    delete it->second;
    blocks.erase(it);
    std::cout << "Deactivating " << name << " timer." << std::endl;
  }
  
  void begin(std::string name) {
    it=blocks.find(name);
    if(it == blocks.end()) {
      return;
    }
    it->second->StartCycleCount = __rdtsc();
  }
  
  void end(std::string name) {
    it=blocks.find(name);
    if(it == blocks.end()) {
      return;
    }
    it->second->CycleCount+= __rdtsc() - it->second->StartCycleCount;
    ++it->second->HitCount;
  }
  friend std::ostream& operator<<(std::ostream&, perf&);
};

std::ostream& operator<< (std::ostream& os, perf& prof) {
  os << "Performance Statistics for " << prof.name << std::endl;
  std::map<std::string, cycle_counter*>::iterator it;
  
  for(it = prof.blocks.begin() ; it != prof.blocks.end() ; ++it) {
    if(it->second->HitCount != 0) {
      uint64_t div = 0;
      div = it->second->CycleCount/it->second->HitCount;
      os << it->first << " : " << it->second->CycleCount 
      << " cycles. " << it->second->HitCount << " hits. "
      << div << " cycles/hit." << std::endl;
      std::cout << std::endl;
    }
    return os;
  }
}
#endif  // ADUPROP_AD_HPP_
