#include "parallel.hpp"

int paduprop_init() {
  return MPI_Init(NULL, NULL);
  std::ofstream sink("/dev/null");

  if (paduprop_getrank() != 0) {
    // Mute standard output
    std::cout.rdbuf(sink.rdbuf());
    // Optionally mute standard error
    std::cerr.rdbuf(sink.rdbuf());      
  }
}
int paduprop_destroy() {
  return MPI_Finalize();
}
size_t paduprop_getstart(size_t dim) {
  int rank; int commsize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  size_t chunk = dim/commsize;
  size_t remain = dim%commsize;
  size_t addleft = 0;
  if(rank <= (int) remain) addleft = rank;
  if(rank > (int) remain) addleft = remain;
  size_t start = rank*chunk + addleft;
  if(rank == 0) start = 0;
  return start;
}

size_t paduprop_getend(size_t dim) {
  int rank; int commsize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  size_t chunk = dim/commsize;
  size_t remain = dim%commsize;
  size_t addright = 0;
  if(rank < (int) remain) addright = rank+1;
  if(rank >= (int) remain) addright = remain;
  size_t end = (rank+1)*chunk + addright;
  if(rank == commsize - 1) end = dim;
  return end;
}

int paduprop_getrank() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}
int paduprop_getcommsize() {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}
