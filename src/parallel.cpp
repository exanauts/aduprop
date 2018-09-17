#include "parallel.hpp"

std::ofstream paduprop_sink("/dev/null");
int paduprop_init() {
#ifdef MPI_VERSION
  int ierr = MPI_Init(NULL, NULL);
#else
  int ierr = 0;
#endif

  if (paduprop_getrank() != 0) {
    // Mute standard output
    std::cout.rdbuf(paduprop_sink.rdbuf());
    // Optionally mute standard error
    std::cerr.rdbuf(paduprop_sink.rdbuf());      
  }
  return ierr;
}
int paduprop_destroy() {
#ifdef MPI_VERSION
  return MPI_Finalize();
#else
  return 0;
#endif
}
size_t paduprop_getstart(size_t dim) {
  int commsize = paduprop_getcommsize();
  int rank = paduprop_getrank();
  size_t chunk = dim/commsize;
  size_t remain = dim%commsize;
  size_t addleft = 0;
  if(rank <= (int) remain) addleft = rank;
  if(rank > (int) remain) addleft = remain;
  size_t start = rank*chunk + addleft;
  if(rank == 0) start = 0;
  return start;
}

size_t paduprop_getstart(size_t dim, int rank) {
  int commsize = paduprop_getcommsize();
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
  int commsize = paduprop_getcommsize();
  int rank = paduprop_getrank();
  size_t chunk = dim/commsize;
  size_t remain = dim%commsize;
  size_t addright = 0;
  if(rank < (int) remain) addright = rank+1;
  if(rank >= (int) remain) addright = remain;
  size_t end = (rank+1)*chunk + addright;
  if(rank == commsize - 1) end = dim;
  return end;
}

size_t paduprop_getend(size_t dim, int rank) {
  int commsize = paduprop_getcommsize();
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
#ifdef MPI_VERSION
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  rank = 0;
#endif
  return rank;
}
int paduprop_getcommsize() {
  int size;
#ifdef MPI_VERSION
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
  size = 1;
#endif
  return size;
}
int paduprop_sum(double &v) {
  int ierr = 0;
#ifdef MPI_VERSION
  ierr = MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return ierr;
}
int paduprop_sum(pMatrix<double> &mat) {
  size_t nrows = mat.nrows();
  size_t ncols = mat.ncols();
  size_t size = nrows*ncols;
  double *ptr = mat.get_datap();
  int ierr = 0;
#ifdef MPI_VERSION
  ierr = MPI_Allreduce(MPI_IN_PLACE, ptr, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return ierr;
}

int paduprop_sum(pVector<double> &vec) {
  size_t size = vec.dim();
  double *ptr = vec.get_datap();
  int ierr = 0;
#ifdef MPI_VERSION
  ierr = MPI_Allreduce(MPI_IN_PLACE, ptr, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return ierr;
}

int paduprop_gather(pTensor3<double> &ten) {
  int ierr = 0;
  size_t dim = ten.get_d1();
  size_t start = paduprop_getstart(dim);
  size_t end = paduprop_getend(dim);
  size_t size = (end-start)*dim*dim;
  double *sendbuf = new double[size];
  for(size_t i = 0; i < size; ++i) sendbuf[i] = *(&ten[0][0][start] + i);
  int *recvcounts = new int[paduprop_getcommsize()];
  int *displs = new int[paduprop_getcommsize()];
  int count = 0;
  for(int i = 0; i < paduprop_getcommsize(); ++i) {
    displs[i] = count;
    recvcounts[i] = (paduprop_getend(dim,i) - paduprop_getstart(dim,i))*dim*dim;
    count += recvcounts[i];
  }
#ifdef MPI_VERSION
  ierr = MPI_Allgatherv(sendbuf, size, MPI_DOUBLE,
                  &ten[0][0][0], recvcounts, displs, MPI_DOUBLE,
                  MPI_COMM_WORLD);
#else
  std::memcpy (sendbuf, &ten[0][0][0], recvcounts[0] * sizeof(double));
#endif
  
  delete [] recvcounts;
  delete [] displs;
  delete [] sendbuf;
  return ierr;
}

int paduprop_gather(pMatrix<double> &mat) {
  int ierr = 0;
  size_t dim = mat.nrows();
  size_t start = paduprop_getstart(dim);
  size_t end = paduprop_getend(dim);
  size_t size = (end-start)*dim;
  double *sendbuf = new double[size];
  for(size_t i = 0; i < size; ++i) sendbuf[i] = *(&mat[0][start] + i);
  int *recvcounts = new int[paduprop_getcommsize()];
  int *displs = new int[paduprop_getcommsize()];
  int count = 0;
  for(int i = 0; i < paduprop_getcommsize(); ++i) {
    displs[i] = count;
    recvcounts[i] = (paduprop_getend(dim,i) - paduprop_getstart(dim,i))*dim;
    count += recvcounts[i];
  }
#ifdef MPI_VERSION
  ierr = MPI_Allgatherv(sendbuf, size, MPI_DOUBLE,
                  &mat[0][0], recvcounts, displs, MPI_DOUBLE,
                  MPI_COMM_WORLD);
#else
  std::memcpy (sendbuf, &mat[0][0], recvcounts[0] * sizeof(double));
#endif
  delete [] recvcounts;
  delete [] displs;
  delete [] sendbuf;
  return ierr;
}

int paduprop_gather(pVector<double> &vec) {
  int ierr = 0;
  size_t dim = vec.dim();
  size_t start = paduprop_getstart(dim);
  size_t end = paduprop_getend(dim);
  size_t size = end-start;
  double *sendbuf = new double[size];
  for(size_t i = 0; i < size; ++i) sendbuf[i] = vec[start+i];
  int *recvcounts = new int[paduprop_getcommsize()];
  int *displs = new int[paduprop_getcommsize()];
  int count = 0;
  for(int i = 0; i < paduprop_getcommsize(); ++i) {
    displs[i] = count;
    recvcounts[i] = (paduprop_getend(dim,i) - paduprop_getstart(dim,i));
    count += recvcounts[i];
  }
  #ifdef MPI_VERSION
  ierr = MPI_Allgatherv(sendbuf, size, MPI_DOUBLE,
                  &vec[0], recvcounts, displs, MPI_DOUBLE,
                  MPI_COMM_WORLD);
  #else
  std::memcpy (sendbuf, &vec[0], recvcounts[0] * sizeof(double));
  #endif
  delete [] recvcounts;
  delete [] displs;
  delete [] sendbuf;
  return ierr;
}
