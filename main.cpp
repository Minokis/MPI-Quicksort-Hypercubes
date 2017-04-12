#include "stdafx.h"
#include "mpi.h"
#include "Header.h"


int read_file(std::ifstream& input, int *result);
int partition(int *A, int lo, int hi);
void quicksort(int *A, int lo, int hi);
int* ProcessInitialization(std::ifstream& input, int &ProcDataSize); 
int lomuto_partition(int * A, int lo, int hi, int pivot);
void PivotDistribution(int *pProcData, int ProcDataSize, int Dim, int Mask, int Iter, int *pPivot);
void DataMerge(int *pMergeData, int MergeDataSize, int *pData, int DataSize, int *pRecvData, int RecvDataSize);
void ParallelHyperQuickSort(int *&ProcData, int *ProcDataSize);
void ProcessTermination(std::ofstream &output, int *ProcData, int ProcDataSize);


int main(int argc, char *argv[]) {
	if (argc != 3) { std::cout << "Usage: enter input and output filenames to perform sorting.\n"; return 1; }
	std::ifstream input(argv[1]);
	if (!input) { std::cout << "Cannot open input file.\n"; return 1; }
	auto start = std::chrono::high_resolution_clock::now();

	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	
	// Reading and sending parts to every process into ProcData
	int ProcDataSize;
	int *ProcData = ProcessInitialization(input, ProcDataSize);
	
	// Sorting
	ParallelHyperQuickSort(ProcData, &ProcDataSize);

	auto finish = std::chrono::high_resolution_clock::now();


	std::string outputName = argv[2];
	std::ofstream output(outputName);
	if (!output) { std::cout << "Cannot open output file\n"; return 1; }

	quicksort(ProcData, 0, ProcDataSize - 1);
	
	ProcessTermination(output, ProcData, ProcDataSize);
	
	MPI_Finalize();
	return 0;
}