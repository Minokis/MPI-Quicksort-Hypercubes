#ifndef QSORTMPI
#define QSORTMPI

#include <iostream>
#include <math.h>
#include <chrono>
#include <fstream>
#include <string>

// I don't know how many numbers there are in my input file. Let's assume there will be less than 500 000.
const int MAX_SIZE = 500000;
int ProcNum, ProcRank;
/*
			Function      read_file
			   Input      file(ifstream object), int array of MAX_SIZE
			  Output      quantity of numbers in file
			     Job      copies numbers from the file to array, returns number of working elements of the array
*/

int read_file(std::ifstream& input, int *result) {
	int i = 0;
	while (input >> result[i] && i < MAX_SIZE)
	{
		i++;
	}
	return i;
}

/*
			Function      partition
			   Input      array of int numbers, indices of sorting interval
			  Output	  index of element which divides the array
			     Job      aids to quicksort func. Pivot is picked as the last element.

*/
int partition(int *A, int lo, int hi) {
	int pivot = A[lo];
	int i = lo - 1;
	int j = hi + 1;
	while (true) {
		do
			i = i + 1;
		while (A[i] < pivot);

		do
			j = j - 1;
		while (A[j] > pivot);
		if (i >= j)
			return j;

		int temp = A[i];
		A[i] = A[j];
		A[j] = temp;
	}
}

/*
			Function      mean
			   Input      array of int numbers, indices of start and end of array
			  Output	  average of numbers from start to end
			     Job      calculates average. Used for finding a nice pivot.

*/

int mean(int * A, int lo, int hi) {
	long total = 0;
	for (int i = lo; i <= hi; i++)
		total = total + A[i];
	return (int)(total / (hi - lo + 1));
}

/*
		Function      lomuto_partition
		   Input      array of int numbers, indices of sorting interval, pivot
		  Output	  index of pivot, which divides the array
		     Job      divides the array, so first part has numbers only < = pivot; moves pivot to its final place
*/

int lomuto_partition(int * A, int lo, int hi, int pivot) {
	int i = lo - 1;
	for (int j = lo; j < hi; j++) {
		if (A[j] <= pivot) {
			i++;
			int temp = A[i];
			A[i] = A[j];
			A[j] = temp;
		}
	}
	int temp = A[i + 1];
	A[i + 1] = A[hi];
	A[hi] = temp;
	return (A[i + 1] > pivot) ? i : (i + 1);
}

/*
			Function      quicksort
			   Input      array of int numbers, indices of sorting interval
			  Output
			     Job      performs sequiential sorting
*/

void quicksort(int *A, int lo, int hi) {
	if (lo < hi) {
		int p = partition(A, lo, hi);
		quicksort(A, lo, p);
		quicksort(A, p + 1, hi);
	}
}


/*
			Function      ProcessInitialization
			   Input      ifstream object for reading, address for writing a size of subarray
			              for each process
			  Output      An equal part (subarray) of initial array for each of the processes
			     Job      Reads file, defines size of subarray, allocates memory for each subarray,
                          sends parts to processes and makes processes receive them.
*/

int* ProcessInitialization(std::ifstream& input, int &ProcDataSize) {
	if (ProcRank == 0) {

		int *Data = (int *)malloc(MAX_SIZE * sizeof(int));
		int DataSize = read_file(input, Data);

		ProcDataSize = DataSize / ProcNum;
		int *ProcData = (int *)malloc(ProcDataSize * sizeof(int));

		for (int i = 0; i < ProcDataSize; i++)
			ProcData[i] = Data[i];

		for (int rank = 1; rank < ProcNum; rank++) {
			MPI_Send(&ProcDataSize, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
			MPI_Send((Data + rank*ProcDataSize), ProcDataSize, MPI_INT, rank, 1, MPI_COMM_WORLD);
		}

		free (Data);
		return ProcData;
	}
	else {
		MPI_Recv(&ProcDataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int *ProcData = (int *)malloc(ProcDataSize * sizeof(int));
		MPI_Recv(ProcData, ProcDataSize, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		return ProcData;
	}
}

// TODO: write descriptions

void  PivotDistribution(int *pProcData, int ProcDataSize, int Dim, int Mask, int Iter, int *pPivot) {
	// УЄУЅУЊУЋУ УАУ УЖУЈУП УЃУАУГУЏУЏУЛ, УЏУЎУЄУЃУАУГУЏУЏУЛ УЈ УЏУЎУЄУЊУЎУЌУЌУ 
	MPI_Group WorldGroup;
	MPI_Group SubcubeGroup; // УЏУЎУЄУЃУЈУЏУЅУАУЊУГУЁ 
	MPI_Comm  SubcubeComm;  // УУЎУЌУЌУГУ­УЈУЊУ УВУЎУА УЏУЎУЄУЃУЈУЏУЅУАУЊУГУЁУ  


	// У УЏУЅУАУЂУЎУЉ УЈУВУЅУАУ УЖУЈУЈ (i=2) УЃУЈУЏУЅУАУЊУГУЁ = 4 УЏУАУЎУЖУ , УЏУЎУВУЎУЌ 2, УЏУЎУВУЎУЌ 1
	int  GroupNum = ProcNum / (int)pow(2, Dim - Iter);
	
	// УБУЎУЇУЄУ У­УЈУЅ УЌУ УБУБУЈУЂУ  УЏУАУЎУЖУЅУБУБУЎУАУЎУЂ УЃУЈУЏУЅУАУЊУГУЁУ 
	int *ProcRanks = (int *)malloc(GroupNum * sizeof(int));
	
	// УДУЎУАУЌУЈУАУЎУЂУ У­УЈУЅ УБУЏУЈУБУЊУ  УАУ У­УЃУЎУЂ УЏУАУЎУЖУЅУБУБУЎУЂ УЄУЋУП УЃУЈУЏУЅУАУЊУГУЁУ 
	int StartProc = ProcRank - GroupNum;
	if (StartProc < 0) StartProc = 0;

	int EndProc = ProcRank + GroupNum;
	if (EndProc > ProcNum) EndProc = ProcNum;

	int j = 0;
	for (int proc = StartProc; proc < EndProc; proc++) {
		// УУБУЋУЈ УБУАУ УЂУ­УЈУЂУ УЅУЌУЛУЉ УЁУЈУВ ("УБУВУ УАУИУЈУЉ" УЂ УЃУЈУЏУЅУАУЊУГУЁУЅ) УЎУЄУЈУ­У УЊУЎУЂУЛУЉ
		if ((ProcRank & Mask) >> (Iter) == (proc & Mask) >> (Iter)) {
			ProcRanks[j++] = proc;
		}
	}

	// УУЁУКУЅУЄУЈУ­УЅУ­УЈУЅ УЏУАУЎУЖУЅУБУБУЎУЂ УЏУЎУЄУЃУЈУЏУЅУАУЊУГУЁУ  УЂ УЎУЄУ­УГ УЃУАУГУЏУЏУГ
	MPI_Comm_group(MPI_COMM_WORLD, &WorldGroup);
	// УБУЎУЇУЄУ УЅУВ У­УЎУЂУГУО УЃУАУГУЏУЏУГ
	MPI_Group_incl(WorldGroup, GroupNum, ProcRanks, &SubcubeGroup);
	MPI_Comm_create(MPI_COMM_WORLD, SubcubeGroup, &SubcubeComm);


	// УУЎУЈУБУЊ УЈ УАУ УБУБУЛУЋУЊУ  УЂУЅУЄУГУЙУЅУЃУЎ УНУЋУЅУЌУЅУ­УВУ  УЂУБУЅУЌ УЏУАУЎУЖУЅУБУБУ УЌ УЏУЎУЄУЃУЈУЏУЅУАУЊУГУЁУ 
	if (ProcRank == ProcRanks[0]) {
		*pPivot = mean(pProcData, 0, ProcDataSize - 1);
		std::cout << "Iter is " << Iter << ", pivot is " << *pPivot << std::endl;
	}

	MPI_Bcast(pPivot, 1, MPI_INT, 0, SubcubeComm);
	MPI_Group_free(&SubcubeGroup);
	MPI_Comm_free(&SubcubeComm);
	free(ProcRanks);
}

void DataMerge(int *pMergeData, int MergeDataSize, int *pData, int DataSize, int *pRecvData, int RecvDataSize) {
	for (int i = 0; i < MergeDataSize; i++) {
		if (i < DataSize)
			pMergeData[i] = pData[i];
		else
			pMergeData[i] = pRecvData[i - DataSize];
	}
}

void ParallelHyperQuickSort(int *&ProcData, int *ProcDataSize) {
	MPI_Status status;
	int CommProcRank; // fellow process
	int *pData,     // УЗУ УБУВУМ УЁУЋУЎУЊУ , УЊУЎУВУЎУАУ УП УЎУБУВУ УЅУВУБУП У­У  УЏУАУЎУЖУЅУБУБУЎУАУЅ
		*pSendData;  // УЗУ УБУВУМ УЁУЋУЎУЊУ , УЏУЅУАУЅУЄУ УЂУ УЅУЌУ УП УЏУАУЎУЖУЅУБУБУЎУАУГ fellowProc
					 // УЗУ УБУВУМ УЁУЋУЎУЊУ , УЏУЎУЋУГУЗУ УЅУЌУ УП УЎУВ fellowProc
					 // УЁУЋУЎУЊ УЄУ У­У­УЛУЕ, УЏУЎУЋУГУЗУ УЅУЌУЛУЉ УЏУЎУБУЋУЅ УБУЋУЈУПУ­УЈУП

	int DataSize, SendDataSize, RecvDataSize, MergeDataSize;
	int HypercubeDim = (int)(log(ProcNum) / log(2)); // УАУ УЇУЌУЅУАУ­УЎУБУВУМ УЃУЈУЏУЅУАУЊУГУЁУ 
	int Mask = ProcNum;
	int Pivot;
	
	for (int i = HypercubeDim; i > 0; i--) {
		PivotDistribution(ProcData, *ProcDataSize, HypercubeDim, Mask, i, &Pivot);
		Mask = Mask >> 1;

		// УУЅУАУЂУЎУ­У УЗУ УЋУМУ­У УП УБУЎУАУВУЈУАУЎУЂУЊУ  УЏУЎ pivot + УЏУЎУЋУГУЗУЅУ­УЈУЅ УЏУЎУЇУЈУЖУЈУЈ УЏУЎУБУЋУЅУЄУ­УЅУЃУЎ УЗУЈУБУЋУ  < pivot
		int pos = lomuto_partition(ProcData, 0, *ProcDataSize - 1, Pivot);

		if ((ProcRank & Mask) == 0) {  // УБУВУ УАУИУЈУЉ УЁУЈУВ = 0
			CommProcRank = ProcRank + Mask;

			pSendData = &ProcData[pos + 1]; // УЏУЎУБУЋУЅ УЏУЎУЇУЈУЖУЈУЈ pivot, УВУЅ, УЗУВУЎ УЁУЎУЋУМУИУЅ
			SendDataSize = *ProcDataSize - pos - 1;
			if (SendDataSize < 0) SendDataSize = 0;

			pData = &ProcData[0];
			DataSize = pos + 1;
		}
		else {
			CommProcRank = ProcRank - Mask;
			pSendData = &ProcData[0];
			SendDataSize = pos + 1;
			if (SendDataSize > *ProcDataSize) SendDataSize = pos;
			pData = &ProcData[pos + 1];
			DataSize = *ProcDataSize - pos - 1;
			if (DataSize < 0) DataSize = 0;
		}

		// УУЅУАУЅУБУЛУЋУЊУ  УАУ УЇУЌУЅУАУЎУЂ УЗУ УБУВУЅУЉ УЁУЋУЎУЊУЎУЂ УЄУ У­У­УЛУЕ
		MPI_Sendrecv(&SendDataSize, 1, MPI_INT, CommProcRank, 0, &RecvDataSize, 1, MPI_INT, CommProcRank, 0, MPI_COMM_WORLD, &status);

		// УУЅУАУЅУБУЛУЋУЊУ  УЗУ УБУВУЅУЉ УЁУЋУЎУЊУЎУЂ УЄУ У­У­УЛУЕ 
		int *pRecvData = (int *)malloc(RecvDataSize * sizeof(int));
		MPI_Sendrecv(pSendData, SendDataSize, MPI_INT, CommProcRank, 0, pRecvData, RecvDataSize, MPI_INT, CommProcRank, 0, MPI_COMM_WORLD, &status);

		// УУЋУЈУПУ­УЈУЅ УЗУ УБУВУЅУЉ
		MergeDataSize = DataSize + RecvDataSize;

		int *pMergeData = (int *)malloc(MergeDataSize * sizeof(int));
		DataMerge(pMergeData, MergeDataSize, pData, DataSize, pRecvData, RecvDataSize);
		
		ProcData = pMergeData;
		*ProcDataSize = MergeDataSize;
		
		free(pRecvData);
	}
}

void ProcessTermination(std::ofstream &output, int *ProcData, int ProcDataSize) {

	MPI_Status sts;
	//quicksort(ProcData, 0, *ProcDataSize - 1);


	long DataSize;
	MPI_Reduce(&ProcDataSize, &DataSize, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	int *Data, *rcounts, *displs;

	rcounts = (int*)malloc(ProcNum*sizeof(int));
	Data = (int*)malloc(DataSize*sizeof(int));
	displs = (int*)malloc(ProcNum*sizeof(int));

	MPI_Gather(&ProcDataSize, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (ProcRank == 0) {
		int acc = 0;
		for (int r = 0; r < ProcNum; r++) {
			displs[r] = acc;
			acc += rcounts[r];
		}
	}
	MPI_Gatherv(ProcData, ProcDataSize, MPI_INT, Data, rcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
	if (ProcRank == 0) {
		for (int i = 0; i < DataSize; i++)
			output << Data[i] << " ";
		output.close();
	}

	// УЗУЅУАУЅУЇ 0?
	/*if (ProcRank != 0) {
		MPI_Send(&ProcDataSize, 1, MPI_INT, 0, 26, MPI_COMM_WORLD);	
	} 
	else {
		Data = (int*)malloc(DataSize*sizeof(int));
		rcounts = (int*)malloc(ProcNum*sizeof(int));
		displs = (int*)malloc(ProcNum*sizeof(int));

		for (int r = 0; r < ProcNum; r++) {
			MPI_Recv((rcounts + r), 1, MPI_INT, r, 26, MPI_COMM_WORLD, &sts);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		int acc = 0;
		for (int r = 0; r < ProcNum; r++) {
			displs[r] = acc;
			acc += rcounts[r];
		}
		std::cout << rcounts[0] << " " << rcounts[1] << " " << rcounts[3] << std::endl;
		std::cout << displs[0] << " " << displs[1] << " " << displs[3] << std::endl;
		MPI_Gatherv(ProcData, ProcDataSize, MPI_INT, Data, rcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

		for (int i = 0; i < DataSize; i++)
			output << Data[i] << " ";
		output.close();
	}
	
*/

	//MPI_Sendrecv(&ProcDataSize, 1, MPI_INT, 0, 25, (rcounts+ProcRank), 1, MPI_INT, ProcRank, 25, MPI_COMM_WORLD, &sts);
	
	//MPI_Barrier(MPI_COMM_WORLD);

	free(ProcData);
}


#endif