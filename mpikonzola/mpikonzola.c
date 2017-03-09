// mpikonzola.cpp : Defines the entry point for the console application.
//

#include <mpi.h>
#include "stdafx.h"
#include <iostream>

int main(int argc, char** argv) {

	MPI_Init(&argc, &argv);
	int a[4];
	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	for (int i = 0; i < 4; i++) {
		if (myrank == 0)a[i] = i;
		else a[i] = 0;
	}
/*
#ifdef _DEBUG
	if (myrank == 0)
	{
		system("pause");
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	
	char hostname[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(hostname, &name_len);
	*/
	std::cout << "proc " << myrank << " original a = " << a[0] << "," <<a[1] << "," << a[2] << "," << a[3] << "\n";
	MPI_Bcast(a, 4, MPI_INT, 0, MPI_COMM_WORLD);
	std::cout << "proc " << myrank << " received a = " << a[0] << "," << a[1] << "," << a[2] << "," << a[3] << "\n";

	system("pause");
	MPI_Finalize();
	
	return 0;
}