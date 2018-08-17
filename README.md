# MPI-Quicksort-Hypercubes

The goal:
Realization of quicksort algorithm with means of MPI library. The set of processes is imagined and used as a hypercube. 
(If you don't know, what is a hypercube, please google it, because it is hard to explain in a nutshell without pictures.)

Input data:
Text file with int numbers, separated with space.
(There is a sample file with 100 random int numbers in a repo. You may use it to test the program.)

NOTE:
The code is written mostly in educational purposes by a programming noob.

A simplified scheme of a program:

1. An array of numbers is divided in equal parts and spreaded among all processes.
2. The first process in this hypercube calculates pivot as a mean of numbers in his array.
3. Each process divides its part in two according to pivot. Process with first bit of binary rank equal to 0 sends numbers > pivot to its partner. Partners (first bit = 1) send smaller numbers. Thus 0* processes will have smaller numbers, 1* prosesses - greater numbers.
4. Hypercube is divided in two hypercubes. Now we won't look on the first bit, the criterion of partner-choosing will be the second bit.
5. Steps 2-4 are repeated N times, where 2^N = number of processes.
6. After N iterations the numbers of process 0 <  the numbers of process 1 <  the numbers of process 2 < ... and so on. Processes can now sort numbers sequentially with normal quicksort.
7. Results are gathered and written by the master process.

Credit:
I used ideas suggested in a textbook by V. Gergel, Theory and Practice of Parallel Calculations

How to configure your Visual Studio to use MPI:
There is a good and simple article here: https://blogs.technet.microsoft.com/windowshpc/2015/02/02/how-to-compile-and-run-a-simple-ms-mpi-program/



