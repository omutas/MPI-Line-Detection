/* Student Name: Mahmut Uğur TAŞ
* Student Number: 2011104090
* Compile Status: Compiling
* Program Status: Working
* Notes: Project gives results same as latest output 
* I didn't check whether given processor number divide 200 or not
*/
#include <stdio.h>
#include "mpi.h"
//Reads text file and fills array 
void readMyFile(int arraySize, int myArray[arraySize], char **argv, int processors){
	//open input file
	FILE* file = fopen (argv[1], "r");
	//file entry
	int i = 0;
	//position 
	int j;
	//fill master processor's local array with 0
	for(j=0;j<40000/(processors-1);j++){
    	myArray[j]=0;
	}

	fscanf (file, "%d", &i);  
	while (!feof (file))
	{
		myArray[j] = i;
		j++;
		fscanf (file, "%d", &i);      
	}
	fclose (file);
}
    /////////////////////////************************////////////////////////
    //                          SMOOTHING OPERATION                        //
    /////////////////////////************************////////////////////////
//Multiplies array elements with 1/9 then returns their sum for smoothing operation
int mySmooth(int array1[3][3]){
	int i,j;
	double sum = 0;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			sum = sum + array1[i][j];
		}
	}
	sum = sum*1.0/9.0;
	return sum;
}
//takes 2 array and gives its result
int mulArray(int array1[3][3], int array2[3][3]){
	int i,j;
	int sum = 0;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			sum = sum + array1[i][j]*array2[i][j];
		}
	}
	return sum;
}
    /////////////////////////************************////////////////////////
    //                        CONVOLUTION OPERATION                        //
    /////////////////////////************************////////////////////////
//
int convol(int array1[3][3], int threshold){
	int hArray[3][3] = {  
	   {-1, -1, -1} ,
	   { 2,  2,  2} ,   
	   {-1, -1, -1}   
	};
	int vArray[3][3] = {  
	   {-1,  2, -1} ,
	   {-1,  2, -1} ,   
	   {-1,  2, -1}  
	};
	int cOblArray[3][3] = {  
	   {-1, -1,  2} ,
	   {-1,  2, -1} ,   
	   { 2, -1, -1}  
	};
	int dOblArray[3][3] = {  
	   { 2, -1, -1} ,
	   {-1,  2, -1} ,   
	   {-1, -1,  2}  
	};
	int hValue = mulArray(array1,hArray);
	int vValue = mulArray(array1,vArray);
	int cOblValue = mulArray(array1,cOblArray);
	int dOblValue = mulArray(array1,dOblArray);
	if (threshold<hValue)
	{
		return 255;
	}else if (threshold<vValue)
	{
		return 255;
	}else if (threshold<cOblValue)
	{
		return 255;
	}else if (threshold<dOblValue)
	{
		return 255;
	}
	return 0;
}
    /////////////////////////************************////////////////////////
    //                   GET 3X3 MATRIX FROM GIVEN MATRIX                  //
    /////////////////////////************************////////////////////////
//Turn back 3x3 matrix from given array and given position
void getArray(int rowSize, int colSize, int myArray[rowSize][colSize], int position1,int position2, int newArray[3][3]){
	int i,j;
	for (i = position1-1; i < position1 + 2; i++)
	{
		for (j = position2-1; j < position2 + 2; j++)
		{
			newArray[i-position1+1][j-position2+1] = myArray[i][j];
		}
	}
}

int main(int argc, char* argv[])
{
	//Some values to use later 
    int rank, size, i, j;
	MPI_Status status;

    MPI_Init(&argc, &argv); //Initialize the MPI execution environment 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Determines the rank of the calling process in the communicator 
    MPI_Comm_size(MPI_COMM_WORLD, &size); //Determines the size of the group associated with a communicator
    //Threshold
    int threshold = atoi(argv[3]);
    //takes integer array come from previous processor
    int beforeArray[200];
    //takes integer array come from next processor
    int afterArray[200];
    //This array keeps all data which is writen in input file
    int globalArray[40000+40000/(size-1)];
    //This array keeps processors' data 
	int localArray[40000/(size-1)];
	//This matrix keeps processors' data and wait for before and after arrays
	int localMatrix[200/(size-1)+2][200];
	//This matrix keeps smoothed version of local matrix
	int localSmoothMatrix[200/(size-1)][198];
	//This matrix keeps smoothed version of local matrix with additions
	int localCompleteSmoothMatrix[200/(size-1)+2][198];
	//This array keeps convolution results
	int localConvolArray[200/(size-1)*196];
	//This array keeps convolution results at master
	int masterConvolArray[39200];
    /////////////////////////************************////////////////////////
    //                                READ FILE                            //
    /////////////////////////************************////////////////////////
	//read file with master processor
    if (rank == 0)
    {
    	readMyFile(40000+40000/(size-1),globalArray, argv, size);
    }
    //This barrier can be remove later
    MPI_Barrier(MPI_COMM_WORLD);
    /////////////////////////************************////////////////////////
    //                           DISTRIBUTE DATA                           //
    /////////////////////////************************////////////////////////
    MPI_Scatter(globalArray, 40000/(size-1), MPI_INT, localArray, 40000/(size-1), MPI_INT, 0, MPI_COMM_WORLD);
    //convert received array to matrix	
    for (i = 0; i < 40000/(size-1); i++)
    {
    	localMatrix[(i/200)+1][i%200]=localArray[i];
    }
    /////////////////////////************************////////////////////////
    //                       PASS DATA BETWEEN SLAVES                      //
    /////////////////////////************************////////////////////////
    if (rank%2==0 && rank != 0 && rank!=size-1)
    {
    	//send data to before and after processors
    	MPI_Send(localArray, 200, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
    	//send last row's data to next processor
    	MPI_Send(&localArray[(200/(size-1)-1)*200], 200, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    	//receive data from both
    	MPI_Recv(beforeArray, 200, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
    	MPI_Recv(afterArray, 200, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &status);
    }else if(rank%2==0 && rank==size-1){
    	//send just before
    	MPI_Send(localArray, 200, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
    	//receive just before
    	MPI_Recv(beforeArray, 200, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
    }else if(rank%2==1 && rank != 1 && rank!=size-1){
    	//receive data from both
    	MPI_Recv(beforeArray, 200, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
    	MPI_Recv(afterArray, 200, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &status);
    	//send data to before and after processors
    	MPI_Send(localArray, 200, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
    	//send last row's data to next processor
    	MPI_Send(&localArray[(200/(size-1)-1)*200], 200, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }else if (rank==1)
    {
    	//receive just after
    	MPI_Recv(afterArray, 200, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &status);
    	//send last row's data to next processor
    	MPI_Send(&localArray[(200/(size-1)-1)*200], 200, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }else if (rank%2==1 && rank==size-1)
    {
    	//receive just before
    	MPI_Recv(beforeArray, 200, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
    	//send data to before
    	MPI_Send(localArray, 200, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
    }
    ///////////////////////////************************//////////////////////////
    //before and after arrays received data we need to add them to local matrix//
    ///////////////////////////************************//////////////////////////
    if (rank == 1)
    {
    	//this just have after array
    	for (i = 0; i < 200; i++)
    	{
    		localMatrix[200/(size-1)+1][i]=afterArray[i];
    	}
    }else if (rank == size-1)
    {
    	//this just have before array
    	for (i = 0; i < 200; i++)
    	{
    		localMatrix[0][i]=beforeArray[i];
    	}
    }else if (rank!=0)
    {
    	//rest has both array
    	//this for before array
    	for (i = 0; i < 200; i++)
    	{
    		localMatrix[0][i]=beforeArray[i];
    	}
    	//this for after array
    	for (i = 0; i < 200; i++)
    	{
    		localMatrix[200/(size-1)+1][i]=afterArray[i];
    	}
    }
    MPI_Barrier(MPI_COMM_WORLD); //Blocks until all processes in the communicator have reached this routine.
    /////////////////////////************************////////////////////////
    //                         SMOOTHING OPERATION                         //
    /////////////////////////************************////////////////////////
    int array2[3][3];
    if (rank == 1)//If processor is first slave we don't smooth its first line
    {
    	for (i = 2; i < 200/(size-1)+1; i++)
		{
			for (j = 1; j < 199; j++)
			{
				getArray(200/(size-1)+2,200,localMatrix,i,j,array2);
				localSmoothMatrix[i-1][j-1] = mySmooth(array2);
			}
		}
    }else if ( rank == size-1)//If processor is last slave we don't smooth its last line
    {
    	for (i = 1; i < 200/(size-1); i++)
		{
			for (j = 1; j < 199; j++)
			{
				getArray(200/(size-1)+2,200,localMatrix,i,j,array2);
				localSmoothMatrix[i-1][j-1] = mySmooth(array2);
			}
		}
    }else if ( rank != 0)//Smooth all lines of these processors' matrix
    {
    	for (i = 1; i < 200/(size-1)+1; i++)
		{
			for (j = 1; j < 199; j++)
			{
				getArray(200/(size-1)+2,200,localMatrix,i,j,array2);
				localSmoothMatrix[i-1][j-1] = mySmooth(array2);
			}
		}
    }
    MPI_Barrier(MPI_COMM_WORLD); //Wait until all slaves complete smoothing operation
    /////////////////////////************************////////////////////////
    //           PASS REQUIRED SMOOTHED DATA FROM SLAVE TO SLAVE           //
    /////////////////////////************************////////////////////////
    if (rank%2==0 && rank != 0 && rank!=size-1)
    {
    	//send data to before and after processors
    	MPI_Send(localSmoothMatrix, 198, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
    	//send last row's data to next processor
    	MPI_Send(&localSmoothMatrix[(200/(size-1)-1)][0], 198, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    	//receive data from both
    	MPI_Recv(beforeArray, 198, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
    	MPI_Recv(afterArray, 198, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &status);
    }else if(rank%2==0 && rank==size-1){
    	//send just before
    	MPI_Send(localSmoothMatrix, 198, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
    	//receive just before
    	MPI_Recv(beforeArray, 198, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
    }else if(rank%2==1 && rank != 1 && rank!=size-1){
    	//receive data from both
    	MPI_Recv(beforeArray, 198, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
    	MPI_Recv(afterArray, 198, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &status);
    	//send data to before and after processors
    	MPI_Send(localSmoothMatrix, 198, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
    	//send last row's data to next processor
    	MPI_Send(&localSmoothMatrix[(200/(size-1)-1)][0], 198, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }else if (rank==1)
    {
    	//receive just after
    	MPI_Recv(afterArray, 198, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &status);
    	//send last row's data to next processor
    	MPI_Send(&localSmoothMatrix[(200/(size-1)-1)][0], 198, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }else if (rank%2==1 && rank==size-1)
    {
    	//receive just before
    	MPI_Recv(beforeArray, 198, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
    	//send data to before
    	MPI_Send(localSmoothMatrix, 198, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
    }
    /////////////////////////************************////////////////////////
    //                ADD TAKEN DATA TO COMPLETE SMOOTHED MATRIX           //
    /////////////////////////************************////////////////////////
    if (rank == 1)
    {
    	//First line of first slave's matrix should be filled with 0 it will not be used
    	for (i = 0; i < 199; i++)
    	{
    		localCompleteSmoothMatrix[0][i]= 0;
    	}
    	//this just have after array
    	for (j = 1; j < 200/(size-1)+1; j++){	
	    	for (i = 0; i < 199; i++)
	    	{
	    		localCompleteSmoothMatrix[j][i]=localSmoothMatrix[j-1][i];
	    	}
	    }
	    for (i = 0; i < 199; i++)
    	{
    		localCompleteSmoothMatrix[200/(size-1)+1][i]= afterArray[i];
    	}
    }else if (rank == size-1)
    {
    	//This one just have before array
    	for (i = 0; i < 199; i++)
    	{
    		localCompleteSmoothMatrix[0][i]= beforeArray[i];
    	}
    	for (j = 1; j < 200/(size-1)+1; j++){	
	    	for (i = 0; i < 199; i++)
	    	{
	    		localCompleteSmoothMatrix[j][i]=localSmoothMatrix[j-1][i];
	    	}
	    }
	    //Last line of last slave's matrix should be filled with 0 it will not be used
	    for (i = 0; i < 199; i++)
    	{
    		localCompleteSmoothMatrix[200/(size-1)+1][i]= 0;
    	}
    }else if (rank!=0)
    {
    	for (i = 0; i < 199; i++)
    	{
    		localCompleteSmoothMatrix[0][i]= beforeArray[i];
    	}
    	for (j = 1; j < 200/(size-1)+1; j++){	
	    	for (i = 0; i < 199; i++)
	    	{
	    		localCompleteSmoothMatrix[j][i]=localSmoothMatrix[j-1][i];
	    	}
	    }
	    for (i = 0; i < 199; i++)
    	{
    		localCompleteSmoothMatrix[200/(size-1)+1][i]= afterArray[i];
    	}
    }
    /////////////////////////************************////////////////////////
    //                        CONVOLUTION OPERATION                        //
    /////////////////////////************************////////////////////////
    if (rank == 1 && size != 101)//First two lines of first slave's matrix will not be evaluated 
    {
    	for (i = 3; i < 200/(size-1)+1; i++)
		{
			for (j = 1; j < 197; j++)
			{
				getArray(200/(size-1)+2,198,localCompleteSmoothMatrix,i,j,array2);
				localConvolArray[(i-1)*196+j-1] = convol(array2, threshold);
			}
		}
    }else if ( rank == size-1 && size != 101)//Last two lines of last slave's matrix will not be evaluated
    {
    	for (i = 1; i < 200/(size-1); i++)
		{
			for (j = 1; j < 197; j++)
			{
				getArray(200/(size-1)+2,198,localCompleteSmoothMatrix,i,j,array2);
				localConvolArray[(i-1)*196+j-1] = convol(array2, threshold);
			}
		}
    }else if ( rank != 0 && rank!=1 && rank!=size-1)
    {
    	for (i = 1; i < 200/(size-1)+1; i++)
		{
			for (j = 1; j < 197; j++)
			{
				getArray(200/(size-1)+2,198,localCompleteSmoothMatrix,i,j,array2);
				localConvolArray[(i-1)*196+j-1] = convol(array2, threshold);
			}
		}
    }
    /////////////////////////************************////////////////////////
    //                GATHER FINALIZED ARRAYS FROM SLAVES                  //
    /////////////////////////************************////////////////////////    
    //If processor is slave send its data to master
    if (rank != 0)
    {
    	MPI_Send(localConvolArray, 200/(size-1)*196, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    //If processor is master receive data from slaves and print it to given file
    if (rank==0)
    {
    	for (i = 1; i < size; i++)
    	{
    		MPI_Recv(&masterConvolArray[200/(size-1)*196*(i-1)], 200/(size-1)*196, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
    	}
    	FILE *outFile = fopen(argv[2], "w");
    	for (i = 2; i < 198; i++)
		{
			for (j = 0; j < 196; j++)
			{
				fprintf(outFile, "%d ", masterConvolArray[i*196+j]);
			}
			fprintf(outFile, "\n");
		}
		fclose(outFile);
    }

    MPI_Barrier(MPI_COMM_WORLD); //Blocks until all processes in the communicator have reached this routine. 
    MPI_Finalize();//Terminates MPI execution environment 

    return 0;
}