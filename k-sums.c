#include <stdio.h>
#include <stdlib.h>
#include "FileHandling.h"
#include <math.h>
#include <string.h>
#include <time.h>
#include "k-sumsH.h"


int main(int argc, char *argv[])
{

    int dim;         // Dimensions of Elements
    int n;           // n Elements
    int i,j,d;           // i counter for n, j counter for Neighborhood Elements, d counter for dimensions
    clock_t start, end;  // Variables for counting execution time
    char filename[40];
    int kPoints;          // Minimum points for a cluster to be created
    int h,k;              // h counter for every other Element, k counter for files
    int newN;        // Holder of new elements
    int oldN;           //temp Holder
    int newKp;      //new K-points
    int choise4 = 1;   //Switch var
    int choise3 = 1;   //Switch var
    int choise2 = 1;   //Switch var
    int choise = 1;    //Switch var
    FILE* AddedFileMerge; //File var for merging
    FILE* OrderFileMerge;  //File var for merging
    char RemovedFile[40];
    sprintf(RemovedFile,"NewlyAdded.txt");
    FILE* AddedCheck;  //File var for checking
    char Data2[40];
    int RecentlyAdded; //temp Holder
    /* -------------------------------------------------------------------------- */
    /*Initializing values */
    if(dim != 0) dim = 0;
    if(n != 0) n = 0;
    if(newN != 0) newN = 0;
    if(newKp != 0) newKp = 0;
    if(RecentlyAdded != 0) RecentlyAdded = 0;
    /* -------------------------------------------------------------------------- */
    // Reading Dataset, counting n elements, assigning to X Array

    FILE* Dataset;

    printf("\n Give the DataSet file:");
    scanf("%s", filename);
    printf("\n");

    Dataset = fopen(filename, "r");
    if (!Dataset)
    {
        printf("\n There is something wrong with the Dataset file! \n\n");
        return -1;
    }
    /* -------------------------------------------------------------------------- */
    dim = getColumns(Dataset);

    n = getRows(Dataset);      // Getting the Elements of the Dataset
    rewind(Dataset);

    printf("\n Elements:%d", n-1);
    printf("\n Dimensions:%d\n\n", dim);

    n--;

    printf("Give the amount of K Points: ");   // kPoints
    scanf("%d",&kPoints );
    printf("\n");
    /* -------------------------------------------------------------------------- */
// All the necessary memory allocation

    float *X;  // Array of Elements
    X = (float*)calloc(n*dim, sizeof(float));



    float *distance;   // Array for holding Distances for each ELement with each other Element
    distance = (float*)calloc(n*n, sizeof(float));


    float *kDistance;   // Array for holding k-Distances for the First Core point of a cluster
    kDistance = (float*)calloc(n, sizeof(float));


    float *SumOfDistances;   // Array for holding Distances for each ELement with each other Element
    SumOfDistances =(float*) calloc(n, sizeof(float));


    unsigned int *Neighborhood; //Array for holding which element belong to which Neighborhood
    Neighborhood =(int*) calloc(n*n, sizeof(int));


    unsigned int *NeighborhoodSize; //Array for holding the size of each Neighborhood
    NeighborhoodSize = (int*)calloc(n,sizeof(int));

    float *Mean;   // Array for holding Distances for each ELement with each other Element
    Mean =(float*)calloc(n, sizeof(float));

    float *tempDistance;   // Array for holding Distances for each ELement with each other Element
    tempDistance =(float*) calloc(n*n, sizeof(float));

    float *OrderedList;
    OrderedList =(float*) calloc(n*dim,sizeof(float));


    float *tmp2;
    tmp2 =(float*) calloc(n*dim,sizeof(float));



    /* -------------------------------------------------------------------------- */
    // Passing elements to Array X[n][dim]

    for(unsigned int i = 0; i < n; ++i)
    {
        for(unsigned int d = 0; d < dim; ++d)
            fscanf(Dataset,"%f, ",&X[i*dim + d]);
    }





    for(unsigned int i = 0; i < n; ++i)
    {
        for(unsigned int d = 0; d < dim; ++d)
        {
            OrderedList[i*dim + d] = X[i*dim + d];
        }
    }
    fclose(Dataset);



    /* -------------------------------------------------------------------------- */

    start = clock();
    /*-----------------------Starting Algorithm--------------------------------------*/
    /* -------------------------------------------------------------------------- */
    //STEP 1

    /*T TIME COST EFFECTIVE STEP */

    distance = GetDistance(X,0,n,dim,distance);
    tempDistance = GetEqual(distance,0,n,tempDistance);
    tempDistance = qSort2D(tempDistance,n,0);

    // Finding the k-Distance of each element in the dataset
    kDistance = GetK_Distance(tempDistance,0,n,kPoints,kDistance);

    /* -------------------------------------------------------------------------- */

    //STEP 2

    //Finding the Neighborhood of each element
    Neighborhood = GetNeighborhood(0,n,kDistance,distance,Neighborhood);

    /* -------------------------------------------------------------------------- */
    //STEP 3

    //Getting the Neighborhood size of each element
    NeighborhoodSize = GetNeighborhoodSize(Neighborhood,0,n,NeighborhoodSize);
    /* -------------------------------------------------------------------------- */
    //STEP 4
    //Calculating the sum of Neighborhood Elements distances for each element

    SumOfDistances = GetNeighborhoodSum(distance,0,n,Neighborhood,SumOfDistances);

    /* -------------------------------------------------------------------------- */
    //STEP 5


//Calculating Mean for each element's Neighborhood

    Mean = GetMean(SumOfDistances,NeighborhoodSize,0,n,Mean);




    /* -------------------------------------------------------------------------- */

//STEP 6
//Ordering the data according to Mean from max to min
    for(unsigned int i = 0; i < n; ++i)
    {
        for(unsigned int h = 0; h < n; ++h)
        {
            if(Mean[h] <= Mean[i])
            {
                double tmp = Mean[i];
                Mean[i] = Mean[h];
                Mean[h] = tmp;



                for(unsigned int d = 0; d < dim; ++d)
                {
                    tmp2[i*dim + d] = OrderedList[i*dim + d];
                    OrderedList[i*dim + d] = OrderedList[h*dim + d];
                    OrderedList[h*dim + d] = tmp2[i*dim + d];
                }
            }
        }
    }


    /*---------------------------------------------------------------------------*/
    //End of Algorithm
    end = clock();
    /*---------------------------------------------------------------------------*/
    double total_time = ((double) (end - start)) / CLOCKS_PER_SEC;

    FILE* OrderFile;
    OrderFile = fopen("OrderFile.txt","w");



    for(unsigned int i = 0; i < n; ++i)
    {
        fprintf(OrderFile, "%0.4f ",1/Mean[i]);
        fprintf(OrderFile, "%0.4f ",Mean[i]);
        for(unsigned int d = 0; d < dim; ++d)
        {
            fprintf(OrderFile,"%0.4f ",OrderedList[i*dim + d] );
        }
        fprintf(OrderFile, "\n");
    }

    fclose(OrderFile);



    FILE* K_File;
    char kname[40];
    sprintf(kname,"output%d.txt",kPoints);
    K_File = fopen(kname,"w");

    for(unsigned int i = 0; i < n; ++i)
    {
        fprintf(K_File, "%0.4f ",1/Mean[i]);
        fprintf(K_File, "%0.4f ",Mean[i]);
        for(unsigned int d = 0; d < dim; ++d)
        {
            fprintf(K_File,"%0.4f ",OrderedList[i*dim + d] );
        }
        fprintf(K_File, "\n");
    }

    fclose(K_File);

    printf("\n");
    printf("\n Time of Algorithm Execution: %lf \n\n",total_time);


//Giving oldN the current amout of Elements
    oldN = n;
    //Starting the menu
    while(choise != 0)
    {
        printf("\n");
        printf("1 : Add new input or change K-Points\n");
        printf("2 : Rerun with the most recent Data and K-Points\n");
        printf("3 : Merge NewlyAdded with OrderFile\n");
        printf("0 : Exit\n");
        printf("Choise : ");
        scanf("%d",&choise);
        printf("\n");

//switch number 1
        switch (choise) {
        case 0:
            /*---------------------------------------------------------------------------*/
//Exit Menu

            if((AddedCheck = fopen("NewlyAdded.txt","r")) != NULL)
            {

                while(choise4 != 0)
                {
                    printf("\n NewlyAdded.txt has been detected!\nDo you want to be Merged with OrderFile or Rerun?\n");
                    printf("1 : Merge and Exit!\n");
                    printf("2 : Rerun and Exit!\n");
                    printf("0 : Just Exit(It will get deleted)!\n");
                    printf("Choise : ");
                    scanf("%d",&choise4);
                    printf("\n");

                    switch (choise4) {
                    case 0:
                        remove(RemovedFile);

                        break;
                    case 1:

                        AddedFileMerge = fopen("NewlyAdded.txt","r");
                        OrderFileMerge = fopen("OrderFile.txt","r");

                        unsigned int dimLocal;
                        unsigned int nAdd;
                        unsigned int nOrder;
                        unsigned int nTotal;
                        dimLocal = getColumns(OrderFileMerge);
                        rewind(OrderFileMerge);
                        nAdd = getRows(AddedFileMerge);
                        rewind(AddedFileMerge);
                        nAdd--;
                        nOrder = getRows(OrderFileMerge);
                        rewind(OrderFileMerge);
                        nOrder--;
                        nTotal =  nAdd + nOrder;

                        float *Merged;
                        Merged = (float *)calloc(nTotal*dimLocal,sizeof(float));

                        float *tmp4;
                        tmp4 =(float*) calloc(nTotal*dimLocal,sizeof(float));

                        start = clock();

                        for(unsigned int i = 0; i < nAdd; ++i)
                        {
                            for(unsigned int d = 0; d < dimLocal; ++d)
                            {
                                fscanf(AddedFileMerge,"%f ",&Merged[i*dimLocal + d]);
                            }
                        }

                        for(unsigned int i = nAdd; i < nTotal; ++i)
                        {
                            for(unsigned int d = 0; d < dimLocal; ++d)
                            {
                                fscanf(OrderFileMerge,"%f ",&Merged[i*dimLocal + d]);
                            }
                        }





                        for(unsigned int i = 0; i < nTotal; ++i)
                        {
                            for(unsigned int h = 0; h < nTotal; ++h)
                            {
                                if(Merged[h*dimLocal + 1] <= Merged[i*dimLocal + 1])
                                {

                                    for(unsigned int d = 0; d < dimLocal; ++d)
                                    {
                                        tmp4[i*dimLocal + d] = Merged[i*dimLocal + d];
                                        Merged[i*dimLocal + d] = Merged[h*dimLocal + d];
                                        Merged[h*dimLocal + d] = tmp4[i*dimLocal + d];
                                    }
                                }
                            }
                        }




                        OrderFile = fopen("OrderFile.txt","w");

                        for(unsigned int i = 0; i < nTotal; ++i)
                        {
                            for(unsigned int d = 0; d < dimLocal; ++d)
                            {
                                fprintf(OrderFile,"%0.4f ",Merged[i*dimLocal + d] );
                            }
                            fprintf(OrderFile, "\n");
                        }

                        fclose(OrderFile);

                        end = clock();
                        /*---------------------------------------------------------------------------*/
                        double total_time = ((double) (end - start)) / CLOCKS_PER_SEC;
                        printf("\n");
                        printf("\n Time of Merging: %lf \n\n",total_time);


                        remove(RemovedFile);
                        fclose(OrderFileMerge);
                        fclose(AddedFileMerge);
                        free(Merged);
                        free(tmp4);
                        choise4 = 0;
                        break;

                    case 2:
                        start = clock();






                        if(newKp == 0)
                            newKp = kPoints;

                        if(newN == 0)
                            newN = n;

                        for(unsigned int i = 0; i < newN; ++i)
                        {
                            for(unsigned int d = 0; d < dim; ++d)
                                OrderedList[i*dim + d] = X[i*dim + d];
                        }





                        distance = GetDistance(X,0,newN,dim,distance);
                        tempDistance = GetEqual(distance,0,newN,tempDistance);
                        tempDistance = qSort2D(tempDistance,newN,0);

                        // Finding the k-Distance of each element in the dataset
                        kDistance = GetK_Distance(tempDistance,0,newN,newKp,kDistance);

                        /* -------------------------------------------------------------------------- */

                        //STEP 2

                        //Finding the Neighborhood of each element
                        Neighborhood = GetNeighborhood(0,newN,kDistance,distance,Neighborhood);

                        /* -------------------------------------------------------------------------- */
                        //STEP 3

                        //Getting the Neighborhood size of each element
                        NeighborhoodSize = GetNeighborhoodSize(Neighborhood,0,newN,NeighborhoodSize);
                        /* -------------------------------------------------------------------------- */
                        //STEP 4
                        //Calculating the sum of Neighborhood Elements distances for each element

                        SumOfDistances = GetNeighborhoodSum(distance,0,newN,Neighborhood,SumOfDistances);

                        /* -------------------------------------------------------------------------- */
                        //STEP 5


                        //Calculating Mean for each element's Neighborhood

                        Mean = GetMean(SumOfDistances,NeighborhoodSize,0,newN,Mean);






                        /* -------------------------------------------------------------------------- */

                        //STEP 6
                        //Ordering the data according to Mean from max to min

                        for(unsigned int i = 0; i < newN; ++i)
                        {
                            for(unsigned int h = 0; h < newN; ++h)
                            {
                                if(Mean[h] <= Mean[i])
                                {
                                    double tmp = Mean[i];
                                    Mean[i] = Mean[h];
                                    Mean[h] = tmp;

                                    for(unsigned int d = 0; d < dim; ++d)
                                    {
                                        tmp2[i*dim + d] = OrderedList[i*dim + d];
                                        OrderedList[i*dim + d] = OrderedList[h*dim + d];
                                        OrderedList[h*dim + d] = tmp2[i*dim + d];
                                    }
                                }
                            }
                        }


                        /*---------------------------------------------------------------------------*/
                        //End of Algorithm
                        end = clock();
                        /*---------------------------------------------------------------------------*/
                        double total_time2 = ((double) (end - start)) / CLOCKS_PER_SEC;



                        OrderFile = fopen("OrderFile.txt","w");



                        for(unsigned int i = 0; i < newN; ++i)
                        {
                            fprintf(OrderFile, "%0.4f ",1/Mean[i]);
                            fprintf(OrderFile, "%0.4f ",Mean[i]);
                            for(unsigned int d = 0; d < dim; ++d)
                            {
                                fprintf(OrderFile,"%0.4f ",OrderedList[i*dim + d] );
                            }
                            fprintf(OrderFile, "\n");
                        }

                        fclose(OrderFile);

                        sprintf(kname,"output%d.txt",newKp);
                        K_File = fopen(kname,"w");

                        for(unsigned int i = 0; i < newN; ++i)
                        {
                            fprintf(K_File, "%0.4f ",1/Mean[i]);
                            fprintf(K_File, "%0.4f ",Mean[i]);
                            for(unsigned int d = 0; d < dim; ++d)
                            {
                                fprintf(K_File,"%0.4f ",OrderedList[i*dim + d] );
                            }
                            fprintf(K_File, "\n");
                        }

                        fclose(K_File);


                        remove(RemovedFile);


                        printf("\n");
                        printf("\n Time of Algorithm Execution: %lf \n\n",total_time2);
                        choise4 = 0;
                        break;

                    default:
                        printf("\n Invalid Option(4) \n\n");
                    }

                }


            }

            printf("\nExited Successfuly\n");
            printf("\n");
            /*---------------------------------------------------------------------------*/
            break;
        /*---------------------------------------------------------------------------*/
        case 1:
            /*---------------------------------------------------------------------------*/
            while(choise2 != 0)
            {
                printf("1 : Update Current Data File \n");
                printf("2 : Add New Data File\n");
                printf("3 : Change K\n");
                printf("0 : Return\n");
                printf("Choise : ");
                scanf("%d",&choise2);
                printf("\n");

                switch (choise2)
                {
                case 0:
                    /*---------------------------------------------------------------------------*/
                    printf("\nReturned Successfuly\n" );
                    printf("\n");
                    /*---------------------------------------------------------------------------*/
                    break;
                /*---------------------------------------------------------------------------*/
                case 1:

                    /*---------------------------------------------------------------------------*/
                    if(newKp == 0)
                        newKp = kPoints;

                    if(newN != 0)
                        oldN = newN;

                    unsigned int NowAdded = 0;

                    if(RecentlyAdded == 0)
                    {
                        RecentlyAdded = n;
                    }

                    unsigned int localRows = 0;
                    int w;

                    Dataset = fopen(filename, "r");
                    localRows = getRows(Dataset);      // Getting the Elements of the Dataset
                    rewind(Dataset);
                    localRows--;

                    NowAdded = localRows - RecentlyAdded;
                    RecentlyAdded = localRows;
                    newN = NowAdded + oldN;

                    printf("\n Elements:%d \n", newN);

                    printf("\n Elements Added : %d\n",NowAdded);

                    if(NowAdded != 0)
                    {
                        /*----------------------------------------------------------------------------*/
//Memory Reallocation
                        X = realloc(X,newN*sizeof(*X)*dim);



                        distance = realloc(distance,newN*newN*sizeof(*distance));

                        kDistance = realloc(kDistance,newN*sizeof(*kDistance));

                        SumOfDistances = realloc(SumOfDistances,newN*sizeof(*SumOfDistances));

                        Neighborhood = realloc(Neighborhood,newN*newN*sizeof(*Neighborhood));

                        NeighborhoodSize = realloc(NeighborhoodSize,newN*sizeof(*NeighborhoodSize));

                        Mean = realloc(Mean,newN*sizeof(*Mean));

                        tempDistance = realloc(tempDistance,newN*newN*sizeof(*tempDistance));

                        OrderedList = realloc(OrderedList,newN*sizeof(*OrderedList)*dim);

                        tmp2 = realloc(tmp2,newN*dim*sizeof(*tmp2));

                        int *Added;
                        Added = (int *)calloc(newN,sizeof(int));

                        float *tempX;
                        tempX = (float *)calloc(localRows*dim,sizeof(float));

                        for(unsigned int i = 0; i < localRows; ++i)
                        {
                            for(unsigned int d = 0; d < dim; ++d)
                                fscanf(Dataset,"%f, ",&tempX[i*dim + d]);
                        }


                        fclose(Dataset);

                        w = oldN;

                        for(unsigned int i = localRows-NowAdded; i < localRows; ++i)
                        {
                            for(unsigned int d = 0; d < dim; ++d)
                                X[w*dim + d] = tempX[i*dim + d];

                            ++w;
                        }

                        for(unsigned int i = oldN; i < newN; ++i)
                        {
                            for(unsigned int d = 0; d < dim; ++d)
                            {
                                OrderedList[i*dim + d] = X[i*dim + d];
                            }
                        }

                        for(unsigned int i = 0; i < oldN; ++i)
                            Added[i] = 0;

                        for(unsigned int i = oldN; i < newN; ++i)
                            Added[i] = 1;



                        /* -------------------------------------------------------------------------- */

                        start = clock();
                        /*-----------------------Starting Algorithm--------------------------------------*/
                        /* -------------------------------------------------------------------------- */
                        distance = GetDistance(X,oldN,newN,dim,distance);
                        tempDistance = GetEqual(distance,oldN,newN,tempDistance);
                        tempDistance = qSort2D(tempDistance,newN,oldN);

                        // Finding the k-Distance of each element in the dataset
                        kDistance = GetK_Distance(tempDistance,oldN,newN,newKp,kDistance);

                        /* -------------------------------------------------------------------------- */

                        //STEP 2

                        //Finding the Neighborhood of each element
                        Neighborhood = GetNeighborhood(oldN,newN,kDistance,distance,Neighborhood);

                        /* -------------------------------------------------------------------------- */
                        //STEP 3

                        //Getting the Neighborhood size of each element
                        NeighborhoodSize = GetNeighborhoodSize(Neighborhood,oldN,newN,NeighborhoodSize);
                        /* -------------------------------------------------------------------------- */
                        //STEP 4
                        //Calculating the sum of Neighborhood Elements distances for each element

                        SumOfDistances = GetNeighborhoodSum(distance,oldN,newN,Neighborhood,SumOfDistances);

                        /* -------------------------------------------------------------------------- */
                        //STEP 5


                        //Calculating Mean for each element's Neighborhood

                        Mean = GetMean(SumOfDistances,NeighborhoodSize,oldN,newN,Mean);




                        /* -------------------------------------------------------------------------- */

//STEP 6
//Ordering the data according to Mean from max to min

                        for(unsigned int i = 0; i < newN; ++i)
                        {
                            for(unsigned int h = 0; h < newN; ++h)
                            {
                                if(Mean[h] <= Mean[i])
                                {
                                    double tmp = Mean[i];
                                    Mean[i] = Mean[h];
                                    Mean[h] = tmp;

                                    int tmp3 = Added[i];
                                    Added[i] = Added[h];
                                    Added[h] = tmp3;

                                    for(unsigned int d = 0; d < dim; ++d)
                                    {
                                        tmp2[i*dim + d] = OrderedList[i*dim + d];
                                        OrderedList[i*dim + d] = OrderedList[h*dim + d];
                                        OrderedList[h*dim + d] = tmp2[i*dim + d];
                                    }
                                }
                            }
                        }


                        /*---------------------------------------------------------------------------*/
//End of Algorithm
                        end = clock();
                        /*---------------------------------------------------------------------------*/
                        double total_time = ((double) (end - start)) / CLOCKS_PER_SEC;

                        FILE* AddedFile;
                        if((AddedFile = fopen("NewlyAdded.txt","r")) != NULL)
                        {
                            AddedFile = fopen("NewlyAdded.txt","a");
                        }
                        else
                        {
                            AddedFile = fopen("NewlyAdded.txt","w");
                        }



                        for(unsigned int i = 0; i < newN; ++i)
                        {
                            if(Added[i] == 1)
                            {
                                fprintf(AddedFile, "%0.4f ",1/Mean[i]);
                                fprintf(AddedFile, "%0.4f ",Mean[i]);
                                for(unsigned int d = 0; d < dim; ++d)
                                {
                                    fprintf(AddedFile,"%0.4f ",OrderedList[i*dim + d] );
                                }
                                fprintf(AddedFile, "\n");
                            }
                        }

                        fclose(AddedFile);

                        printf("\n");
                        printf("\n Time of Algorithm Execution: %lf \n\n",total_time);
                        free(tempX);
                        free(Added);
                    } else
                    {
                        printf("\nNo new Elements added\n\n" );
                    }
                    /*---------------------------------------------------------------------------*/
                    break;
                /*---------------------------------------------------------------------------*/
                case 2:
                    /*---------------------------------------------------------------------------*/



                    printf("\nLoad new File: ");
                    scanf("%s",Data2);
                    printf("\n\n");

                    FILE *NewDataset;
                    NewDataset = fopen(Data2,"r");
                    if(NewDataset == NULL)
                    {
                        printf("\nThere is something wrong with the file!\n\n");
                        break;
                    }

                    if(newKp == 0)
                        newKp = kPoints;

                    unsigned int localDim = 0;
                    unsigned int localN = 0;

                    if(newN != 0)
                        oldN = newN;

                    localDim = getColumns(NewDataset);
                    if(localDim != dim)
                    {
                        printf("\nUnequal Data\n\n");
                        break;
                    }
                    localN = getRows(NewDataset);      // Getting the Elements of the Dataset
                    rewind(Dataset);
                    localN--;

                    newN = oldN + localN;

                    printf("\n Elements:%d \n", newN);
                    printf("\n Elements Added : %d\n",localN);


                    /*----------------------------------------------------------------------------*/
//Memory Reallocation
                    X = realloc(X,newN*sizeof(*X)*dim);



                    distance = realloc(distance,newN*newN*sizeof(*distance));

                    kDistance = realloc(kDistance,newN*sizeof(*kDistance));

                    SumOfDistances = realloc(SumOfDistances,newN*sizeof(*SumOfDistances));

                    Neighborhood = realloc(Neighborhood,newN*newN*sizeof(*Neighborhood));

                    NeighborhoodSize = realloc(NeighborhoodSize,newN*sizeof(*NeighborhoodSize));

                    Mean = realloc(Mean,newN*sizeof(*Mean));

                    tempDistance = realloc(tempDistance,newN*newN*sizeof(*tempDistance));

                    OrderedList = realloc(OrderedList,newN*sizeof(*OrderedList)*dim);

                    tmp2 = realloc(tmp2,newN*dim*sizeof(*tmp2));

                    int *Added;
                    Added = (int *)calloc(newN,sizeof(int));

                    for(unsigned int i = oldN; i < newN; ++i)
                    {
                        for(unsigned int d = 0; d < dim; ++d)
                            fscanf(NewDataset,"%f, ",&X[i*dim + d]);
                    }


                    fclose(NewDataset);

                    for(unsigned int i = oldN; i < newN; ++i)
                    {
                        for(unsigned int d = 0; d < dim; ++d)
                        {
                            OrderedList[i*dim + d] = X[i*dim + d];
                        }
                    }

                    for(unsigned int i = 0; i < oldN; ++i)
                        Added[i] = 0;

                    for(unsigned int i = oldN; i < newN; ++i)
                        Added[i] = 1;



                    /* -------------------------------------------------------------------------- */

                    start = clock();
                    /*-----------------------Starting Algorithm--------------------------------------*/
                    /* -------------------------------------------------------------------------- */
                    //STEP 1

//**WARNING,MOST TIME COST EFFECTIVE STEP **//

                    // Finding the k-Distance of each element in the dataset
                    distance = GetDistance(X,oldN,newN,dim,distance);
                    tempDistance = GetEqual(distance,oldN,newN,tempDistance);
                    tempDistance = qSort2D(tempDistance,newN,oldN);

                    // Finding the k-Distance of each element in the dataset
                    kDistance = GetK_Distance(tempDistance,oldN,newN,newKp,kDistance);

                    /* -------------------------------------------------------------------------- */

                    //STEP 2

                    //Finding the Neighborhood of each element
                    Neighborhood = GetNeighborhood(oldN,newN,kDistance,distance,Neighborhood);

                    /* -------------------------------------------------------------------------- */
                    //STEP 3

                    //Getting the Neighborhood size of each element
                    NeighborhoodSize = GetNeighborhoodSize(Neighborhood,oldN,newN,NeighborhoodSize);
                    /* -------------------------------------------------------------------------- */
                    //STEP 4
                    //Calculating the sum of Neighborhood Elements distances for each element

                    SumOfDistances = GetNeighborhoodSum(distance,oldN,newN,Neighborhood,SumOfDistances);

                    /* -------------------------------------------------------------------------- */
                    //STEP 5


                    //Calculating Mean for each element's Neighborhood

                    Mean = GetMean(SumOfDistances,NeighborhoodSize,oldN,newN,Mean);


                    /* -------------------------------------------------------------------------- */

//STEP 6
//Ordering the data according to Mean from max to min

                    for(unsigned int i = 0; i < newN; ++i)
                    {
                        for(unsigned int h = 0; h < newN; ++h)
                        {
                            if(Mean[h] <= Mean[i])
                            {
                                double tmp = Mean[i];
                                Mean[i] = Mean[h];
                                Mean[h] = tmp;

                                int tmp3 = Added[i];
                                Added[i] = Added[h];
                                Added[h] = tmp3;

                                for(unsigned int d = 0; d < dim; ++d)
                                {
                                    tmp2[i*dim + d] = OrderedList[i*dim + d];
                                    OrderedList[i*dim + d] = OrderedList[h*dim + d];
                                    OrderedList[h*dim + d] = tmp2[i*dim + d];
                                }
                            }
                        }
                    }


                    /*---------------------------------------------------------------------------*/
//End of Algorithm
                    end = clock();
                    /*---------------------------------------------------------------------------*/
                    double total_time = ((double) (end - start)) / CLOCKS_PER_SEC;

                    FILE* AddedFile;
                    if((AddedFile = fopen("NewlyAdded.txt","r")) != NULL)
                    {
                        AddedFile = fopen("NewlyAdded.txt","a");
                    }
                    else
                    {
                        AddedFile = fopen("NewlyAdded.txt","w");
                    }



                    for(unsigned int i = 0; i < newN; ++i)
                    {
                        if(Added[i] == 1)
                        {
                            fprintf(AddedFile, "%0.4f ",1/Mean[i]);
                            fprintf(AddedFile, "%0.4f ",Mean[i]);
                            for(d = 0; d < dim; ++d)
                            {
                                fprintf(AddedFile,"%0.4f ",OrderedList[i*dim + d] );
                            }
                            fprintf(AddedFile, "\n");
                        }
                    }

                    fclose(AddedFile);

                    printf("\n");
                    printf("\n Time of Algorithm Execution: %lf \n\n",total_time);

                    free(Added);

                    /*---------------------------------------------------------------------------*/
                    break;
                /*---------------------------------------------------------------------------*/
                case 3:
                    /*---------------------------------------------------------------------------*/
                    if(newN == 0)
                        newN = n;

                    printf("\n");
                    printf("Pick a new K: ");   // kPoints
                    scanf("%d",&newKp );
                    printf("\n\n");

                    if(newKp != 0)
                    {

                        while(choise3 != 0)
                        {
                            printf("\n");
                            printf("Want to Rerun now?\n");
                            printf("1 : Yes\n");
                            printf("0 : No\n");
                            printf("Choise : ");
                            scanf("%d",&choise3);
                            printf("\n");

                            switch (choise3) {
                            /*---------------------------------------------------------------------------*/
                            case 0:
                                /*---------------------------------------------------------------------------*/
                                printf("\nReturned Successfuly\n" );
                                printf("\n");
                                /*---------------------------------------------------------------------------*/
                                break;
                            /*---------------------------------------------------------------------------*/
                            case 1:
                                /*---------------------------------------------------------------------------*/
                                start = clock();

                                if(newN == 0)
                                    newN = n;

                                for(unsigned int i = 0; i < newN; ++i)
                                {
                                    for(unsigned int d = 0; d < dim; ++d)
                                        OrderedList[i*dim + d] = X[i*dim + d];
                                }


                                //STEP 1

                                //**WARNING,MOST TIME COST EFFECTIVE STEP **//

                                // Finding the k-Distance of each element in the dataset
                                distance = GetDistance(X,0,newN,dim,distance);
                                tempDistance = GetEqual(distance,0,newN,tempDistance);
                                tempDistance = qSort2D(tempDistance,newN,0);

                                // Finding the k-Distance of each element in the dataset
                                kDistance = GetK_Distance(tempDistance,0,newN,newKp,kDistance);

                                /* -------------------------------------------------------------------------- */

                                //STEP 2

                                //Finding the Neighborhood of each element
                                Neighborhood = GetNeighborhood(0,newN,kDistance,distance,Neighborhood);

                                /* -------------------------------------------------------------------------- */
                                //STEP 3

                                //Getting the Neighborhood size of each element
                                NeighborhoodSize = GetNeighborhoodSize(Neighborhood,0,newN,NeighborhoodSize);
                                /* -------------------------------------------------------------------------- */
                                //STEP 4
                                //Calculating the sum of Neighborhood Elements distances for each element

                                SumOfDistances = GetNeighborhoodSum(distance,0,newN,Neighborhood,SumOfDistances);

                                /* -------------------------------------------------------------------------- */
                                //STEP 5


                                //Calculating Mean for each element's Neighborhood

                                Mean = GetMean(SumOfDistances,NeighborhoodSize,0,newN,Mean);


                                /* -------------------------------------------------------------------------- */

//STEP 6
//Ordering the data according to Mean from max to min

                                for(unsigned int i = 0; i < newN; ++i)
                                {
                                    for(unsigned int h = 0; h < newN; ++h)
                                    {
                                        if(Mean[h] <= Mean[i])
                                        {
                                            double tmp = Mean[i];
                                            Mean[i] = Mean[h];
                                            Mean[h] = tmp;

                                            for(unsigned int d = 0; d < dim; ++d)
                                            {
                                                tmp2[i*dim + d] = OrderedList[i*dim + d];
                                                OrderedList[i*dim + d] = OrderedList[h*dim + d];
                                                OrderedList[h*dim + d] = tmp2[i*dim + d];
                                            }
                                        }
                                    }
                                }


                                /*---------------------------------------------------------------------------*/
//End of Algorithm
                                end = clock();
                                /*---------------------------------------------------------------------------*/
                                double total_time2 = ((double) (end - start)) / CLOCKS_PER_SEC;



                                OrderFile = fopen("OrderFile.txt","w");



                                for(unsigned int i = 0; i < newN; ++i)
                                {
                                    fprintf(OrderFile, "%0.4f ",1/Mean[i]);
                                    fprintf(OrderFile, "%0.4f ",Mean[i]);
                                    for(unsigned int d = 0; d < dim; ++d)
                                    {
                                        fprintf(OrderFile,"%0.4f ",OrderedList[i*dim + d] );
                                    }
                                    fprintf(OrderFile, "\n");
                                }

                                fclose(OrderFile);

                                remove(RemovedFile);

                                sprintf(kname,"output%d.txt",newKp);
                                K_File = fopen(kname,"w");

                                for(unsigned int i = 0; i < newN; ++i)
                                {
                                    fprintf(K_File, "%0.4f ",1/Mean[i]);
                                    fprintf(K_File, "%0.4f ",Mean[i]);
                                    for(unsigned int d = 0; d < dim; ++d)
                                    {
                                        fprintf(K_File,"%0.4f ",OrderedList[i*dim + d] );
                                    }
                                    fprintf(K_File, "\n");
                                }

                                fclose(K_File);


                                printf("\n");
                                printf("\n Time of Algorithm Execution: %lf \n\n",total_time2);
                                /*---------------------------------------------------------------------------*/
                                choise3 = 0;
                                break;
                            /*---------------------------------------------------------------------------*/

                            default:
                                printf("\n Invalid Option (3)\n\n" );
                            }

                        }
                    } else
                    {
                        printf("\n No new K has been chosen!\n\n");
                    }
                    choise3 = 1;
                    /*---------------------------------------------------------------------------*/
                    break;
                /*---------------------------------------------------------------------------*/
                default:
                    printf("\n\n Wrong Input, Choose one of valid options (2)\n\n");
                }
            }
            choise2 = 1;
            /*---------------------------------------------------------------------------*/
            break;
        /*---------------------------------------------------------------------------*/
        case 2:
            /*---------------------------------------------------------------------------*/
            start = clock();
            if(newKp == 0)
                newKp = kPoints;

            if(newN == 0)
                newN = n;

            for(unsigned int i = 0; i < newN; ++i)
            {
                for(unsigned int d = 0; d < dim; ++d)
                    OrderedList[i*dim + d] = X[i*dim + d];
            }


            distance = GetDistance(X,0,newN,dim,distance);
            tempDistance = GetEqual(distance,0,newN,tempDistance);
            tempDistance = qSort2D(tempDistance,newN,0);

            // Finding the k-Distance of each element in the dataset
            kDistance = GetK_Distance(tempDistance,0,newN,newKp,kDistance);

            /* -------------------------------------------------------------------------- */

            //STEP 2

            //Finding the Neighborhood of each element
            Neighborhood = GetNeighborhood(0,newN,kDistance,distance,Neighborhood);

            /* -------------------------------------------------------------------------- */
            //STEP 3

            //Getting the Neighborhood size of each element
            NeighborhoodSize = GetNeighborhoodSize(Neighborhood,0,newN,NeighborhoodSize);
            /* -------------------------------------------------------------------------- */
            //STEP 4
            //Calculating the sum of Neighborhood Elements distances for each element

            SumOfDistances = GetNeighborhoodSum(distance,0,newN,Neighborhood,SumOfDistances);

            /* -------------------------------------------------------------------------- */
            //STEP 5


            //Calculating Mean for each element's Neighborhood

            Mean = GetMean(SumOfDistances,NeighborhoodSize,0,newN,Mean);


            /* -------------------------------------------------------------------------- */

//STEP 6
//Ordering the data according to Mean from max to min

            for(unsigned int i = 0; i < newN; ++i)
            {
                for(unsigned int h = 0; h < newN; ++h)
                {
                    if(Mean[h] <= Mean[i])
                    {
                        double tmp = Mean[i];
                        Mean[i] = Mean[h];
                        Mean[h] = tmp;

                        for(unsigned int d = 0; d < dim; ++d)
                        {
                            tmp2[i*dim + d] = OrderedList[i*dim + d];
                            OrderedList[i*dim + d] = OrderedList[h*dim + d];
                            OrderedList[h*dim + d] = tmp2[i*dim + d];
                        }
                    }
                }
            }


            /*---------------------------------------------------------------------------*/
//End of Algorithm
            end = clock();
            /*---------------------------------------------------------------------------*/
            double total_time2 = ((double) (end - start)) / CLOCKS_PER_SEC;



            OrderFile = fopen("OrderFile.txt","w");



            for(unsigned int i = 0; i < newN; ++i)
            {
                fprintf(OrderFile, "%0.4f ",1/Mean[i]);
                fprintf(OrderFile, "%0.4f ",Mean[i]);
                for(unsigned int d = 0; d < dim; ++d)
                {
                    fprintf(OrderFile,"%0.4f ",OrderedList[i*dim + d] );
                }
                fprintf(OrderFile, "\n");
            }

            fclose(OrderFile);

            sprintf(kname,"output%d.txt",newKp);
            K_File = fopen(kname,"w");

            for(unsigned int i = 0; i < newN; ++i)
            {
                fprintf(K_File, "%0.4f ",1/Mean[i]);
                fprintf(K_File, "%0.4f ",Mean[i]);
                for(unsigned int d = 0; d < dim; ++d)
                {
                    fprintf(K_File,"%0.4f ",OrderedList[i*dim + d] );
                }
                fprintf(K_File, "\n");
            }

            fclose(K_File);


            remove(RemovedFile);


            printf("\n");
            printf("\n Time of Algorithm Execution: %lf \n\n",total_time2);



            /*---------------------------------------------------------------------------*/
            break;
        /*---------------------------------------------------------------------------*/
        case 3:
            /*---------------------------------------------------------------------------*/


            if((AddedFileMerge = fopen("NewlyAdded.txt","r")) != NULL)
            {
                AddedFileMerge = fopen("NewlyAdded.txt","r");
            } else
            {
                printf("\nNewlyAdded.txt Do not Exist\n");
                break;
            }

            OrderFileMerge = fopen("OrderFile.txt","r");

            unsigned int dimLocal;
            unsigned int nAdd;
            unsigned int nOrder;
            unsigned int nTotal;
            dimLocal = getColumns(OrderFileMerge);
            rewind(OrderFileMerge);
            nAdd = getRows(AddedFileMerge);
            rewind(AddedFileMerge);
            nAdd--;
            nOrder = getRows(OrderFileMerge);
            rewind(OrderFileMerge);
            nOrder--;
            nTotal =  nAdd + nOrder;

            float *Merged;
            Merged = (float *)calloc(nTotal*dimLocal,sizeof(float));

            float *tmp4;
            tmp4 =(float*) calloc(nTotal*dimLocal,sizeof(float));

            start = clock();

            for(unsigned int i = 0; i < nAdd; ++i)
            {
                for(unsigned int d = 0; d < dimLocal; ++d)
                {
                    fscanf(AddedFileMerge,"%f ",&Merged[i*dimLocal + d]);
                }
            }

            for(unsigned int i = nAdd; i < nTotal; ++i)
            {
                for(unsigned int d = 0; d < dimLocal; ++d)
                {
                    fscanf(OrderFileMerge,"%f ",&Merged[i*dimLocal + d]);
                }
            }





            for(unsigned int i = 0; i < nTotal; ++i)
            {
                for(unsigned int h = 0; h < nTotal; ++h)
                {
                    if(Merged[h*dimLocal + 1] <= Merged[i*dimLocal + 1])
                    {

                        for(unsigned int d = 0; d < dimLocal; ++d)
                        {
                            tmp4[i*dimLocal + d] = Merged[i*dimLocal + d];
                            Merged[i*dimLocal + d] = Merged[h*dimLocal + d];
                            Merged[h*dimLocal + d] = tmp4[i*dimLocal + d];
                        }
                    }
                }
            }




            OrderFile = fopen("OrderFile.txt","w");

            for(unsigned int i = 0; i < nTotal; ++i)
            {
                for(unsigned int d = 0; d < dimLocal; ++d)
                {
                    fprintf(OrderFile,"%0.4f ",Merged[i*dimLocal + d] );
                }
                fprintf(OrderFile, "\n");
            }

            fclose(OrderFile);

            end = clock();
            /*---------------------------------------------------------------------------*/
            double total_time = ((double) (end - start)) / CLOCKS_PER_SEC;
            printf("\n");
            printf("\n Time of Merging: %lf \n\n",total_time);


            remove(RemovedFile);
            fclose(OrderFileMerge);
            fclose(AddedFileMerge);
            free(Merged);
            free(tmp4);
            /*---------------------------------------------------------------------------*/
            break;
        /*---------------------------------------------------------------------------*/

        /*---------------------------------------------------------------------------*/
        default:
            printf("\n\n Wrong Input, Choose one of valid options (1)\n\n");
        }

    }



    free(kDistance);
    free(distance);
    free(Neighborhood);
    free(NeighborhoodSize);
    free(SumOfDistances);
    free(X);
    free(Mean);
    free(tempDistance);
    free(OrderedList);
    free(tmp2);



    return 0;
}
