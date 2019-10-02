#include <math.h>

float *GetDistance(float *,int ,int,int,float *);
float *GetEqual(float *,int , int , float *);
float * qSort1D(float *,int);
float * qSort2D(float *,int, int);
float *GetK_Distance(float *,int,int,int,float *);
int *GetNeighborhood(int, int, float *,float *,int *);
int *GetNeighborhoodSize(int *,int , int , int *);
float *GetNeighborhoodSum(float *,int ,int ,int *,float *);
float* GetMean(float *,int *,int ,int ,float *);







float *GetDistance(float *X,int start,int end,int dim,float *distance)
{

    for(unsigned int i = start; i < end; ++i) //for each Element
    {
        distance[i*end + i] =  9999;
        for(unsigned int h = 0; h < end; ++h)
        {
            if(h != i)
            {
                /* Calculating the distance of each element with i element
                and then we store it to an array so that we can order it later,
                we set the distance of i to max value so that we can avoid it at
                the ordering */
                if(distance[i*end + h] != 0) distance[i*end + h] =  0;

                for(unsigned int d = 0; d < dim; d++)
                    distance[i*end + h] += (X[h*dim + d] - X[i*dim + d])*(X[h*dim + d] - X[i*dim + d]);

                distance[i*end + h] = sqrt(distance[i*end +h]);

            }
        }
    }

    return distance;

}

float *GetEqual(float *distance,int start, int end, float *tempDistance)
{


    for(unsigned int i = start; i < end; ++i) //for each Element
    {
        tempDistance[i*end + i] = 99999;
        for(unsigned int h = 0; h < end; ++h)
        {
            if(h != i)
            {
                /* Calculating the distance of each element with i element
                and then we store it to an array so that we can order it later,
                we set the distance of i to max value so that we can avoid it at
                the ordering */

                tempDistance[i*end + h] = distance[i*end + h];


            }
        }
    }

    return tempDistance;

}


float *qSort1D(float *arr, int elements) {

#define  MAX_LEVELS  300

    float  piv, beg[MAX_LEVELS], end[MAX_LEVELS],  swap ;
    int i=0, L, R;
    beg[0]=0;
    end[0]=elements;
    while (i>=0) {
        L=beg[i];
        R=end[i]-1;
        if (L<R) {
            piv=arr[L];
            while (L<R) {
                while (arr[R]>=piv && L<R) R--;
                if (L<R) arr[L++]=arr[R];
                while (arr[L]<=piv && L<R) L++;
                if (L<R) arr[R--]=arr[L];
            }
            arr[L]=piv;
            beg[i+1]=L+1;
            end[i+1]=end[i];
            end[i++]=L;
            if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
                swap=beg[i];
                beg[i]=beg[i-1];
                beg[i-1]=swap;
                swap=end[i];
                end[i]=end[i-1];
                end[i-1]=swap;
            }
        }
        else {
            i--;
        }
    }
    return arr;
}


float *qSort2D(float *Array,int elements,int start)
{
#define  MAX_LEVELS  300
    int i;

    for(i = start; i < elements; ++i) //for each Element
    {



        /* -------------------------------------------------------------------------- */
        float  piv, beg[MAX_LEVELS], end[MAX_LEVELS],  swap ;
        int k=0, L, R;
        beg[0]=0;
        end[0]=elements;
        while (k>=0) {
            L=beg[k];
            R=end[k]-1;
            if (L<R) {
                piv=Array[i*elements + L];
                while (L<R) {
                    while (Array[i*elements + R]>=piv && L<R) R--;
                    if (L<R) Array[i*elements + L++]=Array[i*elements + R];
                    while (Array[i*elements + L]<=piv && L<R) L++;
                    if (L<R) Array[i*elements + R--]=Array[i*elements + L];
                }
                Array[i*elements + L]=piv;
                beg[k+1]=L+1;
                end[k+1]=end[k];
                end[k++]=L;
                if (end[k]-beg[k]>end[k-1]-beg[k-1]) {
                    swap=beg[k];
                    beg[k]=beg[k-1];
                    beg[k-1]=swap;
                    swap=end[k];
                    end[k]=end[k-1];
                    end[k-1]=swap;
                }
            }
            else {
                k--;
            }
        }
    }

    return Array;

}

float *GetK_Distance(float *tempDistance,int start, int end,int kPoints,float *kDistance)
{

    kPoints--;
    for(unsigned int i = start; i < end; ++i) //for each Element
    {
        kDistance[i] = tempDistance[i*end + (kPoints)];

    }


    return kDistance;
}

int *GetNeighborhood(int start, int end, float *kDistance,float *distance,int *Neighborhood)
{

    for(unsigned int i = start; i < end; ++i)//for each Element
    {
        if(Neighborhood[i*end + i] != 0) Neighborhood[i*end + i] = 0;
        for(unsigned int h = 0; h < end; ++h)
        {
            if(h != i)
            {
                /* Calculating the distance of each element with i element
                and then we store it to an array so that we can determine if that
                element belongs to the Neighborhood of i element */
                if(Neighborhood[i*end + h] != 0) Neighborhood[i*end + h] = 0;


                /* Comparing the distance of each potential Neighborhood element
                with the k-distance of i element, if found less, set
                Neighborhood[master][currentElement] to 1, that means that it belongs
                to i element's Neighborhood */
                if(distance[i*end + h] <= kDistance[i])
                {
                    Neighborhood[i*end + h] = 1;
                }
            }
        }
    }

    return Neighborhood;
}

int *GetNeighborhoodSize(int *Neighborhood,int start, int end, int *NeighborhoodSize)
{


    for(unsigned int i = start; i < end; ++i)//for each Element
    {
        unsigned int counter;
        if(counter != 0) counter = 0;
        for(unsigned int h = 0; h < end; ++h)
        {
            if(h != i)
            {

                if(Neighborhood[i*end + h] != 0)
                {
                    ++counter;
                }
            }
        }
        NeighborhoodSize[i] = counter;

        // printf("%d Neighborhood: %d\n",i,NeighborhoodSize[i] );
    }

    return NeighborhoodSize;
}

float *GetNeighborhoodSum(float *distance,int start,int end,int *Neighborhood,float *SumOfDistances)
{

    for(unsigned int i = start; i < end; ++i)//for each Element
    {
        if(SumOfDistances[i] != 0) SumOfDistances[i] = 0;
        for(unsigned int h = 0; h < end; ++h)
        {
            if(h != i)
            {
                if(Neighborhood[i*end + h] != 0)
                {
                    SumOfDistances[i] += distance[i*end + h];
                }
            }
        }
    }

    return SumOfDistances;
}

float *GetMean(float *SumOfDistances,int *NeighborhoodSize,int start,int end,float *Mean)
{

    for(unsigned int i = start; i < end; ++i)
    {
        Mean[i] = SumOfDistances[i]/NeighborhoodSize[i];
    }
    return Mean;
}
