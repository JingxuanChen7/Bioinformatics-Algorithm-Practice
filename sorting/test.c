#include <stdio.h>
#include <string.h>
void QuickSort(int array[], int min, int max);
void swap(int *x, int *y);
int Partition(int array[], int min, int max);

int main()
{
    int array[] = {8,3,734,23,8,44,985,653};
    //int len = sizeof(array) / sizeof(array[0]), i;
    int len = 8, i;
    QuickSort(array, 0, len);
    for (i = 1; i <= len; i++)
    {
        printf("%d\n",array[i]);
    }
    return 0;
}

void QuickSort(int array[], int min, int max)
{
    if(min<max)
    {
        int pos = Partition(array, min, max);
        QuickSort(array, min , pos-1);
        QuickSort(array, pos+1, max);
    }
}

int Partition(int array[], int min, int max)
{
    int pivot = array[max], x = min, y = max-1;
    while (1)
    {
        while (pivot > array[x])
            x++;
        while (pivot < array[y])
            y--;
        if(x<y) 
        {
            swap(&array[x],&array[y]);
        }
        else
        {
            swap(&array[x],&array[max]);
            return x;
        }       
    }
}

void swap(int *x, int *y)
{
    int tmp = *y;
    *y = *x;
    *x = tmp;
}