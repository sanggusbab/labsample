#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#define num_steps 10000000
#define NUM_THREADS 4

int main(int args, char* argv[])
{
	double sum[NUM_THREADS], step, x, pi;
	int i, id;
	clock_t start, end;

	step = 1./(double)num_steps;

	omp_set_num_threads(NUM_THREADS);

	printf("# of steps: %d\n\n", num_steps);

	for(i = 0; i<NUM_THREADS; i++)
	{
		sum[i] = 0;
	}
	start = clock();
#pragma omp parallel private(id)
	{
		id = omp_get_thread_num();
#pragma omp for private(x)
		for(i = id; i<num_steps; i++)
		{
			x = (i +0.5)*step;
			sum[id] += 4.0/(1.0+x*x);
		}
	}
	for(i = 0, pi = 0.0; i<NUM_THREADS;i++)
	{
		pi += sum[i]*step;
	}
	end = clock();
	printf("multi threading(cores: %d) : %.3lf seconds\n", NUM_THREADS,(end-start)/(double)1000000);
        printf(" numerical pi = %.15f \n", pi);
        printf(" analytical pi = %.15f \n", acos(-1.0));
        printf(" Error = %E \n", fabs(acos(-1.0)-pi));

	printf("\n ----------------------- \n\n");
	pi = 0.0;
	start = clock();
	for(i=0; i<num_steps;i++)
	{
		x = (i+0.5)*step;
		pi += 4.0/(1.0+x*x);
	}
	end = clock();
	pi *= step;
	printf("single threading : %.3lf seconds\n",(end-start)/(double)1000000);
	printf(" numerical pi = %.15f \n", pi);
	printf(" analytical pi = %.15f \n", acos(-1.0));
	printf(" Error = %E \n", fabs(acos(-1.0)-pi));

        printf("\n ----------------------- \n\n");

	printf("\n all variables new allocated\n\n");

	double x2, step2;
	step2 = 1./(double)num_steps;
	int j;
        double pi2 = 0.0;
        start = clock();
        for(j=0; j<num_steps;j++)
        {
                x2 = (j+0.5)*step2;
                pi2 += 4.0/(1.0 + x2 * x2);
        }
        end = clock();
        pi2 *= step2;
        printf("single threading : %.3lf seconds\n",(end-start)/(double)1000000);
        printf(" numerical pi = %.15f \n", pi2);
        printf(" analytical pi = %.15f \n", acos(-1.0));
        printf(" Error = %E \n", fabs(acos(-1.0)-pi2));

        printf("\n ----------------------- \n\n");

        for(i = 0; i<NUM_THREADS; i++)
        {
                sum[i] = 0;
        }
        start = clock();
#pragma omp parallel private(id)
        {
                id = omp_get_thread_num();
#pragma omp for private(x)
                for(i = id; i<num_steps; i++)
                {
                        x = (i +0.5)*step;
                        sum[id] += 4.0/(1.0+x*x);
                }
        }
        for(i = 0, pi = 0.0; i<NUM_THREADS;i++)
        {
                pi += sum[i]*step;
        }
        end = clock();
        printf("multi threading(cores: %d) : %.3lf seconds\n", NUM_THREADS,(end-start)/(double)1000000);
        printf(" numerical pi = %.15f \n", pi);
        printf(" analytical pi = %.15f \n", acos(-1.0));
        printf(" Error = %E \n", fabs(acos(-1.0)-pi));


	return 0;
}
