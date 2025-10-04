#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>     //adds uint_32
#include <inttypes.h>   //for formatting uint_32 so it can be printed
#include <gmp.h>
//compile with -lgmp
#include <omp.h>        //openMP
#include <time.h>       
#include <math.h>

void split_exponent(uint32_t n);
void fib_long(int n);
void fib_gmp(int n);

unsigned int matrix_square(int n);
//void fib_mat(int n);


uint32_t exponent[32];

int main(int argc, char *argv[])
{   
    if(argc != 2){
        printf("kein Eingabeparameter!\n");
        return 0;    
    }
    int n = atoi(argv[1]);
    uint32_t a = (unsigned)atoi(argv[1]);

    int actual = pow(2,n);
    printf("%i th fib: ", actual);
    matrix_square(n);

    //split_exponent(a);
    //fib_long(n);
    //fib_gmp(n);
    //fib_mat(n);
    
    //printf("n: %d\n", n);
}

void fib_long(int n)
{
    unsigned long t1, t2, next;
    t1 = 1;
    t2 = 0;
    next = 1;
    int i;

    for(i = 1; i<=n; i++)
    {   
        next = t1 + t2;
        printf("n: %d result: %lu\n", i, next);
        t1 = t2;
        t2 = next;
    }
    printf("result: %lu\n", next);
}

void fib_gmp(int n)
{
    mpz_t g1, g2, gnext;
    mpz_init_set_ui(g1,1);
    mpz_init_set_ui(g2,0);
    mpz_init(gnext);

    clock_t begin = clock();

    long int i;
    for(i=1; i<=n; i++)
    {   
        mpz_add(gnext, g1, g2);
        mpz_set(g1, g2);
        mpz_set(g2, gnext);
    }

    mpz_clear(g1);
    mpz_clear(g2);
    
    clock_t end = clock();

    gmp_printf("% Zd \n", gnext);
    mpz_clear(gnext);

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("exec_time: %f \n", time_spent);
}


void split_exponent(uint32_t n)
{    
    for(int i=31; i >= 0; i--)
    {   
        if((n >> i)){
            exponent[i] = i;
            n -= 1 << i;
            printf("value: %i \n", i);
        }
    } 
    
    uint32_t a;
    for(int k =0; k<32; k++)
    {   
        a = exponent[k];
        printf("value: %zu \n", a);
    }
}

unsigned int matrix_square(int n)
{
    //    |a b|  = |fn+1    fn|
    //    |c d|    |fn    fn-1|

    //n: 2¹=2  2²=4  2³=8  2⁴=16 
//------------------------
    mpz_t a_x, b_x, d_x;
    mpz_init_set_ui(a_x, 1);
    mpz_init_set_ui(b_x, 1);
    mpz_init_set_ui(d_x, 0);
    
    mpz_t a_y, b_y, d_y;
    mpz_init(a_y);
    mpz_init(b_y);
    mpz_init(d_y);

    mpz_t a_x_temp, b_x_temp, d_x_temp;
    mpz_init(a_x_temp);
    mpz_init(b_x_temp);
    mpz_init(d_x_temp);

    int i;
    for(i=1; i<=n; i++)
    {
        mpz_mul(a_x_temp, a_x, a_x);
        mpz_mul(b_x_temp, b_x, b_x);
        mpz_add(a_y, a_x_temp, b_x_temp);

        mpz_mul(a_x_temp, a_x, b_x);
        mpz_mul(b_x_temp, b_x, d_x);
        mpz_add(b_y, a_x_temp, b_x_temp);
        
        mpz_mul(b_x_temp, b_x, b_x);
        mpz_mul(d_x_temp, d_x, d_x);
        mpz_add(d_y, b_x_temp, d_x_temp);
      
        mpz_set(a_x, a_y);
        mpz_set(b_x, b_y);
        mpz_set(d_x, d_y);
    }

    mpz_clear(a_x);
    mpz_clear(b_x);
    mpz_clear(d_x);
    mpz_clear(a_x_temp);
    mpz_clear(b_x_temp);
    mpz_clear(d_x_temp);
    mpz_clear(a_y);
    mpz_clear(d_y);

    
    gmp_printf("% Zd \n", b_y);
    mpz_clear(b_y);
}


/*
void fib_mat(int n)
{   
    
    //    |a b|  = |fn+1    fn|
    //    |c d|    |fn    fn-1|
    //
    //              |a b|
    //            * |c d| 
    //    |a_n b_n| = (a*a_n + b_n*c)  (a_n*b + b_n*d)
    //    |c_n d_n| = (c_n*a + d_n*c)  (c_n*b + d_n*d)
    

    unsigned long a = 1;
    unsigned long b = 1;
    unsigned long d = 0;

    unsigned long a_temp;
    unsigned long b_temp;
    unsigned long d_temp;
    int i;
    
    //clock_t begin = clock();
    


    //matrix mult for fib, a = fn+1 -> do n-1 loop iterations and return a to get value of fn
    a_temp = a + b;
    b_temp = a;
    d_temp = b;

    a = a_temp;
    b = b_temp;
    d = d_temp;
    //printf("i: %d result: %lu\n",i,a);
    //printf("result: %ld\n",a);
        
    
    //clock_t end = clock();
    
    printf("result: %lu\n",a);
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("exec_time: %f \n", time_spent);
}
*/


/*
int a_x = 1;
    int b_x = 1;
    //int c_x = 1;
    int d_x = 0;
    int a_y, b_y, c_y, d_y;
    
    int i;
    for(i=1; i<=n; i++)
    {
        a_y = a_x*a_x + b_x*b_x;
        b_y = a_x*b_x + b_x*d_x;
        //c_y = c_x*a_x + d_x*c_x;
        c_y = b_y;
        d_y = b_x*b_x + d_x*d_x;

        a_x = a_y;
        b_x = b_y;
        //c_x = c_y;
        d_x = d_y; 
        printf("%dth fib: %u \n", i, b_x);
    } 

*/
