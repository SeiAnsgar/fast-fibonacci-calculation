#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>     //adds uint_32
#include <inttypes.h>   //for formatting uint_32 so it can be printed
#include <gmp.h>        //link with -lgmp
#include <omp.h>        //openMP
#include <time.h>       
#include <math.h>       //link with -lm


void fib_long(int n);
void fib_gmp(int n);

int split_exponent(int n);
int fib_fraction(int n);
void matrix_square(int exp, mpz_t collect_values);
//void fib_mat(int n);

int exponent[32];

int main(int argc, char *argv[])
{   
    if(argc != 2){
        printf("missing input value!\n");
        return 0;    
    }
    int n = atoi(argv[1]);
    uint32_t a = (unsigned)atoi(argv[1]);

    int fraction;
    fraction = fib_fraction(split_exponent(a));
    printf("fraction: %i\n",fraction);

    mpz_t collect;
    mpz_init(collect);
    mpz_init_set_ui(collect, 1);
    mpz_mul_ui(collect, collect, fraction);

    for(int i=0; i<=26; i++)
    {
        if(exponent[i]>0){
            matrix_square(exponent[i], collect);
        }
    }

    printf("%i th fib: ", n);
    gmp_printf("% Zd \n", collect);
    mpz_clear(collect);
}

int fib_fraction(int n)
{
    int t1, t2, next;
    t1 = 1;
    t2 = 0;
    next = 1;
    
    for(int i = 1; i<=n; i++)
    {   
        next = t1 + t2;
        t1 = t2;
        t2 = next;
    }
    return next;
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


int split_exponent(int n)
{    
    int small_fraction = 0;

    for(int i=31; i >= 0; i--)
    {   
        if((n >> i)){
            exponent[31-i] = i; 
            n -= 1 << i;
            //printf("value: %i \n", i);
        }
    }
    
    for(int k = 27; k<=31; k++)
        small_fraction += pow(2,exponent[k]);

    return small_fraction;
}


void matrix_square(int exp, mpz_t collect_values)
{
    //    |a b|  = |fn+1    fn|
    //    |c d|    |fn    fn-1|

    //n: 2¹=2  2²=4  2³=8  2⁴=16 

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
    for(i=1; i<=exp; i++)
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
    mpz_clear(d_y);
    
    mpz_mul(collect_values, collect_values, b_y);
    gmp_printf("b_y: %Zd \n", b_y);
    
    mpz_clear(a_y);
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
    

    //clock_t begin = clock();  
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
