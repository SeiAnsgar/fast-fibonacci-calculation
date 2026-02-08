#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>     //adds uint_32
#include <inttypes.h>   //for formatting uint_32 so it can be printed
#include <gmp.h>        //link with -lgmp
#include <omp.h>        //openMP
#include <time.h>       
#include <math.h>       //link with -lm

typedef struct{
    mpz_t a,b,d;
    //-> fn+1 fn fn-1
}Matrix;
//    |a b|  = |fn+1    fn|
//    |c d|    |fn    fn-1|

int split_exponent(uint32_t n);
void calculate_small_fraction(int i, Matrix *mf);
void matrix_mul(Matrix *m1, Matrix *m2, Matrix *rop); //f: M x M -> M, (Ma, Mb) -> Ma*Mbs
void matrix_mul_square(int exp, Matrix *set_m); //change Matrix  m to *m where *m points to collector matrix

void final_product(int fraction);

int exponent[32];

int main(int argc, char *argv[])
{   
    int fraction;
    uint32_t a;

    if(argc != 2){
        printf("Missing input value!\n");
        return 0;    
    }

    a = (unsigned)atoi(argv[1]);

    
    fraction = split_exponent(a);
    final_product(fraction);
}

void matrix_mul_square(int exp, Matrix *set_m)
{   
    Matrix m_init;
    mpz_init_set_ui(m_init.a, 1);
    mpz_init_set_ui(m_init.b, 1);
    mpz_init_set_ui(m_init.d, 0);

    Matrix m_temp;
    mpz_inits(m_temp.a, m_temp.b, m_temp.d, NULL);

    int edited_most_recent = 0;

    for(int i=1; i<=exp; i++)
    {   
        //override eachother alternately -> no need to initalize m_init with m_temp after eacht iteration
        if(i%2){
            //printf("i mod 2\n");
            matrix_mul(&m_init, &m_init, &m_temp);
            edited_most_recent = 1;
        }
        else{
            matrix_mul(&m_temp, &m_temp, &m_init);
            edited_most_recent = 2;
        }   
    }
    if(edited_most_recent == 1){
        //gmp_printf("matrix_square:%Zd \n", m_temp.b);
        mpz_set(set_m->a, m_temp.a);
        mpz_set(set_m->b, m_temp.b);
        mpz_set(set_m->d, m_temp.d);
    }
    else{
        //gmp_printf("matrix_square:%Zd \n", m_init.b);
        mpz_set(set_m->a, m_init.a);
        mpz_set(set_m->b, m_init.b);
        mpz_set(set_m->d, m_init.d);
    }
        
    mpz_clears(m_init.a, m_init.b, m_init.d, NULL);
    mpz_clears(m_temp.a, m_temp.b, m_temp.d, NULL);
}


void matrix_mul(Matrix *m1, Matrix *m2, Matrix *rop)
{
    //void mpz_addmul (mpz t rop, const mpz t op1, const mpz t op2)
    
    //gmp_printf("m1.a: %Zd \n", m1->a);
    //gmp_printf("m2.a: %Zd \n", m2->a);
    int m1_is_rop = 0;

    if(m1 == rop){

        Matrix tmp;

        m1_is_rop = 1;

        mpz_inits(tmp.a, tmp.b, tmp.d, NULL);
        mpz_set(tmp.a, m1->a);
        mpz_set(tmp.b, m1->b);
        mpz_set(tmp.d, m1->d);

        mpz_mul(rop->a, tmp.a, m2->a);
        mpz_addmul(rop->a, tmp.b, m2->b);

        mpz_mul(rop->b, tmp.a, m2->b);
        mpz_addmul(rop->b, tmp.b, m2->d);

        mpz_mul(rop->d, tmp.b, m2->b);
        mpz_addmul(rop->d, tmp.d, m2->d);

        mpz_clears(tmp.a, tmp.b, tmp.d, NULL);
    }
    else{
        mpz_mul(rop->a, m1->a, m2->a);
        //gmp_printf("rop.a: %Zd \n", rop->a);
        mpz_addmul(rop->a, m1->b, m2->b);
        //gmp_printf("% Zd \n", rop.a);

        mpz_mul(rop->b, m1->a, m2->b);
        mpz_addmul(rop->b, m1->b, m2->d);

        mpz_mul(rop->d, m1->b, m2->b);
        mpz_addmul(rop->d, m1->d, m2->d);
    }
}


int split_exponent(uint32_t n)
{    
    int small_fraction = 0;
    int slot = 0;
    int power_zero = 0;
    
    //check if 2⁰ bit exists
    if(n%2)
        power_zero=1;

    for(int i=31; i >= 0; i--)
    {   
        if((n >> i)){
            exponent[slot] = i; 
            n -= 1 << i;
            slot++;
            //printf("value: %i \n", i);
        }
    }
    //separate small portion of input-number, matrix exponentiation is only worth it for large n
    for(int k = 0; k<=31; k++)
    {   
        if(exponent[k] && exponent[k] < 4){
            small_fraction += pow(2,exponent[k]);
        }
        else if(exponent[k] == 0)
            break;
    }
    small_fraction += power_zero;
        
    printf("small_fraction:%i \n",small_fraction);
    return small_fraction;
}


void calculate_small_fraction(int n, Matrix *mf)
{   
    unsigned long i, t1, t2, next, pre;
    t1 = 1;
    t2 = 0;
    next = 1;

    for(i=1; i<=n+1; i++)
    {
        pre = t1;
        next = t1 + t2;
        t1 = t2;
        t2 = next;
    }
    //printf("small fration result t1:%i\n",t1);
    mpz_set_ui(mf->a,t2);
    mpz_set_ui(mf->b,t1);
    mpz_set_ui(mf->d,pre);
}


void final_product(int fraction)
{   
    clock_t begin = clock();

    Matrix m_collector;
    mpz_inits(m_collector.a, m_collector.b, m_collector.d, NULL);

    Matrix tmp;
    mpz_inits(tmp.a, tmp.b, tmp.d, NULL);

    calculate_small_fraction(fraction, &m_collector);
    /*
    printf("after calculate_small_fraction():\n");
    gmp_printf("mcollector.a result:%Zd \n", m_collector.a);
    gmp_printf("mcollector.b result:%Zd \n", m_collector.b);
    gmp_printf("mcollector.d result:%Zd \n", m_collector.d);
    */
    

    for(int i = 0; i<=31; i++)
    {
        if(exponent[i]>=4){
            //printf("i:%i  exponent[i]: %i\n", i, exponent[i]);

            matrix_mul_square(exponent[i], &tmp);
            /*
            printf("after matrix_mul_square():\n");
            gmp_printf("temp.a result:%Zd \n", tmp.a);
            gmp_printf("temp.b result:%Zd \n", tmp.b);
            gmp_printf("temp.d result:%Zd \n", tmp.d);
            */
            

            matrix_mul(&m_collector, &tmp, &m_collector);
        }
        else
            break;
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    gmp_printf("Final result:%Zd \n", m_collector.b);
    printf("\nexec_time: %f \n", time_spent);

    mpz_clears(m_collector.a, m_collector.b, m_collector.d, NULL);
    mpz_clears(tmp.a, tmp.b, tmp.d, NULL);
}