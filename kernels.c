/********************************************************
 * Kernels to be optimized for the CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"

/* ---------------------------
   Portable RESTRICT macro
   ---------------------------
   Place this right after the #include lines.
*/
#ifndef RESTRICT
# if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#  define RESTRICT restrict
# else
#  define RESTRICT /* nothing */
# endif
#endif


/* 
 * Please fill in the following team struct 
 */
team_t team = {
    "Image Filters Optimization",            

    "Disha Suthar",     /* Member full name */
    "",  

    "",                  
    ""                    
};

/***************
 * ROTATE KERNEL
 ***************/

/******************************************************
 * Your different versions of the rotate kernel go here
 ******************************************************/

/* 
 * naive_rotate - The naive baseline version of rotate 
 */
char naive_rotate_descr[] = "naive_rotate: Naive baseline implementation";
void naive_rotate(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
}

/* 
 * rotate_block - Blocking optimization for rotate
 */
char rotate_block_descr[] = "rotate_block: blocking optimization (B=32)";
void rotate_block(int dim, pixel *src, pixel *dst) 
{
    int i, j, i0, j0;
    const int B = 32; /* block size, try 16 or 32 */

    for (i0 = 0; i0 < dim; i0 += B) {
        for (j0 = 0; j0 < dim; j0 += B) {
            for (i = i0; i < i0 + B; i++) {
                for (j = j0; j < j0 + B; j++) {
                    dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
                }
            }
        }
    }
}

/* 
 * rotate_block_unrolled - Blocking + loop unrolling optimization
 */
char rotate_block_unrolled_descr[] = "rotate_block_unrolled: blocking + unrolling (B=32, unroll=4)";
void rotate_block_unrolled(int dim, pixel *src, pixel *dst) 
{
    int i0, j0, i, j;
    const int B = 32;

    for (i0 = 0; i0 < dim; i0 += B) {
        for (j0 = 0; j0 < dim; j0 += B) {
            for (i = i0; i < i0 + B; i++) {
                int src_idx = i * dim + j0;
                int dst_idx = (dim - 1 - j0) * dim + i;

                for (j = j0; j < j0 + B; j += 4) {
                    /* Unroll by 4 */
                    dst[dst_idx]            = src[src_idx];
                    dst[dst_idx - dim]      = src[src_idx + 1];
                    dst[dst_idx - 2*dim]    = src[src_idx + 2];
                    dst[dst_idx - 3*dim]    = src[src_idx + 3];

                    src_idx += 4;
                    dst_idx -= 4*dim;
                }
            }
        }
    }
}



/* helper to stringify B */
#define xstr(s) str(s)
#define str(s) #s

/* B = 16 */
#undef B
#define B 16
char rotate_block16_descr[] = "rotate_block: blocking optimization (B=16)";
void rotate_block16(int dim, pixel *src, pixel *dst) {
    int i, j, ii, jj;
    for (i = 0; i < dim; i += B) {
        for (j = 0; j < dim; j += B) {
            for (ii = i; ii < i + B && ii < dim; ii++) {
                for (jj = j; jj < j + B && jj < dim; jj++) {
                    dst[RIDX(dim-1-jj, ii, dim)] = src[RIDX(ii, jj, dim)];
                }
            }
        }
    }
}

/* B = 32 */
#undef B
#define B 32
char rotate_block32_descr[] = "rotate_block: blocking optimization (B=32)";
void rotate_block32(int dim, pixel *src, pixel *dst) {
    int i, j, ii, jj;
    for (i = 0; i < dim; i += B) {
        for (j = 0; j < dim; j += B) {
            for (ii = i; ii < i + B && ii < dim; ii++) {
                for (jj = j; jj < j + B && jj < dim; jj++) {
                    dst[RIDX(dim-1-jj, ii, dim)] = src[RIDX(ii, jj, dim)];
                }
            }
        }
    }
}

/* B = 64 */
#undef B
#define B 64
char rotate_block64_descr[] = "rotate_block: blocking optimization (B=64)";
void rotate_block64(int dim, pixel *src, pixel *dst) {
    int i, j, ii, jj;
    for (i = 0; i < dim; i += B) {
        for (j = 0; j < dim; j += B) {
            for (ii = i; ii < i + B && ii < dim; ii++) {
                for (jj = j; jj < j + B && jj < dim; jj++) {
                    dst[RIDX(dim-1-jj, ii, dim)] = src[RIDX(ii, jj, dim)];
                }
            }
        }
    }
}

/* Optimized rotate: 32x32 blocking + unroll 8 + index hoisting + restrict */
char rotate_block32_unroll8_descr[] =
    "rotate_block32_unroll8: 32x32 block, unroll 8, code-motion, restrict";

void rotate_block32_unroll8(int dim, pixel * RESTRICT src, pixel * RESTRICT dst)
{
    const int BLK = 32;  /* block size */
    int i0, j0, i, j;

    for (i0 = 0; i0 < dim; i0 += BLK) {
        for (j0 = 0; j0 < dim; j0 += BLK) {
            for (i = i0; i < i0 + BLK; i++) {
                int src_row = i * dim + j0;
                int dst_col = i;

                for (j = 0; j < BLK; j += 8) {
                    int jj = j0 + j;
                    int dst_base = (dim - 1 - jj) * dim + dst_col;

                    dst[dst_base]           = src[src_row + j];
                    dst[dst_base -     dim] = src[src_row + j + 1];
                    dst[dst_base - 2 * dim] = src[src_row + j + 2];
                    dst[dst_base - 3 * dim] = src[src_row + j + 3];
                    dst[dst_base - 4 * dim] = src[src_row + j + 4];
                    dst[dst_base - 5 * dim] = src[src_row + j + 5];
                    dst[dst_base - 6 * dim] = src[src_row + j + 6];
                    dst[dst_base - 7 * dim] = src[src_row + j + 7];
                }
            }
        }
    }
}

/* Optimized rotate: 16x16 block + unroll 8 + code motion */
char rotate_block16_unroll8_descr[] =
    "rotate_block16_unroll8: 16x16 block, unroll 8, code-motion";

void rotate_block16_unroll8(int dim, pixel * RESTRICT src, pixel * RESTRICT dst)
{
    const int BLK = 16;  /* block size */
    int i0, j0, i, j;

    for (i0 = 0; i0 < dim; i0 += BLK) {
        for (j0 = 0; j0 < dim; j0 += BLK) {
            for (i = i0; i < i0 + BLK; i++) {
                int src_row = i * dim + j0;
                int dst_col = i;

                for (j = 0; j < BLK; j += 8) {
                    int jj = j0 + j;
                    int dst_base = (dim - 1 - jj) * dim + dst_col;

                    dst[dst_base]           = src[src_row + j];
                    dst[dst_base -     dim] = src[src_row + j + 1];
                    dst[dst_base - 2 * dim] = src[src_row + j + 2];
                    dst[dst_base - 3 * dim] = src[src_row + j + 3];
                    dst[dst_base - 4 * dim] = src[src_row + j + 4];
                    dst[dst_base - 5 * dim] = src[src_row + j + 5];
                    dst[dst_base - 6 * dim] = src[src_row + j + 6];
                    dst[dst_base - 7 * dim] = src[src_row + j + 7];
                }
            }
        }
    }
}


/* Rotate via transpose + reverse rows (uses a temporary buffer) */
char rotate_transpose_reverse_descr[] =
    "rotate_transpose_reverse: 2-step (transpose + reverse rows)";

void rotate_transpose_reverse(int dim, pixel * RESTRICT src, pixel * RESTRICT dst)
{
    /* allocate temporary buffer on heap */
    pixel *tmp = (pixel *) malloc(dim * dim * sizeof(pixel));
    if (!tmp) return; /* fallback if malloc fails */

    int i, j;

    /* Step 1: transpose src -> tmp */
    for (i = 0; i < dim; i++) {
        int src_row = i * dim;
        for (j = 0; j < dim; j++) {
            tmp[j * dim + i] = src[src_row + j];
        }
    }

    /* Step 2: reverse rows of tmp -> dst */
    for (i = 0; i < dim; i++) {
        int dst_row = i * dim;
        int tmp_row = (dim - 1 - i) * dim;
        for (j = 0; j < dim; j++) {
            dst[dst_row + j] = tmp[tmp_row + j];
        }
    }

    free(tmp);
}

char rotate_block32_strength_descr[] =
    "rotate_block32_strength: blocking + pointer strength reduction";

void rotate_block32_strength(int dim, pixel * RESTRICT src, pixel * RESTRICT dst)
{
    const int BLK = 32;
    int i0, j0, i, j;

    for (i0 = 0; i0 < dim; i0 += BLK) {
        for (j0 = 0; j0 < dim; j0 += BLK) {
            for (i = i0; i < i0 + BLK; i++) {
                pixel *s = src + i * dim + j0;                   /* start of src row */
                pixel *d = dst + (dim - 1 - j0) * dim + i;       /* start of dst col */

                for (j = 0; j < BLK; j++) {
                    *d = *s;     /* copy pixel */
                    s++;         /* next src pixel (row-major) */
                    d -= dim;    /* next dst pixel (col-major downwards) */
                }
            }
        }
    }
}

/* 
 * rotate - Your current working version of rotate
 * IMPORTANT: This is the version you will be graded on
 */
char rotate_descr[] = "rotate: Current working version";
void rotate(int dim, pixel *src, pixel *dst) 
{
    naive_rotate(dim, src, dst);
}

/*********************************************************************
 * register_rotate_functions - Register all of your different versions
 *     of the rotate kernel with the driver by calling the
 *     add_rotate_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/
void register_rotate_functions() 
{
    add_rotate_function(&naive_rotate, naive_rotate_descr);   //naive approach
    add_rotate_function(&rotate, rotate_descr);   
    add_rotate_function(&rotate_block16, rotate_block16_descr);
    add_rotate_function(&rotate_block32, rotate_block32_descr);
    add_rotate_function(&rotate_block64, rotate_block64_descr);
    add_rotate_function(&rotate_block32_unroll8, rotate_block32_unroll8_descr);
    add_rotate_function(&rotate_block16_unroll8, rotate_block16_unroll8_descr);
    add_rotate_function(&rotate_transpose_reverse, rotate_transpose_reverse_descr);
    add_rotate_function(&rotate_block32_strength, rotate_block32_strength_descr);



}



/***************
 * SMOOTH KERNEL
 **************/

/***************************************************************
 * Various typedefs and helper functions for the smooth function
 * You may modify these any way you like.
 **************************************************************/

/* A struct used to compute averaged pixel value */
typedef struct {
    int red;
    int green;
    int blue;
    int num;
} pixel_sum;

/* Compute min and max of two integers, respectively */
static int min(int a, int b) { return (a < b ? a : b); }
static int max(int a, int b) { return (a > b ? a : b); }

/* 
 * initialize_pixel_sum - Initializes all fields of sum to 0 
 */
static void initialize_pixel_sum(pixel_sum *sum) 
{
    sum->red = sum->green = sum->blue = 0;
    sum->num = 0;
    return;
}

/* 
 * accumulate_sum - Accumulates field values of p in corresponding 
 * fields of sum 
 */
static void accumulate_sum(pixel_sum *sum, pixel p) 
{
    sum->red += (int) p.red;
    sum->green += (int) p.green;
    sum->blue += (int) p.blue;
    sum->num++;
    return;
}

/* 
 * assign_sum_to_pixel - Computes averaged pixel value in current_pixel 
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum) 
{
    current_pixel->red = (unsigned short) (sum.red/sum.num);
    current_pixel->green = (unsigned short) (sum.green/sum.num);
    current_pixel->blue = (unsigned short) (sum.blue/sum.num);
    return;
}

/* 
 * avg - Returns averaged pixel value at (i,j) 
 */
static pixel avg(int dim, int i, int j, pixel *src) 
{
    int ii, jj;
    pixel_sum sum;
    pixel current_pixel;

    initialize_pixel_sum(&sum);
    for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) 
	for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++) 
	    accumulate_sum(&sum, src[RIDX(ii, jj, dim)]);

    assign_sum_to_pixel(&current_pixel, sum);
    return current_pixel;
}

/******************************************************
 * Your different versions of the smooth kernel go here
 ******************************************************/

/*
 * naive_smooth - The naive baseline version of smooth 
 */
/* Better smooth: inline avg, no function calls, clamped loops */
char smooth_better_descr[] =
    "smooth_better: inlined avg with clamped loops";

void smooth_better(int dim, pixel * RESTRICT src, pixel * RESTRICT dst)
{
    int i, j, ii, jj;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            int r=0, g=0, b=0, count=0;

            /* 3x3 neighborhood with clamped indices */
            for (ii = (i > 0 ? i-1 : i); ii <= (i < dim-1 ? i+1 : i); ii++) {
                for (jj = (j > 0 ? j-1 : j); jj <= (j < dim-1 ? j+1 : j); jj++) {
                    pixel p = src[RIDX(ii, jj, dim)];
                    r += p.red;
                    g += p.green;
                    b += p.blue;
                    count++;
                }
            }

            dst[RIDX(i,j,dim)].red   = r / count;
            dst[RIDX(i,j,dim)].green = g / count;
            dst[RIDX(i,j,dim)].blue  = b / count;
        }
    }
}

char smooth_optimized_descr[] = "smooth_optimized: case-split (corners/edges/interior) + unrolled sums + strength reduction";
void smooth_optimized(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    /* ---------- Corners ---------- */
    {
        /* Top-left corner (0,0) */
        int r = src[0].red + src[1].red + src[dim].red + src[dim+1].red;
        int g = src[0].green + src[1].green + src[dim].green + src[dim+1].green;
        int b = src[0].blue + src[1].blue + src[dim].blue + src[dim+1].blue;
        dst[0].red   = r >> 2;
        dst[0].green = g >> 2;
        dst[0].blue  = b >> 2;

        /* Top-right corner (0,dim-1) */
        int idx = dim-1;
        r = src[idx].red + src[idx-1].red + src[idx+dim].red + src[idx+dim-1].red;
        g = src[idx].green + src[idx-1].green + src[idx+dim].green + src[idx+dim-1].green;
        b = src[idx].blue + src[idx-1].blue + src[idx+dim].blue + src[idx+dim-1].blue;
        dst[idx].red   = r >> 2;
        dst[idx].green = g >> 2;
        dst[idx].blue  = b >> 2;

        /* Bottom-left corner (dim-1,0) */
        idx = (dim-1)*dim;
        r = src[idx].red + src[idx+1].red + src[idx-dim].red + src[idx-dim+1].red;
        g = src[idx].green + src[idx+1].green + src[idx-dim].green + src[idx-dim+1].green;
        b = src[idx].blue + src[idx+1].blue + src[idx-dim].blue + src[idx-dim+1].blue;
        dst[idx].red   = r >> 2;
        dst[idx].green = g >> 2;
        dst[idx].blue  = b >> 2;

        /* Bottom-right corner (dim-1,dim-1) */
        idx = dim*dim - 1;
        r = src[idx].red + src[idx-1].red + src[idx-dim].red + src[idx-dim-1].red;
        g = src[idx].green + src[idx-1].green + src[idx-dim].green + src[idx-dim-1].green;
        b = src[idx].blue + src[idx-1].blue + src[idx-dim].blue + src[idx-dim-1].blue;
        dst[idx].red   = r >> 2;
        dst[idx].green = g >> 2;
        dst[idx].blue  = b >> 2;
    }

    /* ---------- Top and Bottom edges ---------- */
    for (j = 1; j < dim-1; j++) {
        /* Top row (i=0) */
        int r = src[j-1].red + src[j].red + src[j+1].red +
                src[j-1+dim].red + src[j+dim].red + src[j+1+dim].red;
        int g = src[j-1].green + src[j].green + src[j+1].green +
                src[j-1+dim].green + src[j+dim].green + src[j+1+dim].green;
        int b = src[j-1].blue + src[j].blue + src[j+1].blue +
                src[j-1+dim].blue + src[j+dim].blue + src[j+1+dim].blue;
        dst[j].red   = r / 6;
        dst[j].green = g / 6;
        dst[j].blue  = b / 6;

        /* Bottom row (i=dim-1) */
        int idx = (dim-1)*dim + j;
        r = src[idx-1].red + src[idx].red + src[idx+1].red +
            src[idx-1-dim].red + src[idx-dim].red + src[idx+1-dim].red;
        g = src[idx-1].green + src[idx].green + src[idx+1].green +
            src[idx-1-dim].green + src[idx-dim].green + src[idx+1-dim].green;
        b = src[idx-1].blue + src[idx].blue + src[idx+1].blue +
            src[idx-1-dim].blue + src[idx-dim].blue + src[idx+1-dim].blue;
        dst[idx].red   = r / 6;
        dst[idx].green = g / 6;
        dst[idx].blue  = b / 6;
    }

    /* ---------- Left and Right edges ---------- */
    for (i = 1; i < dim-1; i++) {
        int row = i * dim;

        /* Left column (j=0) */
        int r = src[row].red + src[row+1].red + 
                src[row-dim].red + src[row-dim+1].red +
                src[row+dim].red + src[row+dim+1].red;
        int g = src[row].green + src[row+1].green + 
                src[row-dim].green + src[row-dim+1].green +
                src[row+dim].green + src[row+dim+1].green;
        int b = src[row].blue + src[row+1].blue + 
                src[row-dim].blue + src[row-dim+1].blue +
                src[row+dim].blue + src[row+dim+1].blue;
        dst[row].red   = r / 6;
        dst[row].green = g / 6;
        dst[row].blue  = b / 6;

        /* Right column (j=dim-1) */
        int idx = row + dim - 1;
        r = src[idx].red + src[idx-1].red + 
            src[idx-dim].red + src[idx-dim-1].red +
            src[idx+dim].red + src[idx+dim-1].red;
        g = src[idx].green + src[idx-1].green + 
            src[idx-dim].green + src[idx-dim-1].green +
            src[idx+dim].green + src[idx+dim-1].green;
        b = src[idx].blue + src[idx-1].blue + 
            src[idx-dim].blue + src[idx-dim-1].blue +
            src[idx+dim].blue + src[idx+dim-1].blue;
        dst[idx].red   = r / 6;
        dst[idx].green = g / 6;
        dst[idx].blue  = b / 6;
    }

    /* ---------- Interior ---------- */
    for (i = 1; i < dim-1; i++) {
        int row = i * dim;
        for (j = 1; j < dim-1; j++) {
            int idx = row + j;

            int r = src[idx].red + src[idx-1].red + src[idx+1].red +
                    src[idx-dim].red + src[idx-dim-1].red + src[idx-dim+1].red +
                    src[idx+dim].red + src[idx+dim-1].red + src[idx+dim+1].red;
            int g = src[idx].green + src[idx-1].green + src[idx+1].green +
                    src[idx-dim].green + src[idx-dim-1].green + src[idx-dim+1].green +
                    src[idx+dim].green + src[idx+dim-1].green + src[idx+dim+1].green;
            int b = src[idx].blue + src[idx-1].blue + src[idx+1].blue +
                    src[idx-dim].blue + src[idx-dim-1].blue + src[idx-dim+1].blue +
                    src[idx+dim].blue + src[idx+dim-1].blue + src[idx+dim+1].blue;

            dst[idx].red   = r / 9;
            dst[idx].green = g / 9;
            dst[idx].blue  = b / 9;
        }
    }
}
/*
 * smooth_hybrid - Case-split for boundaries + sliding window for interior
 */
char smooth_hybrid_descr[] = "smooth_hybrid: case-split corners/edges + sliding window interior";

void smooth_hybrid(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    /* ---- Handle corners ---- */
    // Top-left corner (0,0)
    {
        int sum_r=0, sum_g=0, sum_b=0;
        int count=0;
        for (int ii=0; ii<=1; ii++)
            for (int jj=0; jj<=1; jj++) {
                pixel p = src[RIDX(ii,jj,dim)];
                sum_r+=p.red; sum_g+=p.green; sum_b+=p.blue; count++;
            }
        dst[RIDX(0,0,dim)].red   = sum_r/count;
        dst[RIDX(0,0,dim)].green = sum_g/count;
        dst[RIDX(0,0,dim)].blue  = sum_b/count;
    }

    // Top-right corner (0,dim-1)
    {
        int sum_r=0, sum_g=0, sum_b=0;
        int count=0;
        for (int ii=0; ii<=1; ii++)
            for (int jj=dim-2; jj<dim; jj++) {
                pixel p = src[RIDX(ii,jj,dim)];
                sum_r+=p.red; sum_g+=p.green; sum_b+=p.blue; count++;
            }
        dst[RIDX(0,dim-1,dim)].red   = sum_r/count;
        dst[RIDX(0,dim-1,dim)].green = sum_g/count;
        dst[RIDX(0,dim-1,dim)].blue  = sum_b/count;
    }

    // Bottom-left corner (dim-1,0)
    {
        int sum_r=0, sum_g=0, sum_b=0;
        int count=0;
        for (int ii=dim-2; ii<dim; ii++)
            for (int jj=0; jj<=1; jj++) {
                pixel p = src[RIDX(ii,jj,dim)];
                sum_r+=p.red; sum_g+=p.green; sum_b+=p.blue; count++;
            }
        dst[RIDX(dim-1,0,dim)].red   = sum_r/count;
        dst[RIDX(dim-1,0,dim)].green = sum_g/count;
        dst[RIDX(dim-1,0,dim)].blue  = sum_b/count;
    }

    // Bottom-right corner (dim-1,dim-1)
    {
        int sum_r=0, sum_g=0, sum_b=0;
        int count=0;
        for (int ii=dim-2; ii<dim; ii++)
            for (int jj=dim-2; jj<dim; jj++) {
                pixel p = src[RIDX(ii,jj,dim)];
                sum_r+=p.red; sum_g+=p.green; sum_b+=p.blue; count++;
            }
        dst[RIDX(dim-1,dim-1,dim)].red   = sum_r/count;
        dst[RIDX(dim-1,dim-1,dim)].green = sum_g/count;
        dst[RIDX(dim-1,dim-1,dim)].blue  = sum_b/count;
    }

    /* ---- Handle edges (without corners) ---- */
    // Top & Bottom rows
    for (j=1; j<dim-1; j++) {
        // Top row (0,j)
        {
            int sum_r=0, sum_g=0, sum_b=0;
            int count=0;
            for (int ii=0; ii<=1; ii++)
                for (int jj=j-1; jj<=j+1; jj++) {
                    pixel p = src[RIDX(ii,jj,dim)];
                    sum_r+=p.red; sum_g+=p.green; sum_b+=p.blue; count++;
                }
            dst[RIDX(0,j,dim)].red   = sum_r/count;
            dst[RIDX(0,j,dim)].green = sum_g/count;
            dst[RIDX(0,j,dim)].blue  = sum_b/count;
        }
        // Bottom row (dim-1,j)
        {
            int sum_r=0, sum_g=0, sum_b=0;
            int count=0;
            for (int ii=dim-2; ii<dim; ii++)
                for (int jj=j-1; jj<=j+1; jj++) {
                    pixel p = src[RIDX(ii,jj,dim)];
                    sum_r+=p.red; sum_g+=p.green; sum_b+=p.blue; count++;
                }
            dst[RIDX(dim-1,j,dim)].red   = sum_r/count;
            dst[RIDX(dim-1,j,dim)].green = sum_g/count;
            dst[RIDX(dim-1,j,dim)].blue  = sum_b/count;
        }
    }

    // Left & Right columns
    for (i=1; i<dim-1; i++) {
        // Left col (i,0)
        {
            int sum_r=0, sum_g=0, sum_b=0;
            int count=0;
            for (int ii=i-1; ii<=i+1; ii++)
                for (int jj=0; jj<=1; jj++) {
                    pixel p = src[RIDX(ii,jj,dim)];
                    sum_r+=p.red; sum_g+=p.green; sum_b+=p.blue; count++;
                }
            dst[RIDX(i,0,dim)].red   = sum_r/count;
            dst[RIDX(i,0,dim)].green = sum_g/count;
            dst[RIDX(i,0,dim)].blue  = sum_b/count;
        }
        // Right col (i,dim-1)
        {
            int sum_r=0, sum_g=0, sum_b=0;
            int count=0;
            for (int ii=i-1; ii<=i+1; ii++)
                for (int jj=dim-2; jj<dim; jj++) {
                    pixel p = src[RIDX(ii,jj,dim)];
                    sum_r+=p.red; sum_g+=p.green; sum_b+=p.blue; count++;
                }
            dst[RIDX(i,dim-1,dim)].red   = sum_r/count;
            dst[RIDX(i,dim-1,dim)].green = sum_g/count;
            dst[RIDX(i,dim-1,dim)].blue  = sum_b/count;
        }
    }

    /* ---- Handle interior with sliding window ---- */
    for (i=1; i<dim-1; i++) {
        int sum_r=0, sum_g=0, sum_b=0;

        /* Initialize first window (i-1..i+1 , 0..2) */
        for (int ii=i-1; ii<=i+1; ii++)
            for (int jj=0; jj<=2; jj++) {
                pixel p = src[RIDX(ii,jj,dim)];
                sum_r+=p.red; sum_g+=p.green; sum_b+=p.blue;
            }

        dst[RIDX(i,1,dim)].red   = sum_r/9;
        dst[RIDX(i,1,dim)].green = sum_g/9;
        dst[RIDX(i,1,dim)].blue  = sum_b/9;

        /* Slide window across row */
        for (j=2; j<dim-1; j++) {
            // remove left col
            for (int ii=i-1; ii<=i+1; ii++) {
                pixel p = src[RIDX(ii,j-2,dim)];
                sum_r-=p.red; sum_g-=p.green; sum_b-=p.blue;
            }
            // add right col
            for (int ii=i-1; ii<=i+1; ii++) {
                pixel p = src[RIDX(ii,j+1,dim)];
                sum_r+=p.red; sum_g+=p.green; sum_b+=p.blue;
            }
            dst[RIDX(i,j,dim)].red   = sum_r/9;
            dst[RIDX(i,j,dim)].green = sum_g/9;
            dst[RIDX(i,j,dim)].blue  = sum_b/9;
        }
    }
}


char naive_smooth_descr[] = "naive_smooth: Naive baseline implementation";
void naive_smooth(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(i, j, dim)] = avg(dim, i, j, src);
}

/*
 * smooth - Your current working version of smooth. 
 * IMPORTANT: This is the version you will be graded on
 */
char smooth_descr[] = "smooth: Current working version";
void smooth(int dim, pixel *src, pixel *dst) 
{
    naive_smooth(dim, src, dst);
}


/********************************************************************* 
 * register_smooth_functions - Register all of your different versions
 *     of the smooth kernel with the driver by calling the
 *     add_smooth_function() for each test function.  When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/


void register_smooth_functions() {
    add_smooth_function(&smooth, smooth_descr);
    add_smooth_function(&naive_smooth, naive_smooth_descr);
    add_smooth_function(&smooth_better, smooth_better_descr); 
    add_smooth_function(&smooth_optimized, smooth_optimized_descr);
    add_smooth_function(&smooth_hybrid, smooth_hybrid_descr);


    /* ... Register additional test functions here */
}

