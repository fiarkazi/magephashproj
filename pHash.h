#ifndef _PHASH_H
#define _PHASH_H


#include <pHash-config.h>
#include <limits.h>
#include <math.h>
#include <dirent.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define __STDC_CONSTANT_MACROS

#include <stdint.h>

#if defined(HAVE_IMAGE_HASH)
#define cimg_debug 0
#define cimg_display 0
#include "CImg.h"
using namespace cimg_library;
#endif


using namespace std;

#define SQRT_TWO 1.4142135623730950488016887242097

#ifndef ULLONG_MAX
#define ULLONG_MAX 18446744073709551615ULL
#endif

#define ROUNDING_FACTOR(x) (((x) >= 0) ? 0.5 : -0.5) 

#if defined( _MSC_VER) || defined(_BORLANDC_)
typedef unsigned _uint64 ulong64;
typedef signed _int64 long64;
#else
typedef unsigned long long ulong64;
typedef signed long long long64;
#endif

#ifdef __cplusplus
extern "C" {
#endif

const int MaxFileSize = (1<<30); /* 1GB file size limit (for mvp files) */
const off_t HeaderSize = 64;     /* header size for mvp file */


typedef struct ph_file_offset {
    off_t offset;
    uint8_t fileno;
} FileIndex;


/* structure for a single hash */
typedef struct ph_datapoint {
    char *id;
    void *hash;
    float *path;
    uint32_t hash_length;
    uint8_t hash_type;
}DP;

typedef struct ph_slice {
    DP **hash_p;
    int n;
    void *hash_params;
} slice;

struct BinHash 
{
	uint8_t *hash;
	uint32_t bytelength;
	uint32_t byteidx; // used by addbit()
	uint8_t bitmask;  // used by addbit()

	/*
	 * add a single bit to hash. the bits are 
	 * written from left to right.
	 */
	int addbit(uint8_t bit)
	{
		if (bitmask == 0) 
		{
			bitmask = 128; // reset bitmask to "10000000"
			byteidx++;     // jump to next byte in array
		}

		if (byteidx >= bytelength) return -1;
		
		if (bit == 1) *(hash + byteidx) |= bitmask;
		bitmask >>=1;
		return 0;
	}	
};

BinHash* _ph_bmb_new(uint32_t bytelength);
void ph_bmb_free(BinHash *binHash);

/*! /brief Radon Projection info
 */
#ifdef HAVE_IMAGE_HASH
typedef struct ph_projections {
    CImg<uint8_t> *R;           //contains projections of image of angled lines through center
    int *nb_pix_perline;        //the head of int array denoting the number of pixels of each line
    int size;                   //the size of nb_pix_perline
}Projections;
#endif

/*! /brief feature vector info
 */
typedef struct ph_feature_vector {
    double *features;           //the head of the feature array of double's
    int size;                   //the size of the feature array
}Features;

/*! /brief Digest info
 */
typedef struct ph_digest {
    char *id;                   //hash id
    uint8_t *coeffs;            //the head of the digest integer coefficient array
    int size;                   //the size of the coeff array
} Digest;





/* /brief alloc a single data point
 *  allocates path array, does nto set id or path
 */
 DP* ph_malloc_datapoint(int hashtype);

/** /brief free a datapoint and its path
 *
 */
void ph_free_datapoint(DP *dp);


/*! /brief radon function
 *  Find radon projections of N lines running through the image center for lines angled 0
 *  to 180 degrees from horizontal.
 *  /param img - CImg src image
 *  /param  N  - int number of angled lines to consider.
 *  /param  projs - (out) Projections struct 
 *  /return int value - less than 0 for error
 */
#ifdef HAVE_IMAGE_HASH
int ph_radon_projections(const CImg<uint8_t> &img,int N,Projections &projs);

/*! /brief feature vector
 *         compute the feature vector from a radon projection map.
 *  /param  projs - Projections struct
 *  /param  fv    - (out) Features struct
 *  /return int value - less than 0 for error
*/
int ph_feature_vector(const Projections &projs,Features &fv);

/*! /brief dct 
 *  Compute the dct of a given vector
 *  /param R - vector of input series
 *  /param D - (out) the dct of R
 *  /return  int value - less than 0 for error
*/
int ph_dct(const Features &fv, Digest &digest);

/*! /brief cross correlation for 2 series
 *  Compute the cross correlation of two series vectors
 *  /param x - Digest struct
 *  /param y - Digest struct
 *  /param pcc - double value the peak of cross correlation
 *  /param threshold - double value for the threshold value for which 2 images
 *                     are considered the same or different.
 *  /return - int value - 1 (true) for same, 0 (false) for different, < 0 for error
 */

int ph_crosscorr(const Digest &x,const Digest &y,double &pcc, double threshold = 0.90);

/*! /brief image digest
 *  Compute the image digest for an image given the input image
 *  /param img - CImg object representing an input image
 *  /param sigma - double value for the deviation for a gaussian filter function 
 *  /param gamma - double value for gamma correction on the input image
 *  /param digest - (out) Digest struct
 *  /param N      - int value for the number of angles to consider. 
 *  /return       - less than 0 for error
 */
int _ph_image_digest(const CImg<uint8_t> &img,double sigma, double gamma,Digest &digest,int N=180);

/*! /brief image digest
 *  Compute the image digest given the file name.
 *  /param file - string value for file name of input image.
 *  /param sigma - double value for the deviation for gaussian filter
 *  /param gamma - double value for gamma correction on the input image.
 *  /param digest - Digest struct
 *  /param N      - int value for number of angles to consider
 */
int ph_image_digest(const char *file, double sigma, double gamma, Digest &digest,int N=180);


/*! /brief compare 2 images
 *  /param imA - CImg object of first image 
 *  /param imB - CImg object of second image
 *  /param pcc   - (out) double value for peak of cross correlation
 *  /param sigma - double value for the deviation of gaussian filter
 *  /param gamma - double value for gamma correction of images
 *  /param N     - int number for the number of angles of radon projections
 *  /param theshold - double value for the threshold
 *  /return int 0 (false) for different images, 1 (true) for same image, less than 0 for error
 */
int _ph_compare_images(const CImg<uint8_t> &imA,const CImg<uint8_t> &imB,double &pcc, double sigma = 3.5, double gamma = 1.0,int N=180,double threshold=0.90);

/*! /brief compare 2 images
 *  Compare 2 images given the file names
 *  /param file1 - char string of first image file
 *  /param file2 - char string of second image file
 *  /param pcc   - (out) double value for peak of cross correlation
 *  /param sigma - double value for deviation of gaussian filter
 *  /param gamma - double value for gamma correction of images
 *  /param N     - int number for number of angles
 *  /return int 0 (false) for different image, 1 (true) for same images, less than 0 for error
 */
int ph_compare_images(const char *file1, const char *file2,double &pcc, double sigma = 3.5, double gamma=1.0, int N=180,double threshold=0.90);

/*! /brief return dct matrix, C
 *  Return DCT matrix of sqare size, N
 *  /param N - int denoting the size of the square matrix to create.
 *  /return CImg<double> size NxN containing the dct matrix
 */
static CImg<float>* ph_dct_matrix(const int N);

/*! /brief compute dct robust image hash
 *  /param file string variable for name of file
 *  /param hash of type ulong64 (must be 64-bit variable)
 *  /return int value - -1 for failure, 1 for success
 */
int ph_dct_imagehash(const char* file,ulong64 &hash);

int ph_bmb_imagehash(const char *file, uint8_t method, BinHash **ret_hash);
#endif


#ifdef HAVE_IMAGE_HASH
int ph_hamming_distance(const ulong64 hash1,const ulong64 hash2);

/** /brief create a list of datapoint's directly from a directory of image files
 *  /param dirname - path and name of directory containg all image file names
 *  /param capacity - int value for upper limit on number of hashes
 *  /param count - number of hashes created (out param)
 *  /return pointer to a list of DP pointers (NULL for error)
 */

DP** ph_read_imagehashes(const char *dirname,int capacity, int &count);

/** /brief create MH image hash for filename image
*   /param filename - string name of image file
*   /param N - (out) int value for length of image hash returned
*   /param alpha - int scale factor for marr wavelet (default=2)
*   /param lvl   - int level of scale factor (default = 1)
*   /return uint8_t array
**/
uint8_t* ph_mh_imagehash(const char *filename, int &N, float alpha=2.0f, float lvl = 1.0f);
#endif
/** /brief count number bits set in given byte
*   /param val - uint8_t byte value
*   /return int value for number of bits set
**/
int ph_bitcount8(uint8_t val);

/** /brief compute hamming distance between two byte arrays
 *  /param hashA - byte array for first hash
 *  /param lenA - int length of hashA 
 *  /param hashB - byte array for second hash
 *  /param lenB - int length of hashB
 *  /return double value for normalized hamming distance
 **/
double ph_hammingdistance2(uint8_t *hashA, int lenA, uint8_t *hashB, int lenB);

/** /brief get all the filenames in specified directory
 *  /param dirname - string value for path and filename
 *  /param cap - int value for upper limit to number of files
 *  /param count - int value for number of file names returned
 *  /return array of pointers to string file names (NULL for error)
 **/

char** ph_readfilenames(const char *dirname,int &count);


/** /brief textual hash for file
 *  /param filename - char* name of file
 *  /param nbpoints - int length of array of return value (out)
 *  /return TxtHashPoint* array of hash points with respective index into file.
 **/

#ifdef __cplusplus
}
#endif

#endif
