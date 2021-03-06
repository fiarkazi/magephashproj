#include "pHash.h"
#include <stdio.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <algorithm>

#ifdef HAVE_PTHREAD
#include <pthread.h>

int ph_num_threads()
{
	int numCPU = 1;
#ifdef __GLIBC__
		numCPU = sysconf( _SC_NPROCESSORS_ONLN );
#else
		int mib[2];
		size_t len; 

		mib[0] = CTL_HW;
		mib[1] = HW_AVAILCPU;

		sysctl(mib, 2, &numCPU, &len, NULL, 0);

		if( numCPU < 1 ) 
		{
     			mib[1] = HW_NCPU;
     			sysctl( mib, 2, &numCPU, &len, NULL, 0 );

		     	if( numCPU < 1 )
     			{
          			numCPU = 1;
     			}
		}

#endif
	return numCPU;
}
#endif

#ifdef HAVE_IMAGE_HASH
int ph_radon_projections(const CImg<uint8_t> &img,int N,Projections &projs){

    int width = img.width();
    int height = img.height();
    int D = (width > height)?width:height;
    float x_center = (float)width/2;
    float y_center = (float)height/2;
    int x_off = (int)std::floor(x_center + ROUNDING_FACTOR(x_center));
    int y_off = (int)std::floor(y_center + ROUNDING_FACTOR(y_center));

    projs.R = new CImg<uint8_t>(N,D,1,1,0);
    projs.nb_pix_perline = (int*)calloc(N,sizeof(int));

    if (!projs.R || !projs.nb_pix_perline)
	return EXIT_FAILURE;

    projs.size = N;

    CImg<uint8_t> *ptr_radon_map = projs.R;
    int *nb_per_line = projs.nb_pix_perline;

    for (int k=0;k<N/4+1;k++){
        double theta = k*cimg::PI/N;
        double alpha = std::tan(theta);
        for (int x=0;x < D;x++){
	    double y = alpha*(x-x_off);
            int yd = (int)std::floor(y + ROUNDING_FACTOR(y));
            if ((yd + y_off >= 0)&&(yd + y_off < height) && (x < width)){
		*ptr_radon_map->data(k,x) = img(x,yd + y_off);
                nb_per_line[k] += 1;
	    }
            if ((yd + x_off >= 0) && (yd + x_off < width) && (k != N/4) && (x < height)){
		*ptr_radon_map->data(N/2-k,x) = img(yd + x_off,x);
                nb_per_line[N/2-k] += 1;
	    }
	}
    }
    int j= 0;
    for (int k=3*N/4;k<N;k++){
		double theta = k*cimg::PI/N;
			double alpha = std::tan(theta);
			for (int x=0;x < D;x++){
			double y = alpha*(x-x_off);
				int yd = (int)std::floor(y + ROUNDING_FACTOR(y));
				if ((yd + y_off >= 0)&&(yd + y_off < height) && (x < width))
				{
			*ptr_radon_map->data(k,x) = img(x,yd + y_off);
					nb_per_line[k] += 1;
				}
				if ((y_off - yd >= 0)&&(y_off - yd<width)&&(2*y_off-x>=0)&&(2*y_off-x<height)&&(k!=3*N/4))
				{
			*ptr_radon_map->data(k-j,x) = img(-yd+y_off,-(x-y_off)+y_off);
					nb_per_line[k-j] += 1;
				}
            
	}
        j += 2;
    }

    return EXIT_SUCCESS;

}
int ph_feature_vector(const Projections &projs, Features &fv)
{

    CImg<uint8_t> *ptr_map = projs.R;
    CImg<uint8_t> projection_map = *ptr_map;
    int *nb_perline = projs.nb_pix_perline;
    int N = projs.size;
    int D = projection_map.height();

    fv.features = (double*)malloc(N*sizeof(double));
    fv.size = N;
    if (!fv.features)
	return EXIT_FAILURE;

    double *feat_v = fv.features;
    double sum = 0.0;
    double sum_sqd = 0.0;
    for (int k=0; k < N; k++){
	double line_sum = 0.0;
        double line_sum_sqd = 0.0;
        int nb_pixels = nb_perline[k];
	for (int i=0;i<D;i++){
	    line_sum += projection_map(k,i);
    	    line_sum_sqd += projection_map(k,i)*projection_map(k,i);
	}
	feat_v[k] = (line_sum_sqd/nb_pixels) - (line_sum*line_sum)/(nb_pixels*nb_pixels);
        sum += feat_v[k];
        sum_sqd += feat_v[k]*feat_v[k];
    }
    double mean = sum/N;
    double var  = sqrt((sum_sqd/N) - (sum*sum)/(N*N));

    for (int i=0;i<N;i++){
    	feat_v[i] = (feat_v[i] - mean)/var;
    }

    return EXIT_SUCCESS;
} 
int ph_dct(const Features &fv,Digest &digest)
{
    int N = fv.size;
    const int nb_coeffs = 40;

    digest.coeffs = (uint8_t*)malloc(nb_coeffs*sizeof(uint8_t));
    if (!digest.coeffs)
	return EXIT_FAILURE;

    digest.size = nb_coeffs;

    double *R = fv.features;

    uint8_t *D = digest.coeffs;

    double D_temp[nb_coeffs];
    double max = 0.0;
    double min = 0.0;
    for (int k = 0;k<nb_coeffs;k++){
	double sum = 0.0;
        for (int n=0;n<N;n++){
	    double temp = R[n]*cos((cimg::PI*(2*n+1)*k)/(2*N));
            sum += temp;
	}
        if (k == 0)
	    D_temp[k] = sum/sqrt((double)N);
        else
            D_temp[k] = sum*SQRT_TWO/sqrt((double)N);
        if (D_temp[k] > max)
            max = D_temp[k];
        if (D_temp[k] < min)
            min = D_temp[k];
    }
       
    for (int i=0;i<nb_coeffs;i++){

	D[i] = (uint8_t)(UCHAR_MAX*(D_temp[i] - min)/(max - min));

    }
    
    return EXIT_SUCCESS;
}

int ph_crosscorr(const Digest &x,const Digest &y,double &pcc,double threshold){

    int N = y.size;
    int result = 0;

    uint8_t *x_coeffs = x.coeffs;
    uint8_t *y_coeffs = y.coeffs;

    double *r = new double[N];
    double sumx = 0.0;
    double sumy = 0.0;
    for (int i=0;i < N;i++){
	sumx += x_coeffs[i];
        sumy += y_coeffs[i];
    }
    double meanx = sumx/N;
    double meany = sumy/N;
    double max = 0;
    for (int d=0;d<N;d++){
        double num = 0.0;
        double denx = 0.0;
        double deny = 0.0;
	for (int i=0;i<N;i++){
	    num  += (x_coeffs[i]-meanx)*(y_coeffs[(N+i-d)%N]-meany);
            denx += pow((x_coeffs[i]-meanx),2);
            deny += pow((y_coeffs[(N+i-d)%N]-meany),2);
	}
        r[d] = num/sqrt(denx*deny);
        if (r[d] > max)
	    max = r[d];
    }
    delete[] r;
    pcc = max;
    if (max > threshold)
	    result = 1;

    return result;
}

#ifdef max
#undef max
#endif

int _ph_image_digest(const CImg<uint8_t> &img,double sigma, double gamma,Digest &digest, int N){
    
    int result = EXIT_FAILURE;
    CImg<uint8_t> graysc;
    if (img.spectrum() >= 3){
	graysc = img.get_RGBtoYCbCr().channel(0);
    }
    else if (img.spectrum() == 1){
	graysc = img;
    }
    else {
	return result;
    }
	
 
    graysc.blur((float)sigma);
 
    (graysc/graysc.max()).pow(gamma);
     
    Projections projs;
    if (ph_radon_projections(graysc,N,projs) < 0)
	goto cleanup;
 
    Features features;
    if (ph_feature_vector(projs,features) < 0)
	goto cleanup;
    
    if (ph_dct(features,digest) < 0)
        goto cleanup;
 
    result = EXIT_SUCCESS;

cleanup:
    free(projs.nb_pix_perline);
    free(features.features);

    delete projs.R;
    return result;
}

#define max(a,b) (((a)>(b))?(a):(b))

int ph_image_digest(const char *file, double sigma, double gamma, Digest &digest, int N){
    
    CImg<uint8_t> *src = new CImg<uint8_t>(file);
	int res = -1;
	if(src)
	{
    		int result = _ph_image_digest(*src,sigma,gamma,digest,N);
		delete src;
    		res = result;
	}
	return res;
}

int _ph_compare_images(const CImg<uint8_t> &imA,const CImg<uint8_t> &imB,double &pcc, double sigma, double gamma,int N,double threshold){

    int result = 0;
    Digest digestA;
    if (_ph_image_digest(imA,sigma,gamma,digestA,N) < 0)
	goto cleanup;

    Digest digestB;
    if (_ph_image_digest(imB,sigma,gamma,digestB,N) < 0)
	goto cleanup;

    if (ph_crosscorr(digestA,digestB,pcc,threshold) < 0)
	goto cleanup;

    if  (pcc  > threshold)
        result = 1;

cleanup:

    free(digestA.coeffs);
    free(digestB.coeffs);
    return result;
}

int ph_compare_images(const char *file1, const char *file2,double &pcc, double sigma, double gamma, int N,double threshold){

    CImg<uint8_t> *imA = new CImg<uint8_t>(file1);
    CImg<uint8_t> *imB = new CImg<uint8_t>(file2);
    
    int res = _ph_compare_images(*imA,*imB,pcc,sigma,gamma,N,threshold);

    delete imA;
    delete imB;
    return res;
}

CImg<float>* ph_dct_matrix(const int N){
    CImg<float> *ptr_matrix = new CImg<float>(N,N,1,1,1/sqrt((float)N));
    const float c1 = sqrt(2.0/N); 
    for (int x=0;x<N;x++){
	for (int y=1;y<N;y++){
	    *ptr_matrix->data(x,y) = c1*cos((cimg::PI/2/N)*y*(2*x+1));
	}
    }
    return ptr_matrix;
}

int ph_dct_imagehash(const char* file,ulong64 &hash){

    if (!file){
	return -1;
    }
    CImg<uint8_t> src;
    try {
	src.load(file);
    } catch (CImgIOException ex){
	return -1;
    }
    CImg<float> meanfilter(7,7,1,1,1);
    CImg<float> img;
    if (src.spectrum() == 3){
        img = src.RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else if (src.spectrum() == 4){
	int width = img.width();
        int height = img.height();
        int depth = img.depth();
	img = src.crop(0,0,0,0,width-1,height-1,depth-1,2).RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else {
	img = src.channel(0).get_convolve(meanfilter);
    }

    img.resize(32,32);
    CImg<float> *C  = ph_dct_matrix(32);
    CImg<float> Ctransp = C->get_transpose();

    CImg<float> dctImage = (*C)*img*Ctransp;

    CImg<float> subsec = dctImage.crop(1,1,8,8).unroll('x');;
   
    float median = subsec.median();
    ulong64 one = 0x0000000000000001;
    hash = 0x0000000000000000;
    for (int i=0;i< 64;i++){
	float current = subsec(i);
        if (current > median)
	    hash |= one;
	one = one << 1;
    }
  
    delete C;

    return 0;
}

#ifdef HAVE_PTHREAD
void *ph_image_thread(void *p)
{
        slice *s = (slice *)p;
        for(int i = 0; i < s->n; ++i)
        {
                DP *dp = (DP *)s->hash_p[i];
		ulong64 hash;
		int ret = ph_dct_imagehash(dp->id, hash);
                dp->hash = (ulong64*)malloc(sizeof(hash));
		memcpy(dp->hash, &hash, sizeof(hash));
                dp->hash_length = 1;
        }
}

DP** ph_dct_image_hashes(char *files[], int count, int threads)
{
       	if(!files || count <= 0)
                return NULL;

        int num_threads;
        if(threads > count)
        {
                num_threads = count;
        }
	else if(threads > 0)
        {
                num_threads = threads;
        }
	else
	{
                num_threads = ph_num_threads();
        }

	DP **hashes = (DP**)malloc(count*sizeof(DP*));

        for(int i = 0; i < count; ++i)
        {
                hashes[i] = (DP *)malloc(sizeof(DP));
                hashes[i]->id = strdup(files[i]);
	}

	pthread_t thds[num_threads];

        int rem = count % num_threads;
        int start = 0;
        int off = 0;
        slice *s = new slice[num_threads];
        for(int n = 0; n < num_threads; ++n)
        {
                off = (int)floor((count/(float)num_threads) + (rem>0?num_threads-(count % num_threads):0));

                s[n].hash_p = &hashes[start];
                s[n].n = off;
                s[n].hash_params = NULL;
                start += off;
                --rem;
                pthread_create(&thds[n], NULL, ph_image_thread, &s[n]);
        }
	for(int i = 0; i < num_threads; ++i)
        {
                pthread_join(thds[i], NULL);
        }
	delete[] s;

        return hashes;

}
#endif


#endif


int ph_hamming_distance(const ulong64 hash1,const ulong64 hash2){
    ulong64 x = hash1^hash2;
    const ulong64 m1  = 0x5555555555555555ULL;
    const ulong64 m2  = 0x3333333333333333ULL;
    const ulong64 h01 = 0x0101010101010101ULL;
    const ulong64 m4  = 0x0f0f0f0f0f0f0f0fULL;
    x -= (x >> 1) & m1;
    x = (x & m2) + ((x >> 2) & m2);
    x = (x + (x >> 4)) & m4;
    return (x * h01)>>56;
}

DP* ph_malloc_datapoint(int hashtype){
    DP* dp = (DP*)malloc(sizeof(DP));
    dp->hash = NULL;
    dp->id = NULL;
    dp->hash_type = hashtype;
    return dp;
}
void ph_free_datapoint(DP *dp){
    if (!dp){
	return;
    }
    free(dp);
    dp=NULL;
    return;
}


#ifdef HAVE_IMAGE_HASH
DP** ph_read_imagehashes(const char *dirname,int pathlength, int &count){

    count = 0;
    struct dirent *dir_entry;
    DIR *dir = opendir(dirname);
    if (!dir)
        return NULL;

    while ((dir_entry = readdir(dir)) != 0){
	if (strcmp(dir_entry->d_name,".")&& strcmp(dir_entry->d_name,"..")){
	    count++;
	}
    }
  
    DP **hashlist = (DP**)malloc(count*sizeof(DP**));
    if (!hashlist)
    {
        closedir(dir);
        return NULL;
    }

    DP *dp = NULL;
    int index = 0;
    errno = 0;
    ulong64 tmphash = 0;
    char path[100];
    path[0] = '\0';
    rewinddir(dir);
    while ((dir_entry = readdir(dir)) != 0){
	if (strcmp(dir_entry->d_name,".") && strcmp(dir_entry->d_name,"..")){
	    strcat(path, dirname);
	    strcat(path, "/");
	    strcat(path, dir_entry->d_name);
	    if (ph_dct_imagehash(path, tmphash) < 0)  //calculate the hash
		continue;
	    dp = ph_malloc_datapoint(UINT64ARRAY);
	    dp->id = strdup(path);
	    dp->hash = (void*)&tmphash;
	    dp->hash_length = 1;
	    hashlist[index++] = dp;
	}
	errno = 0;
    path[0]='\0';
    }
        
    closedir(dir);
    return hashlist;

}

CImg<float>* GetMHKernel(float alpha, float level){
    int sigma = (int)4*pow((float)alpha,(float)level);
    static CImg<float> *pkernel = NULL;
    float xpos, ypos, A;
    if (!pkernel){
	pkernel = new CImg<float>(2*sigma+1,2*sigma+1,1,1,0);
        cimg_forXY(*pkernel,X,Y){
	    xpos = pow(alpha,-level)*(X-sigma);
            ypos = pow(alpha,-level)*(Y-sigma);
            A = xpos*xpos + ypos*ypos;
            pkernel->atXY(X,Y) = (2-A)*exp(-A/2);
	}
    }
    return pkernel;
}

uint8_t* ph_mh_imagehash(const char *filename, int &N,float alpha, float lvl){
    if (filename == NULL){
	return NULL;
    }
    uint8_t *hash = (unsigned char*)malloc(72*sizeof(uint8_t));
    N = 72;

    CImg<uint8_t> src(filename);
    CImg<uint8_t> img;

    if (src.spectrum() == 3){
	img = src.get_RGBtoYCbCr().channel(0).blur(1.0).resize(512,512,1,1,5).get_equalize(256);
    } else{
	img = src.channel(0).get_blur(1.0).resize(512,512,1,1,5).get_equalize(256);
    }
    src.clear();

    CImg<float> *pkernel = GetMHKernel(alpha,lvl);
    CImg<float> fresp =  img.get_correlate(*pkernel);
    img.clear();
    fresp.normalize(0,1.0);
    CImg<float> blocks(31,31,1,1,0);
    for (int rindex=0;rindex < 31;rindex++){
		for (int cindex=0;cindex < 31;cindex++){
			blocks(rindex,cindex) = fresp.get_crop(rindex*16,cindex*16,rindex*16+16-1,cindex*16+16-1).sum();
		}
    }
    int hash_index;
    int nb_ones = 0, nb_zeros = 0;
    int bit_index = 0;
    unsigned char hashbyte = 0;
    for (int rindex=0;rindex < 31-2;rindex+=4){
		CImg<float> subsec;
		for (int cindex=0;cindex < 31-2;cindex+=4){
			subsec = blocks.get_crop(cindex,rindex, cindex+2, rindex+2).unroll('x');
			float ave = subsec.mean();
			cimg_forX(subsec, I){
				hashbyte <<= 1;
				if (subsec(I) > ave){
					hashbyte |= 0x01;
					nb_ones++;
				} else {
					nb_zeros++;
				}
				bit_index++;
				if ((bit_index%8) == 0){
					hash_index = (int)(bit_index/8) - 1; 
					hash[hash_index] = hashbyte;
					hashbyte = 0x00;
				}
			}
		}
	}

    return hash;
}
#endif

char** ph_readfilenames(const char *dirname,int &count){
    count = 0;
    struct dirent *dir_entry;
    DIR *dir = opendir(dirname);
    if (!dir)
        return NULL;

    
    while ((dir_entry = readdir(dir)) != NULL){
	if (strcmp(dir_entry->d_name, ".") && strcmp(dir_entry->d_name,".."))
	    count++;
    }
    
    
    char **files = (char**)malloc(count*sizeof(*files));
    if (!files)
	return NULL;

    errno = 0;
    int index = 0;
    char path[1024];
    path[0] = '\0';
    rewinddir(dir);
    while ((dir_entry = readdir(dir)) != 0){
	if (strcmp(dir_entry->d_name,".") && strcmp(dir_entry->d_name,"..")){
	    strcat(path, dirname);
	    strcat(path, "/");
	    strcat(path, dir_entry->d_name);
	    files[index++] = strdup(path);
	}
        path[0]='\0';
    }
    if (errno)
	return NULL;
    closedir(dir);
    return files;
}


int ph_bitcount8(uint8_t val){
    int num = 0;
    while (val){
	++num;
	val &= val - 1;
    }
    return num;
}



double ph_hammingdistance2(uint8_t *hashA, int lenA, uint8_t *hashB, int lenB){
    if (lenA != lenB){
	return -1.0;
    }
    if ((hashA == NULL) || (hashB == NULL) || (lenA <= 0)){
	return -1.0;
    }
    double dist = 0;
    uint8_t D = 0;
    for (int i=0;i<lenA;i++){
	D = hashA[i]^hashB[i];
	dist = dist + (double)ph_bitcount8(D);
    }
    double bits = (double)lenA*8;
    return dist/bits;

}

/*
	thx stackoverflow :^)
*/
using namespace std;

#define TRUE 1
#define FALSE 0

//структура для хэша и айди иоз
struct ph_imagepoint{
    ulong64 hash;
    char *id;
};

//aux function to create imagepoint data struct
struct ph_imagepoint* ph_malloc_imagepoint(){

    return (struct ph_imagepoint*)malloc(sizeof(struct ph_imagepoint));

}

//Дополнительная функция для сортировки списка хэшей
bool cmp_lt_imp(struct ph_imagepoint dpa, struct ph_imagepoint dpb){
    int result = strcmp(dpa.id, dpb.id);
    if (result < 0)
	return TRUE;
    return FALSE;
}
\
int main(int argc, char **argv){


    if (argc < 3){
	printf("no input args\n");
	printf("expected: \"test_imagephash [dir name] [dir_name]\"\n");
	exit(1);
    }
    const char *dir_name = argv[1];
    const char *dir_name2 = argv[2];
    struct dirent *dir_entry;
    vector<ph_imagepoint> hashlist1; //for hashes in first directory
    vector<ph_imagepoint> hashlist2; //for hashes in second directory
    ph_imagepoint *dp = NULL;

    //first directory
    DIR *dir = opendir(dir_name);
    if (!dir){
	printf("unable to open directory\n");
	exit(1);
    }
    errno = 0;
    int i = 0;
    ulong64 tmphash;
    char path[100];
    path[0] = '\0';
    while ((dir_entry = readdir(dir)) != 0){
	if (strcmp(dir_entry->d_name,".") && strcmp(dir_entry->d_name,"..")){
	    strcat(path, dir_name);
	    strcat(path, "/");
	    strcat(path, dir_entry->d_name);
	    if (ph_dct_imagehash(path, tmphash) < 0)  //calculate the hash
		continue;
            dp = ph_malloc_imagepoint();              //store in structure with file name
	    dp->id = dir_entry->d_name;
	    dp->hash = tmphash;
	    hashlist1.push_back(*dp);
	    i++;
	}
	errno = 0;
        path[0]='\0';
    }

    if (errno){
	printf("error reading directory\n");
	exit(1);
    }

    sort(hashlist1.begin(),hashlist1.end(),cmp_lt_imp);

    //second directory
    dir_entry = NULL;
    DIR *dir2 = opendir(dir_name2);
    if (!dir){
	printf("unable to open directory\n");
	exit(1);
    }
    errno = 0;
    path[0] = '\0';
    i=0;
    while ((dir_entry = readdir(dir2)) != 0){
	if (strcmp(dir_entry->d_name,".") && strcmp(dir_entry->d_name,"..")){
	    strcat(path,dir_name2);
	    strcat(path,"/");
	    strcat(path,dir_entry->d_name);
	    if (ph_dct_imagehash(path,tmphash) < 0)    //calculate the hash
		continue;
    	    dp = ph_malloc_imagepoint();               //store in structure with filename
	    dp->id = dir_entry->d_name;
	    dp->hash = tmphash;
	    hashlist2.push_back(*dp);
	    i++;
	}
	errno = 0;
	path[0] = '\0';
    }

    if (errno){
	printf("error reading directory\n");
	exit(1);
    }

    sort(hashlist2.begin(),hashlist2.end(),cmp_lt_imp);

    int nbfiles1 = hashlist1.size();
    int nbfiles2 = hashlist2.size();
    int nbfiles = nbfiles1;
    if (nbfiles1 != nbfiles2){
	nbfiles = (nbfiles2 > nbfiles1) ? nbfiles2:nbfiles1;
    }
    
    int distance = -1;
    printf("**************************\n");
    printf("similar images distance comparisons\n");
    printf("**************************\n");
    for (i=0;i<nbfiles;i++){
	printf(" %d %s %s ",i,hashlist1[i].id, hashlist2[i].id);

	//calculate distance
	distance = ph_hamming_distance(hashlist1[i].hash,hashlist2[i].hash);

	printf(" dist = %d\n",distance);
    }


    printf("**************************\n");
    printf("different distance comparisons\n");
    printf("**************************\n");
    for (i=0;i<nbfiles1;i++){
	for (int j=i+1;j<nbfiles1;j++){
	    printf(" %s %s ", hashlist1[i].id, hashlist1[j].id);

	    //calculate distance
	    distance = ph_hamming_distance(hashlist1[i].hash,hashlist1[j].hash);

	    printf(" dist = %d\n",distance);
	}
    }
    printf("**************************\n");


    return 0;
}
