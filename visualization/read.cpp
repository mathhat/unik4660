#include "fstream"      /*writes to file*/
#include "hdf5.h"       /* reads from .h5 files*/
#include "iostream"     /* cout, endl  */
#include <stdlib.h>     /* srand, rand */
#include <string>        //stringsss
#include <string.h>      //memcpy
#include <math.h>       /* sqrt */
#define FILE            "isabel_2d.h5"
#define DATASET         "Velocity/X-comp"
#define DATASET2        "Velocity/Y-comp"
#define DIM0            500
#define DIM1            500
using namespace std;

void pgm(const unsigned char *noise, int width, int height,int L);
void makenoise(unsigned char *noise, int width, int height);
void writer(const double* rdata,const double* rdata2,int noise_scale);
void streamwriter(const double* linex,const double* liney, int L);
void init(hid_t file, hid_t dset,hid_t dset2, herr_t status,herr_t status2,double* rdata, double* rdata2);
void euler(const double* rdata,const double* rdata2,double *streamlinex,double *streamliney, int x0, int y0,int L,int width, int height);
void rk(const double* rdata,const double* rdata2,double *streamlinex,double *streamliney, int x0, int y0,int L,int width, int height);
void interpolate(const double* rdata,const double* rdata2, double* x, double* y,int width, int height);
double box_convolution(const double *streamlinex,const double *streamliney, const unsigned char *noise,double p, int L,int width);

int
main (void)
{
    bool write = true;
    int i, j;
    hid_t       file, space, dset,dset2;
    herr_t      status,status2;

    int noise_scale = 0; //how many times do you wanna interpolate?
    noise_scale = pow(2,noise_scale); //each interpolation quadruples the pixels
    int width = noise_scale*DIM0;
    int height = noise_scale*DIM1;


    double *x;
    x = new double[height*width];    
    double *y;
    y = new double[height*width];


    double *rdata;
    rdata = new double[DIM0*DIM1];
    double *rdata2;
    rdata2 = new double[DIM0*DIM1];

    
    unsigned char *noise;
    noise = new unsigned char[height*width];
    unsigned char *image;
    image = new unsigned char[height*width];
    
    makenoise(noise, width, height);
    
    //HERE WE GATHER VECTOR FIELD DATA INTO RDATA and RDATA2 (and x and y)

    //special case of the isabel set

    if (FILE == "isabel_2d.h5"){

        init(file,dset,dset2,status,status2,y,x); //YOU INVERTED THEM
        for ( i = 0; i <DIM0; i++){
            for ( j = 0; j<DIM1 ; j++) {
                rdata[i+j*DIM0] = x[i+j*DIM0];

            }
        }        
        for ( i = 0; i<DIM0 ; i++) {
            for ( j = 0; j<DIM1 ; j++) {
                rdata2[i+j*DIM0] = y[i+j*DIM0];

            }  
        }

    }
    else{
        init(file,dset,dset2,status,status2,rdata,rdata2); 
    }
    //interpolation was a bit tricky since I'm used to flexible array dimensions
    int interp = 0;
    if (noise_scale > 1){
        interpolate(rdata,rdata2,x,y,2*DIM0,2*DIM1);
        for(interp = 2; interp < noise_scale;){
            delete[] rdata;
            delete[] rdata2;
            rdata = new double[interp*interp*DIM0*DIM1];
            rdata2 = new double[interp*interp*DIM0*DIM1];
            memcpy(rdata,x,interp*interp*DIM0*DIM1*sizeof(double));
            memcpy(rdata2,y,interp*interp*DIM0*DIM1*sizeof(double));
            delete[] x;
            delete[] y;
            interp*=2;
            x = new double[interp*interp*DIM0*DIM1];
            y = new double[interp*interp*DIM0*DIM1];
            interpolate(rdata,rdata2,x,y,interp*DIM0,interp*DIM1);

        }
    }


    //HERE WE INTIALIZE AND INTEGRATE THE STREAMLINE
    int x0,y0;
    int L = 30;        //half of the streamline's length
    double streamlinex[2*L+1]; //the streamline's x coordinates
    double streamliney[2*L+1]; //the streamline's y coordinates
    double p;
    double tol = 0.0002;
    int index;

    
    for (j = 0; j < height; j+=1){
        for (i = 0; i < width; i+=1){
        p = 0;

        x0 = i;
        y0 = j;
        index = x0+y0*width;
        if (sqrt(x[index]*x[index]+y[index]*y[index]) > tol){
            euler(x,y,streamlinex,streamliney,x0,y0,L,width,height);
            image[i + j*width] = static_cast<unsigned char>(box_convolution(streamlinex,streamliney,noise,p,L,width));
        
        }
        else{
            image[i + j*width] = 0;
        }
        }
    }
    
    

    //streamline writer  
    /*
    ofstream sfile;
    sfile.open ("lines.txt");
     
    for (j = 0; j < height; j+=height/100){
        for (i = 0; i < width; i+=width/100){
        p = 0;
        x0 = i;
        y0 = j;
        index = x0+y0*width;
        euler(x,y,streamlinex,streamliney,x0,y0,L,width,height);
        streamwriter(streamlinex,streamliney,L);
        }
    }
    
    sfile.close();
    */
    pgm(image, width, height,L);
    
    
    /*if (write){
        writer(x ,y,noise_scale);
    
    }
    */
    delete[] rdata;
    delete[] rdata2;
    delete[] noise;
    delete[] image;
    delete[] x;
    delete[] y;
    return 0;
}

//  END of main

double 
box_convolution(const double *streamlinex,const double *streamliney, const unsigned char *noise,double p, int L,int width){
    for(int i = 0; i<2*L+1;++i){
        p += noise[int(round(streamlinex[i])) +int(round(streamliney[i]))*width];
    }

    return p/((2*L+1));
}

void 
interpolate(const double* rdata,const double* rdata2, double* x, double* y,int width, int height){
    for (int i = 0; i < width; i+=2){
        for (int j = 0; j < height; j+=2){
            x[i + j*width] = rdata[i/2 + j/2*width/2];
            y[i + j*width] = rdata2[i/2 + j/2*width/2];
        }
    }
    for (int i = 0; i < width-1; i+=2){ //might need -2 or -0 in bool instead of -1
        for (int j = 0; j < height-1; j+=2){
            x[i + 1 + (j)*width]  = (x[i + (j)*width] + x[i+2 + (j)*width])/2;
            x[i + (j+1)*width]    = (x[i + (j)*width] + x[i + (j+2)*width])/2;
            y[i + 1 + (j)*width]  = (y[i + (j)*width] + y[i+2 + (j)*width])/2;
            y[i + (j+1)*width]    = (y[i + (j)*width] + y[i + (j+2)*width])/2;
        }
    }

    for (int i = 0; i < width; i+=2){
        for (int j = 0; j < height; j+=2){
            x[i + 1 + (j+1)*width]= (x[i + 1 + (j)*width] + x[i + 1 + (j+2)*width])/2;
            y[i + 1 + (j+1)*width]= (y[i + 1 + (j)*width] + y[i + 1 + (j+2)*width])/2;
        }
    }
}

void 
makenoise(unsigned char *noise, int width, int height)
{
    int maxColorValue = 255;
    for(int i=0;i<height*width;++i)
        noise[i] = round(float(rand())/RAND_MAX *maxColorValue);
    //if(wannaFlush)
    //    f << std::flush;
} // block scope closes file, which flushes anyway.


void 
pgm(const unsigned char *image, int width, int height,int L)
{
    string l = "isabel_linelength_"+to_string(L)+".pgm";
    cout <<"watup"<<endl;
    ofstream f(l,    ios_base::out
                              |ios_base::binary
                              |ios_base::trunc
                   );

    int maxColorValue = 255;
    f << "P5\n" << width << " " << height << "\n" << maxColorValue << "\n";

    for(int j=0;j<height;++j)
    {
        for (int i = 0; i<width;++i)
        {    
        f<< image[i + j*width] ;
        }
    }
    f.close();

} // block scope closes file, which flushes anyway.


void 
writer(const double* rdata,const double* rdata2,int noise_scale){
    ofstream xfile;
    ofstream yfile;
    xfile.open ("x.txt");
    yfile.open ("y.txt");
    for (int j = 0; j < noise_scale*DIM1;++j){
        for (int i = 0; i < noise_scale*DIM0; ++i){
            yfile << rdata2[i + j*noise_scale*DIM0]<<endl;
            xfile << rdata[i + j*noise_scale*DIM0]<<endl;
        }
    } 
    xfile.close();
    yfile.close();

}

void 
streamwriter(const double* linex,const double* liney, int L){
    ofstream sfile;
    int length = 2*L+1;
    sfile.open ("lines.txt",::fstream::app);
    for (int j = 0; j < length;++j){
        sfile << linex[j] <<" " <<liney[j]<<endl;
        } 
}

void 
init(hid_t file, hid_t dset,hid_t dset2, herr_t status,herr_t status2,double* rdata,double* rdata2){
        file = H5Fopen (FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
        dset = H5Dopen (file, DATASET, H5P_DEFAULT);
        dset2 = H5Dopen (file, DATASET2, H5P_DEFAULT);

        status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    rdata);
        status2 = H5Dread (dset2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    rdata2);
    }

void 
euler(const double* rdata,const double* rdata2,double *streamlinex,double *streamliney, int x0, int y0,int L,int width, int height){
    double tol = 0.0002;
    double timestep = 1e-0;
    streamlinex[L] = x0;
    streamliney[L] = y0;
    int xprev = x0;
    int yprev = y0;
    double x = x0;
    double y = y0;
    int index;
    double v_mag;
    for (int i = 1;i < L+1;++i){
        xprev = round(x);
        yprev = round(y);
        index = xprev + yprev*width;
        if (index > width*height){
            cout << "out of bounds forward, index:"<< yprev<<" "<< xprev<<endl;
            streamlinex[L+i] = streamlinex[L+i-1];
            streamliney[L+i] = streamliney[L+i-1];
            continue;
            }
        v_mag = sqrt(rdata[index]*rdata[index] 
                            + rdata2[index]*rdata2[index]);

        if (((xprev < width-1) && (xprev >0)) && ((yprev < height-1) && (yprev >0)) &&(v_mag > tol)){

            x = streamlinex[L+i-1] + (rdata[index]/v_mag)*timestep;
            y = streamliney[L+i-1] + (rdata2[index]/v_mag)*timestep;
            streamlinex[L+i] = x;
            streamliney[L+i] = y;}
        else{

            streamlinex[L+i] = streamlinex[L+i-1];
            streamliney[L+i] = streamliney[L+i-1];}            
        }

    x = x0;
    y = y0;
    
    for (int i = L;i > 0;--i){
        xprev = round(x);
        yprev = round(y);

        index = xprev + yprev*width;
        if (index > width*height){
            cout << "out of bounds back, index:"<< yprev<<" "<< xprev<<endl;
            streamlinex[i-1] = streamlinex[i];
            streamliney[i-1] = streamliney[i];
            continue;
            }
        v_mag = sqrt(rdata[index]*rdata[index] 
                            + rdata2[index]*rdata2[index]);

        if (((xprev < width-1) && (xprev >0)) && ((yprev < height-1) && (yprev >0)) &&(v_mag > tol)){

            x = streamlinex[i] - rdata[index]/v_mag*timestep;
            y = streamliney[i] - rdata2[index]/v_mag*timestep;

            streamlinex[i-1] = x;
            streamliney[i-1] = y;
        }
        else{

            streamlinex[i-1] = streamlinex[i];
            streamliney[i-1] = streamliney[i];}            
    }
}


void rk(const double* rdata,const double* rdata2,double *streamlinex,double *streamliney, int x0, int y0,int L,int width, int height){
    streamlinex[L] = x0;
    streamliney[L] = y0;

    double h = 1e-0;
    double tol = 0.0002;
    int xprev = x0;
    int yprev = y0;

    double x = x0;
    double y = y0;
    double x_;
    double y_; 
    double k1x;
    double k1y;
    double k2x;
    double k2y;
    double k3x;
    double k3y;
    double k4x;
    double k4y;
    double kx;
    double ky;
    
    double v_mag_k1;
    double v_mag_k2;
    double v_mag_k3;
    double v_mag_k4;
    int index;
    for (int i = 1;i < L+1;++i){
        xprev = round(x);
        yprev = round(y);
        k1x = rdata[xprev + yprev*width];
        k1y = rdata2[xprev + yprev*width];
        v_mag_k1 = sqrt(k1x*k1x+k1y*k1y);
        if (((xprev < width) && (xprev >0)) && ((yprev < height) && (yprev >0)) && (v_mag_k1 > tol)){
            k1x /= v_mag_k1;
            k1y /= v_mag_k1;
            x_ = x + k1x*h/2;
            y_ = y + k1y*h/2;
            xprev = round(x_); //rounding off instead of triangulating velocities
            yprev = round(y_); //cause anything else sounds horrible


            k2x = rdata[xprev + yprev*width];
            k2y = rdata2[xprev + yprev*width];
            v_mag_k2 = sqrt(k2x*k2x+k2y*k2y); 

            if (v_mag_k2 > tol){
                k2x /= v_mag_k2;
                k2y /= v_mag_k2;}
            else{
                streamlinex[L+i] = streamlinex[L+i-1];
                streamliney[L+i] = streamliney[L+i-1];  

                continue;
            }
            x_ = x + k2x*h/2;
            y_ = y + k2y*h/2;
            xprev = round(x_);
            yprev = round(y_); 

            k3x = rdata[xprev + yprev*width];
            k3y = rdata2[xprev + yprev*width];

            v_mag_k3 = sqrt(k3x*k3x+k3y*k3y); 
            if (v_mag_k3 > tol){
                k3x /= v_mag_k3;
                k3y /= v_mag_k3;}
            else{
                streamlinex[L+i] = streamlinex[L+i-1];
                streamliney[L+i] = streamliney[L+i-1];  
                continue;
            }

            x_ = x + k3x*h;
            y_ = y + k3y*h;
            xprev = round(x_);
            yprev = round(y_); 
            k4x = rdata[xprev + yprev*width];
            k4y = rdata2[xprev + yprev*width];
            v_mag_k4 = sqrt(k4x*k4x+k4y*k4y); 
            if (v_mag_k4 > tol){
                k4x /= v_mag_k4;
                k4y /= v_mag_k4;}
            else{
                streamlinex[L+i] = streamlinex[L+i-1];
                streamliney[L+i] = streamliney[L+i-1];  
                continue;
            }
            
            kx = double(1)/6*(k1x+2*k2x+2*k3x+k4x);
            ky = double(1)/6*(k1y+2*k2y+2*k3y+k4y);

            x +=  kx*h;
            y +=  ky*h;

            xprev = round(x);

            yprev = round(y);

            streamlinex[L+i] = x;
            streamliney[L+i] = y;}
        else{
            streamlinex[L+i] = streamlinex[L+i-1];
            streamliney[L+i] = streamliney[L+i-1];}            
    }



    xprev = x0; x = x0;
    yprev = y0; y = y0;
    //backward integrate?
    for (int i = L;i > 0;--i){
        xprev = round(x);
        yprev = round(y);

        k1x = rdata[xprev + yprev*width];
        k1y = rdata2[xprev + yprev*width];
        v_mag_k1 = sqrt(k1x*k1x+k1y*k1y); 

        if (((xprev < width) && (xprev >0)) && ((yprev < height) && (yprev >0)) && (v_mag_k1 > tol)){
                                
            k1x /= v_mag_k1;
            k1y /= v_mag_k1;

            x_ = x - k1x*h/2;
            y_ = y - k1y*h/2;

            xprev = round(x_); //rounding off instead of triangulating velocities
            yprev = round(y_); //cause anything else sounds horrible
            k2x = rdata[xprev + yprev*width];
            k2y = rdata2[xprev + yprev*width];
            v_mag_k2 = sqrt(k2x*k2x+k2y*k2y); 
            if (v_mag_k2 > tol){
                k2x /= v_mag_k2;
                k2y /= v_mag_k2;}
            else{
                streamlinex[i-1] = streamlinex[i];
                streamliney[i-1] = streamliney[i];  
                continue;
            }

            x_ = x - k2x*h/2;
            y_ = y - k2y*h/2;
            xprev = round(x_);
            yprev = round(y_); 
            k3x = rdata[xprev + yprev*width];
            k3y = rdata2[xprev + yprev*width];
            v_mag_k3 = sqrt(k3x*k3x+k3y*k3y); 
            if (v_mag_k3 > tol){
                k3x /= v_mag_k3;
                k3y /= v_mag_k3;}
            else{
                streamlinex[i-1] = streamlinex[i];
                streamliney[i-1] = streamliney[i];  
                continue;
            }

            x_ = x - k3x*h;
            y_ = y - k3y*h;
            xprev = round(x_);
            yprev = round(y_); 
            k4x = rdata[xprev + yprev*width];
            k4y = rdata2[xprev + yprev*width];
            v_mag_k4 = sqrt(k4x*k4x+k4y*k4y);
 
            if (v_mag_k4 > tol){
                k4x /= v_mag_k4;
                k4y /= v_mag_k4;}
            else{
                streamlinex[i-1] = streamlinex[i];
                streamliney[i-1] = streamliney[i];  
                continue;
            }

            kx = double(1)/6*(k1x+2*k2x+2*k3x+k4x);
            ky = double(1)/6*(k1y+2*k2y+2*k3y+k4y);

            x -=  kx*h;
            y -=  ky*h;

            xprev = round(x);
            yprev = round(y);
            streamlinex[i-1] = x;
            streamliney[i-1] = y;

            }

        else{

            streamlinex[i-1] = streamlinex[i];
            streamliney[i-1] = streamliney[i];}            
        }
}