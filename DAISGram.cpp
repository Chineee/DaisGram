#include <iostream>
#include <string>

#include "dais_exc.h"
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

using namespace std;

DAISGram::DAISGram() {

}

DAISGram::~DAISGram() {

}

void DAISGram::load_image(string filename){
    BmpImg img = BmpImg();

    img.read(filename.c_str());

    const int h = img.get_height();
    const int w = img.get_width();

    data = Tensor(h, w, 3, 0.0);

    for(int i=0;i<img.get_height();i++){
        for(int j=0;j<img.get_width();j++){ 
            data(i,j,0) = (float) img.red_at(j,i);
            data(i,j,1) = (float) img.green_at(j,i);    
            data(i,j,2) = (float) img.blue_at(j,i);   
        }                
    }
}

void DAISGram::save_image(string filename){

    data.clamp(0,255);

    BmpImg img = BmpImg(getCols(), getRows());

    img.init(getCols(), getRows());

    for(int i=0;i<getRows();i++){
        for(int j=0;j<getCols();j++){
            img.set_pixel(j,i,(unsigned char) data(i,j,0),(unsigned char) data(i,j,1),(unsigned char) data(i,j,2));                   
        }                
    }

    img.write(filename);
}

int DAISGram::getRows() const{
    return this->data.rows();
}

int DAISGram::getCols() const{
    return this->data.cols();
}

int DAISGram::getDepth() const{
    return this->data.depth();
}

void DAISGram::generate_random(int h, int w, int d){
    data = Tensor(h,w,d,0.0);
    data.init_random(128,50);
    data.rescale(255);
}

DAISGram DAISGram::brighten(float bright) {
    DAISGram new_img;

    new_img.data = this->data;
    new_img.data = new_img.data + bright;
    new_img.data.clamp(0, 255);
    
    return new_img;
}

DAISGram DAISGram::grayscale() {
    DAISGram new_img;
    
    new_img.data = this->data;
    float avg {0};

    for (int i = 0; i < getRows(); i++) {
        //I just accessed the row

        for(int j = 0; j < getCols(); j++) {
            //I just accessed the column

            avg = 0; //resetting the average

            for(int k = 0; k < getDepth(); k++) {
                //I just accessed the depth and now I'm going to sum the values
                avg += this->data.operator()(i, j, k);
            }

            //Calculating the average
            avg /= getDepth();

            //Setting the new value in new_img
            for(int k = 0; k < getDepth(); k++) {
                new_img.data.operator()(i, j, k) = avg;
            }
        }
    }

    return new_img;
}

DAISGram DAISGram::blend(const DAISGram & rhs, float alpha) {
    if (this->getRows() != rhs.getRows() || this->getCols() != rhs.getCols() || this->getDepth() != rhs.getDepth()) {
        throw (dimension_mismatch());
    } else {
        DAISGram new_img;
        new_img.data = Tensor(getRows(), getCols(), getDepth());
        new_img.data = this->data * alpha + rhs.data*(1-alpha);
        return new_img;
    }
}

DAISGram DAISGram::emboss() {
 
    DAISGram new_img;
    DAISGram filter;

    new_img.data = Tensor(this->data.rows(), this->data.cols(), this->data.depth());
    filter.data = Tensor(3, 3, 3); 
    
    for (int i = 0; i < filter.getRows(); i++) {
        for (int j = 0; j < filter.getCols(); j++) {
            if (i == 0 && j == 0) {
                filter.data.operator()(i, j, 0) = -2;
                filter.data.operator()(i, j, 1) = -2;
                filter.data.operator()(i, j, 2) = -2;
            }
            else if (i == 2 && j == 2) {
                filter.data.operator()(i, j, 0) = 2;
                filter.data.operator()(i, j, 1) = 2;
                filter.data.operator()(i, j, 2) = 2;
            }
            else if ( (i == 1 && j == 0) || (i == 0 && j == 1)) {
                filter.data.operator()(i, j, 0) = -1;
                filter.data.operator()(i, j, 1) = -1;
                filter.data.operator()(i, j, 2) = -1;
            }
            else if ( (i == 2 && j == 1) || (i == 1 && j == 2) || (i == j)) {
                filter.data.operator()(i, j, 0) = 1;
                filter.data.operator()(i, j, 1) = 1;
                filter.data.operator()(i, j, 2) = 1;
            }
            else {
                filter.data.operator()(i, j, 0) = 0;
                filter.data.operator()(i, j, 1) = 0;
                filter.data.operator()(i, j, 2) = 0;
            }
        }
    }

    new_img.data = this->data.convolve(filter.data);
    new_img.data.clamp(0, 255);

    return new_img;

}

DAISGram DAISGram:: smooth(int h){
  
    DAISGram new_image;
    DAISGram filter;

    new_image.data = Tensor (data.rows(), data.cols(), data.depth());
    filter.data = Tensor (h, h, 3);

    float c = 1.0f / (float)(h * h);

    for (int i = 0; i < filter.getRows(); i++) {
        for (int j = 0; j < filter.getCols(); j++) {
            filter.data.operator()(i, j, 0) = c;
            filter.data.operator()(i, j, 1) = c;
            filter.data.operator()(i, j, 2) = c;
        }
    }

    new_image.data = this->data.convolve(filter.data);

    return new_image;
}

DAISGram DAISGram::warhol(){
    DAISGram new_img;
    new_img.data = Tensor(this->data);                  //t1 top-left
    Tensor t2(data.rows(),data.cols(),data.depth());    //t2 top-right
    Tensor t3(data.rows(),data.cols(),data.depth());    //t3 bot-left
    Tensor t4(data.rows(),data.cols(),data.depth());    //t4 bot-right
    
    for (int z = 0; z < data.depth(); z++)
    {
        for (int i = 0; i < data.rows(); i++)
        {
            for (int j = 0; j < data.cols(); j++)
            {
                if (z==0){
                    t2(i,j,z) = data(i,j,1); //R~~G
                    t3(i,j,z) = data(i,j,z); //R~~R
                    t4(i,j,z) = data(i,j,2); //R~~B
                }else if (z==1){
                    t2(i,j,z) = data(i,j,0); //G~~R
                    t3(i,j,z) = data(i,j,2); //G~~B
                    t4(i,j,z) = data(i,j,z); //G~~G
                }else{
                    t2(i,j,z) = data(i,j,z); //B~~B
                    t3(i,j,z) = data(i,j,1); //B~~G
                    t4(i,j,z) = data(i,j,0); //B~~R
                }
            }   
        }
    }
    new_img.data = new_img.data.concat(t2,1);
    t3 = t3.concat(t4,1);
    new_img.data = new_img.data.concat(t3,0);

    return new_img; 
}

DAISGram DAISGram::sharpen(){
    /*
    *  filter[3][3]
    *  0  -1  0
    *  -1  5 -1
    *  0  -1  0
    */

    DAISGram new_img;
    DAISGram filter;

    new_img.data = Tensor(data.rows(), data.cols(), data.depth());
    filter.data = Tensor(3, 3, 3); 
    
    for (int i = 0; i < filter.getRows(); i++) {
        for (int j = 0; j < filter.getCols(); j++) {
            if (j % 2 == 0 && i % 2 == 0) {  //colonna e riga pari
                filter.data(i, j, 0) = 0;
                filter.data(i, j, 1) = 0;
                filter.data(i, j, 2) = 0;
            }
            else if (i == 1 && j == 1) {    //cella centrale
                filter.data(i, j, 0) = 5;
                filter.data(i, j, 1) = 5;
                filter.data(i, j, 2) = 5;
            }
            else {
                filter.data(i, j, 0) = -1;
                filter.data(i, j, 1) = -1;
                filter.data(i, j, 2) = -1; //altri
            }
        }
    }

    new_img.data = this->data.convolve(filter.data);
    new_img.data.clamp(0, 255);

    return new_img;
    
}

DAISGram DAISGram::edge(){
    
   
    DAISGram new_img;
    DAISGram gray;
    DAISGram filter;

    new_img.data = Tensor(this->data.rows(), this->data.cols(), this->data.depth());
    filter.data = Tensor(3, 3, 3); 
    
    for (int i = 0; i < filter.getRows(); i++) {
        for (int j = 0; j < filter.getCols(); j++) {
            if (j == 1 && i == 1) {  //cella centrale
                filter.data(i, j, 0) = 8;
                filter.data(i, j, 1) = 8;
                filter.data(i, j, 2) = 8;
            }
            else {
                filter.data(i, j, 0) = -1; 
                filter.data(i, j, 1) = -1; 
                filter.data(i, j, 2) = -1; 
            }//altri
        }
    }

    gray = this->grayscale(); //convert the image to grayscale before running the convolution.
    new_img.data = gray.data.convolve(filter.data);
    new_img.data.clamp(0, 255);

    return new_img;
    
}

/*
 *
 *      PARTI OPZIONALI
 *    
*/


DAISGram DAISGram::equalize() {

    DAISGram new_img;
    new_img.data.init(this->getRows(), this->getCols(), 3);

    float M = (float) this->getRows();
    float N = (float) this->getCols();

    float TOT_SIZE = (float) M*N;

    for (int k = 0; k < 3; k++) {
        float istogramma[256] = {};
        float cdf_min = data(0,0,0);

        for (int i = 0; i < getRows(); i++) {
            for (int j = 0; j < getCols(); j++) {
                float pos = this->data(i, j, k);
                istogramma[(int)pos] += 1.0f;
                if (pos < cdf_min)
                    cdf_min = this->data(i,j,k);
            }
        }
       
        cdf_min = istogramma[(int)cdf_min];
        float last_position = 0.0f;    
        for (unsigned int i = 0; i < 256; i++) {
            if (istogramma[i] != 0) {
                istogramma[i] += last_position; 
                last_position = istogramma[i];
                istogramma[i] =  roundf(( (istogramma[i] - cdf_min) / (TOT_SIZE - cdf_min) ) * (float) 255 ); // casting ad int come soluzione alternativa
                
            }
        }

        for (int i = 0; i < getRows(); i++) {
            for (int j = 0; j < getCols(); j++) {
                float pos = this->data(i, j, k);
                new_img.data(i,j,k) = istogramma[(int)pos];
            }
        }
    }

    return new_img;
}

DAISGram DAISGram::greenscreen(DAISGram & bkg, int rgb[], float threshold[]) {

    bool apply_green = true;

    // (rgb - threshold) <= pixel <= (rgb + threshold)
    DAISGram new_img;
    new_img.data = this->data;

    for (int i = 0; i < getRows(); i++) {
        for (int j = 0; j < getCols(); j++) {
            for (int k = 0; k < getDepth(); k++) {
                if ( (data(i, j, k) < (rgb[k] - threshold[k])) || (data(i, j, k) > (rgb[k] + threshold[k]))) {
                    apply_green = false;
                }
            }
            if (apply_green) {
                new_img.data(i, j, 0) = bkg.data(i, j, 0);
                new_img.data(i, j, 1) = bkg.data(i, j, 1);
                new_img.data(i, j, 2) = bkg.data(i, j, 2);
            }
            apply_green = true;
        }
    }
    return new_img;
}