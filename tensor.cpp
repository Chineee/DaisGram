#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */
#define EPSILON 0.000001f  /* the rounding precision for comparing floats */

using namespace std;

Tensor::Tensor() {
    data = nullptr;
    r = 0, c = 0, d = 0;
}

Tensor::Tensor(int r, int c, int d, float v) {

    this->init(r, c, d, v);
}


Tensor::Tensor(const Tensor &that){
    *this = that; 
}

Tensor::~Tensor() {
    // deallocco data[i][j] (che punta all'array che corrisponde alla profondità)

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            delete[] data[i][j];
        }
    }

    //ora deallocco data[i] (che punta all'array che corrisponde alle colonne)

    for (int i = 0; i < r; i++) {
        delete[] data[i];
    }

    //infine dealloco data (che corrisponde all'array delle righe)

    delete[] data;
    this->r = 0;
    this->c = 0;
    this->d = 0;
}

bool Tensor::operator==(const Tensor& rhs) const{
    if (r != rhs.r || c != rhs.c || d != rhs.d)
        throw dimension_mismatch();
    
    bool out = true;
    int z = 0, i, j;

    while (z < d && out){           //Per ogni Canale
        i = 0;
        while (i < r && out){       //Per ogni Riga
            j = 0;
            while (j < c && out){   //Per ogni Colonna
                out &= std::abs( (*this)(i,j,z) - rhs(i,j,z) ) < EPSILON;
                 if (!out) cout << rhs(i,j,z) << " " << (*this)(i,j,z) << " IN POSIZIONE " << i << " " << j << " " << z << endl;
                j++; 
            }
            i++;  
        }
        z++;
    }
    return out;
}

Tensor Tensor::operator-(const Tensor &rhs){
    if(r == rhs.r && c == rhs.c && d == rhs.d) {
        Tensor new_Tensor (r, c, d);

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < d; k++) {
                    new_Tensor.data[i][j][k] = this->data[i][j][k] - rhs.data[i][j][k];
                }
            }
        }

        return new_Tensor;
    }else
        throw dimension_mismatch();
}

Tensor Tensor::operator+(const Tensor &rhs){
    if(r == rhs.r && c == rhs.c && d == rhs.d) {
        Tensor new_Tensor (r, c, d);
        
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < d; k++) {
                    new_Tensor.data[i][j][k] = this->data[i][j][k] + rhs.data[i][j][k];
                }
            }
        }

        return new_Tensor;
    }else
        throw dimension_mismatch();
}

Tensor Tensor::operator*(const Tensor& rhs) {

    if (this->r == rhs.r && this->c == rhs.c && this->d == rhs.d) {
        // Creo il tensore da ritornare in modo da non modificare l'oggetto this
        Tensor new_Tensor(this->r, this->c, this->d);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < d; k++) {
                    new_Tensor.data[i][j][k] = this->data[i][j][k] * rhs.data[i][j][k];
                }
            }
        }
        // ritorno il tensore nuovo
        return new_Tensor;
    } else {
        throw(dimension_mismatch());
    }
}

Tensor Tensor::operator/(const Tensor& rhs) {
    if (this->r == rhs.r && this->c == rhs.c && this->d == rhs.d) {
        // Creo il tensore da ritornare in modo da non modificare l'oggetto this
        Tensor new_Tensor(this->r, this->c, this->d);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < d; k++) {
                    if (rhs.data[i][j][k] == 0) throw(unknown_exception());
                    new_Tensor.data[i][j][k] = this->data[i][j][k] / rhs.data[i][j][k];
                }
            }
        }
        // ritorno il tensore nuovo
        return new_Tensor;

    } else {
        throw(dimension_mismatch());
    }
}

Tensor Tensor::operator-(const float& rhs) {

    Tensor new_Tensor(r, c, d);

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                new_Tensor.data[i][j][k] = data[i][j][k] - rhs;
            }
        }
    }
    
    return new_Tensor;
}

Tensor Tensor::operator+(const float& rhs) {

    Tensor new_Tensor(r, c, d);

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                new_Tensor.data[i][j][k] = data[i][j][k] + rhs;
            }
        }
    }
    
    return new_Tensor;
}

Tensor Tensor::operator*(const float& rhs) const{

    Tensor new_Tensor(this->r, this->c, this->d);

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                new_Tensor.data[i][j][k] = this->data[i][j][k] * rhs;
            }
        }
    }
    
    return new_Tensor;
}

Tensor Tensor::operator/(const float& rhs) {
    Tensor new_Tensor(this->r, this->c, this->d);

    if (rhs == 0) throw (unknown_exception());

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                new_Tensor.data[i][j][k] = this->data[i][j][k] / rhs;
            }
        }
    }
    
    return new_Tensor;
}


const float Tensor::operator()(int i, int j, int k) const{

    if  ( i < 0 || j < 0 || k < 0 ||
          i >= r || j >= c || k >= d)
        throw index_out_of_bound();

    return data[i][j][k];

}

float& Tensor::operator()(int i, int j, int k){

    if  ( i < 0 || j < 0 || k < 0 ||
          i >= r || j >= c || k >= d)
        throw index_out_of_bound();

    return data[i][j][k];

}

Tensor &Tensor::operator=(const Tensor &sother){

    if (this == &sother) //Check Self Assignment
        return *this;
    
    if (this->r != sother.r || this->c != sother.c || this->d != sother.d) {
        //Delete "this"
        if (r != 0 && c != 0 && d != 0)
            this->~Tensor();
        //Create new
        this->init(sother.r,sother.c,sother.d);
    }

    //Assignment
    for (int z = 0; z < d; z++)         //canale
        for (int i = 0; i < r; i++)     //coord Row
            for (int j = 0; j < c; j++) //coord Column
                (*this)(i,j,z) = sother.data[i][j][z];

    return *this;   
}

void Tensor::init_random(float mean, float std){
    if(data){

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean,std);

        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    this->operator()(i,j,k)= distribution(generator);
                }
            }
        }    

    }else{
        throw(tensor_not_initialized());
    }
}

void Tensor::init(int r, int c, int d, float v){

    this->r = r;
    this->c = c;
    this->d = d;

    if (r == 0 || c == 0 || d == 0)
        throw(dimension_mismatch());
        
    data = new float**[r];

    // create the third dimension matrix

    for (int i = 0; i < r; i++) {
        data[i] = new float*[c];
    }

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            data[i][j] = new float[d];
        }
    }

    // initialized the value to the matrix

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                data[i][j][k] = v;
            }
        }
    }
}

void Tensor::clamp(float low, float high){
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            for (int k = 0; k < d; k++) {
                if (operator()(i, j, k) < low)
                    operator()(i, j, k) = low;
                if (operator()(i, j, k) > high)
                    operator()(i, j, k) = high;
            }
        }
    }
}

void Tensor::rescale(float new_max){
    //Controllo Errore
    if (!data)
        throw tensor_not_initialized();
    
    for (int z = 0; z < d; z++){
        float min = getMin(z);
        float max = getMax(z);
        for (int i = 0; i < r; i++){
            for (int j = 0; j < c; j++){
                if(min == max)
                    (*this)(i,j,z) = new_max;
                else
                    (*this)(i,j,z) = ( ((*this)(i,j,z) - min) / (max - min) ) * new_max;
            }
        }
    }
      
}

Tensor Tensor::padding(int pad_h, int pad_w){  
    int rows = this->r + 2*pad_h;
    int cols = this->c + 2*pad_w;

    Tensor new_Tensor(rows, cols, this->d);

    for (int i = pad_h; i < new_Tensor.rows()-pad_h; i++) {
        for (int j = pad_w; j < new_Tensor.cols()-pad_w; j++) {
            for (int k = 0; k < new_Tensor.depth(); k++) {
                new_Tensor.data[i][j][k] = this->data[i-pad_h][j-pad_w][k];
            }
        }
    }
    return new_Tensor;
}

Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end){
    Tensor subsetted_tensor (row_end - row_start, col_end - col_start, depth_end - depth_start);

    for (unsigned int i = row_start; i < row_end; i++) {
        for (unsigned int j = col_start; j < col_end; j++) {
            for (unsigned int k = depth_start; j < depth_end; j++) {
                subsetted_tensor.operator()(i, j, k) = this -> data[i][j][k];
            }
        }
    }

    return subsetted_tensor;
}

Tensor Tensor::concat(const Tensor &rhs, int axis){
    if ((axis == 0 && (c != rhs.c || d != rhs.d)) ||
        (axis == 1 && (r != rhs.r || d != rhs.d)) ||
        (axis == 2 && (r != rhs.r || c != rhs.c)) ||
        (axis > 2))
        throw concat_wrong_dimension();

    Tensor newt;
    int r_, c_, d_;
    r_ = r + (rhs.r * (axis == 0));
    c_ = c + (rhs.c * (axis == 1));
    d_ = d + (rhs.d * (axis == 2));
    newt.init(r_, c_, d_);

    for (int z = 0; z < d_; z++){
        for (int i = 0; i < r_; i++){
            for (int j = 0; j < c_; j++){
                if (i >= r && axis == 0)
                    newt(i,j,z) = rhs(i-r,j,z);

                else if (j >= c && axis == 1)
                    newt(i,j,z) = rhs(i,j-c,z);

                else if (z >= d && axis == 2)
                    newt(i,j,z) = rhs(i,j,z-d);

                else
                    newt(i,j,z) = (*this)(i,j,z);
            }
        }     
    } 
    return newt;
}

Tensor Tensor::convolve(const Tensor &f){
    
    if (f.r % 2 != 1 || f.c % 2 != 1 || f.r != f.c || f.d != this->depth()) {
        throw (filter_odd_dimensions());
    } else {
    
        int pad_h = (f.rows() - 1) / 2;
        int pad_w = (f.cols() - 1) / 2;
        
        Tensor padded_img = this->padding(pad_h, pad_w);  //l'immagine paddata invece ha dimensioni maggiori della immagine da modificare, bordi posti a 0 e resto uguale
        
        Tensor convolved_img(this->r, this->c, this->d); //l'immagine convolved ha le stesse dimensioni della immagine da modificare

        int number_of_trasl_x = padded_img.rows() - f.rows();
        int number_of_trasl_y = padded_img.cols() - f.cols();

        for (int depth = 0; depth < padded_img.depth(); depth++) {
            for (int i = 0; i <= number_of_trasl_x; i++) {
                for (int j = 0; j <= number_of_trasl_y; j++) {
                    for (int k = i; k < i + f.rows(); k++) {
                        for (int x = j; x < j + f.cols(); x++) {
                            convolved_img(i,j,depth) += padded_img(k,x,depth) * f(k-i, x-j, depth);
                        }
                    }
                }
            }
        }
    
    return convolved_img;

    }
}

/* UTILITY */

int Tensor::rows() const{
    return r;
}

int Tensor::cols() const{
    return c;
}

int Tensor::depth() const{
    return d;
}

float Tensor::getMin(int k){
    float min = FLT_MAX;

    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            if (data[i][j][k] < min)
                min = data[i][j][k];

    return min;
}

float Tensor::getMax(int k){
    float max = 0; 
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            if (data[i][j][k] > max)
                max = data[i][j][k];

    return max;
}

void Tensor::showSize(){
    /* The format is the following:
    "rows" x "colums" x "depth" */
    cout << r << " x " << c << " x " << d;
}

/* IOSTREAM */

ostream& operator<< (ostream& stream, const Tensor& obj){
    //decladed as friend in the header file

    //Print Verticale
    for (int z = 0; z < obj.d; z++){
        stream <<std::endl<< "Depth: " << z << std::endl;
        for (int i = 0; i < obj.r; i++){
            stream << " | ";
            for (int j = 0; j < obj.c; j++){
                if (obj(i,j,z) < 10) cout << "00";
                else if (obj(i,j,z) < 100) cout << "0";
                stream << obj(i,j,z)<<" | ";
            }
            stream << std::endl;
        }
    }
    return stream;
}

void Tensor::read_file(string filename){
    ifstream ifile;
    string s{""};

    ifile.open(filename);
    if(!ifile) 
        throw unable_to_read_file();

    int r_, c_, d_;
    ifile >> r_ >> c_ >> d_;
    if(r != r_ || c != c_ || d != d_){
        if(data)
            this -> ~Tensor();
        this -> init(r_, c_, d_);
    }  

    for (int k = 0; k < this ->  d; k++) {
        for (int i = 0; i < this -> r; i++) {
            for (int j = 0; j < this -> c; j++) {
                ifile >> s;
                this -> data[i][j][k] = stof(s);
            }
        }
    }

    ifile.close();
}

void Tensor::write_file(string filename){
    ofstream myFile;
    myFile.open(filename);
    myFile << this->r << "\n" << this->c << "\n" << this->d << "\n" << flush;
    for (int k = 0; k < d; k++) {
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                myFile << this->data[i][j][k] << "\n" << flush;
            }
        }
    }
    myFile.close();
}
